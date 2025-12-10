#!/usr/bin/env python3
"""
Metadata Processing and Quality Control for CellProfiler Morphology Data

This script performs:
1. Metadata recoding and standardization
2. DMSO replicate selection via Wasserstein distance medoids per cell line
3. Optional one-hot encoding for downstream ML tasks

Editor: NZ
Date: 12/10/25
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import wasserstein_distance

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('metadata_processing.log')
    ]
)
logger = logging.getLogger(__name__)

def validate_required_columns(df: pd.DataFrame, required: List[str]) -> None:
    """Validate that required columns exist in dataframe."""
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
def fix_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply metadata corrections and standardization.
    
    Corrections:
    - Recode DMSO as NC (Negative Control) with concentration 0
    - Recode MS/MT (Mechanistic Standards) as PC (Positive Control)
    - Recode Clavulanate as NC
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with metadata columns
        
    Returns
    -------
    pd.DataFrame
        Dataframe with corrected metadata
    """
    df = df.copy()
    
    # Validate required columns
    required = ['Metadata_COND', 'Metadata_CONC', 'Metadata_CMPD']
    validate_required_columns(df, required)
    
    original_cond = df['Metadata_COND'].value_counts()
    logger.info(f"Original Metadata_COND distribution:\n{original_cond}")
    
    # 1. Recode DMSO
    dmso_mask = df['Metadata_CMPD'].str.upper() == 'DMSO'
    n_dmso = dmso_mask.sum()
    df.loc[dmso_mask, 'Metadata_COND'] = 'NC'
    df.loc[dmso_mask, 'Metadata_CONC'] = 0.0
    logger.info(f"Recoded {n_dmso} DMSO samples as NC with concentration 0")
    
    # 2. Recode MS/MT as PC
    ms_mt_mask = df['Metadata_COND'].isin(['MS', 'MT'])
    n_ms_mt = ms_mt_mask.sum()
    df.loc[ms_mt_mask, 'Metadata_COND'] = 'PC'
    logger.info(f"Recoded {n_ms_mt} MS/MT samples as PC")
    
    # 3. Recode Clavulanate as NC
    clav_mask = df['Metadata_CMPD'].str.contains('Clavulanate', na=False)
    n_clav = clav_mask.sum()
    df.loc[clav_mask, 'Metadata_COND'] = 'NC'
    logger.info(f"Recoded {n_clav} Clavulanate samples as NC")
    
    final_cond = df['Metadata_COND'].value_counts()
    logger.info(f"Final Metadata_COND distribution:\n{final_cond}")
    
    return df

def add_condition_encoding(df: pd.DataFrame, 
                           drop_original: bool = False) -> pd.DataFrame:
    """
    Add binary encoding for Metadata_COND for ML classification.
    
    Creates a single binary column 'Metadata_COND_binary' where:
    - 1 = PC (Positive Control)
    - 0 = NC (Negative Control) or other conditions
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    drop_original : bool, optional
        Whether to drop original Metadata_COND column
        
    Returns
    -------
    pd.DataFrame
        Dataframe with binary encoded condition column
    """
    df = df.copy()
    validate_required_columns(df, ['Metadata_COND'])
    
    # Create single binary column: 1 for PC, 0 for everything else NC
    df['Metadata_COND_binary'] = (df['Metadata_COND'] == 'PC').astype(int)
    
    # Insert after original column
    cond_idx = df.columns.get_loc('Metadata_COND')
    binary_col = df.pop('Metadata_COND_binary')
    df.insert(cond_idx + 1, 'Metadata_COND_binary', binary_col)
    
    if drop_original:
        df = df.drop(columns=['Metadata_COND'])
        logger.info(f"Added binary encoding (1=PC, 0=other) and dropped...")
    else:
        logger.info(f"Added binary encoding: Metadata_COND_binary (1=PC, 0=other)")
    
    return df

def get_feature_columns(df: pd.DataFrame) -> List[str]:
    """
    Extract feature columns (non-metadata numeric columns).
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
        
    Returns
    -------
    List[str]
        List of feature column names
    """
    # Identify metadata columns
    meta_pattern = 'Metadata|Location|Center|Object_Number'
    meta_cols = df.columns[df.columns.str.contains(meta_pattern, case=False, regex=True)]
    
    # Get numeric columns that aren't metadata
    feature_cols = (df.drop(columns=meta_cols)
                    .select_dtypes(include=[np.number])
                    .columns
                    .tolist())
    
    logger.info(f"Identified {len(feature_cols)} feature columns")
    return feature_cols

def compute_wasserstein_distances(profiles: np.ndarray) -> np.ndarray:
    """
    Compute pairwise Wasserstein distances between sample profiles.
    
    Uses vectorized operations for efficiency.
    
    Parameters
    ----------
    profiles : np.ndarray
        Array of shape (n_samples, n_features)
        
    Returns
    -------
    np.ndarray
        Distance matrix of shape (n_samples, n_samples)
    """
    n_samples = profiles.shape[0]
    distances = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            d = wasserstein_distance(profiles[i], profiles[j])
            distances[i, j] = d
            distances[j, i] = d
    
    return distances

def select_representative_medoids(df_subset: pd.DataFrame, 
                                   feature_cols: List[str],
                                   k: int = 3) -> pd.DataFrame:
    """
    Select k most representative samples using Wasserstein distance medoids.
    
    Selects samples that are most centrally located in the feature space,
    i.e., have minimal median distance to all other samples.
    
    Missing values are imputed using column-wise medians. Infinite values
    are replaced with the column min/max (for -inf/+inf respectively).
    
    Parameters
    ----------
    df_subset : pd.DataFrame
        Subset of samples to select from (e.g., DMSOs from one donor)
    feature_cols : List[str]
        Feature column names
    k : int, optional
        Number of medoids to select (default: 3)
        
    Returns
    -------
    pd.DataFrame
        Selected medoid samples with metadata
    """
    if len(df_subset) < k:
        logger.warning(
            f"Only {len(df_subset)} samples available, returning all instead of {k}"
        )
        return df_subset
    
    # Extract feature matrix
    X = df_subset[feature_cols].values.copy()  # Copy to avoid modifying original
    
    # Handle infinite values first (replace with column min/max)
    if np.any(np.isinf(X)):
        n_inf = np.isinf(X).sum()
        logger.warning(f"Detected {n_inf} infinite values, replacing with column min/max")
        
        for col_idx in range(X.shape[1]):
            col_data = X[:, col_idx]
            
            # Replace +inf with column max (excluding inf)
            if np.any(np.isposinf(col_data)):
                finite_max = np.max(col_data[np.isfinite(col_data)])
                X[np.isposinf(X[:, col_idx]), col_idx] = finite_max
            
            # Replace -inf with column min (excluding inf)
            if np.any(np.isneginf(col_data)):
                finite_min = np.min(col_data[np.isfinite(col_data)])
                X[np.isneginf(X[:, col_idx]), col_idx] = finite_min
    
    # Handle missing values with median imputation
    if np.any(np.isnan(X)):
        n_missing = np.isnan(X).sum()
        logger.warning(f"Detected {n_missing} missing values, performing median imputation")
        
        # Compute column-wise medians (ignoring NaN)
        col_medians = np.nanmedian(X, axis=0)
        
        # Find indices of NaN values and replace with column median
        nan_indices = np.where(np.isnan(X))
        X[nan_indices] = np.take(col_medians, nan_indices[1])
        
        logger.info(f"Imputed {n_missing} values using column medians")
    
    # Verify no remaining invalid values
    if np.any(np.isnan(X)) or np.any(np.isinf(X)):
        logger.error("Failed to clean all NaN/Inf values, cannot compute distances")
        return df_subset.iloc[:k]  # Return first k as fallback
    
    # Compute pairwise Wasserstein distances
    distances = compute_wasserstein_distances(X)
    
    # Find medoids (samples with minimum median distance to others)
    median_distances = np.median(distances, axis=1)
    medoid_indices = np.argsort(median_distances)[:k]
    
    selected = df_subset.iloc[medoid_indices].copy()
    selected['_median_wasserstein_distance'] = median_distances[medoid_indices]
    
    logger.info(
        f"Selected {k} medoids with median distances: "
        f"{median_distances[medoid_indices]}"
    )
    
    return selected

def filter_dmso_replicates(df: pd.DataFrame, 
                           feature_cols: List[str],
                           n_replicates: int = 3) -> pd.DataFrame:
    """
    Filter DMSO samples to retain top n representative replicates per donor.
    
    Uses Wasserstein distance medoids to select the most concordant replicates
    for each cell line/donor.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe (must have metadata fixes applied)
    feature_cols : List[str]
        Feature column names
    n_replicates : int, optional
        Number of replicates to retain per donor (default: 3)
        
    Returns
    -------
    pd.DataFrame
        Dataframe with filtered DMSO samples and all non-DMSO samples
    """
    validate_required_columns(df, ['Metadata_COND', 'Metadata_Donor'])
    
    # Separate DMSO (NC) and non-DMSO samples
    is_dmso = df['Metadata_CMPD'].str.upper() == 'DMSO'  # ✅ Only DMSO
    dmso_samples = df[is_dmso].copy()
    other_samples = df[~is_dmso].copy()
    
    logger.info(f"Processing {len(dmso_samples)} DMSO samples across donors")
    logger.info(f"Retaining all {len(other_samples)} non-DMSO samples")
    
    # Select representative DMSOs per donor
    selected_dmso_list = []
    for donor, group in dmso_samples.groupby('Metadata_Donor'):
        logger.info(f"Donor {donor}: {len(group)} DMSO samples")
        selected = select_representative_medoids(group, feature_cols, k=n_replicates)
        selected_dmso_list.append(selected)
    
    # Combine selected DMSOs with all other samples
    selected_dmso = pd.concat(selected_dmso_list, ignore_index=True)
    result = pd.concat([selected_dmso, other_samples], ignore_index=True)
    
    logger.info(
        f"Final dataset: {len(selected_dmso)} DMSO + "
        f"{len(other_samples)} other = {len(result)} total samples"
    )
    
    return result

def main():
    parser = argparse.ArgumentParser(
        description="Process metadata and perform QC filtering for morphology data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python 02_process_metadata_and_qc.py input.csv -o output.csv
  
  # With one-hot encoding
  python 02_process_metadata_and_qc.py input.csv -o output.csv --encode-condition
  
  # Custom DMSO replicate count
  python 02_process_metadata_and_qc.py input.csv -o output.csv --n-dmso-replicates 5
        """
    )
    
    parser.add_argument(
        'input',
        type=Path,
        help='Input CSV file (e.g., pct_processed.csv)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output CSV file path'
    )
    
    parser.add_argument(
        '--encode-condition',
        action='store_true',
        help='Add one-hot encoding for Metadata_COND'
    )
    
    parser.add_argument(
        '--n-dmso-replicates',
        type=int,
        default=3,
        help='Number of DMSO replicates to retain per donor (default: 3)'
    )
    
    parser.add_argument(
        '--skip-dmso-filter',
        action='store_true',
        help='Skip DMSO filtering step'
    )
    
    args = parser.parse_args()
    
    # Validate input
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Reading input from: {args.input}")
    df = pd.read_csv(args.input)
    logger.info(f"Loaded {len(df)} rows, {len(df.columns)} columns")
    
    # Step 1: Fix metadata
    logger.info("=" * 60)
    logger.info("STEP 1: Fixing metadata")
    logger.info("=" * 60)
    df = fix_metadata(df)
    
    # Step 2: Optional one-hot encoding
    if args.encode_condition:
        logger.info("=" * 60)
        logger.info("STEP 2: Adding condition encoding")
        logger.info("=" * 60)
        df = add_condition_encoding(df)
    
    # Step 3: DMSO filtering
    if not args.skip_dmso_filter:
        logger.info("=" * 60)
        logger.info("STEP 3: Filtering DMSO replicates")
        logger.info("=" * 60)
        feature_cols = get_feature_columns(df)
        df = filter_dmso_replicates(df, feature_cols, n_replicates=args.n_dmso_replicates)
    
    # Save output
    logger.info("=" * 60)
    logger.info(f"Saving output to: {args.output}")
    df.to_csv(args.output, index=False)
    logger.info(f"Final dataset: {len(df)} rows, {len(df.columns)} columns")
    logger.info("Processing complete!")


if __name__ == "__main__":
    main()