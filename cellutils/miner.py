import os
import re

import numpy as np
import pandas as pd

import pycytominer as pct
from pycytominer.cyto_utils.cells import SingleCells

from .utils import get_data_cols

NEG_CONTROL_QUERY = 'Metadata_Compound == DMSO"'

def plate_level(df, plate_id='Metadata_PlateID'):
    # plate aggregation and plate outlier correction
    meta_cols, data_cols = get_data_cols(df=df)

    # check if nan fills 
    
    aggregated = pct.aggregated(
        population_df=df,
        strata=meta_cols,
        features=data_cols,
        operation='median'        
    )

    dfs = []
    for pid in aggregated[plate_id].unique():
        temp = aggregated.loc[aggregated[plate_id]==pid]
        norm_plate = pct.normalize(
            profilers=temp,
            features=data_cols,
            meta_features=meta_cols,
            samples=NEG_CONTROL_QUERY # to change later for metadata reasons
        )
        dfs.append(norm_plate)
    plate_norm = pd.concat(dfs)
    return plate_norm


def feature_selection(df, data_cols, meta_cols, 
                      na_cutoff=0.05, freq_cut=0.05, unique_cut=0.01, outlier_cut=500,
                      corr_thres=0.90, corr_method="spearman", 
                      noise_perturb_group="Metadata_CMPD", nosise_std_cut=1.0):
    
    print("Running Feature Selection...")
    # print("Blocklist")
    # fs = pct.feature_select(
    #     profiles=df,
    #     features=data_cols,
    #     operation="blocklist"
    # )

    _, selected_features = get_data_cols(df)

    print('Missing Values')
    fs = pct.feature_select(
        profiles=df,
        features=selected_features,
        operation="drop_na_columns",
        na_cutoff=na_cutoff
    )

    _, selected_features = get_data_cols(fs)

    print("Low Variance")
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="variance_threshold",
        freq_cut=freq_cut,
        unique_cut=unique_cut
    )

    _, selected_features = get_data_cols(fs)

    print('Outlier-prone Features')
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="drop_outliers",
        outlier_cutoff=outlier_cut
    )

    _, selected_features = get_data_cols(fs)

    print("Highly Correlated Features")
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="correlation_threshold",
        corr_threshold=corr_thres,
        corr_method=corr_method,
    )
    _, selected_features = get_data_cols(fs)

    print("Noisy Features")
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="noise_removal",
        samples=NEG_CONTROL_QUERY,
        noise_removal_perturb_groups=noise_perturb_group,
        noise_removal_stdev_cutoff=nosise_std_cut
    )
    _, selected_features = get_data_cols(fs)

    sphere = pct.normalize(
        profiles=fs,
        feature_selection=selected_features,
        meta_cols=[c for c in fs.columns if c not in selected_features]
        samples=NEG_CONTROL_QUERY,
        method="spherize"
    )

    return sphere, selected_features

