import os
import re

import numpy as np
import pandas as pd

import pycytominer as pct
from pycytominer.cyto_utils.cells import SingleCells

from glob import glob
from tqdm import tqdm

from .utils import get_data_cols

import click

NEG_CONTROL_QUERY = 'Metadata_Compound == DMSO"'

def aggregate_robust(df, fname, data_cols, strata, neg_control_query='Metadata_Compound == "DMSO"', out_dest='interim' ):
    """_summary_

    Args:
        df (_type_): _description_
        fname (_type_): _description_
        data_cols (_type_): _description_
        strata (_type_): _description_
        neg_control_query (str, optional): _description_. Defaults to 'Metadata_Compound == "DMSO"'.
        out_dest (str, optional): _description_. Defaults to 'interim'.

    Returns:
        _type_: _description_
    """
    interim_name = os.path.join(out_dest, "robust_aggregate_"+os.path.basename(fname))
    df.fillna(value=0, inplace=True) # subject to change

    aggregated = pct.aggregated(
        population_df=df,
        strata=strata,
        features=data_cols,
        operation='median'        
    )
    norm_plate = pct.normalize(
        profiles=aggregated,
        features=data_cols,
        strata=strata,
        samples=neg_control_query,
        method='mad_robustize'
    )
    norm_plate.to_csv(interim_name, index=False)
    return norm_plate


def feature_selection(df, data_cols, meta_cols, 
                      na_cutoff=0.05, freq_cut=0.05, unique_cut=0.01, outlier_cut=500,
                      corr_thres=0.90, corr_method="spearman", 
                      noise_perturb_group="Metadata_CMPD", nosise_std_cut=1.0,
                      interim_dir='intrim', interim_name='feature_selected.csv'):
    
    print("Running Feature Selection...")
    print("Blocklist")
    fs = pct.feature_select(
        profiles=df,
        features=data_cols,
        operation="blocklist"
    )
    print(fs.shape) # save to log later

    _, selected_features = get_data_cols(df)

    print('Missing Values')
    fs = pct.feature_select(
        profiles=df,
        features=selected_features,
        operation="drop_na_columns",
        na_cutoff=na_cutoff
    )
    print(fs.shape)

    _, selected_features = get_data_cols(fs)

    print("Low Variance")
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="variance_threshold",
        freq_cut=freq_cut,
        unique_cut=unique_cut
    )
    print(fs.shape)

    _, selected_features = get_data_cols(fs)

    print('Outlier-prone Features')
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="drop_outliers",
        outlier_cutoff=outlier_cut
    )
    print(fs.shape)

    _, selected_features = get_data_cols(fs)

    print("Highly Correlated Features")
    fs = pct.feature_select(
        profiles=fs,
        features=selected_features,
        operation="correlation_threshold",
        corr_threshold=corr_thres,
        corr_method=corr_method,
    )
    print(fs.shape)
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
    print(fs.shape)
    _, selected_features = get_data_cols(fs)

    fs.to_csv(os.path.join(interim_dir, interim_name), index=False)

    return fs, selected_features

def spherize(df, data_cols, neg_control_query='Metadata_COMPOUND == "DMSO"'):
    
    sphere = pct.normalize(
        profiles=df,
        feature_selection=data_cols,
        meta_cols=[c for c in df.columns if c not in data_cols],
        samples=NEG_CONTROL_QUERY,
        method="spherize"
    )
    
    return sphere

def DMSO_QC(df):
    print("DO THE DMSO CORRECTION THING HERE")
    return df

@click.command()
@click.argument('src')
@click.argument('dest')
def miner(src, dest):
    assert os.path.exists(src), "Source Directory not found"
    if not os.path.exists(dest):
        os.mkdir(dest)
    
    neg_contorl_query = 'Metadata_CMPD  == "DMSO"'
    interim_dest = os.path.join(dest, "interim")
    processed_dest = os.path.join(dest, "processed")
    os.mkdir(interim_dest)
    os.mkdir(processed_dest)

    files = glob(os.path.join(src, '*.csv'))
    dfs = []
    for f in tqdm(files, desc="Plate Aggregation and Outlier Correction"):
        df = pd.read_csv(f)
        fname = os.path.basename(f)
        meta_cols, data_cols = get_data_cols(df)
        dfs.append(aggregate_robust(df, fname, data_cols, meta_cols, interim_dest))
    df = pd.concat(dfs, ignore_index=True)
    meta_cols, data_cols = get_data_cols(df) # reset data_cols and meta_cols since we need those
    fs, reduced_data_cols = feature_selection(
        df, data_cols=data_cols,
        interim_dir=interim_dest, interim_name="feature_selected.csv"
    )
    sphere = spherize(fs, reduced_data_cols, neg_control_query=neg_contorl_query)
    sphere.to_csv(os.path.join(interim_dest, "spherized_data.csv"), index=False)
    out_df = DMSO_QC(sphere)
    out_df.to_csv(os.path.join(processed_dest, "morphology_profiles_QCfiltered.csv"), index=False)
