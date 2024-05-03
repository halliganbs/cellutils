import os
import re
import pickle
import numpy as np
import pandas as pd

def zpad(df:pd.DataFrame, wellid:str):
    """
    Fix zpadding in well ids e.g. A1 => A01
    Args:
        df (pd.DataFrame): Source Dataframe might change to just list of well ids
        wellid (str): Well id column name

    Returns:
        np.Array: fixed well ids in np.Array
    """
    wids = df[wellid].values
    return np.array([f"{w[0]+str(w[1:]).zfill(2)}" for w in wids])

def split_metadata():
    print("split metadata of file name")
    
def zprime(pcs, ncs):
    """
    Zprime calculation for assay quality, value >= 0.5 are best
    Args:
        pcs (numpy ndarray): positive control (plateau) values (score, viability, etc)
        ncs (numpy ndarray): negative control (baseline) values (score, viability, etc)
    """
    a = (3*pcs.std())+(3*ncs.std())
    b = np.abs(pcs.mean()-ncs.mean())
    return 1-(a/b)

def make_well(df:pd.DataFrame, meta_cols, data_cols, well_id='Image_Metadata_WellID', score='score'):
    """
    Convert Cell level scored data into well level means

    Args:
        df (pd.DataFrame): cell lvel data
        meta_cols (list like): List-like of metadata columns
        data_cols (list like): List-like of data columns used in scoring
        id (str, optional): Well ID. Defaults to 'Image_Metadata_WellID'.
        score (str, optional): Score column name. Defaults to 'score'.

    Returns:
        pd.DataFrame: well level data
    """
    data = {x:[] for x in meta_cols+data_cols}
    data[score] = []
    data['count'] = []
    data['index'] = []
    for i, wid in enumerate(df[well_id].unique().tolist()):
        temp = df.loc[df[well_id]==wid]
        for c in meta_cols:
            data[c].append(temp[c].unique()[0])
        for c in data_cols:
            data[c].append(temp[c].mean())
        data[score].append(temp[score].mean())
        data['count'].append(len(temp))
        data['index'].append(i)

    well = pd.DataFrame(data=data)
    print(well.shape)
    return well

def get_data_cols(df:pd.DataFrame, extra=[""], save=False, fname='data_cols'):
    """_summary_

    Args:
        df (pd.DataFrame): _description_
        extra (list, optional): _description_. Defaults to [""].
        save (bool, optional): _description_. Defaults to False.
        fname (str, optional): _description_. Defaults to 'data_cols'.

    Returns:
        _type_: _description_
    """
    LOC_COLS = "Location|Center|Children|Parent"
    pattern = "|".join([LOC_COLS, extra])
    meta_cols = df.columns[df.columns.str.contains(pat=pattern)].tolist()
    data_cols = df.drop(columns=meta_cols).select_dtypes(include='float64').columns.tolist()
    if save:
        with open(fname, 'wb') as f:
            pickle.dump(data_cols, f)
    return data_cols