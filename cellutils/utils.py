import os
import re

import numpy as np
import pandas as pd

def zpad():
    print("PAD ZEROS INTO WELL ID")

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

def make_well(df:pd.DataFrame, meta_cols, data_cols, id='Image_Metadata_WellID', score='score'):
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
    data['index'] = []
    for i, wid in enumerate(df[id].unique().tolist()):
        temp = df.loc[df[id]==wid]
        for c in meta_cols:
            data[c].append(temp[c].unique()[0])
        for c in data_cols:
            data[c].append(temp[c].mean())
        data[score].append(temp[score].mean())
        data['index'].append(i)

    well = pd.DataFrame(data=data)
    print(well.shape)
    return well