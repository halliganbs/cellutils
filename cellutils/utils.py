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
        df (pd.DataFrame): cell level data
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

def char_range(c1, c2):
    """
    Provide a range over letters, e.g. A P gives letters ABCDEFGHIJKLMNOP
    Args:
        c1 (str): start letter
        c2 (str): end letter
    """
    for c in range(ord(c1), ord(c2)+1):
        yield(chr(c))

def add_controls(df:pd.DataFrame, rows=('A', 'P'), nc_cols=(1, 2), pc_cols=(23,24), pc_val='Test', nc_val='DMSO', well_id='wellid', cmpd_col='Compound', extra_cols={}, test=False):
    """
    Adds Rows for controls, typically when HPD300 is used for controls and echo for compounds
    Args:
        df (pd.DataFrame): Source dataframe
        rows (tuple, optional): Letter row ids . Defaults to ('A', 'P').
        pc_cols (tuple, optional): start and end cols . Defaults to (1, 2).
        nc_cols (tuple, optional): _description_. Defaults to (23,24).
        pc_val (str, optional): _description_. Defaults to 'Test'.
        nc_val (str, optional): _description_. Defaults to 'DMSO'.
        well_id (str, optional): _description_. Defaults to 'wellid'.
        cmpd_col (str, optional): _description_. Defaults to 'Compound'.
        extra_cols (dict, optional): _description_. Defaults to {}.

    Returns:
        _type_: _description_
    """
    _, scr_cols = df.shape
    data = {c:[] for c in df.columns.tolist()}
    for c in char_range(rows[0], rows[1]):
        for i in range(pc_cols[0], pc_cols[1]+1):
            name = c+str(i).zfill(2)
            data[well_id].append(name)
            data[cmpd_col].append(pc_val)
            for k in extra_cols.keys():
                data[k].append(extra_cols[k])
        for i in range(nc_cols[0], nc_cols[1]+1):
            name = c+str(i).zfill(2)
            data[well_id].append(name)
            data[cmpd_col].append(nc_val)
            for k in extra_cols.keys():
                data[k].append(extra_cols[k])
    if test:
        for k in data.keys():
            print(k, len(data[k]))
    dt = pd.DataFrame(data=data)
    _, out_cols = dt.shape
    assert scr_cols==out_cols, "Number of Columns do not match"
    return pd.concat([df, dt])