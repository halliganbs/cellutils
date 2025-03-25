import os
import re
import pickle
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

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

def synergy_convert(df, drug2, drugs, make_dmso=False):
    """_summary_

    Args:
        df (pandas DataFrame): Dataframe with columns WellID, Cell_Line, Compound, Concentration, Cells, Viability
        drug2 (String): Drug 2 to compare against
        drugs (list of Strings): Drug1 drug lists
        make_dmso (bool, optional): Adds zero-zero concentration rows. Defaults to False.

    Returns:
        DataFrame: New synergy dataframe
    """
    df['Concentration'] = round(df['Concentration'], 3)
    df['score'] = 100 - df['Viability']
    
    a = df.loc[df['Compound']!=drug2]
    b = df.loc[df['Compound']==drug2]
    
    a['Drug1'] = a['Compound']
    a['Conc1'] = a['Concentration']
    b['Drug2'] = b['Compound']
    b['Conc2'] = b['Concentration']

    a = a[['WellID', 'Cell_Line', 'Cells','Viability', 'score', 'Drug1', 'Conc1']]
    b = b[['WellID', 'Cell_Line', 'Cells','Viability', 'score', 'Drug2', 'Conc2']]
    
    dt = pd.merge(a,b, on=['WellID', 'Cell_Line', 'Cells','Viability', 'score'], how='outer')
    
    dt['Drug2'] = drug2
    dt['Conc1'].fillna(value=0.0, inplace=True)
    dt['Conc2'].fillna(value=0.0, inplace=True)
    
    idx = dt[dt['Drug1'].isna()].sort_values(by=['Conc1', 'Conc2']).index.tolist()
    step = 0
    for i in idx:
        dt['Drug1'].iloc[i]=drugs[step]
        step+=1
        if step >= 3:
            step=0
    
    drugs_list = []
    for i in range(3):
        drugs_list+=drugs
        
    if make_dmso:
        dmso = {
            'Viability':[100,100,100,100,100,100,100,100,100],
            'score':[0,0,0,0,0,0,0,0,0],
            'Drug1':drugs_list,
            'Conc1':[0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0],
            'Drug2':[drug2 for _ in range(9)],
            'Conc2':[0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0]
        }
        dmso = pd.DataFrame(data=dmso)
        dt = pd.concat([dt,dmso])

    return dt

def find_number_components(df, data_cols, variance_threshold=0.9):
    """
    Finds the number of components needed for greater than variance threshold

    Args:
        df (pd:DataFrame): Dataframe of data
        data_cols (list<str>): Relevant data columns (excluding metadata, location, etc.)
        variance_threshold (float, optional): Standard variance threshold. Defaults to 0.9.

    Returns:
        int: Optimal number of Principal Compents 
    """
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df[data_cols])
    pca = PCA()
    _ = pca.fit_transform(df_scaled)
    num_components = np.where(np.cumsum(pca.explained_variance_ratio_)>=variance_threshold)[0][0]
    print(f"{variance_threshold} found after {num_components} Components...")
    return num_components