import os
import re
import glob
import click
import sqlite3
import pickle
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from xgboost  import XGBRegressor

from cellutils.utils import make_well
from .logger import logger

@click.command()
@click.argument('db')
@click.argument('colname')
@click.option('--table', '-t', default='Per_Object', required=True, help='Database table name')
@click.option('--fname', '-f', default='Per_Object', required=True, help='Ouput file name')
@click.option('--log', '-l', default='score_log', required=False, help='Logfile output name')
def score_plate(db, colname,  modelname, table, fname, log):
    # windows check
    if os.name == 'nt':
        os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
    
    # check if files exist
    assert os.path.exists(db), 'sqlite file does not exist'
    assert os.path.exists(colname), 'Data column file does not exist'
    assert os.path.exists(modelname), "Model file does not exist"
    
    logger(log, locals())
    
    with open(colname, 'rb') as f:
        data_cols = pickle.load(f)
    
    model = XGBRegressor()
    model.load_model(modelname)
    
    print(f"Reading table: {table}... ")
    con = sqlite3.connect(db)
    query = f"SELECT * FROM {table}"
    dfs = []
    for chunk in pd.read_sql_query(query, con, chunksize=int(1e6)):
        print(f"Read chunk with size of: {len(chunk)} rows")
        dfs.append(chunk)
        del chunk
    df = pd.concat(dfs)
    
    print("Scoring...")
    data = StandardScaler().fit_transform(df[data_cols])
    df['score'] = model.predict(data)
    df.to_csv(fname+".csv", index=False)
    