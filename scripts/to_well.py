from email.policy import default
import os
import re

import glob
import pickle

import numpy as np
import pandas as pd

from cellutils.utils import make_well
from .logger import logger

from progress.bar import Bar
import click

def mass_well(indir, outdir, data, meta, exten):
    if os.name == 'nt':
        os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
    assert os.path.exists(indir), 'Input directory not found...'
    if not os.path.exists(outdir):
        print("Creating out directory: ", outdir)
        os.mkdir(outdir)
    
    assert os.path.exists(data), 'Data column file not found...'
    assert os.path.exists(meta), 'Metadata column file not found...'
    
    files = glob.glob(indir+"/*")
    assert len(files) > 0, 'No files found...'
    
    with open(data, 'rb') as f:
        data_cols = pickle.load(f)
    with open(meta, 'rb') as f:
        meta_cols = pickle.load(f)
    
    logger(outdir, locals())
    
    bar = Bar('Reducing to well...', max=len(files))
    for f in files:
        fname = os.path.join(
            outdir,
            os.path.basename(f).split('.')[0]+"_well.csv")
        df = pd.read_csv(f)
        well = make_well(df, meta_cols=meta_cols, data_cols=data_cols)
        well.to_csv(fname)
        bar.next()
    bar.finish()