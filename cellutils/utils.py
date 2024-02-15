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