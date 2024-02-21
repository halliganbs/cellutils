import os
import re

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from scipy.spatial import distance

def euclidean(df:pd.DataFrame, data_cols, n_pcas=30):
    """
    Calculate Euclidean distance of n number of PCA columns
    Args:
        df (pd.DataFrame): Source data
        data_cols (list-like): list of relevant measurement columns
        n_pcas (int, optional): Number of PCAs. Defaults to 30.

    Returns:
        list-like: new euclidean distance columns
    """
    data = StandardScaler().fit_transform(df[data_cols])
    ar_pcas = PCA(n_components=n_pcas).fit_transform(data)
    rows, _ = ar_pcas.shape
    x0 = [0 for i in range(n_pcas)]
    euclidean_dist = []
    for i in range(rows):
        euclidean_dist.append(distance.euclidean(x0, ar_pcas[i]))
    return euclidean_dist
        

def mahalanobis(df:pd.DataFrame, data_cols, n_pcas=30):
    """
    Calculate Mahalonbis distance of n number of PCA columns
    Args:
        df (pd.DataFrame): Souce data
        data_cols (list-like): list of relevant measurement columns
        n_pcas (int, optional): number of PCAs. Defaults to 30.
    """
    data = StandardScaler().fit_transform(df[data_cols])
    pcas = PCA(n_components=n_pcas).fit_transform(data)
    rows, _ = pcas.shape
    x0 = [0 for i in range(n_pcas)]
    mahala_dist = []
    for i in range(rows):
        cm = np.cov(pcas[i])
        icm = np.linalg.inv(cm)
        # todo fix this ^^^
        mahala_dist.append(distance.mahalanobis(x0, pcas[i],icm))
    
    return mahala_dist