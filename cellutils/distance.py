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
        

def mahalanobis(df:pd.DataFrame, data_cols, n_pcas=30, cov=None):
    """
   https://www.geeksforgeeks.org/how-to-calculate-mahalanobis-distance-in-python/
    Calculate Mahalonbis distance of n number of PCA columns
    Args:
        df (pd.DataFrame): Souce data
        data_cols (list-like): list of relevant measurement columns
        n_pcas (int, optional): number of PCAs. Defaults to 30.
    """
    # pcas = PCA(n_components=n_pcas).fit_transform(df[data_cols])
    # rows, _ = pcas.shape
    # y_mu = df - np.mean(df[data_cols]) 
    # if not cov:
    #     cov = np.cov(df[data_cols].values.T) 
    # inv_covmat = np.linalg.inv(cov)
    # left = np.dot(y_mu, inv_covmat)
    # mahal = np.dot(left, y_mu.T)
    # return mahal.diagonal()
    rows = df[data_cols].shape
    x0 = np.zeros(len(data_cols))
    covariance = np.array(df[data_cols].values)
    mahal_dist = []
    for _, row  in df.iterrows():
        mahal = distance.mahalanobis(x0, row[data_cols], np.linalg.inv(covariance))
        mahal_dist.append(mahal)
    return mahal_dist