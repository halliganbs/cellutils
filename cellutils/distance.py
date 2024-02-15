import os
import re

import numpy as numpy
import pandas as pd

from skelarn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from scipy.spatial import distance

def euclidean():
    """
    Generate  euclidean distance score using PCAs of cell measurements
    returns new column of distances
    """
    print("this is a function")

def mahalanobis():
    """
    Generate mahalanobis distnace, including inverse covarience matrix
    returns new column of distances
    """
    print("do somethings")