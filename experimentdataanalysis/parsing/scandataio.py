# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 18:54:51 2016

@author: Michael
"""

import pickle

# %% NEEDS TEST, SPHINX DOCUMENTATION
def save_scandata(scandata, filepath):
    """
    Saves scandata to file. Currently implements using pickle, but
    this should be changed to JSON to improve stability and external
    compatibility. This means custom conversion and parsing to standard
    JSON-compatible types.
    """
    with open(filepath, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(scandata, f)


# %% NEEDS TEST, SPHINX DOCUMENTATION
def load_scandata(filepath):
    """
    Loads scandata from file. Currently implements using pickle, but
    this should be changed to JSON to improve stability and external
    compatibility. This means custom conversion and parsing to standard
    JSON-compatible types.
    """
    with open(filepath, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        scandata = pickle.load(f)
    return scandata


# %% NEEDS TEST, SPHINX DOCUMENTATION
def save_scandataset(scandataset, filepath):
    """
    Saves scandata to file. Currently implements using pickle, but
    this should be changed to JSON to improve stability and external
    compatibility. This means custom conversion and parsing to standard
    JSON-compatible types - and more complicated than for scandata!
    """
    with open(filepath, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(scandataset, f)


# %% NEEDS TEST, SPHINX DOCUMENTATION
def load_scandataset(filepath):
    """
    Loads scandata from file. Currently implements using pickle, but
    this should be changed to JSON to improve stability and external
    compatibility. This means custom conversion and parsing to standard
    JSON-compatible types - and more complicated than for scandata!
    """
    with open(filepath, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        scandataset = pickle.load(f)
    return scandataset

