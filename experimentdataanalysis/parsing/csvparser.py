# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 12:29:56 2016

@author: Michael
"""

from collections import namedtuple
import csv
import os

import experimentdataanalysis.guis.guistarter as guistarter


# %%
def is_numeric(s):
    """For a quick check if first row is numeric values or strings"""
    try:
        float(s)
        return True
    except ValueError:
        return False


# %%
def parse_csv(filepath, delimiter=','):
    """
    Returns a tuple of lists which contain the data of each column.
    If the columns have a header line with names, the data object will
    be a namedtuple from the collections module.

    returns filepath (including filename) (String), data (DataStruct)

    param 1: filepath to file to parse (including filename)
    param 2 (optional): delimiter to use, default ','
    """
    with open(filepath, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delimiter)
        firstrow = next(csvreader)
        if all(not is_numeric(value) for value in firstrow):  # row is names
            DataStruct = namedtuple("DataStruct", firstrow)
            data = DataStruct(*([] for value in firstrow))
        else:  # row is data
            data = tuple([float(value)] for value in firstrow)
        for row in csvreader:
            for value, field in zip(row, data):
                field.append(float(value))
        return filepath, data


# %%
def parse_csv_gui(defaultpath=None, delimiter=','):
    """
    Displays a GUI to choose a file, then sends it to
    parse_csv() and returns the results.

    param 1: default filepath to browse from
    param 2 (optional): delimiter to use, default ','
    """
    extfilter = "CSVs (*.csv *.dat);;TSVs (*.tsv *.dat)"
    filepath = guistarter.get_file_dialog(defaultpath=defaultpath,
                                          extensionfilter=extfilter)
    return parse_csv(filepath, delimiter=delimiter)


# %%
def parse_csv_directory(directorypath, delimiter=','):
    """
    Returns an iterator whose elements are tuples corresponding to each
    csv file in directory and subdirectories. Subdirectories' contents are
    guaranteed to be returned back-to-back.
    Returned tuple contains the following:
    1st element: filepath (including filename)
    2nd element: result of parse_csv(filepath))

    param 1: filepath to directory
    param 2 (optional): delimiter to use, default ','
    """
    subdirs = (x[0] for x in os.walk(directorypath))  # includes dir itself
    dirfiles = (x[2] for x in os.walk(directorypath))
    for subdir, files in zip(subdirs, dirfiles):
        if (len(files) > 0):
            for file in files:
                if file.lower().endswith(('.csv', '.tsv', '.dat')):
                    # TODO: add support for multiple formats
                    yield parse_csv(subdir + "\\" + file,
                                    delimiter=delimiter)


# %%
def parse_csv_directory_gui(defaultpath=None, delimiter=','):
    """
    Displays a GUI to choose a directory, then sends it to
    parse_csv_directory() and returns the results.

    param 1: default filepath to browse from
    param 2 (optional): delimiter to use, default ','
    """
    directory = guistarter.get_dir_dialog(defaultpath=defaultpath)
    return parse_csv_directory(directory, delimiter=delimiter)
