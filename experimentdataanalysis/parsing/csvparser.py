# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 12:29:56 2016

@author: Michael
"""

from collections import namedtuple
import csv
import os

import experimentdataanalysis.guis.guistarter as guistarter


# returned from csv parsing: first field is a tuple containing field names,
# second field is a nested tuple, where each tuple holds each column's data.
RawData = namedtuple("RawData", ["columnheaders", "columndata"])


# %%
def parse_csv(filepath, delimiter='\t'):
    """
    Returns a tuple of lists which contain the data of each column.
    If the columns have a header line with names, the data object will
    be a namedtuple from the collections module.

    returns filepath (including filename) (String), data (RawData)

    param 1: filepath to file to parse (including filename)
    param 2 (optional): delimiter to use, default '\t'
    """
    with open(filepath, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delimiter)
        firstrow = next(csvreader)
        ncols = len(firstrow)
        if all(not is_numeric(value) for value in firstrow):  # row is names
            columnheaders = firstrow
            columndata = [[] for col in range(ncols)]
        else:  # row is data
            columnheaders = ["Column {}".format(x + 1) for x in range(ncols)]
            columndata = [[float(value)] for value in firstrow]
        for row in csvreader:
            if len(row) >= ncols:  # ignore empty crap rows, such as at end
                for colindex, value in enumerate(row[:ncols]):  # or extra ""s
                    try:
                        columndata[colindex].append(float(value))
                    except ValueError:
                        columndata[colindex].append(value)
        csvrawdata = RawData(tuple(columnheaders),
                             tuple(tuple(col) for col in columndata))
        return filepath, csvrawdata


# %%
def parse_csv_gui(defaultpath=None, delimiter='\t'):
    """
    Displays a GUI to choose a file, then sends it to
    parse_csv() and returns the results.

    param 1: default filepath to browse from
    param 2 (optional): delimiter to use, default '\t'
    """
    extfilter = "CSVs (*.csv *.dat);;TSVs (*.tsv *.dat)"
    filepath = guistarter.get_file_dialog(defaultpath=defaultpath,
                                          extensionfilter=extfilter)
    return parse_csv(filepath, delimiter=delimiter)


# %%
def parse_csv_directory(directorypath, delimiter='\t'):
    """
    Returns an iterator whose elements are tuples corresponding to each
    csv file in directory and subdirectories. Subdirectories' contents are
    guaranteed to be returned back-to-back.
    
    Operates lazily; only reads files when iterated.
    
    Returned tuple contains the following:
    1st element: filepath (including filename)
    2nd element: result of parse_csv(filepath))

    param 1: filepath to directory
    param 2 (optional): delimiter to use, default '\t'
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
def parse_csv_directory_gui(defaultpath=None, delimiter='\t'):
    """
    Displays a GUI to choose a directory, then sends it to
    parse_csv_directory() and returns the results.

    param 1: default filepath to browse from
    param 2 (optional): delimiter to use, default '\t'
    """
    directory = guistarter.get_dir_dialog(defaultpath=defaultpath)
    return parse_csv_directory(directory, delimiter=delimiter)


# %%
def is_numeric(s):
    """For a quick check if first row is numeric values or strings"""
    try:
        float(s)
        return True
    except ValueError:
        return False
