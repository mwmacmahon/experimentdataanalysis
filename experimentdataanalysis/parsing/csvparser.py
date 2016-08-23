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
        csvrows = [row for row in csvreader]  # read file into memory
    csvrows = [[col.strip() for col in row]  # strip whitespace and newlines
               for row in csvrows]
    csvrows = [row for row in csvrows  # excise any rows with no values
               if any(row)]

    ncols = 1  # find number of rows as # of cols in last non-footer line
    table_start = 0
    table_end = 0
    for rownum, row in reversed(list(enumerate(csvrows))):  # note: 2x RAM used
        if ncols == 1:
            ncols = len(row)
            table_end = rownum + 1
        elif ncols == len(row):
            table_start = rownum
        else:
            break
    headerfooterrows = csvrows[:table_start] + csvrows[table_end:]
    headerfooterstr = "\n".join(" ".join(row)
                                for row in headerfooterrows)

    if not all(is_numeric(col) for col in csvrows[table_start]):  # header row
        colheaders = tuple(csvrows[table_start])
        table_start += 1
    else:
        colheaders = tuple("Column " + str(colnum + 1)
                           for colnum in range(ncols))
    coltuples = tuple(tuple(float(row[colnum])
                            for row in csvrows[table_start:table_end])
                      for colnum in range(ncols))

    csvrawdata = RawData(colheaders, coltuples)
    return filepath, headerfooterstr, csvrawdata


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
