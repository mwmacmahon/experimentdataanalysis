# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:07:53 2016

@author: Michael
"""

import os.path
import time

import experimentdataanalysis.parsing.csvparser as csvparser
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries


# %%
def fetch_dir_as_unfit_scandata_iterator(directorypath=None,
                                         attribute='lockin2x'):
    """
    Takes a directory and returns an iterator, which upon each call gives
    the contents of a csv file in that directory or its subdirectories in
    DataSeries form, along with the filename and filepath.

    Yield format: (String filepath, String filename, DataSeries data)

    Keyword arguments:
    directorypath -- if given, the target directory. If not given, a
        GUI is used to select the directory. (default None)
    attribute -- the csv column to be used for the DataSeries values.
        (default 'lockin2x')
    """
    # Grab raw data
    if directorypath is not None:
        csvfiles = csvparser.parse_csv_directory(directorypath, delimiter='\t')
    else:
        csvfiles = csvparser.parse_csv_directory_gui('C:\\Data\\',
                                                     delimiter='\t')
    for filepath, rawdata in csvfiles:
        scandata = parsed_csv_to_unfit_scandata(filepath, rawdata, attribute)
        yield scandata


# %%
def fetch_csv_as_unfit_scandata(filepath=None,
                                attribute='lockin2x'):
    # Grab raw data
    if filepath is not None:
        filepath, rawdata = csvparser.parse_csv(filepath, delimiter='\t')
    else:
        filepath, rawdata = csvparser.parse_csv_gui('C:\\Data\\',
                                                    delimiter='\t')
    scandata = parsed_csv_to_unfit_scandata(filepath, rawdata, attribute)
    return scandata


# %%
def parsed_csv_to_unfit_scandata(filepath, rawdata, attribute='lockin2x'):
    colnames, coldata = rawdata
    try:
        attribute_index = colnames.index(attribute)
    except ValueError:  # invalid attribute, so no matching index found
        print("Warning: Attribute {} not found in ".format(attribute) +
              "csv file, order from csv file will be kept.")
        attribute_index = 0
    # Work out field ordering based on given attribute - attribute to front
    rawdata_indices = list(range(len(colnames)))
    rawdata_indices.remove(attribute_index)
    rawdata_indices = [attribute_index] + rawdata_indices

    # Assemble ScanData
    scaninfo = analyze_scan_filepath(filepath)
    fields = tuple(colnames[ind] for ind in rawdata_indices)
    dataseries = tuple(DataSeries(zip(coldata[0], coldata[ind]))
                       for ind in rawdata_indices)
    fitdata = tuple(None for ind in rawdata_indices)
    scandata = ScanData(filepath, scaninfo, fields, dataseries, fitdata)
    return scandata


# %%
def analyze_scan_filepath(filepath):
    """
    Scans the filepath (including filename) of a scan, and returns
    a dict containing all the info and tags it can match.

    Custom keyword lists can be passed by the keywordlists keyword
    argument, where None for a keywordlist type means
    """
    scaninfo = {'Filepath': filepath}
    scaninfo['File Last Modified'] = time.ctime(os.path.getmtime(filepath))
    # SEARCH TERMS:
    # 1. If first string found, register second string as
    #    tag containing True
    this_element_keyword_list = [("Channel1", "Channel 1"),
                                 ("Channel2", "Channel 2"),
                                 ("Channel3", "Channel 3"),
                                 ("Channel4", "Channel 4")]
    # 2. Grab next element(s) if this one CONTAINS first string,
    #    tag next element(s) as second string(s)
    next_element_keyword_list = [("Voltage", "Voltage"),
                                 ("Ind", "FastScanIndex"),
                                 ("2Dscan", ["MiddleScanType",
                                             "FastScanType"])]
    # 3. Grab this element if it CONTAINS first string,
    #    tag remainder as second string
    inside_this_element_keyword_list = [("K", "SetTemperature"),
                                        ("nm", "Wavelength"),
                                        ("x.dat", "MiddleScanCoord")]
    for segment in filepath.split("\\"):
        # get rid of idiosyncratic delimiters by swapping with _
        segment = segment.replace(" ", "_")
        next_element_tags = []
        for element in segment.split("_"):
            if len(next_element_tags) > 0:
                try:
                    value = float(element)
                except ValueError:
                    value = element
                scaninfo[next_element_tags.pop(0)] = value
            for matchstr, tag in this_element_keyword_list:
                if element == matchstr:
                    scaninfo[tag] = True
            for matchstr, tags in next_element_keyword_list:
                if element == matchstr:
                    if isinstance(tags, str):  # only one string
                        next_element_tags = [tags]
                    else:
                        next_element_tags = tags
            for matchstr, tag in inside_this_element_keyword_list:
                if matchstr in element:
                    value = element.replace(matchstr, "")
                    try:
                        value = float(value)
                    except ValueError:
                        value = value
                    scaninfo[tag] = value
    return scaninfo
