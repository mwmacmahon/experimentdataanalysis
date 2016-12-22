# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:07:53 2016

@author: Michael
"""

import os.path
import time

import numpy as np

import experimentdataanalysis.parsing.csvparser as csvparser
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData


# %%
def fetch_dir_as_unfit_scandata_iterator(directorypath=None,
                                         key_field=None,
                                         delimiter=None,
                                         parser='tab-delimited',
                                         key_field_error_val=None,
                                         parsing_keywordlists=None):
    """
    Takes a directory and returns an iterator, which upon each call gives
    the contents of a csv file in that directory or its subdirectories in
    ScanData form. Other file types can be read via other parsers.

    Important note: parsers should return ITERATORS, NOT LISTS.
    This function should parse that iterator one element at a time
    and yield a corresponding scandata, so that in theory one could sort
    through thousands of files with only one invocation of this function.
        read file, convert to scandata, yield scandata, repeat...
    See csvparser.py and the csv parsing "if" branch for an example.

    Yield format: (String filepath, String filename, DataSeries data)

    Keyword arguments:
    directorypath -- if given, the target directory. If not given, a
        GUI is used to select the directory. (default None)
    key_field -- if found, the given field will be moved to index 0
        in the list (default 'lockin2x')
    """
    if parser == 'tab-delimited':  # default case
        if not delimiter:
            delimiter = '\t'
        if directorypath is not None:
            csvfiles = csvparser.parse_csv_directory(directorypath,
                                                     delimiter)
        else:
            csvfiles = csvparser.parse_csv_directory_gui('C:\\Data\\',
                                                         delimiter)
        # NOTE: csvfiles is an iterator, only reads one
        # file on demand each time through loop:
        for filepath, headerfooterstr, tabledata in csvfiles:
            scandata = tabledata_to_unfit_scandata(filepath, headerfooterstr,
                                                   tabledata, key_field,
                                                   key_field_error_val,
                                                   parsing_keywordlists)
            yield scandata
#    elif parser == '???':
#        if not delimiter:
#            delimiter = 'vOv'
#       ...
#       yield scandata
    else:
        raise AttributeError("fetch_dir_as_unfit_scandata_iterator: parsing " +
                             "type invalid or not given.")


# %%
def fetch_csv_as_unfit_scandata(filepath=None,
                                key_field=None,
                                delimiter=None,
                                parser='tab-delimited',
                                key_field_error_val=None,
                                parsing_keywordlists=None):
    if parser == 'tab-delimited':  # default case
        if not delimiter:
            delimiter = '\t'
        if filepath is not None:
            filepath, headerfooterstr, tabledata = \
                csvparser.parse_csv(filepath, delimiter)
        else:
            filepath, headerfooterstr, tabledata = \
                csvparser.parse_csv_gui('C:\\Data\\', delimiter)
        scandata = tabledata_to_unfit_scandata(filepath, headerfooterstr,
                                               tabledata, key_field,
                                               key_field_error_val,
                                               parsing_keywordlists)
        return scandata
#    elif parser == '???':
#        if not delimiter:
#            delimiter = 'vOv'
#       ...
#       yield scandata
    else:
        raise AttributeError("fetch_csv_as_unfit_scandata: parsing " +
                             "type invalid or not given.")


# %%
def tabledata_to_unfit_scandata(filepath, headerfooterstr,
                                rawdata, key_field=None,
                                key_field_error_val=None,
                                parsing_keywordlists=None):
    colnames, coldata = rawdata
    field_names = list(colnames)
    field_arrays = [np.array(column) for column in coldata]
    scaninfo = analyze_scan_filepath(filepath,
                                     parsing_keywordlists)
    scaninfo = analyze_string_for_dict_pairs(headerfooterstr, scaninfo)

    # if key_field given, attempt to set error array in scaninfo
    if key_field and key_field_error_val:
        if key_field in field_names:
            error_array = key_field_error_val * np.ones(len(field_arrays[0]))
            field_names.append(key_field + '_error')
            field_arrays.append(error_array)

    # Assemble ScanData
    scandata = ScanData(field_names,
                        field_arrays,
                        scaninfo,
                        x_field_name=None,  # first column defaults to x anyway
                        y_field_name=key_field)

    return scandata


# %%
def analyze_scan_filepath(filepath, scaninfo={}, keywordlists=None):
    """
    Scans the filepath (including filename) of a scan, and returns
    a dict containing all the info and tags it can match.

    Custom keyword lists can be passed by the keywordlists keyword
    argument, where None for a keywordlist type means use defaults
    """
    scaninfo['Filepath'] = filepath
    scaninfo['File Last Modified'] = time.ctime(os.path.getmtime(filepath))
    if keywordlists:
        this_element_keyword_list, \
            next_element_keyword_list, \
            inside_this_element_keyword_list = keywordlists
    else:
        # DEFAULT SEARCH TERMS AND SEARCH RULES:
        # 1. If first string found, register second string as
        #    tag containing third string/value
        #        e.g. if keyword_list contains ("warmup", "Warmup?", "Yes"):
        #             "...warmup..." -> {"Warmup?": "Yes"}
        this_element_keyword_list = []
        # 2. Grab next element(s) if this one CONTAINS first string,
        #    tag next element(s) as second string(s)
        #        e.g. "..._Ind_3_..." -> {"FastScanIndex": 3}
        #        e.g. "..._2Dscan_MirrorY_MirrorZ_..."
        #                 -> {"MiddleScanType": "MirrorY",
        #                     "FastScanType": "MirrorZ"}
        next_element_keyword_list = [("Voltage", "Voltage"),
                                     ("Ind", "FastScanIndex"),
                                     ("2Dscan", ["MiddleScanType",
                                                 "FastScanType"])]
        # 3. Grab this element if it CONTAINS first string,
        #    tag remainder as second string
        #        e.g. "..._30K_..." -> {"SetTemperature": 30}
        inside_this_element_keyword_list = [("Channel", "Channel"),
                                            ("Vcm", "Electric Field"),
                                            ("mT", "Magnetic Field (mT)"),
                                            ("K", "SetTemperature"),
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
                except ValueError:  # if not a numeric value:
                    # ignore trailing keywords, e.g. units
                    value = element
                    for matchstr, _ in inside_this_element_keyword_list:
                        if matchstr in element:
                            value = element.replace(matchstr, "")
                    try:
                        value = float(value)
                    except ValueError:  # still non-numeric
                        pass
                scaninfo[next_element_tags.pop(0)] = value
            else:
                for matchstr, tags in next_element_keyword_list:
                    if element == matchstr:
                        if isinstance(tags, str):  # only one string
                            next_element_tags = [tags]
                        else:
                            next_element_tags = tags
            for matchstr, tag, value in this_element_keyword_list:
                if element == matchstr:
                    scaninfo[tag] = value
            for matchstr, tag in inside_this_element_keyword_list:
                if matchstr in element:
                    value = element.replace(matchstr, "")
                    try:
                        value = float(value)
                    except ValueError:
                        value = value
                    scaninfo[tag] = value
    return scaninfo


# %%
def analyze_string_for_dict_pairs(infostr, scaninfo={}):
    """
    Currently looks for key, value pairs in form "key: value" in the
    provided strings and adds them to the dict given (or otherwise
    creates a new dict).
    """
    strrows = infostr.splitlines()
    for row in strrows:
        key, value = "", ""
        try:
            key, value = row.split(":")
        except ValueError:
            pass
        if key:
            key = key.strip()
            value = value.strip()
            scaninfo[key] = value
    return scaninfo
