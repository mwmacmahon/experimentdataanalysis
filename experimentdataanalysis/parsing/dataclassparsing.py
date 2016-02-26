# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:07:53 2016

@author: Michael
"""

import os.path
import time

import experimentdataanalysis.parsing.csvparser as csvparser
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, TimeSeries


# %%
def fetch_dir_as_unfit_scandata_iterator(directorypath=None,
                                         attribute='lockin2x',
                                         scale_by_laserpower=False):
    """
    Takes a directory and returns an iterator, which upon each call gives
    the contents of a csv file in that directory or its subdirectories in
    TimeSeries form, along with the filename and filepath.

    Yield format: (String filepath, String filename, TimeSeries data)

    Keyword arguments:
    directorypath -- if given, the target directory. If not given, a
        GUI is used to select the directory. (default None)
    attribute -- the csv column to be used for the TimeSeries values.
        (default 'lockin2x')
    scale_by_laserpower -- if True, the TimeSeries values consist of the
        attribute value divided by the laserpower value at each time point.
        (default False)
    """

    if directorypath is not None:
        csvfiles = csvparser.parse_csv_directory(directorypath, delimiter='\t')
    else:
        csvfiles = csvparser.parse_csv_directory_gui('C:\\Data\\',
                                                     delimiter='\t')
    for csvfile in csvfiles:
        filepath = csvfile[0]
        rawdata = csvfile[1]
        timeseries = unpack_raw_csv_data_as_timeseries(rawdata, attribute,
                                                       scale_by_laserpower)
        scaninfo = analyze_scan_filepath(filepath)
        scaninfo['filepath'] = filepath
        scaninfo['file last modified'] = time.ctime(os.path.getmtime(filepath))
        scaninfo['attribute'] = attribute
        scaninfo['csv raw data'] = rawdata
        scandata = ScanData(filepath, scaninfo, timeseries, None)
        yield scandata


# %%
def fetch_csv_as_unfit_scandata(filepath=None,
                                attribute='lockin2x',
                                scale_by_laserpower=False):
    if filepath is not None:
        filepath, rawdata = csvparser.parse_csv(filepath, delimiter='\t')
    else:
        filepath, rawdata = csvparser.parse_csv_gui('C:\\Data\\',
                                                    delimiter='\t')
    timeseries = unpack_raw_csv_data_as_timeseries(rawdata, attribute,
                                                   scale_by_laserpower)
    scaninfo = analyze_scan_filepath(filepath)
    scaninfo['Filepath'] = filepath
    scaninfo['File Last Modified'] = time.ctime(os.path.getmtime(filepath))
    scaninfo['Attribute'] = attribute
    scaninfo['CSV Raw Data'] = rawdata
    scandata = ScanData(filepath, scaninfo, timeseries, None)
    return scandata


# %%
def unpack_raw_csv_data_as_timeseries(rawdata, attribute,
                                      scale_by_laserpower):
    tdata = rawdata.scancoord
    zdataraw = rawdata.__getattribute__(attribute)
    if scale_by_laserpower:
        laserpower = rawdata.laserpower
        zdata = [zdatapt/laserpowerpt
                 for zdatapt, laserpowerpt in zip(zdataraw, laserpower)]
        timeseries = TimeSeries(zip(tdata, zdata))
    else:
        timeseries = TimeSeries(zip(tdata, zdataraw))
    return timeseries


# %%
def analyze_scan_filepath(filepath):
    """
    Scans the filepath (including filename) of a scan, and returns
    a dict containing all the info and tags it can match.

    Custom keyword lists can be passed by the keywordlists keyword
    argument, where None for a keywordlist type means
    """
    scaninfo = {}
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
