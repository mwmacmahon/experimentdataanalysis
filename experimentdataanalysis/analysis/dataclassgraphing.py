# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:27:59 2016

@author: Michael
"""

import matplotlib.pyplot as plt

import experimentdataanalysis.analysis.dataclassfitting as dcfitting
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
import experimentdataanalysis.parsing.dataclassparsing as dcparsing


# %%
def graph_fitted_csv(filepath=None,
                     attribute='lockin2x',
                     scale_by_laserpower=False,
                     time_offset=False,
                     fit_drift=False,
                     dataseriesfitfunction=None):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after dataseries and returns a tuple
    (drift-corrected DataSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    scandata = dcparsing.fetch_csv_as_unfit_scandata(
                            filepath, attribute, scale_by_laserpower)
    if time_offset:
        scandata = dcfitting.get_time_offset_scandata(scandata)
    scandata = dcfitting.fit_scandata(scandata,
                                      dataseriesfitfunction, fit_drift)
    scandata = graph_scandata(scandata)
    return scandata


def graph_fitted_csv_directory(directorypath=None,
                               attribute='lockin2x',
                               scale_by_laserpower=False,
                               dataseriesfitfunction=None,
                               time_offset=False,
                               fit_drift=False,
                               multiprocessing=False):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after dataseries and returns a tuple
    (drift-corrected DataSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    scandata_iterator = dcparsing.fetch_dir_as_unfit_scandata_iterator(
                            directorypath, attribute, scale_by_laserpower)
    if time_offset:
        scandata_iterable = [dcfitting.get_time_offset_scandata(scandata)
                             for scandata in scandata_iterator]
    else:
        scandata_iterable = list(scandata_iterator)
    scandata_iterable = dcfitting.fit_scandata_iterable(
                            scandata_iterable,
                            dataseriesfitfunction,
                            fit_drift)
    graph_scandata_iterable(scandata_iterable, dataseriesfitfunction,
                            fit_drift, multiprocessing)
    return scandata_iterable


# %%
def graph_scandata(scandata):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after dataseries and returns a tuple
    (drift-corrected DataSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    fig, ax = plt.subplots()
    if scandata.fitdata is not None:
        prefix = "Fitted_"
    else:
        prefix = "Plot_"
    suffix = "_" + scandata.scaninfo['Attribute']
    file_name = collapse_filepath_one_dir_up(scandata.filepath,
                                             filename_prefix=prefix,
                                             filename_suffix=suffix,
                                             newextension=".png")
    try:
        plt.clf()
        plot_title = ""
        plot_dataseries(scandata.dataseries,
                        plot_title, ax, scandata.fitdata)
        plt.savefig(file_name)
    except RuntimeError:
        print("Error generating {}, cancelling...".format(file_name))
    return scandata


def graph_scandata_iterable(scandata_iterable):
    fig, ax = plt.subplots()
    for scandata in scandata_iterable:
        if scandata.fitdata is not None:
            prefix = "Fitted_"
        else:
            prefix = "Plot_"
        suffix = "_" + scandata.scaninfo['Attribute']
        file_name = collapse_filepath_one_dir_up(scandata.filepath,
                                                 filename_prefix=prefix,
                                                 filename_suffix=suffix,
                                                 newextension=".png")
        try:
            plt.clf()
            plot_title = ""
            plot_dataseries(scandata.dataseries,
                            plot_title, ax, scandata.fitdata)
            plt.savefig(file_name)
        except RuntimeError:
            print("Error generating {}".format(file_name))
            print("skipping file...")
    plt.close()
    return scandata_iterable


# %%
def plot_scandata(scandata, title="", axes=None):
    plot_dataseries(scandata.dataseries, title=title,
                    axes=axes, fitdata=scandata.fitdata)


def plot_dataseries(dataseries, title="", axes=None, fitdata=None):
    if axes is None:
        fig, axes = plt.subplots()
    plt.plot(dataseries.xvals(unfiltered=True),
             dataseries.yvals(unfiltered=True), 'b.')
    if fitdata is not None:
        fitdataseries = fitdata.fitdataseries
        plt.plot(fitdataseries.xvals(unfiltered=True),
                 fitdataseries.yvals(unfiltered=True), 'r.')
        # text box with parameters of fit
        props = dict(boxstyle='round', facecolor='palegreen',
                     alpha=0.5)
        textstr = fitdata.fitparamstring
        plt.text(0.95, 0.95, textstr, transform=axes.transAxes,
                 fontsize=14, verticalalignment='top',
                 horizontalalignment='right', multialignment='left',
                 bbox=props)
        plt.title(title)
        plt.draw()


# %%
def collapse_filepath_one_dir_up(filepath,
                                 filename_prefix=None,
                                 filename_suffix=None,
                                 newextension=None):
    segments = filepath.split("\\")
    if filename_prefix is None:
        filepath_minus_last_segment = "\\".join(segments[:-1])
    else:
        filepath_except_last_two = "\\".join(segments[:-2])
        filepath_minus_last_segment = filepath_except_last_two + \
            "\\" + filename_prefix + segments[-2]
    filepath = filepath_minus_last_segment + "_" + segments[-1]
    if filename_suffix is not None:
        filepath = filepath[:-4] + filename_suffix + filepath[-4:]
    if newextension is not None:
        filepath = filepath[:-4] + newextension
    return filepath
