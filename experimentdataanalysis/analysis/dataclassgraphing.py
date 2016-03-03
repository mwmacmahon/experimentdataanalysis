# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:27:59 2016

@author: Michael
"""

import matplotlib.pyplot as plt

import experimentdataanalysis.analysis.dataclassfitting as dcfitting
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, TimeSeries
import experimentdataanalysis.parsing.dataclassparsing as dcparsing


# %%
def graph_fitted_csv(filepath=None,
                     attribute='lockin2x',
                     scale_by_laserpower=False,
                     fit_drift=False,
                     timeseriesfitfunction=None):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after timeseries and returns a tuple
    (drift-corrected TimeSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    scandata = dcparsing.fetch_csv_as_unfit_scandata(
                            filepath, attribute, scale_by_laserpower)
    scandata = dcfitting.get_time_offset_scandata(scandata)
    scandata = dcfitting.fit_scandata(scandata,
                                      timeseriesfitfunction, fit_drift)
    scandata = graph_scandata(scandata)
    return scandata


def graph_fitted_csv_directory(directorypath=None,
                               attribute='lockin2x',
                               scale_by_laserpower=False,
                               timeseriesfitfunction=None,
                               fit_drift=False,
                               multiprocessing=False):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after timeseries and returns a tuple
    (drift-corrected TimeSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    scandata_iterator = dcparsing.fetch_dir_as_unfit_scandata_iterator(
                            directorypath, attribute, scale_by_laserpower)
    scandata_iterable = [dcfitting.get_time_offset_scandata(scandata)
                         for scandata in scandata_iterator]
    scandata_iterable = dcfitting.fit_scandata_iterable(
                            scandata_iterable,
                            timeseriesfitfunction,
                            fit_drift)
    graph_scandata_iterable(scandata_iterable, timeseriesfitfunction,
                            fit_drift, multiprocessing)
    return scandata_iterable


# %%
def graph_scandata(scandata):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after timeseries and returns a tuple
    (drift-corrected TimeSeries, FitData) instead of just FitData
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
        plot_timeseries(scandata.timeseries,
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
            plot_timeseries(scandata.timeseries,
                            plot_title, ax, scandata.fitdata)
            plt.savefig(file_name)
        except RuntimeError:
            print("Error generating {}".format(file_name))
            print("skipping file...")
    plt.close()
    return scandata_iterable


# %%
def plot_scandata(scandata, title="", axes=None):
    plot_timeseries(scandata.timeseries, title=title,
                    axes=axes, fitdata=scandata.fitdata)


def plot_timeseries(timeseries, title="", axes=None, fitdata=None):
    if axes is None:
        fig, axes = plt.subplots()
    plt.plot(timeseries.times(unfiltered=True),
             timeseries.values(unfiltered=True), 'b.')
    if fitdata is not None:
        fittimeseries = fitdata.fittimeseries
        plt.plot(fittimeseries.times(unfiltered=True),
                 fittimeseries.values(unfiltered=True), 'r.')
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
