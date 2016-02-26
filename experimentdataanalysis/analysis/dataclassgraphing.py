# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:27:59 2016

@author: Michael
"""

import matplotlib.pyplot as plt

import experimentdataanalysis.analysis.curvefitting as curvefitting
import experimentdataanalysis.analysis.dataclassfitting as dataclassfitting
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, TimeSeries
import experimentdataanalysis.parsing.dataclassparsing as dataparsing


# %%
def graph_csv(filepath=None,
              attribute='lockin2x',
              scale_by_laserpower=False,
              timeseriesfitfunction=None):
    unfit_scandata = dataparsing.fetch_csv_as_unfit_scandata(filepath)


def graph_csv_directory(directorypath=None,
                        attribute='lockin2x',
                        scale_by_laserpower=False,
                        timeseriesfitfunction=None,
                        multiprocess_fit=False):
    scandata_iterator = dataparsing.fetch_dir_as_unfit_scandata_iterator(
        directorypath, attribute, scale_by_laserpower=False)
    graph_scandata_iterator(timeseriesfitfunction, multiprocess_fit)


# %%
def graph_scandata(scandata, timeseriesfitfunction=None):
    offset_scandata = dataclassfitting.get_time_offset_scandata(scandata)
    if timeseriesfitfunction is not None:
        timeseries, fitdata = timeseriesfitfunction(timeseries)
    else:
        fitdata = None
    fig, ax = plt.subplots()
    plot_title = ""
    plot_timeseries(timeseries, plot_title, ax, fitdata)


def graph_scandata_iterator(scandata_iterator,
                            timeseriesfitfunction=None,
                            multiprocess_fit=False):
    offset_scandata_iterable = \
        [dataclassfitting.get_time_offset_scandata(scandata)
         for scandata in scandata_iterator]

    if timeseriesfitfunction is not None:
        if multiprocess_fit:
            fitoutput_iterator = multiprocessable_map(timeseriesfitfunction,
                                                    timeseries_iterator)
        else:
            fitoutput_iterator = (timeseriesfitfunction(timeseries)
                                  for timeseries in timeseries_iterator)
    fig, ax = plt.subplots()
    for filepath, filename in zip(filepaths, filenames):
        # TODO: extract voltage, etc. fit data, and save results to file
        if timeseriesfitfunction is not None:
            timeseries, fitdata = next(fitoutput_iterator)
        else:
            timeseries = next(timeseries_iterator)
            fitdata = None
        try:
            plt.clf()
            plot_title = filepath + "\\" + filename
            plot_timeseries(timeseries, plot_title, ax, fitdata)
            plt.savefig(filepath + "_" + filename + ".png")
        except RuntimeError:
            print("Error plotting {}".format(filepath + "_" + filename))
            print("skipping file...")
    plt.close()


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


