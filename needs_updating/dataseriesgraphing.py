# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:27:59 2016

NOTE: this library is not currently used by anything, and may not follow
the current practices of the rest of the library.

It's mostly here for useful code in saving and plotting data,
or in case anyone wants to overhaul it and put it back into use.

@author: Michael
"""

import matplotlib.pyplot as plt


# %%
def plot_and_save_scandata(scandata, field_index=0, plot_options={}):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after dataseries and returns a tuple
    (drift-corrected DataSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    fig, ax = plt.subplots()
    if scandata.fitdata_list[field_index] is not None:
        prefix = "Fitted_"
    else:
        prefix = "Plot_"
    suffix = "_" + scandata.fields[field_index]
    filepath = scandata.scaninfo_list[field_index]['Filepath']
    file_name = collapse_filepath_one_dir_up(filepath,
                                             filename_prefix=prefix,
                                             filename_suffix=suffix,
                                             newextension=".png")
    try:
        plt.clf()
        plot_title = ""
        plot_scandata(scandata, field_index, plot_title, ax, plot_options)
        plt.savefig(file_name)
    except RuntimeError:
        print("Error generating {}, cancelling...".format(file_name))
    return scandata


def plot_and_save_scandata_iterable(scandata_iterable,
                                    field_index=0, plot_options={}):
    fig, ax = plt.subplots()
    for scandata in scandata_iterable:
        if scandata.fitdata_list[field_index] is not None:
            prefix = "Fitted_"
        else:
            prefix = "Plot_"
        suffix = "_" + scandata.fields[field_index]
        filepath = scandata.scaninfo_list[field_index]['Filepath']
        file_name = collapse_filepath_one_dir_up(filepath,
                                                 filename_prefix=prefix,
                                                 filename_suffix=suffix,
                                                 newextension=".png")
        try:
            plt.clf()
            plot_title = ""
            plot_scandata(scandata, field_index, plot_title, ax, plot_options)
            plt.savefig(file_name)
        except RuntimeError:
            print("Error generating {}".format(file_name))
            print("skipping file...")
    plt.close()
    return scandata_iterable


# %%
def plot_scandata(scandata, field_index=0, title=None, datatype=None,
                  axes=None, plot_options={}):
    try:
        xlabel = scandata.scaninfo_list[field_index]['FastScanType']
    except KeyError:
        xlabel = None
    plot_dataseries(scandata.dataseries_list[field_index],
                    scandata.error_dataseries_list[field_index],
                    title=title, xlabel=xlabel,
                    ylabel=scandata.fields[field_index], axes=axes,
                    fitdata=scandata.fitdata_list[field_index],
                    plot_options=plot_options)


def plot_dataseries(dataseries, error_dataseries=None,
                    title=None, xlabel=None, ylabel=None,
                    axes=None, fitdata=None, plot_options={}):
    if axes is None:
        fig, axes = plt.subplots()
    if error_dataseries is not None:
        axes.errorbar(dataseries.xvals(),
                      dataseries.yvals(),
                      error_dataseries.yvals(), fmt='b.')
    else:
        axes.plot(dataseries.xvals(),
                  dataseries.yvals(), 'b.')
    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)
    if title:
        axes.set_title(title)
    if fitdata is not None:
        axes.hold(True)
        fitdataseries = fitdata.fitdataseries
        axes.plot(fitdataseries.xvals(),
                  fitdataseries.yvals(), 'r-')
        if not plot_options.get('suppress_legend'):
            # text box with parameters of fit
            props = dict(boxstyle='round', facecolor='palegreen',
                         alpha=0.5)
            textstr = fitdata.fitparamstring
            axes.text(0.95, 0.95, textstr, transform=axes.transAxes,
                      fontsize=14, verticalalignment='top',
                      horizontalalignment='right', multialignment='left',
                      bbox=props)


def plot_additional_dataseries(dataseries, error_dataseries=None,
                               axes=None, plot_options={}):
    if axes is None:
        fig, axes = plt.subplots()
    axes.hold(True)
    if error_dataseries is not None:
        axes.errorbar(dataseries.xvals(),
                      dataseries.yvals(),
                      error_dataseries.yvals(), fmt='b.')
    else:
        axes.plot(dataseries.xvals(),
                  dataseries.yvals(), 'b.')


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
