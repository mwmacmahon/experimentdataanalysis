# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:14:54 2016

@author: Michael
"""

import experimentdataanalysis.analysis.curvefitting as curvefitting
from experimentdataanalysis.analysis.generalutilities \
    import multiprocessable_map
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries


# %%
def fit_scandata_iterable(scandata_iterable,
                          dataseriesfitfunction=None,
                          fit_drift=False,
                          multiprocessing=False):
    """
    To use drift fitting, must send a function that accepts a fit_drift
    flag as 2nd argument after dataseries and returns a tuple
    (drift-corrected DataSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    if not multiprocessing:
        for scandata in scandata_iterable:
            yield fit_scandata(scandata, dataseriesfitfunction, fit_drift)
    else:  # more complicated, must break down
        filepaths, scaninfos, dataserieses, _ = zip(*scandata_iterable)
        # need this multiple times: it MUST be a list comp, NOT genexp:
        dataserieses = [get_time_offset_dataseries(dataseries)
                        for dataseries in dataserieses]
        if fit_drift:
            fitfunction_args_iter = ([dataseries, True]  # fit_drift flag
                                     for dataseries in dataserieses)
        else:
            fitfunction_args_iter = ([dataseries]
                                     for dataseries in dataserieses)
        fitoutput_iterator = multiprocessable_map(dataseriesfitfunction,
                                                  fitfunction_args_iter,
                                                  multiprocessing=True)
        for filepath, scaninfo, dataseries, fitoutput in zip(
                filepaths, scaninfos, dataserieses, fitoutput_iterator):
            if fit_drift:
                newscaninfo = scaninfo.copy()  # shallow dict copy!
                newscaninfo['pre-fit_drift dataseries'] = dataseries
                yield ScanData(filepath, newscaninfo, *fitoutput)
            else:
                yield ScanData(filepath, scaninfo, dataseries, fitoutput)


# %%
def fit_scandata(scandata, dataseriesfitfunction=None, fit_drift=False):
    """
    Fits a ScanData object with the given function and returns a new
    ScanData object with a new time offset and the fitted result.
    If fit_drift is used, given function must have fit_drift flag as
    2nd argument after dataseries and return
    (drift-corrected DataSeries, FitData) instead of just FitData
    when flag is set to True.
    """
    if dataseriesfitfunction is not None:
        scandata = get_time_offset_scandata(scandata)
        dataseries = scandata.dataseries
        if fit_drift:
            newscaninfo = scandata.scaninfo.copy()  # shallow dict copy!
            newscaninfo['pre-fit_drift dataseries'] = dataseries
            return ScanData(scandata.filepath,
                            newscaninfo,
                            *dataseriesfitfunction(dataseries, True))
        else:
            return ScanData(scandata.filepath,
                            scandata.scaninfo,
                            dataseries,
                            dataseriesfitfunction(dataseries))
    else:
        return ScanData(scandata.filepath,
                        scandata.scaninfo,
                        scandata.dataseries,
                        scandata.fitdata)


# %%
def fit_dataseries_with_one_decaying_cos(dataseries,
                                         fit_drift=False):
    """
    Takes a DataSeries and fits as a function of a decaying cosine
    with variable phase, plus a sharp exponential and a flat offset
    using fit_dataseries__fit_driftting(dataseries, fitfunction)

    if fit_drift is False:
        returns FitData(dataseries)
    if fit_drift is True:
        returns (no_bkgd_dataseries, FitData(dataseries))
            where no_bkgd_dataseries is dataseries with background
            drift subtracted out.

    Note: fit_drift must be default off, so that it may be used in
    functions that are unaware of drift fitting capabilities
    """
    if fit_drift:
        return fit_dataseries_plus_drift_fitting(
                dataseries,
                curvefitting.fit_sorted_series_to_decaying_cos)
    else:
        return curvefitting.fit_sorted_series_to_decaying_cos(dataseries)


def fit_dataseries_with_two_decaying_cos(dataseries,
                                         fit_drift=False):
    """
    Takes a DataSeries and fits as a function of two decaying cosines,
    plus a sharp exponential and a flat offset using
    fit_dataseries_plus_drift_fitting(dataseries, fitfunction)

    if fit_drift is False:
        returns FitData(dataseries)
    if fit_drift is True:
        returns (no_bkgd_dataseries, FitData(dataseries))
            where no_bkgd_dataseries is dataseries with background
            drift subtracted out.
    """
    if fit_drift:
        return fit_dataseries_plus_drift_fitting(
                dataseries,
                curvefitting.fit_sorted_series_to_two_decaying_cos)
    else:
        return curvefitting.fit_sorted_series_to_two_decaying_cos(dataseries)


# %%
def fit_dataseries_plus_drift_fitting(dataseries, fitfunction):
    """
    Takes a DataSeries and fits as a function of two decaying cosines,
    plus a sharp exponential and a flat offset.
    After the first fit attempt, fits a 5th order polynomial to the

    background, subtracts the polynomial, then repeats. The third fit
    (after two background correction steps) is returned, along with
    a background-corrected version of the DataSeries input.

    Return format: (DataSeries time_offset_background_corrected_data,
                        FitData fit_results)

    Positional arguments:
    dataseries -- DataSeries object containing data to fit.
    function (DataSeries->FitData): function mapping dataseries to FitData
    """
    # "dataseries" needs to be a container object, not a pure iterator,
    # since we need to iterate over it several times
    if iter(dataseries) is iter(dataseries):  # if pure iterator
        raise TypeError("fit_dataseries_plus_drift_fitting: dataseries " +
                        "must be a container object, not an iterator.")

    firstfit = fitfunction(dataseries)
    try:
        data_bkgd = dataseries - firstfit.fitdataseries
        firstbkgdfit = \
            curvefitting.fit_unsorted_series_to_polynomial(data_bkgd, 5)
        data_minusbkgd = dataseries - firstbkgdfit.fitdataseries
        secondfit = fitfunction(data_minusbkgd)
    except AttributeError:  # fit failed
        print("Warning: drift fit failed, aborting...")
        return dataseries, firstfit

    try:
        data_bkgd2 = dataseries - secondfit.fitdataseries
        secondbkgdfit = \
            curvefitting.fit_unsorted_series_to_polynomial(data_bkgd2, 5)
        data_minusbkgd2 = dataseries - secondbkgdfit.fitdataseries
        thirdfit = fitfunction(data_minusbkgd2)
    except AttributeError:  # fit failed
        print("Warning: drift fit failed, aborting...")
        return data_minusbkgd, secondfit

    return data_minusbkgd2, thirdfit


# %%
def get_time_offset_scandata(scandata):
    olddataseries = scandata.dataseries
    newscaninfo = scandata.scaninfo.copy()  # shallow dict copy!
    newscaninfo['pre-time_offset dataseries'] = olddataseries
    return ScanData(scandata.filepath,
                    newscaninfo,
                    get_time_offset_dataseries(olddataseries),
                    scandata.fitdata)


def get_time_offset_dataseries(dataseries):
    time_at_datamax, datamax = \
        curvefitting.get_abs_max_value_data_pt(dataseries)
    excluded_times = []
    for time, value in dataseries.datatuples():  # ensure one continuous seg.
        if not excluded_times:  # if no times yet
            if abs(value) > datamax/2:
                excluded_times.append(time)
        else:  # if already started getting times
            if abs(value) >= datamax/2:
                excluded_times.append(time)
            else:
                break
#    excluded_times = [time for time, value in dataseries.datatuples()
#                      if abs(value) >= datamax/2]
    # insert sanity check if too many times excluded?
    if len(excluded_times) > len(dataseries)/10:
        excluded_times = excluded_times[0]
    time_offset = excluded_times[-1]
    exclusion_start = excluded_times[0] - time_offset
    return DataSeries(((time - time_offset, value)
                       for time, value in
                       dataseries.datatuples(unsorted=True)),
                      excluded_intervals=[(exclusion_start, 0)])
