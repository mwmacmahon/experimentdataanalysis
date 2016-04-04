# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 00:35:09 2016

@author: Michael
"""

import numpy as np
from scipy.optimize import curve_fit

from experimentdataanalysis.analysis.dataclasses \
    import FitData, DataSeries

LASER_REPRATE = 13000


# %%
def two_decays(t, b, c, d, e):
    """
    note: does NOT handle "wrapping around" t values less than zero.
    see wrap_negative_times()
    """
    def single_pulse_fcn(t2):
        return b*np.exp(-t2/c) + d*np.exp(-t2/e)

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE
    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return 0 + thispulse + lastpulse + lastpulse2
#    return a + thispulse + lastpulse + lastpulse2


# %%
def two_decays_sinusoidal(t, a, b, c, d, e, f, g):
    """
    note: does NOT handle "wrapping around" t values less than zero.
    see wrap_negative_times()
    """
    def single_pulse_fcn(t2):
        return b*np.exp(-t2/c) + d*np.exp(-t2/e)*np.cos(f*t2 + g)

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE
    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return a + thispulse + lastpulse + lastpulse2


# %%
def two_decays_two_cosines(t, a, b, c, d, e, f, g, h, k, m, n):
    """
    note: does NOT handle "wrapping around" t values less than zero.
    see wrap_negative_times()
    """
    def single_pulse_fcn(t2):
        return (b*np.exp(-t2/c) +
                d*np.exp(-t2/e)*np.cos(f*t2 + m) +
                g*np.exp(-t2/h)*np.cos(k*t2 + n))

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE
    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return a + thispulse + lastpulse + lastpulse2


# %%
def wrap_negative_times(times):
    """
    note: assumes a subscriptable input, e.g. a list or array
    """
    newtimes = times[:]
    for i, t in enumerate(times):
        if t < 0:
            newtimes[i] = t + LASER_REPRATE
    return newtimes


# %%
def fit_nonlinear_fcn_dataseries_sorted(dataseries, fitfcn, initialparams,
                                        *, maxfev=20000):
    """
    note: assumes a fitfcn that handles arrays of elements
    """
    # note: filter times in advance to get rid of wraparound for t < 0
    # if you filter in fcn instead, expect BIG slowdown in curve fit
    # as every fcn call must evaluate t wraparound element-by-element...
    times = wrap_negative_times(dataseries.xvals(unsorted=False))
    tdata = np.array(times)
    zdata = np.array(list(dataseries.yvals(unsorted=False)))

    with np.errstate(all='ignore'):
        fitparams, covariances = curve_fit(fitfcn, tdata, zdata,
                                           p0=initialparams, maxfev=maxfev)
        fitstds = np.sqrt(np.diag(covariances))
    return fitparams, fitstds


# %%
def get_max_value_data_pt(dataseries):
    time_at_valuemax = None
    valuemax = 0
    for time, value in dataseries.datatuples():
        if value > valuemax:  # NOT abs()
            time_at_valuemax = time
            valuemax = value
    return time_at_valuemax, valuemax


# %%
def get_abs_max_value_data_pt(dataseries):
    time_at_valuemax = None
    valuemax = 0
    for time, value in dataseries.datatuples():
        if abs(value) > valuemax:
            time_at_valuemax = time
            valuemax = abs(value)
    return time_at_valuemax, valuemax


# %%
def fit_sorted_series_to_two_decays(dataseries):
    Tshortguess = 100
    Tlongguess = 1000
    time_at_datamax, datamax = get_max_value_data_pt(dataseries)

    initguess = [datamax*0.15, Tshortguess, datamax*0.1, Tlongguess]
    try:
        fitparams, fitstds = fit_nonlinear_fcn_dataseries_sorted(
                                 dataseries, two_decays, initguess)
    except TypeError as err:
        print("Warning: TypeError during " +
              "fit_sorted_series_to_two_decays:")
        print(err)
        return None
    except RuntimeError as err:
        print("Warning: RuntimeError during " +
              "fit_sorted_series_to_two_decays:")
        print(err)
        return None
    else:
        times = dataseries.xvals(unsorted=True, unfiltered=True)
        fittimes = wrap_negative_times(times)
        fitvalues = []
        for time in fittimes:
            if time in wrap_negative_times(dataseries.xvals(unfiltered=False)):
                fitvalues.append(two_decays(time, *fitparams))
            else:
                # print("excluding time {}".format(time))  # for bug testing
                fitvalues.append(0)
        fitdataseries = DataSeries(
                            zip(times, fitvalues),
                            excluded_intervals=dataseries.excluded_intervals())
        fitparamstring = (
            "Fit Properties\n" +
            "Short Amplitude: {:.5g}\n".format(fitparams[0]) +
            "Short Lifetime: {:.5g}\n".format(fitparams[1]) +
            "Long Amplitude: {:.5g}\n".format(fitparams[2]) +
            "Long Lifetime: {:.5g}\n".format(fitparams[3]) +
            "     +-{:.5g}\n".format((fitstds[3])))
        return FitData(fitparams, fitstds, fitparamstring, fitdataseries)


# %%
def fit_sorted_series_to_decaying_cos(dataseries):
    Tshortguess = 100
    Tlongguess = 1000
    time_at_datamax, datamax = get_max_value_data_pt(dataseries)

    initguess = [0, datamax*0.15, Tshortguess,
                 datamax*0.1, Tlongguess, 2*np.pi/800, 0]
    try:
        fitparams, fitstds = fit_nonlinear_fcn_dataseries_sorted(
                                 dataseries, two_decays_sinusoidal, initguess)
    except TypeError as err:
        print("Warning: TypeError during " +
              "fit_sorted_series_to_decaying_cos:")
        print(err)
        return None
    except RuntimeError as err:
        print("Warning: RuntimeError during " +
              "fit_sorted_series_to_decaying_cos:")
        print(err)
        return None
    else:
        times = dataseries.xvals(unsorted=True, unfiltered=True)
        fittimes = wrap_negative_times(times)
        fitvalues = []
        for time in fittimes:
            if time in wrap_negative_times(dataseries.xvals(unfiltered=False)):
                fitvalues.append(two_decays_sinusoidal(time, *fitparams))
            else:
                # print("excluding time {}".format(time))  # for bug testing
                fitvalues.append(0)
        fitdataseries = DataSeries(
                            zip(times, fitvalues),
                            excluded_intervals=dataseries.excluded_intervals())
        fitparamstring = (
            "Fit Properties\n" +
            "Offset: {:.5g}\n".format(fitparams[0]) +
            "Short Amplitude: {:.5g}\n".format(fitparams[1]) +
            "Short Lifetime: {:.5g}\n".format(fitparams[2]) +
            "Long Amplitude: {:.5g}\n".format(fitparams[3]) +
            "Long Lifetime: {:.5g}\n".format(fitparams[4]) +
            "     +-{:.5g}\n".format((fitstds[4])) +
            "Long Frequency: {:.5g}\n".format(fitparams[5]) +
            "     +-{:.5g}\n".format((fitstds[5])) +
            "Long Cosine Phase: {:.5g}\n".format(fitparams[6]))
        return FitData(fitparams, fitstds, fitparamstring, fitdataseries)


# %%
def fit_sorted_series_to_two_decaying_cos(dataseries):
    Tshortguess = 100
    Tlongguess = 1000
    Tlongerguess = 10000
    time_at_datamax, datamax = get_max_value_data_pt(dataseries)

    initguess = [0, datamax*0.15, Tshortguess,
                 datamax*0.1, Tlongguess, 2*np.pi/800,
                 -datamax*0.05, Tlongerguess, 2*np.pi/800,
                 0, 0]
    try:
        fitparams, fitstds = fit_nonlinear_fcn_dataseries_sorted(
            dataseries, two_decays_two_cosines, initguess)
    except TypeError as err:
        print("Warning: TypeError during " +
              "fit_sorted_series_to_two_decaying_cos:")
        print(err)
        return None
    except RuntimeError as err:
        print("Warning: RuntimeError during " +
              "fit_sorted_series_to_two_decaying_cos:")
        print(err)
        return None
    else:
        times = dataseries.xvals(unsorted=True, unfiltered=True)
        fittimes = wrap_negative_times(times)
        fitvalues = []
        for time in fittimes:
            if time in wrap_negative_times(dataseries.xvals(unfiltered=False)):
                fitvalues.append(two_decays_two_cosines(time, *fitparams))
            else:
                # print("excluding time {}".format(time))  # for bug testing
                fitvalues.append(0)
        fitdataseries = DataSeries(
                            zip(times, fitvalues),
                            excluded_intervals=dataseries.excluded_intervals())
        fitparamstring = (
            "Fit Properties\n" +
            "Offset: {:.5g}\n".format(fitparams[0]) +
            "Short Amplitude: {:.5g}\n".format(fitparams[1]) +
            "Short Lifetime: {:.5g}\n".format(fitparams[2]) +
            "Long Amplitude 1: {:.5g}\n".format(fitparams[3]) +
            "Long Amplitude 2: {:.5g}\n".format(fitparams[6]) +
            "Long Lifetime 1: {:.5g}\n".format(fitparams[4]) +
            "Long Lifetime 2: {:.5g}\n".format(fitparams[7]) +
            "Long Frequency 1: {:.5g}\n".format(fitparams[5]) +
            "Long Frequency 2: {:.5g}\n".format(fitparams[8]) +
            "Long Phase 1: {:.5g}\n".format(fitparams[9]) +
            "Long Phase 2: {:.5g}\n".format(fitparams[10]))
        return FitData(fitparams, fitstds, fitparamstring, fitdataseries)


# %%
def fit_unsorted_series_to_polynomial(dataseries, power):
    zdata = np.array(list(dataseries.yvals(unsorted=True)))
    indices = np.array(list(range(len(zdata))))
    polycoeffs = np.polyfit(indices, zdata, power)
    fitfunc = np.poly1d(polycoeffs)  # function of filtered index #, NOT times

    # should only be nonzero within the filtered range!
    # need to map from filtered, unsorted indices to unfiltered, unsorted times
    fittimes = dataseries.xvals(unfiltered=True, unsorted=True)
    fitvalues = []
    fitindex = 0
    for time in fittimes:
        # no distinguishing before/after t=0, since not fcn of delay time
        if time in wrap_negative_times(dataseries.xvals(unfiltered=False)):
            fitvalues.append(fitfunc(fitindex))
            fitindex += 1
        else:
            fitvalues.append(0)
    fitdata = DataSeries(zip(fittimes, fitvalues),
                         excluded_intervals=dataseries.excluded_intervals())
    return FitData(polycoeffs, None, None, fitdata)
