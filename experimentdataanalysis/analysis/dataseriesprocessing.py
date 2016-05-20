# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:35:06 2016

@author: vsih-lab
"""

import inspect

import numpy as np
from scipy.optimize import curve_fit

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries


LASER_REPRATE = 13100  # ps period


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_dataseries_fit(scandata, dataseries_index, fitfunction,
                            free_params, initial_params, param_bounds,
                            weights_dataseries=None, max_fcn_evals=20000):
    """
    """
    new_fitdata = dataseries_fit(scandata.dataseries[dataseries_index],
                                 fitfunction, free_params,
                                 initial_params, param_bounds,
                                 weights_dataseries, max_fcn_evals)
    new_scandata_fitdatalist = list(scandata.fitdata)
    new_scandata_fitdatalist[dataseries_index] = new_fitdata
    return ScanData(scandata.filepath,
                    scandata.scaninfo.copy(),
                    scandata.fields,
                    scandata.dataseries,
                    tuple(new_scandata_fitdatalist))


# %% NEEDS TEST, SPHINX DOCUMENTATION
def dataseries_fit(dataseries, fitfunction,
                   free_params, initial_params, param_bounds,
                   weights_dataseries=None, max_fcn_evals=20000):
    """
    Takes a DataSeries and fits as a function of an arbitrary single-
    valued scalar function whose first parameter is assumed to correspond
    to "x" values and whose output is assumed to correspond to "y" values.

    Positional arguments:
    :dataseries: DataSeries object containing data points to fit.
    :function (x,...->y): function mapping 1+ scalars to a scalar output
    :free_params: list/tuple describing whether each non-x parameter should
    be free, should have True/False for each parameter.
    :initial_params: list/tuple of parameter starting guesses, should have
    a value for all non-x parameters
    :param_bounds: list/tuple containing parameter upper and lower bounds,
    needs a 2-tuple for all free non-x parameters (and anything for
    fixed parameters)

    Return type:
    :rtype: fitparams, fitparamstds - lists of fitted parameters and their
    uncertainties, where fixed parameters are given an uncertainty (std) of 0
    """
    # "dataseries" needs to be a true dataseries object or at least a
    # container object, since we need to iterate over it several times
    if iter(dataseries) is iter(dataseries):  # if pure iterator
        raise TypeError("fit_function_to_dataseries: first argument " +
                        "is an iterator, must be a dataseries object.")

    # should check validiy of all arguments, inc. nonzero free params, etc.
    # also convert "lower bound = upper bound" params to fixed.
    fcn_sig = inspect.signature(fitfunction)
    num_nonx_args = len(fcn_sig.parameters) - 1
    free_param_indices = \
        check_fit_parameter_consistency(num_nonx_args, free_params,
                                        initial_params, param_bounds)

    # should create a partial function with _only_ free parameters
    # and a function to get full parameter list from only partial list
    def get_full_param_list(free_arg_list, alternate_fixed_arg_list=None):
        if len(free_arg_list) != len(free_param_indices):
            raise TypeError("fit_function_to_dataseries: non-matching "
                            "length of parameters given.")
        free_ind = 0
        all_params = []
        for ind in range(num_nonx_args):
            if ind in free_param_indices:
                all_params.append(free_arg_list[free_ind])
                free_ind += 1
            else:
                if alternate_fixed_arg_list:
                    all_params.append(alternate_fixed_arg_list[ind])
                else:
                    all_params.append(initial_params[ind])
        return all_params

    def partialfcn(x, *free_args):
        all_args = [x] + get_full_param_list(free_args)
        return fitfunction(*all_args)

    # convert initial_params, param_bounds to versions omitting fixed params
    # for use in partial function
    partial_initial_params = [val for ind, val in enumerate(initial_params)
                              if ind in free_param_indices]
    partial_lower_bound = [val for ind, (val, _) in enumerate(param_bounds)
                           if ind in free_param_indices]
    partial_upper_bound = [val for ind, (_, val) in enumerate(param_bounds)
                           if ind in free_param_indices]
    partial_bounds = (partial_lower_bound, partial_upper_bound)

    # fit partial function with curve_fit, use sorted (but filtered) data
    xvals = dataseries.xvals(unsorted=False)
    yvals = dataseries.yvals(unsorted=False)
    if weights_dataseries:  # first ensure same xvals range
        weights_xvals = weights_dataseries.xvals(unsorted=False)
        if weights_xvals == xvals:
            weights = weights_dataseries.yvals(unsorted=False)
            use_weights = True
        else:
            weights = None
            use_weights = False
    else:
        weights = None
        use_weights = False
    with np.errstate(all='ignore'):
        rawfitparams, rawcovariances = curve_fit(partialfcn, xvals, yvals,
                                                 p0=partial_initial_params,
                                                 sigma=weights,
                                                 absolute_sigma=use_weights,
                                                 bounds=partial_bounds,
                                                 max_nfev=max_fcn_evals)
        rawfitstds = np.sqrt(np.diag(rawcovariances))

    # convert fit outputs to reintroduce fixed parameters
    fitparams = get_full_param_list(rawfitparams)
    fitparamstds = get_full_param_list(rawfitstds, [0]*num_nonx_args)

    # convert results to FitData object and return it
    # note fitdataseries should have same sorting and filters
    fitdataseries = DataSeries([(xval, fitfunction(xval, *fitparams))
                                for xval in dataseries.xvals(raw=True)],
                               excluded_intervals=\
                                   dataseries.excluded_intervals())
    meansquarederror = (1./len(dataseries))* \
                        sum((y_fit - y_real)**2 for y_fit, y_real in
                            zip(fitdataseries.yvals(), dataseries.yvals()))
    return FitData(fitparams, fitparamstds,
                   "fitparamstring", fitdataseries, meansquarederror)


# %% NEEDS TEST, SPHINX DOCUMENTATION
def check_fit_parameter_consistency(num_nonx_args, free_params,
                                    initial_params, param_bounds):
    """
    Checks the validity of the parameters sent to dataseries_fit.
    Any problems cause an exception to be thrown, and any jury-rigged
    fixed parameters (lower bound = upper bound = starting guess) are removed
    from the free parameter list. This filtered list is returned.

    Positional arguments:
    :free_params: list/tuple describing whether each non-x parameter should
    be free, should have True/False for each parameter.
    :initial_params: list/tuple of parameter starting guesses, should have
    a value for all non-x parameters
    :param_bounds: list/tuple containing parameter upper and lower bounds,
    needs a 2-tuple for all free non-x parameters (and anything for
    fixed parameters)

    Return type:
    :rtype: list containing indices of free parameters of fit function
    """
    if (len(free_params) != num_nonx_args or
        len(initial_params) != num_nonx_args or
            len(param_bounds) != num_nonx_args):
                raise TypeError("fit_function_to_dataseries: inconsistent " +
                                "parameter counts, recall first argument " +
                                "is x, not a parameter.")
    free_params_bool = [bool(x) for x in free_params]
    free_param_indices = [i for i, x in enumerate(free_params_bool)
                          if x is True]
    if len(free_param_indices) == 0:
        raise TypeError("fit_function_to_dataseries: not enough free " +
                        "parameters, recall first argument is x, " +
                        "not a parameter.")
    final_free_param_indices = free_param_indices[:]
    for index in free_param_indices:
        tupl = param_bounds[index]
        try:
            tupl_len = len(tupl)
        except TypeError:
            raise TypeError("fit_function_to_dataseries: param_bounds " +
                            "must contain only 2-tuples for free parameters")
        else:
            if tupl_len != 2:
                raise TypeError("fit_function_to_dataseries: param_bounds " +
                                "must contain only 2-tuples for free " +
                                "parameters")
            if tupl[0] > tupl[1]:
                raise TypeError("fit_function_to_dataseries: param_bounds " +
                                "contains a lower bound above the " +
                                "corresponding upper bound")
            param_guess = initial_params[index]
            if tupl[0] > param_guess or param_guess > tupl[1]:
                raise TypeError("fit_function_to_dataseries: initial_params " +
                                "contains a starting estimate inconsistent " +
                                "with the corresponding param_bounds limits")
            if tupl[0] == param_guess and param_guess == tupl[1]:
                final_free_param_indices.remove(index)  # NOT pop()!
    return final_free_param_indices


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_positive_time_delay_dataseries(dataseries, zero_delay_offset=0):
    """
    """
    oldxvals = dataseries.xvals(raw=True)
    oldyvals = dataseries.yvals(raw=True)
    newxvals = oldxvals[:]
    for i, x in enumerate(oldxvals):
        if x < zero_delay_offset:
            newxvals[i] = x + LASER_REPRATE
    return DataSeries(zip(newxvals, oldyvals),
                      excluded_intervals=dataseries.excluded_intervals())


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_max_yval_xypair(dataseries):
    """
    """
    xval_at_yvalmax = dataseries.xvals()[0]
    yvalmax = dataseries.yvals()[0]
    for xval, yval in dataseries.datatuples():
        if yval > yvalmax:
            xval_at_yvalmax = xval
            yvalmax = yval
    return xval_at_yvalmax, yvalmax


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_x_offset_dataseries_TRKRstyle(dataseries, otherseries=[]):
    """
    """
    xval_at_datamax, datamax = get_max_yval_xypair(dataseries)
    excluded_xvals = []
    for xval, yval in dataseries.datatuples():  # ensure one continuous seg.
        if not excluded_xvals:  # if no xvals yet
            if  yval > datamax*2/3:
                excluded_xvals.append(xval)
        else:  # if already started getting xvals
            if yval >= datamax*2/3:
                excluded_xvals.append(xval)
            else:
                break
#    excluded_xvals = [xval for xval, yval in dataseries.datatuples()
#                      if abs(yval) >= datamax/2]
    # insert sanity check if too many xvals excluded?
    if len(excluded_xvals) > len(dataseries)/4:
        excluded_xvals = [excluded_xvals[0]]
    xval_offset = excluded_xvals[-1]
    exclusion_start = excluded_xvals[0] - xval_offset
    new_dataseries = DataSeries([(xval - xval_offset, yval)
                                 for xval, yval in
                                 dataseries.datatuples(raw=True)],
                                excluded_intervals=[(exclusion_start, 0)])
    print(xval_offset)
    if otherseries:
        new_otherseries_list = []
        for series in otherseries:
            new_otherseries_list.append(
                    DataSeries([(xval - xval_offset, yval)
                                 for xval, yval in
                                 series.datatuples(raw=True)],
                                excluded_intervals=[(exclusion_start, 0)]))
        return new_dataseries, new_otherseries_list
    else:
        return new_dataseries