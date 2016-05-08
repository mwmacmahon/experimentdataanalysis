# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:35:06 2016

@author: vsih-lab
"""

import inspect

from experimentdataanalysis.analysis.generalutilities \
    import multiprocessable_map
from experimentdataanalysis.analysis.dataclasses \
    import FitData, FitFunc, ScanData, DataSeries


# %% NEEDS TEST, SPHINX DOCUMENTATION
def dataseries_fit(dataseries, fitfunction,
                   free_params, initial_params, param_bounds):
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
    :rtype: FitData
    """
    # "dataseries" needs to be a true dataseries object or at least a
    # container object, since we need to iterate over it several times
    if iter(dataseries) is iter(dataseries):  # if pure iterator
        raise TypeError("fit_function_to_dataseries: first argument" +
                        "is an iterator, must be a dataseries object.")

    # should check validiy of all arguments, inc. nonzero free params, etc.
    fcn_sig = inspect.signature(fitfunction)
    num_nonx_args = len(fcn_sig.parameters) - 1
    free_param_indices = \
        check_fit_parameter_consistency(num_nonx_args, free_params,
                                        initial_params, param_bounds)

    # should create a partial function with _only_ free parameters
    # convert initial_params, param_bounds to versions omitting fixed params
    def partialfcn(x, *free_args):
        if len(free_args) != len(free_param_indices):
            raise TypeError("fit_function_to_dataseries: )
        free_ind = 0
        all_args = [x]
        for ind in range(num_nonx_args):
            if ind in free_param_indices:
                all_args.append(free_args[free_ind])
                free_ind += 1
            else:
                all_args.append(initial_params[ind])
        return fitfunction(*all_args)

    # fit partial function with curve_fit

    # convert fit outputs to reintroduce fixed parameters

    # convert results to FitData object and return it


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
def get_max_yval_xypair(dataseries):
    xval_at_yvalmax = dataseries.xvals()[0]
    yvalmax = dataseries.yvals()[0]
    for xval, yval in dataseries.datatuples():
        if yval > yvalmax:
            xval_at_yvalmax = xval
            yvalmax = yval
    return xval_at_yvalmax, yvalmax


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_x_offset_dataseries_TRKRstyle(dataseries):
    xval_at_datamax, datamax = get_max_yval_xypair(dataseries)
    excluded_xvals = []
    for xval, yval in dataseries.datatuples():  # ensure one continuous seg.
        if not excluded_xvals:  # if no xvals yet
            if abs(yval) > datamax/2:
                excluded_xvals.append(xval)
        else:  # if already started getting xvals
            if abs(yval) >= datamax/2:
                excluded_xvals.append(xval)
            else:
                break
#    excluded_xvals = [xval for xval, yval in dataseries.datatuples()
#                      if abs(yval) >= datamax/2]
    # insert sanity check if too many xvals excluded?
    if len(excluded_xvals) > len(dataseries)/10:
        excluded_xvals = [excluded_xvals[0]]
    xval_offset = excluded_xvals[-1]
    exclusion_start = excluded_xvals[0] - xval_offset
    return DataSeries(((xval - xval_offset, yval)
                       for xval, yval in
                       dataseries.datatuples(unsorted=True)),
                      excluded_intervals=[(exclusion_start, 0)])
