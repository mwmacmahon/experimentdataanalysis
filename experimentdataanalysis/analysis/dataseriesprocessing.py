# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:35:06 2016

The module that handles actual data fitting. Makes heavy use of the
FitData construct from dataclasses.py

@author: vsih-lab
"""

import inspect

import numpy as np
from scipy.optimize import curve_fit
import scipy.ndimage.filters as filters

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries


LASER_REPRATE = 13160  # ps period


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_dataseries_fit(scandata, field_index, fitfunction,
                            free_params, initial_params, param_bounds,
                            max_fcn_evals=20000, excluded_intervals=None):
    """
    """
    new_fitdata = dataseries_fit(scandata.dataseries_list[field_index],
                                 fitfunction, free_params,
                                 initial_params, param_bounds,
                                 scandata.error_dataseries_list[field_index],
                                 max_fcn_evals=max_fcn_evals,
                                 excluded_intervals=excluded_intervals)
    new_scandata_fitdata_list = list(scandata.fitdata_list)
    new_scandata_fitdata_list[field_index] = new_fitdata
    return ScanData(scandata.fields,
                    [scaninfo.copy() for scaninfo in scandata.scaninfo_list],
                    scandata.dataseries_list,
                    scandata.error_dataseries_list,
                    new_scandata_fitdata_list)


# %% NEEDS TEST, SPHINX DOCUMENTATION
def dataseries_fit(dataseries, fitfunction,
                   free_params, initial_params, param_bounds,
                   weights_dataseries=None, max_fcn_evals=20000,
                   excluded_intervals=None):
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

    Optional arguments:
    :excluded_intervals: list/tuple containing pairs (x_start, x_end)
        corresponding to x-intervals in which data should not be used for fit

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

    # fit partial function with curve_fit, use sorted data
    xvals, yvals = dataseries.data(unsorted=False)
    if weights_dataseries:  # first ensure same xvals range
        weights_xvals = weights_dataseries.xvals(unsorted=False)
        if np.array_equal(weights_xvals, xvals):
            weights = weights_dataseries.yvals(unsorted=False)
            use_weights = True
        else:
            weights = None
            use_weights = False
    else:
        weights = None
        use_weights = False
    if excluded_intervals is not None:
        try:
            excluded_indices = []
            for x_start, x_end in excluded_intervals:
                is_excluded = np.logical_and(xvals > x_start, xvals < x_end)
                new_excluded_indices = list(np.argwhere(is_excluded).flatten())
                excluded_indices += new_excluded_indices
            included_indices = [index for index in range(len(xvals))
                                if index not in excluded_indices]
            xvals = xvals[included_indices]
            yvals = yvals[included_indices]
            if use_weights:
                weights = weights[included_indices]
        except ValueError:
            print("Warning: curve_fit given invalid exclusion bounds. " +
                  "Exclusion bounds must be an iterable containing 2-tuples " +
                  "of form (x_start, x_end)")
    with np.errstate(all='ignore'):
        try:
            rawfitparams, rawcovariances = curve_fit(partialfcn, xvals, yvals,
                                                 p0=partial_initial_params,
                                                 sigma=weights,
                                                 absolute_sigma=use_weights,
                                                 bounds=partial_bounds,
                                                 method='trf',
                                                 max_nfev=max_fcn_evals)
            rawfitstds = np.sqrt(np.diag(rawcovariances))
        except RuntimeError:
            print("Warning: curve_fit failed to converge...")
            return None

    # convert fit outputs to reintroduce fixed parameters
    fitparams = get_full_param_list(rawfitparams)
    fitparamstds = get_full_param_list(rawfitstds, [0]*num_nonx_args)

    # convert results to FitData object and return it
    # note fitdataseries should have same sorting and filters
    xvals = dataseries.xvals(raw=True)
    fitdataseries = DataSeries(xvals, fitfunction(xvals, *fitparams))
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
#==============================================================================
# def get_positive_time_delay_scandata(scandata, zero_delay_offset=0):
#     new_dataseries_list = []
#     new_error_dataseries_list = []
#     for dataseries in scandata.dataseries_list:
#         if dataseries is not None:
#             new_dataseries = \
#                 get_positive_time_delay_dataseries(dataseries,
#                                                    zero_delay_offset)
#             new_dataseries_list.append(new_dataseries)
#         else:
#             new_dataseries_list.append(dataseries)
#     for error_dataseries in scandata.error_dataseries_list:
#         if error_dataseries is not None:
#             new_error_dataseries = \
#                 get_positive_time_delay_dataseries(error_dataseries,
#                                                    zero_delay_offset)
#             new_error_dataseries_list.append(new_error_dataseries)
#         else:
#             new_error_dataseries_list.append(error_dataseries)
#     return ScanData(scandata.fields,
#                     [scaninfo.copy() for scaninfo in scandata.scaninfo_list],
#                     new_dataseries_list,
#                     new_error_dataseries_list,
#                     [None for field in scandata.fields])  # invalidates fits
# 
# 
# def get_positive_time_delay_dataseries(dataseries, zero_delay_offset=0):
#     oldxvals, yvals = dataseries.data(raw=True)
#     newxvals = oldxvals.copy()
#     for i, x in enumerate(newxvals):
#         if x < zero_delay_offset:
#             newxvals[i] = x + LASER_REPRATE
#     return DataSeries(newxvals, yvals)
#==============================================================================


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
#==============================================================================
# def get_x_offset_dataseries_TRKRstyle(dataseries, otherseries=[]):
#     """
#     """
#     xval_at_datamax, datamax = get_max_yval_xypair(dataseries)
#     excluded_xvals = []
#     for xval, yval in dataseries.datatuples():  # ensure one continuous seg.
#         if not excluded_xvals:  # if no xvals yet
#             if  yval > datamax*2/3:
#                 excluded_xvals.append(xval)
#         else:  # if already started getting xvals
#             if yval >= datamax*2/3:
#                 excluded_xvals.append(xval)
#             else:
#                 break
# #    excluded_xvals = [xval for xval, yval in dataseries.datatuples()
# #                      if abs(yval) >= datamax/2]
#     # insert sanity check if too many xvals excluded?
#     if len(excluded_xvals) > len(dataseries)/4:
#         excluded_xvals = [excluded_xvals[0]]
#     xval_offset = excluded_xvals[-1]
#     xvals, yvals = dataseries.data(raw=True)
#     new_dataseries = DataSeries(xvals - xval_offset, yvals)
#     print(xval_offset)
# 
#     if otherseries:
#         new_otherseries_list = []
#         for series in otherseries:
#             xvals, yvals = series.data(raw=True)
#         return new_dataseries, new_otherseries_list
#     else:
#         return new_dataseries
#==============================================================================


# %% ----------DATASERIES PROCESSING LIBRARY----------

# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_max_yval_xypair(dataseries):
    """
    """
    xvals, yvals = dataseries.data()
    xval_at_yvalmax = xvals[0]
    yvalmax = yvals[0]
    for xval, yval in dataseries.datatuples():
        if yval > yvalmax:
            xval_at_yvalmax = xval
            yvalmax = yval
    return xval_at_yvalmax, yvalmax


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def process_dataseries_and_error_in_scandata(scandata,
                                             process_fcn, process_fcn_params,
                                             only_process_these_indices=None):
    new_dataseries_list = []
    new_error_dataseries_list = []
    for index, (dataseries, error_dataseries) in enumerate(zip(
            scandata.dataseries_list, scandata.error_dataseries_list)):
        if only_process_these_indices is None or \
                                        index in only_process_these_indices:
            new_dataseries, new_error_dataseries = \
                process_fcn(dataseries, error_dataseries, *process_fcn_params)
            new_dataseries_list.append(new_dataseries)
            new_error_dataseries_list.append(new_error_dataseries)
        else:
            new_dataseries_list.append(dataseries)
            new_error_dataseries_list.append(error_dataseries)
            
    return ScanData(scandata.fields,
                    [scaninfo.copy() for scaninfo in scandata.scaninfo_list],
                    new_dataseries_list,
                    new_error_dataseries_list,
                    [None for field in scandata.fields])  # invalidates fits


# %%
def get_positive_time_delay_scandata(scandata, zero_delay_offset=0,
                                     neg_time_weight_multiplier=1.0):
    def get_positive_time_delay_dataseries(dataseries, error_dataseries,
                                           dataseries_zero_delay_offset,
                                           dataseries_weight_multiplier):
        oldxvals, yvals = dataseries.data(raw=True)
        neg_time_indices = []
        newxvals = oldxvals.copy()
        for i, x in enumerate(newxvals):
            if x < dataseries_zero_delay_offset:
                newxvals[i] = x + LASER_REPRATE
                neg_time_indices.append(i)

        new_dataseries = DataSeries(newxvals, yvals)
        if error_dataseries is not None:
            old_error_yvals = error_dataseries.yvals(raw=True)
            new_error_yvals = old_error_yvals.copy()
            for i in neg_time_indices:  # recall actual weight is error**-2
                new_error_yvals[i] *= np.sqrt(1.0/dataseries_weight_multiplier)
            new_error_dataseries = DataSeries(newxvals, new_error_yvals)
        else:
            new_error_dataseries = None

        return new_dataseries, new_error_dataseries

    process_fcn = get_positive_time_delay_dataseries
    process_fcn_params = [zero_delay_offset, neg_time_weight_multiplier]
    return process_dataseries_and_error_in_scandata(scandata,
                                   process_fcn, process_fcn_params)


# %%
def get_continuous_phase_scandata(scandata):
    def get_continuous_phase_dataseries(dataseries, error_dataseries):
        sortedyvals = dataseries.yvals(unsorted=False)
        sortedyvals_offset = np.unwrap(sortedyvals) - sortedyvals
        new_dataseries = dataseries + sortedyvals_offset  # acts on sorted
        new_error_dataseries = error_dataseries
        return new_dataseries, new_error_dataseries

    process_fcn = get_continuous_phase_dataseries
    process_fcn_params = []
    indices_to_process = [index
                          for index, field in enumerate(scandata.fields)
                          if 'phase' in field]
    return process_dataseries_and_error_in_scandata(scandata, process_fcn,
                                                    process_fcn_params,
                                                    indices_to_process)
                                        

# %%
def get_gaussian_smoothed_scandata(scandata, field_indices_to_process,
                                   gaussian_width, edge_handling='reflect',
                                   subtract_smoothed_data_from_original=False):
    """
    returns data smoothed over with a gaussian integral, unless the flag
    subtract_smoothed_gaussian is set to True, in which case it returns the
    original data with the smoothed gaussian subtracted

    warning: expects and assumes evenly spaced xvalues, not necessarily in
    order though. gets scaling of data from first two xvalues. this means
    DO NOT USE AFTER ADDING 13NS TO NEGATIVE DELAY TIMES.
    """
    def get_gaussian_smoothed_dataseries(dataseries, error_dataseries,
                                         gaussian_width, edge_handling,
                                         subtract_smoothed_data):
        sortedxvals, sortedyvals = dataseries.data(unsorted=False)
        xdiffs = np.abs(sortedxvals[1:] - sortedxvals[:-1])
        x_scale_factor = np.mean(xdiffs)
        if max(xdiffs) > 1.5*x_scale_factor:
            raise ValueError("get_gaussian_smoothed_dataseries: requires a " +
                             "dataseries with roughly even x-spacing " +
                             "(including no discontinuities in spacing)")
        sigma = gaussian_width/x_scale_factor
        smoothedyvals = filters.gaussian_filter1d(sortedyvals, sigma,
                                                  mode=edge_handling)
        # note (DataSeries +- numpy array) acts on _sorted_ DataSeries
        if subtract_smoothed_data:
            new_dataseries = dataseries - smoothedyvals
        else:
            new_dataseries = 0*dataseries + smoothedyvals
        new_error_dataseries = error_dataseries
        return new_dataseries, new_error_dataseries

    process_fcn = get_gaussian_smoothed_dataseries
    process_fcn_params = [gaussian_width, edge_handling,
                          subtract_smoothed_data_from_original]
    return process_dataseries_and_error_in_scandata(scandata, process_fcn,
                                                    process_fcn_params,
                                                    field_indices_to_process)


# %% high pass filter?
def get_high_pass_filtered_scandata(scandata, min_freq_cutoff=0):
    """
    warning: assumes data is in sorted order, will resort otherwise
    TODO: fix this
    """
    print("Warning: get_high_pass_filtered_scandata currently does not " +
          "respect dataseries ordering, consider using " +
          "get_gaussian_smoothed_scandata with a large gaussian width instead")
    def get_high_pass_filtered_dataseries(dataseries, error_dataseries,
                                          dataseries_min_freq_cutoff):
        # extract data and handle error_dataseries & if odd # elements
        xvals, oldyvals = dataseries.data(raw=True)
        if len(xvals) % 2 > 0:
            xvals = xvals[1:]
            oldyvals = oldyvals[1:]
            if error_dataseries is not None:
                errorxvals, erroryvals = error_dataseries.data(raw=True)
                new_error_dataseries = DataSeries(errorxvals[1:],
                                                  erroryvals[1:])
            else:
                new_error_dataseries = None
        else:
            if error_dataseries is not None:
                new_error_dataseries = error_dataseries
            else:
                new_error_dataseries = None

        # now actually calculate fft-filtered yvals
        inverse_sample_rate = xvals[1] - xvals[0]
        f_space_freqs = np.fft.rfftfreq(oldyvals.shape[-1],
                                        d=inverse_sample_rate)
        f_space_yvals = np.fft.rfft(oldyvals)
        f_space_yvals[f_space_freqs < dataseries_min_freq_cutoff] = 0
        newyvals = np.fft.irfft(f_space_yvals)
        new_dataseries = DataSeries(xvals, newyvals)
        return new_dataseries, new_error_dataseries

    process_fcn = get_high_pass_filtered_dataseries
    process_fcn_params = [min_freq_cutoff]
    return process_dataseries_and_error_in_scandata(scandata,
                                   process_fcn, process_fcn_params)

