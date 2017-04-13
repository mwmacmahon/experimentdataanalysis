# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:35:06 2016

The module that handles actual data fitting. Makes heavy use of the
FitData construct from dataclasses.py

@author: vsih-lab
"""

from collections import OrderedDict
from copy import deepcopy
import inspect

import numpy as np
from scipy.optimize import curve_fit
import scipy.ndimage.filters as filters

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData
from experimentdataanalysis.analysis.generalutilities \
    import multiprocessable_map

LASER_REPRATE = 13160  # ps period


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class ScanDataModel():
    """
    Model used to facilitate fitting a function to ScanData.
    
    Attributes designed to be either almost entirely overwritten via keyword
    parameters or the whole class inherited. To avoid long, bloated analysis
    scripts, it is recommended to create a separate module to contain these
    fit models, which can be imported into data analysis scripts.

        e.g.    [in data analysis script:]
                gaussian_model = ScanDataModel(
                                    model_name="Gaussian Model",
                                    =fcn_gaussian_fit,
                                    (etc.),
                                    )  # can leave some attributes as default

        e.g.    [likely in module containing user's fit models:]
                class GaussianModel(ScanDataModel):
                    def __init__(self, **kwargs):
                        # note: all ScanDataModel's attributes must be written!
                        self.model_name = "Gaussian Model"
                        self.yfield = None
                        self.fitfunction = fcn_gaussian_fit
                        (etc.)
                    (other methods are not listed, ScanDataModel's
                     inherited versions will work just fine for most models)

                [in data analysis script:]
                model = my_fit_models_module.GaussianModel()
    """
    def __init__(self, **kwargs):
        self.model_name = "Generic Linear Model"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = lambda x, slope, y0: y0 + slope * x
        self.model_params = \
            OrderedDict([('slope', {'free': True,
                                    'initial value': 0,
                                    'bounds': (-np.inf, np.inf)}),
                         ('y0', {'free': True,
                                 'initial value': 0}),  # param unbounded
                        ])
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        self.ignore_weights = False

        # for ScanDataSets: coord spanning the set
        self.fit_result_scan_coord = "SecondScanCoord"  

        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def get_fit_params_labels(self):
        return np.array(list(self.model_params.keys()))

    def get_fit_params_is_free_param(self):
        is_free = [param_dict.get('free', False)  # default: fixed
                   for param_dict in self.model_params.values()]
        return np.array(is_free, dtype=np.bool)
            
    def get_fit_params_initial_values(self):
        initial_values = [param_dict.get('initial value', 0)  # default: 0
                          for param_dict in self.model_params.values()]
        return np.array(initial_values)
            
    def get_fit_params_bounds(self):
        bounds = [param_dict.get('bounds',
                                 (-np.inf, np.inf))  # default: no bounds
                  for param_dict in self.model_params.values()]
        return np.array(bounds)

    def get_which_fit_params_fit_at_bound(self, fitdata):
        param_fits = fitdata.fitparams
        param_bounds = self.get_fit_params_bounds()
        # TODO: FINISH THIS
        return NotImplementedError('Function not implemented yet!')

    def fitfunction_with_initial_parameters(self, x):
        initial_params = self.get_fit_params_initial_values()
        return self.fitfunction(x, *initial_params)
            
    def copy(self):  # must deep copy all mutable attributes
        kwargs = deepcopy(self.__dict__)
        return self.__class__(**kwargs)
        

# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_model_fit(scandata, model):
    """
    Uses generic_curve_fit to fit a ScanData field with the given model's
    function and parameters.

    The resulting FitData is stored in the ScanData's info dict under the
    key 'fitdata_[yfield]'. 'None' is stored for failed fits.
    The FitData (or None) is also returned.
    """    
    fitdata = scandata_fit(scandata, model.yfield, model.fitfunction,
                           model.get_fit_params_is_free_param(),
                           model.get_fit_params_initial_values(),
                           model.get_fit_params_bounds(), model.max_fcn_evals,
                           model.excluded_intervals, model.ignore_weights)
    if fitdata is not None:
        return FitData(fitdata.fitfunction, fitdata.partialfcn,
                       fitdata.fitparams, fitdata.fitparamstds,
                       model.get_fit_params_labels(), fitdata.fityvals,
                       fitdata.freeparamindices,
                       fitdata.covariancematrix,
                       fitdata.meansquarederror)
    else:
        return None


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_list_model_fit(scandata_list, model, multiprocessing=False):
    """
    Uses generic_curve_fit to fit a ScanData field with the given model's
    function and parameters.

    The resulting FitData is stored in the ScanData's info dict under the
    key 'fitdata_[yfield]'. 'None' is stored for failed fits.
    The FitData (or None) is also returned.
    """
    scandata_list = list(scandata_list)  # ensure an actual list
    if len(scandata_list) == 0:
        return []
    fitdata_list = scandata_list_fit(scandata_list, model.yfield,
                                     model.fitfunction,
                                     model.get_fit_params_is_free_param(),
                                     model.get_fit_params_initial_values(),
                                     model.get_fit_params_bounds(),
                                     model.max_fcn_evals,
                                     model.excluded_intervals,
                                     model.ignore_weights, multiprocessing)
    fixed_fitdata_list = []
    for fitdata in fitdata_list:
        if fitdata is not None:
            fixed_fitdata_list.append(
                        FitData(fitdata.fitfunction, fitdata.partialfcn,
                                fitdata.fitparams, fitdata.fitparamstds,
                                model.get_fit_params_labels(),
                                fitdata.fityvals,
                                fitdata.freeparamindices,
                                fitdata.covariancematrix,
                                fitdata.meansquarederror))
        else:
            fixed_fitdata_list.append(None)
    return fixed_fitdata_list


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_fit(scandata, yfield, fitfunction, free_params,
                 initial_params, param_bounds, max_fcn_evals=20000,
                 excluded_intervals=None, ignore_weights=False):
    """
    Uses generic_curve_fit to fit a ScanData field with the given function
    and parameters.
    The resulting FitData is stored in the ScanData's info dict under the
    key 'fitdata_[yfield]'. 'None' is stored for failed fits.
    The FitData (or None) is also returned.
    """
    if yfield is None:
        yfield = scandata.yfield
    xvals, yvals, yerrvals = scandata.get_field_xyyerr(yfield)
    fitdata = generic_curve_fit(xvals, yvals, yerrvals, fitfunction,
                                free_params, initial_params,
                                param_bounds, max_fcn_evals,
                                excluded_intervals, ignore_weights)
    setattr(scandata, 'fitdata_' + yfield, fitdata)
    return fitdata


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_list_fit(scandata_list, yfield, fitfunction,
                      free_params, initial_params, param_bounds,
                      max_fcn_evals=20000, excluded_intervals=None,
                      ignore_weights=False, multiprocessing=False):
    """
    Takes a list (or other iterable) of ScanData and performs a fit on all
    of them.

    Optionally, each parameter after scandata_list (except multiprocessing)
    can be replaced with a list of equal length to scandata_list. This will
    fit each ScanData with the corresponding parameters from the other lists
    instead of using the same parameters for each fit.

    FitData resulting from fits are stored in the info dict of each ScanData
    under the key 'fitdata_[yfield]'. 'None' is stored for failed fits.
    These FitData (and Nones) are also returned as a list.
    """
    scandata_list = list(scandata_list)  # ensure an actual list
    if len(scandata_list) == 0:
        return []
    try:  # assume first all args are equal length, correspond to ScanData 1:1
        assert len(fitfunction) == len(scandata_list)
        assert len(yfield) == len(scandata_list)
        assert len(free_params) == len(scandata_list)
        assert len(initial_params) == len(scandata_list)
        assert len(param_bounds) == len(scandata_list)
        assert len(max_fcn_evals) == len(scandata_list)
        assert len(excluded_intervals) == len(scandata_list)
        assert len(ignore_weights) == len(scandata_list)
        yfield = [field if field is not None else scandata.yfield
                  for scandata, field in zip(scandata_list, yfield)]
        xyyerr_lists = \
            list(zip(*[scandata.get_field_xyyerr(field)
                       for scandata, field in zip(scandata_list, yfield)]))
        assert len(xyyerr_lists[0]) == len(scandata_list)  # zip error check
        input_args_list = list(zip(*xyyerr_lists,
                                   fitfunction, free_params, initial_params,
                                   param_bounds, max_fcn_evals,
                                   excluded_intervals, ignore_weights))
        assert len(input_args_list) == len(scandata_list)  # zip error check
        yfield_list = yfield
    except (TypeError, AssertionError):  # assume all scandata share fit params
        if yfield is None:
            yfield = scandata_list[0].yfield
        input_args_list = [[*scandata.get_field_xyyerr(yfield),
                            fitfunction, free_params, initial_params,
                            param_bounds, max_fcn_evals,
                            excluded_intervals, ignore_weights]
                           for scandata in scandata_list]
        yfield_list = [yfield] * len(scandata_list)
    fitdata_list = multiprocessable_map(generic_curve_fit, input_args_list,
                                        multiprocessing)
    for scandata, fitdata, field in \
                            zip(scandata_list, fitdata_list, yfield_list):
        setattr(scandata, 'fitdata_' + field, fitdata)
    return fitdata_list


# %% NEEDS TEST, SPHINX DOCUMENTATION
def generic_curve_fit(xvals, yvals, yerrvals, fitfunction, free_params,
                      initial_params, param_bounds, max_fcn_evals=20000,
                      excluded_intervals=None, ignore_weights=False):
    """
    Takes a ScanData and fits one field as a function of an arbitrary single-
    valued scalar function whose first parameter is assumed to correspond
    to "x" values and whose output is assumed to correspond to "y" values.

    Positional arguments:
    :xvals:
    :yvals:
    :yerrvals:
    :function (x,...)->y: function mapping a numpy array x to numpy array y
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
    # don't want to risk changing the arrays in-place!
    original_xvals = xvals.T
    original_yvals = yvals.T
    xvals = np.array(xvals).T
    yvals = np.array(yvals).T
    if yerrvals is not None:
        yerrvals = np.array(yerrvals).T

    # should check validiy of all arguments, inc. nonzero free params, etc.
    # also convert "lower bound = upper bound" params to fixed.
#    if len(xvals) != len(yvals):
#        raise ValueError("len(xvals) != len(yvals)")
    if yerrvals is not None and not ignore_weights:  
        if len(yerrvals) != len(yvals):
            raise ValueError("len(yvals) != len(yerr)")
        else:
            use_errors = True
    else:
        use_errors = False

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
                if alternate_fixed_arg_list is not None:
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
            if use_errors:
                yerrvals = yerrvals[included_indices]
        except ValueError:
            print("Warning: curve_fit given invalid exclusion bounds. " +
                  "Exclusion bounds must be an iterable containing 2-tuples " +
                  "of form (x_start, x_end)")
    with np.errstate(all='ignore'):
        try:
            rawfitparams, rawcovariances = curve_fit(partialfcn, xvals, yvals,
                                                     p0=partial_initial_params,
                                                     sigma=yerrvals,
                                                     absolute_sigma=use_errors,
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
    fityvals = fitfunction(original_xvals, *fitparams)
    meansquarederror = (1./len(original_xvals))* \
                        sum((y_fit - y_real)**2 for y_fit, y_real in
                            zip(fityvals, original_yvals))
    return FitData(fitfunction, partialfcn, fitparams, fitparamstds,
                   len(fitparams) * ['?'], fityvals,
                   free_param_indices, rawcovariances, meansquarederror)


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
                err_msg = ("fit_function_to_dataseries: inconsistent " +
                           "parameter counts, recall first argument " +
                           "is x, not a parameter.")
                err_msg += "\n # Free params: {}".format(len(free_params))
                err_msg += "\n # Non-X params: {}".format(num_nonx_args)
                raise TypeError(err_msg)
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
                                "with the corresponding param_bounds " +
                                "limits: {} vs {}".format(param_guess, tupl))
            if tupl[0] == param_guess and param_guess == tupl[1]:
                final_free_param_indices.remove(index)  # NOT pop()!
    return final_free_param_indices


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_iterable_sort(scandata_iterable,
                           primary_key, secondary_key, numeric_sort=True):
    """
    right now, no mixed numeric/non-numeric keys, because effort
    """
    # Subfunction to use as sort key
    def scandatasortfcn_strings(scandata_2_tuple):
        scandata, _ = scandata_2_tuple
        scaninfo = scandata.info
        try:
            return (str(scaninfo[primary_key]),
                    str(scaninfo[secondary_key]))
        except KeyError:
            try:
                return (str(scaninfo[primary_key]),
                        "")
            except KeyError:
                try:
                    return ("",
                            str(scaninfo[secondary_key]))
                except KeyError:
                    return ("", "")
        except AttributeError:
            print("scandata_iterable_sort: ScanData expected as list element.")
            return ("", "")

    # Subfunction to use as sort key
    def scandatasortfcn_numerical(scandata_2_tuple):
        scandata, _ = scandata_2_tuple
        scaninfo = scandata.info
        try:
            return (float(scaninfo[primary_key]),
                    float(scaninfo[secondary_key]))
        except KeyError:
            try:
                return (float(scaninfo[primary_key]),
                        99999999)
            except KeyError:
                try:
                    return (99999999,
                            float(scaninfo[secondary_key]))
                except KeyError:
                    return (99999999, 99999999)
                except ValueError:
                    print("scandata_iterable_sort: numerical_sort flag on, " +
                          "numerical sort keys only!")
                    return (99999999, 99999999)
            except ValueError:
                print("scandata_iterable_sort: numerical_sort flag on, " +
                      "numerical sort keys only!")
                return (99999999, 99999999)
        except ValueError:
            print("scandata_iterable_sort: numerical_sort flag on, " +
                  "numerical sort keys only!")
            return (99999999, 99999999)
        except AttributeError:
            print("scandata_iterable_sort: ScanData expected as list element.")
            return (99999999, 99999999)

    scandata_list = list(scandata_iterable)
    index_ordering = range(len(scandata_list))
    if numeric_sort:
        key_fcn = scandatasortfcn_numerical
    else:
        key_fcn = scandatasortfcn_strings

    scandata_list, index_ordering = zip(*sorted(
                                        zip(scandata_list, index_ordering),
                                        key=key_fcn))
    return scandata_list, index_ordering


# %% ----------DATASERIES PROCESSING LIBRARY----------

# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_max_yval_xypair(xvals, yvals):
    """
    """
    xval_at_yvalmax = xvals[0]
    yvalmax = yvals[0]
    for xval, yval in zip(xvals, yvals):
        if yval > yvalmax:
            xval_at_yvalmax = xval
            yvalmax = yval
    return xval_at_yvalmax, yvalmax


# %% TODO: NEEDS TEST, SPHINX DOCUMENTATION
def get_processed_scandata(scandata,
                           process_fcn_list, process_fcn_params_list=[],
                           only_process_these_fields=None):
    """
    """
    new_scandata = scandata.copy()  # read from old, replace in new
    for process_fcn, process_fcn_params in \
                            zip(process_fcn_list, process_fcn_params_list):
        process_fcn(new_scandata, *process_fcn_params)
    return new_scandata


# %%
def process_scandata_fields(scandata, xyyerr_fcn,
                            xyyerr_fcn_args=[], xyyerr_fcn_kwargs={},
                            only_process_these_fields=None):
    """
    Processes every field set in scandata of form (xfield, yfield) or
    (xfield, yfield, yfielderror) in place.
    
    Process functions should be of form 
        xyyerr_fcn(x, y, yerr, *args, **kwargs)
        \-> returns (new_x, new_y, new_yerr)
    where these return values will be used to overwrite the old x/y/yerr
    arrays in place if given (a return value of None means don't overwrite).
    Note the keyword argument "field_name" is always provided, so make sure
    the xyyerr_fcn function signature contains (..., **kwargs).

    Note the xfield will be overwritten for each set in sequence! Careful...
    """
    if only_process_these_fields is not None:
        field_list = only_process_these_fields
    else:
        field_list = scandata.fields

    # prune x, yerr fields so that we can use the rest to make x, y, yerr sets
    field_list = [field_name for field_name in field_list
                  if '_error' not in field_name
                  if field_name != scandata.xfield]
    x_copy = np.array(scandata.x)
    for field_name in field_list:
        # run xyyerr_fcn, overwrite attributes with returned values
        # (unless returned values are None)
        y, yerr = scandata.get_field_yyerr(field_name)
        if y is None:
            print("Warning: process_scandata_fields: " +
                  "field {} not found in ScanData, ".format(field_name) +
                  "skipping...")
            continue
        new_x, new_y, new_yerr = \
            xyyerr_fcn(x_copy, y, yerr, *xyyerr_fcn_args,
                       field_name=field_name, **xyyerr_fcn_kwargs)
        if new_x is not None:
            scandata.x = new_x
        if new_y is not None:
            scandata.y = new_y
        if new_yerr is not None:
            scandata.yerr = new_yerr
        # delete old fits, generally now invalidated:
        fitdata_key = 'fitdata_' + field_name
        if fitdata_key in scandata.__dict__:
            setattr(scandata, fitdata_key, None)


# %%
def make_scandata_time_delay_positive(scandata, zero_delay_offset=0,
                                      neg_time_weight_multiplier=1.0):
    def get_positive_time_delay_xyyerr(xvals, yvals, yerrvals,
                                       zero_delay_offset, weight_multiplier,
                                       *args, **kwargs):
        # bump negative time data points up by the laser rep period
        newxvals = np.array(xvals)
        neg_time_indices = []
        for i, x in enumerate(newxvals):
            if x < zero_delay_offset:
                newxvals[i] = x + LASER_REPRATE
                neg_time_indices.append(i)

        # change weights of formerly negative time data points:
        if yerrvals is not None:
            newyerrvals = np.array(yerrvals)  # make copy
            for i in neg_time_indices:  # recall actual weight is error**-2
                newyerrvals[i] *= np.sqrt(1.0/weight_multiplier)
            return newxvals, None, newyerrvals
        else:
            return newxvals, None, None  # only change xvals

    xyyerr_fcn = get_positive_time_delay_xyyerr
    xyyerr_fcn_args = [zero_delay_offset, neg_time_weight_multiplier]
    xyyerr_fcn_kwargs = {}
    process_scandata_fields(scandata, xyyerr_fcn, xyyerr_fcn_args,
                            xyyerr_fcn_kwargs)


# %%
def make_scandata_phase_continuous(scandata):
    def get_continuous_phase_xyyerr(xvals, yvals, yerrvals, *args, **kwargs):
        map_to_sorted = xvals.argsort()
        map_to_unsorted = map_to_sorted.argsort()
        sorted_yvals = np.array(yvals[map_to_sorted])
        sorted_yvals = np.unwrap(sorted_yvals)
        reunsorted_yvals = np.array(sorted_yvals[map_to_unsorted])
        return None, reunsorted_yvals, None  # only change yvals

    xyyerr_fcn = get_continuous_phase_xyyerr
    xyyerr_fcn_args = []
    xyyerr_fcn_kwargs = {}
    fields_to_process = [field_name for field_name in scandata.fields
                         if 'phase' in field_name]
    process_scandata_fields(scandata, xyyerr_fcn, xyyerr_fcn_args,
                            xyyerr_fcn_kwargs, fields_to_process)
                                        

# %%
def gaussian_smooth_scandata(scandata, fields_to_process,
                             gaussian_width, edge_handling='reflect',
                             subtract_smoothed_data_from_original=False):
    """
    Returns data smoothed over with a gaussian integral, unless the flag
    subtract_smoothed_gaussian is set to True, in which case it returns the
    original data with the smoothed gaussian subtracted. Note that the
    scandata is processed in-place; though a reference to the scandata
    is returned identical to the one provided.

    Warning: expects and assumes evenly spaced xvalues, not necessarily in
    order though. Gets scaling of data from average timestep, so this means
    DO NOT USE AFTER ALREADY ADDING 13NS TO NEGATIVE DELAY TIMES.
    """
    def get_gaussian_smoothed_xyyerr(xvals, yvals, yerrvals, gaussian_width,
                                     edge_handling, subtract_smoothed_data,
                                     *args, **kwargs):
        map_to_sorted = xvals.argsort()
        map_to_unsorted = map_to_sorted.argsort()
        sorted_xvals = np.array(xvals[map_to_sorted])
        sorted_yvals = np.array(yvals[map_to_sorted])

        xdiffs = np.abs(sorted_xvals[1:] - sorted_xvals[:-1])
        x_scale_factor = np.mean(xdiffs)
        if max(xdiffs) > 1.5*x_scale_factor:
            raise ValueError("get_gaussian_smoothed_dataseries: requires a " +
                             "dataseries with roughly even x-spacing " +
                             "(including no discontinuities in spacing)")
        sigma = gaussian_width/x_scale_factor
        smoothed_sorted_yvals = filters.gaussian_filter1d(sorted_yvals, sigma,
                                                          mode=edge_handling)
        if subtract_smoothed_data:
            sorted_new_yvals = sorted_yvals - smoothed_sorted_yvals
        else:
            sorted_new_yvals = smoothed_sorted_yvals

        new_yvals = np.array(sorted_new_yvals[map_to_unsorted])
        return None, new_yvals, None  # only change yvals

    xyyerr_fcn = get_gaussian_smoothed_xyyerr
    xyyerr_fcn_args = [gaussian_width, edge_handling,
                       subtract_smoothed_data_from_original]
    xyyerr_fcn_kwargs = {}
    process_scandata_fields(scandata, xyyerr_fcn, xyyerr_fcn_args,
                            xyyerr_fcn_kwargs, fields_to_process)
    return scandata

