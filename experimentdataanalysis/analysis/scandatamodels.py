# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:33:16 2016

@author: Michael
"""

import numpy as np

import experimentdataanalysis.analysis.fitfunctions as fitfcns


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class ScanDataModel:
    """
    Stores the parameters and formulas needed to fit scandata to a
    1D model and extract attributes. Should not store any data, although
    it can be modified after creation (e.g. to increase max_fcn_evals
    attribute or adjust the experimental uncertainty assumed).

    Can be inherited from to change the fit model, may only need to change
    __init__ method if first 3 fit parameters are left the same.

    All models should have fcns all_model_fields and get_model_filter_fcns
    """
    def __init__(self, **kwargs):
        self.model_type = "specific_model_name"
        self.dim2_key = "MiddleScanCoord"  # coord spanning each ScanData_SET_
        self.field_index = 0
        self.fitfunction = None
        # params = DUMMY
        self.free_params = [True]
        self.initial_params = [0]
        self.param_bounds = [(0, 0)]
        self.error_thresholds = [None]
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def get_model_filter_fcns(self):
        field_index = self.field_index
        thresholds = self.error_thresholds

        def filter_fit_error_vs_threshold(scandata):
            errors_ok = True
            errors = scandata.fitdata_list[field_index].fitparamstds
            for error, threshold in zip(errors, thresholds):
                if threshold is not None:
                    if error > threshold:
                        errors_ok = False
            return errors_ok

        return [filter_fit_error_vs_threshold]

    # DUMMY, SHOULD BE OVERWRITTEN BASED ON ATTRIBUTE FUNCTIONS
    # IDEALLY SELF-WRITTEN BASED ON ATTRIBUTE DECORATORS...
    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        fields = ["model_attribute"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.model_attribute(params, param_sigmas)])
        return fields, dataseries_list, uncertainty_dataseries_list

    # DUMMY, SHOULD BE OVERWRITTEN BY ONE OR MORE FUNCTIONS
    def model_attribute(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        attribute = model_param_dataseries_list[0]
        attribute_sigma = model_param_uncertainty_dataseries_list[0]
        return attribute, attribute_sigma

    def __eq__(self, other):
        return all([
            self.model_type == other.model_type,
            self.field_index == other.field_index,
            self.free_params == other.free_params,
            self.initial_params == other.initial_params,
            self.param_bounds == other.param_bounds,
            self.error_thresholds == other.error_thresholds,
            self.max_fcn_evals == other.max_fcn_evals])


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class LinearFitModel(ScanDataModel):
    """
    Stores the parameters and formulas needed to fit scandata to a
    1D model and extract attributes. Should not store any data, although
    it can be modified after creation (e.g. to increase max_fcn_evals
    attribute or adjust the experimental uncertainty assumed).

    Can be inherited from to change the fit model, may only need to change
    __init__ method if first 3 fit parameters are left the same.

    All models should have fcns all_model_fields and get_model_filter_fcns
    """
    def __init__(self, **kwargs):
        self.model_type = "fitfcn_simple_line"
        self.dim2_key = "Voltage"
        self.field_index = 2
        self.fitfunction = fitfcns.fitfcn_simple_line
        # params = slope, offset
        self.free_params = [True, True]
        self.initial_params = [0, 0]
        self.param_bounds = [(-np.inf, np.inf), (-np.inf, np.inf)]
        self.error_thresholds = [None, None]
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        fields = ["slope",
                  "offset"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.slopes(params, param_sigmas),
                  self.offsets(params, param_sigmas),
                  ])
        return fields, dataseries_list, uncertainty_dataseries_list

    def slopes(self, model_param_dataseries_list,
               model_param_uncertainty_dataseries_list):
        slopes = model_param_dataseries_list[0]
        slopes_sigma = model_param_uncertainty_dataseries_list[0]
        return slopes, slopes_sigma

    def offsets(self, model_param_dataseries_list,
                model_param_uncertainty_dataseries_list):
        offsets = model_param_dataseries_list[1]
        offsets_sigma = model_param_uncertainty_dataseries_list[1]
        return offsets, offsets_sigma


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class GaussianModel(ScanDataModel):
    """
    Stores the parameters and formulas needed to fit scandata to a
    1D model and extract attributes. Should not store any data, although
    it can be modified after creation (e.g. to increase max_fcn_evals
    attribute or adjust the experimental uncertainty assumed).

    Can be inherited from to change the fit model, may only need to change
    __init__ method if first 3 fit parameters are left the same.

    All models should have fcns all_model_fields and get_model_filter_fcns
    """
    def __init__(self, **kwargs):
        self.model_type = "1d_gaussian_with_linear_offset"
        self.dim2_key = "MiddleScanCoord"
        self.field_index = 1
        self.fitfunction = fitfcns.fitfcn_1d_gaussian_with_linear_offset
        # params = amplitude, x0, sigma, slope, offset
        self.free_params = [True, True, True, True, True]
        self.initial_params = [0.02, 0, 20, 0, 0]
        self.param_bounds = [(0, 1), (-200, 200),
                             (5, 200), (-0.01, 0.01), (-0.01, 0.01)]
        self.error_thresholds = [0.01, 15, 15, None, None]
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        fields = ["gaussian_amplitudes",
                  "gaussian_widths",
                  "gaussian_centers",
                  "gaussian_areas"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.gaussian_amplitudes(params, param_sigmas),
                  self.gaussian_widths(params, param_sigmas),
                  self.gaussian_centers(params, param_sigmas),
                  self.gaussian_areas(params, param_sigmas),
                  ])
        return fields, dataseries_list, uncertainty_dataseries_list

    def gaussian_amplitudes(self, model_param_dataseries_list,
                            model_param_uncertainty_dataseries_list):
        amplitudes = model_param_dataseries_list[0]
        amplitudes_sigma = model_param_uncertainty_dataseries_list[0]
        return amplitudes, amplitudes_sigma

    def gaussian_widths(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        widths = model_param_dataseries_list[2]
        widths_sigma = model_param_uncertainty_dataseries_list[2]
        return widths, widths_sigma

    def gaussian_centers(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        centers = model_param_dataseries_list[1]
        centers_sigma = model_param_uncertainty_dataseries_list[1]
        return centers, centers_sigma

    def gaussian_areas(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        amplitudes = model_param_dataseries_list[0]
        amplitudes_sigma = model_param_uncertainty_dataseries_list[0]
        widths = model_param_dataseries_list[2]
        widths_sigma = model_param_uncertainty_dataseries_list[2]
        # note: uncertainty calc. assumes width, amplitude independent...
        areas = amplitudes*widths
        areas_sigma = \
            ((amplitudes*widths_sigma)**2 + (widths*amplitudes_sigma)**2)**0.5
        return areas, areas_sigma


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class SpinLifetimeModel(ScanDataModel):
    """
    Stores the parameters and formulas needed to fit scandata to a
    1D model and extract attributes. Should not store any data, although
    it can be modified after creation (e.g. to increase max_fcn_evals
    attribute or adjust the experimental uncertainty assumed).

    Can be inherited from to change the fit model, may only need to change
    __init__ method if first 3 fit parameters are left the same.

    All models should have fcns all_model_fields and get_model_filter_fcns
    """
    def __init__(self, **kwargs):
        self.model_type = "fitfcn_two_exp_decay"
        self.dim2_key = "MiddleScanCoord"  # should be 3rd coord, but unknown!
        self.field_index = 3  # gaussian area
        self.fitfunction = fitfcns.fitfcn_two_exp_decay
        # params = pulse_amp1, lifetime1, pulse_amp2, lifetime2, offset
        self.free_params = [True, True, True, True, False]
        self.initial_params = [0.05, 50, 0.05, 2000, 0]
        self.param_bounds = [(0, 1), (1, 200),
                             (0, 1), (10, 1e6), (-0.01, 0.01)]
        self.error_thresholds = [None, None, None, None, None]
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        fields = ["short_amplitudes",
                  "short_lifetimes",
                  "long_amplitudes",
                  "long_lifetimes"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.short_amplitudes(params, param_sigmas),
                  self.short_lifetimes(params, param_sigmas),
                  self.long_amplitudes(params, param_sigmas),
                  self.long_lifetimes(params, param_sigmas),
                  ])
        return fields, dataseries_list, uncertainty_dataseries_list

    def short_amplitudes(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[0]
        params_sigma = model_param_uncertainty_dataseries_list[0]
        return params, params_sigma

    def short_lifetimes(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[1]
        params_sigma = model_param_uncertainty_dataseries_list[1]
        return params, params_sigma

    def long_amplitudes(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[2]
        params_sigma = model_param_uncertainty_dataseries_list[2]
        return params, params_sigma

    def long_lifetimes(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[3]
        params_sigma = model_param_uncertainty_dataseries_list[3]
        return params, params_sigma


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class SinusoidalSpinLifetimeModel(ScanDataModel):
    """
    Stores the parameters and formulas needed to fit scandata to a
    1D model and extract attributes. Should not store any data, although
    it can be modified after creation (e.g. to increase max_fcn_evals
    attribute or adjust the experimental uncertainty assumed).

    Can be inherited from to change the fit model, may only need to change
    __init__ method if first 3 fit parameters are left the same.

    All models should have fcns all_model_fields and get_model_filter_fcns
    """
    def __init__(self, **kwargs):
        self.model_type = "fitfcn_two_exp_sin_decay"
        self.dim2_key = "Voltage"
        self.field_index = 1  # lockin2x
        self.fitfunction = fitfcns.fitfcn_two_exp_sin_decay
        # params = pulse_amp1, lifetime1, pulse_amp2,
        #          lifetime2, osc_period, phase, offset
        self.free_params = [True, True, True, True,
                            True, True, True, True]
        self.initial_params = [0.05, 50, 0.05, 2000, 800, 0, 0, 0]
        self.param_bounds = [(0, 1), (1, 200),
                             (0, 1), (10, 1e6),
                             (600, 1000), (0, 2*np.pi),
                             (-1e-6, 1e-6), (-0.01, 0.01)]
        self.error_thresholds = [None, None, None, None,
                                 None, None, None, None]
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        # TODO: design a model field decorator fcn to make this automatic
        fields = ["short_amplitudes",
                  "short_lifetimes",
                  "long_amplitudes",
                  "long_lifetimes"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.short_amplitudes(params, param_sigmas),
                  self.short_lifetimes(params, param_sigmas),
                  self.long_amplitudes(params, param_sigmas),
                  self.long_lifetimes(params, param_sigmas),
                  ])
        return fields, dataseries_list, uncertainty_dataseries_list

    def short_amplitudes(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[0]
        params_sigma = model_param_uncertainty_dataseries_list[0]
        return params, params_sigma

    def short_lifetimes(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[1]
        params_sigma = model_param_uncertainty_dataseries_list[1]
        return params, params_sigma

    def long_amplitudes(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[2]
        params_sigma = model_param_uncertainty_dataseries_list[2]
        return params, params_sigma

    def long_lifetimes(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[3]
        params_sigma = model_param_uncertainty_dataseries_list[3]
        return params, params_sigma
