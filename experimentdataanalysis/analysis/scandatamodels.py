# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:33:16 2016

@author: Michael
"""

import numpy as np

import experimentdataanalysis.analysis.fitfunctions as fitfcns


gfactorCONSTANT = 71.44773  # ps*Tesla, = 2*pi*hbar/bohr magneton
GFACTORCONSTANT2 = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar


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
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")
            
        fields = ["model_attribute"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.model_attribute(params, param_sigmas)])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    # DUMMY, SHOULD BE OVERWRITTEN BY ONE OR MORE FUNCTIONS
    def model_attribute(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        attribute = model_param_dataseries_list[0]
        attribute_sigma = model_param_uncertainty_dataseries_list[0]
        return attribute, attribute_sigma

    def copy(self):  # must deep copy all mutable attributes
        kwargs = self.__dict__.copy()
        kwargs['free_params'] = list(self.free_params)
        kwargs['initial_params'] = list(self.initial_params)
        kwargs['param_bounds'] = [tuple(pair)
                                  for pair in self.param_bounds]
        kwargs['error_thresholds'] = list(self.error_thresholds)
        if self.excluded_intervals is not None:
            kwargs['excluded_intervals'] = [tuple(pair)
                                            for pair in self.excluded_intervals]
        return self.__class__(**kwargs)

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
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")
            
        fields = ["slope",
                  "offset"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.slope(params, param_sigmas),
                  self.offset(params, param_sigmas),
                  ])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    def slope(self, model_param_dataseries_list,
               model_param_uncertainty_dataseries_list):
        slope = model_param_dataseries_list[0]
        slope_sigma = model_param_uncertainty_dataseries_list[0]
        return slope, slope_sigma

    def offset(self, model_param_dataseries_list,
                model_param_uncertainty_dataseries_list):
        offset = model_param_dataseries_list[1]
        offset_sigma = model_param_uncertainty_dataseries_list[1]
        return offset, offset_sigma


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
        self.field_index = 0  # key field, probably lockin2x
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
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")
            
        fields = ["gaussian_amplitude",
                  "gaussian_widths",
                  "gaussian_centers",
                  "gaussian_areas"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.gaussian_amplitude(params, param_sigmas),
                  self.gaussian_widths(params, param_sigmas),
                  self.gaussian_centers(params, param_sigmas),
                  self.gaussian_areas(params, param_sigmas),
                  ])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    def gaussian_amplitude(self, model_param_dataseries_list,
                            model_param_uncertainty_dataseries_list):
        amplitude = model_param_dataseries_list[0]
        amplitude_sigma = model_param_uncertainty_dataseries_list[0]
        return amplitude, amplitude_sigma

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
        amplitude = model_param_dataseries_list[0]
        amplitude_sigma = model_param_uncertainty_dataseries_list[0]
        widths = model_param_dataseries_list[2]
        widths_sigma = model_param_uncertainty_dataseries_list[2]
        # note: uncertainty calc. assumes width, amplitude independent...
        areas = amplitude*widths
        areas_sigma = \
            ((amplitude*widths_sigma)**2 + (widths*amplitude_sigma)**2)**0.5
        return areas, areas_sigma


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class RSAFieldScanModel(ScanDataModel):
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
        self.model_type = "fitfcn_rsa_field_scan"
        self.dim2_key = "MiddleScanCoord"  # should be 3rd coord, but unknown!
        self.field_index = 0  # key field, probably lockin2x
        self.fitfunction = fitfcns.fitfcn_rsa_field_scan
        self.model_params = ["pulse_amplitude", "lifetime", "gfactor",
                             "field_offset", "drift_velocity", "phase",
                             "slope", "offset"]
        # params = num_pulses, delay_time,
        #          pulse_amplitude, lifetime, gfactor,
        #          field_offset, drift_velocity, phase, slope, offset
        self.free_params = [False, False,
                            True, True, True,
                            True, False, True, True, True]
        self.initial_params = [40, -160,
                               0.05, 2000, 0.43,
                               0, 0, 0, 0, 0]
        self.param_bounds = [(1, 1000), (-1000, 10000),
                             (0, 1), (10, 1e9), (0.3, 0.6),
                             (-0.1, 0.1), (-1e3, 1e3),
                             (-np.pi, np.pi), (-0.05, 0.05), (-0.01, 0.01)]
        self.error_thresholds = [None, None,
                                 None, None, None,
                                 None, None, None, None, None]
        self.max_fcn_evals = 10000
        self.excluded_intervals = None
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")
            
        fields = ["pulse_amplitude",
                  "lifetime",
                  "gfactor",
                  "osc_period_200mT",
                  "osc_period_300mT",
                  "field_offset",
                  "drift_velocity"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.pulse_amplitude(params, param_sigmas),
                  self.lifetime(params, param_sigmas),
                  self.gfactor(params, param_sigmas),
                  self.osc_period_200mT(params, param_sigmas),
                  self.osc_period_300mT(params, param_sigmas),
                  self.field_offset(params, param_sigmas),
                  self.drift_velocity(params, param_sigmas),
                  ])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    def pulse_amplitude(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[2]
        params_sigma = model_param_uncertainty_dataseries_list[2]
        return params, params_sigma

    def lifetime(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[3]
        params_sigma = model_param_uncertainty_dataseries_list[3]
        return params, params_sigma

    def gfactor(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[4]
        params_sigma = model_param_uncertainty_dataseries_list[4]
        return params, params_sigma

    def osc_period_200mT(self, model_param_dataseries_list,
                          model_param_uncertainty_dataseries_list):
        params = (model_param_dataseries_list[4]*GFACTORCONSTANT2*0.2)**(-1)
        params_sigma = params*model_param_uncertainty_dataseries_list[4]/(
                            model_param_dataseries_list[4])
        return params, params_sigma

    def osc_period_300mT(self, model_param_dataseries_list,
                          model_param_uncertainty_dataseries_list):
        params = (model_param_dataseries_list[4]*GFACTORCONSTANT2*0.3)**(-1)
        params_sigma = params*model_param_uncertainty_dataseries_list[4]/(
                            model_param_dataseries_list[4])
        return params, params_sigma

    def field_offset(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[5]
        params_sigma = model_param_uncertainty_dataseries_list[5]
        return params, params_sigma

    def drift_velocity(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[6]
        params_sigma = model_param_uncertainty_dataseries_list[6]
        return params, params_sigma


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
        self.model_params = ["short_amplitude", "short_lifetime", "long_amplitude",
                             "long_lifetime", "offset"]
        # params = short_amplitude, short_lifetime, long_amplitude, long_lifetime, offset
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
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")
            
        fields = ["short_amplitude",
                  "short_lifetime",
                  "long_amplitude",
                  "long_lifetime"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.short_amplitude(params, param_sigmas),
                  self.short_lifetime(params, param_sigmas),
                  self.long_amplitude(params, param_sigmas),
                  self.long_lifetime(params, param_sigmas),
                  ])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    def short_amplitude(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[0]
        params_sigma = model_param_uncertainty_dataseries_list[0]
        return params, params_sigma

    def short_lifetime(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[1]
        params_sigma = model_param_uncertainty_dataseries_list[1]
        return params, params_sigma

    def long_amplitude(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[2]
        params_sigma = model_param_uncertainty_dataseries_list[2]
        return params, params_sigma

    def long_lifetime(self, model_param_dataseries_list,
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
        self.field_index = 0  # key field, probably lockin2x
        self.fitfunction = fitfcns.fitfcn_two_exp_sin_decay
        self.model_params = ["short_amplitude", "short_lifetime", "long_amplitude",
                             "long_lifetime", "osc_period", "phase", "offset"]
        # params = short_amplitude, short_lifetime, long_amplitude,
        #          long_lifetime, osc_period, phase, offset
        self.free_params = [True, True, True, True,
                            True, True, True, True]
        self.initial_params = [0.05, 50, 0.05, 2000, 800, 0, 0, 0]
        self.param_bounds = [(0, 1), (1, 200),
                             (0, 1), (10, 1e6),
                             (600, 1000), (-np.pi, np.pi),
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
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")

        # TODO: design a model field decorator fcn to make this automatic
        fields = ["short_amplitude",
                  "short_lifetime",
                  "long_amplitude",
                  "long_lifetime"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.short_amplitude(params, param_sigmas),
                  self.short_lifetime(params, param_sigmas),
                  self.long_amplitude(params, param_sigmas),
                  self.long_lifetime(params, param_sigmas),
                  ])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    def short_amplitude(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[0]
        params_sigma = model_param_uncertainty_dataseries_list[0]
        return params, params_sigma

    def short_lifetime(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[1]
        params_sigma = model_param_uncertainty_dataseries_list[1]
        return params, params_sigma

    def long_amplitude(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[2]
        params_sigma = model_param_uncertainty_dataseries_list[2]
        return params, params_sigma

    def long_lifetime(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[3]
        params_sigma = model_param_uncertainty_dataseries_list[3]
        return params, params_sigma


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class IndependentSinusoidalSpinLifetimeModel(ScanDataModel):
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
        self.model_type = "fitfcn_two_indep_exp_sin_decay"
        self.dim2_key = "Voltage"
        self.field_index = 0  # key field, probably lockin2x
        self.fitfunction = fitfcns.fitfcn_two_indep_exp_sin_decay
        self.model_params = ["num_pulses", "short_amplitude", "long_amplitude",
                             "short_lifetime", "long_lifetime", "short_osc_period",
                             "long_osc_period", "short_drift_velocity",
                             "long_drift_velocity", "short_phase", "long_phase",
                             "slope", "offset"]
        # params = num_pulses, short_amplitude, long_amplitude, short_lifetime, long_lifetime,
        #          short_osc_period, long_osc_period, short_drift_velocity, long_drift_velocity,
        #          short_phase, long_phase, slope, offset
        self.free_params = [False, True, True, True, True,
                            True, True, False, False,
                            True, True, True, True]
        self.initial_params = [40, 0.05, 0.05, 50, 2000,
                               800, 800, 0, 0,
                               0, 0, 0, 0]
        self.param_bounds = [(1,1000), (0, 1), (0, 1),
                             (1, 1e9), (1, 1e9),
                             (600, 1000), (600, 1000),
                             (-1e3, 1e3), (-1e3, 1e3),
                             (-np.pi, np.pi), (-np.pi, np.pi),
                             (-1e-6, 1e-6), (-0.01, 0.01)]
        self.error_thresholds = [None, None, None, None, None,
                                 None, None, None, None,
                                 None, None, None, None]
        self.max_fcn_evals = 10000
        self.excluded_intervals = None
        
        # unique to this model!
        self.BField = 0  # in mT

        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def all_model_fields(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list
        param_sigmas = model_param_uncertainty_dataseries_list
        if len(params) == 0 or len(params) != len(param_sigmas):
            raise ValueError("ScanDataModel: tried to extract " +
                             "fit results without successful fit")
            
        # TODO: design a model field decorator fcn to make this automatic
        fields = ["short_amplitude",
                  "long_amplitude",
                  "short_lifetime",
                  "long_lifetime",
                  "short_osc_period",
                  "long_osc_period",
                  "short_gfactor",
                  "long_gfactor",
                  "short_drift_velocity",
                  "long_drift_velocity",
                  "short_phase",
                  "long_phase",
                  "beat_osc_period"]
        dataseries_list, uncertainty_dataseries_list = \
            zip(*[self.short_amplitude(params, param_sigmas),
                  self.long_amplitude(params, param_sigmas),
                  self.short_lifetime(params, param_sigmas),
                  self.long_lifetime(params, param_sigmas),
                  self.short_osc_period(params, param_sigmas),
                  self.long_osc_period(params, param_sigmas),
                  self.short_gfactor(params, param_sigmas),
                  self.long_gfactor(params, param_sigmas),
                  self.short_drift_velocity(params, param_sigmas),
                  self.long_drift_velocity(params, param_sigmas),
                  self.short_phase(params, param_sigmas),
                  self.long_phase(params, param_sigmas),
                  self.beat_osc_period(params, param_sigmas),
                  ])
        dataseries_list = list(dataseries_list)
        uncertainty_dataseries_list = list(uncertainty_dataseries_list)
        return fields, dataseries_list, uncertainty_dataseries_list

    def short_amplitude(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[1]
        params_sigma = model_param_uncertainty_dataseries_list[1]
        return params, params_sigma

    def long_amplitude(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[2]
        params_sigma = model_param_uncertainty_dataseries_list[2]
        return params, params_sigma

    def short_lifetime(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[3]
        params_sigma = model_param_uncertainty_dataseries_list[3]
        return params, params_sigma

    def long_lifetime(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[4]
        params_sigma = model_param_uncertainty_dataseries_list[4]
        return params, params_sigma

    def short_osc_period(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[5]
        params_sigma = model_param_uncertainty_dataseries_list[5]
        return params, params_sigma

    def long_osc_period(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[6]
        params_sigma = model_param_uncertainty_dataseries_list[6]
        return params, params_sigma

    def short_gfactor(self, model_param_dataseries_list,
                             model_param_uncertainty_dataseries_list):
        params = (model_param_dataseries_list[5] * \
                  GFACTORCONSTANT2*self.BField)**(-1)
        params_sigma = params*model_param_uncertainty_dataseries_list[5]/(
                            model_param_dataseries_list[5])
        return params, params_sigma

    def long_gfactor(self, model_param_dataseries_list,
                            model_param_uncertainty_dataseries_list):
        params = (model_param_dataseries_list[6] * \
                  GFACTORCONSTANT2*self.BField)**(-1)
        params_sigma = params*model_param_uncertainty_dataseries_list[6]/(
                            model_param_dataseries_list[6])
        return params, params_sigma

    def short_drift_velocity(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[7]
        params_sigma = model_param_uncertainty_dataseries_list[7]
        return params, params_sigma

    def long_drift_velocity(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[8]
        params_sigma = model_param_uncertainty_dataseries_list[8]
        return params, params_sigma

    def short_phase(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[9]
        params_sigma = model_param_uncertainty_dataseries_list[9]
        return params, params_sigma

    def long_phase(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        params = model_param_dataseries_list[10]
        params_sigma = model_param_uncertainty_dataseries_list[10]
        return params, params_sigma

    def beat_osc_period(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
#        with np.errstate(all='ignore'):
#            params = np.array(pow(np.abs(
#                                pow(model_param_dataseries_list[5],  -1.) -
#                                pow(model_param_dataseries_list[6], -1.)), -1))
#        params_sigma_yvals = np.max(np.vstack([
#                                  model_param_uncertainty_dataseries_list[6].yvals(unsorted=False),
#                                  model_param_uncertainty_dataseries_list[5].yvals(unsorted=False)]),
#                              axis=0)
#        params_sigma = model_param_uncertainty_dataseries_list[5]*0 + params_sigma_yvals
        params = abs(model_param_dataseries_list[6] -
                     model_param_dataseries_list[5])
        params_sigma = model_param_dataseries_list[5]
        return params, params_sigma
