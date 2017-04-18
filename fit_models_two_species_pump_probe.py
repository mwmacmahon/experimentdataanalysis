# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 12:48:01 2017

@author: Michael
"""

from collections import OrderedDict

import numpy as np

from experimentdataanalysis.analysis.scandataprocessing \
    import ScanDataModel


# GLOBAL CONSTANTS
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period


# %% FIT FUNCTIONS
def fitfcn_two_species_twogfactortest_TRKR(t,
                                          config_dict,
                                          num_pulses,
                                          pulse_amplitude, species_amp_ratio,
                                          lifetime1, lifetime2,
                                          osc_period1, osc_period2,
                                          time_offset, phase1, phase2,
                                          drift_velocity1, drift_velocity2,
                                          pump_probe_spot_width,
                                          pump_probe_dist,
                                          slope, offset):
    """
    Expected units:
    (times): ps
    (positions): um
    efield: V/cm
    bfield: mT
    wavelength: nm
    temperature: K
    drift_velocity: um/ps
    slope, offset: varies
    
    """
    pulse_amplitude1 = pulse_amplitude
    pulse_amplitude2 = pulse_amplitude * species_amp_ratio
    osc_ang_freq1 = 2 * np.pi / osc_period1
    osc_ang_freq2 = 2 * np.pi / osc_period2
#    phase1 = phase1 + phase_offset
#    phase2 = phase2 + phase_offset
    sigma_pump = pump_probe_spot_width
    sigma_probe = pump_probe_spot_width
    tvals = np.array(t, copy=False)
    otvals = tvals + time_offset
    null_result = np.zeros(tvals.shape)  # vector/scalar to match input shape

    def single_pulse_fcn(t_pulse):
        xdiffs1 = pump_probe_dist - drift_velocity1 * t_pulse
        xdiffs2 = pump_probe_dist - drift_velocity2 * t_pulse
        drift_signal_factors1 = np.exp(-xdiffs1**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        drift_signal_factors2 = np.exp(-xdiffs2**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        effective_amplitudes1 = pulse_amplitude1 * drift_signal_factors1
        effective_amplitudes2 = pulse_amplitude2 * drift_signal_factors2
        if np.all(drift_signal_factors1 > 1e-6) and \
                np.all(drift_signal_factors2 > 1e-6):  # avoid if unnecessary
            return (
                effective_amplitudes1 * np.exp(-t_pulse / lifetime1)
                    * np.cos(osc_ang_freq1 * t_pulse + phase1) +
                effective_amplitudes2 * np.exp(-t_pulse / lifetime2)
                    * np.cos(osc_ang_freq2 * t_pulse + phase2))
        else:
            return null_result  # such low amplitude it may as well be zero

    pulsesum = sum([single_pulse_fcn(otvals + pulsenum * LASER_REPRATE)
                    for pulsenum in range(int(num_pulses))], null_result)

    # needs to be carefully written to avoid numpy weirdness if tvals scalar: 
    linear_offset = np.array(offset + slope * otvals)
    linear_offset[otvals > 1e4] -= slope * LASER_REPRATE  # in case t > ~13000

    result = linear_offset + pulsesum
    return result


# %% FIT FUNCTIONS
def fitfcn_two_species_exp_sin_decay_TRKR(t,
                                          config_dict,
                                          num_pulses,
                                          pulse_amplitude, species_amp_ratio,
                                          lifetime1, lifetime2, osc_period,
                                          time_offset, phase1, phase2,
                                          drift_velocity1, drift_velocity2,
                                          pump_probe_spot_width,
                                          pump_probe_dist,
                                          slope, offset):
    """
    Expected units:
    (times): ps
    (positions): um
    efield: V/cm
    bfield: mT
    wavelength: nm
    temperature: K
    drift_velocity: um/ps
    slope, offset: varies
    
    """
    pulse_amplitude1 = pulse_amplitude
    pulse_amplitude2 = pulse_amplitude * species_amp_ratio
    osc_ang_freq = 2 * np.pi / osc_period
#    phase1 = phase1 + phase_offset
#    phase2 = phase2 + phase_offset
    sigma_pump = pump_probe_spot_width
    sigma_probe = pump_probe_spot_width
    tvals = np.array(t, copy=False)
    otvals = tvals + time_offset
    null_result = np.zeros(tvals.shape)  # vector/scalar to match input shape

    def single_pulse_fcn(t_pulse):
        xdiffs1 = pump_probe_dist - drift_velocity1 * t_pulse
        xdiffs2 = pump_probe_dist - drift_velocity2 * t_pulse
        drift_signal_factors1 = np.exp(-xdiffs1**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        drift_signal_factors2 = np.exp(-xdiffs2**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        effective_amplitudes1 = pulse_amplitude1 * drift_signal_factors1
        effective_amplitudes2 = pulse_amplitude2 * drift_signal_factors2
        if np.all(drift_signal_factors1 > 1e-6) and \
                np.all(drift_signal_factors2 > 1e-6):  # avoid if unnecessary
            return (
                effective_amplitudes1 * np.exp(-t_pulse / lifetime1)
                    * np.cos(osc_ang_freq * t_pulse + phase1) +
                effective_amplitudes2 * np.exp(-t_pulse / lifetime2)
                    * np.cos(osc_ang_freq * t_pulse + phase2))
        else:
            return null_result  # such low amplitude it may as well be zero

    pulsesum = sum([single_pulse_fcn(otvals + pulsenum * LASER_REPRATE)
                    for pulsenum in range(int(num_pulses))], null_result)

    # needs to be carefully written to avoid numpy weirdness if tvals scalar: 
    linear_offset = np.array(offset + slope * otvals)
    linear_offset[otvals > 1e4] -= slope * LASER_REPRATE  # in case t > ~13000

    result = linear_offset + pulsesum
    return result


def fitfcn_two_species_exp_sin_decay_RSA(bfield,
                                         config_dict,
                                         num_pulses, delay_time,
                                         pulse_amplitude, species_amp_ratio,
                                         lifetime1, lifetime2, gfactor,
                                         phase_mean, phase_diff,
                                         drift_velocity1, drift_velocity2,
                                         pump_probe_spot_width,
                                         pump_probe_dist,
                                         slope, offset):
    """
    Expected units:
    (times): ps
    (positions): um
    efield: V/cm
    bfield: mT
    wavelength: nm
    temperature: K
    drift_velocity: um/ps
    slope, offset: varies
    
    """
    pulse_amplitude1 = pulse_amplitude
    pulse_amplitude2 = pulse_amplitude * species_amp_ratio
    phase1 = phase_mean - 0.5 * phase_diff
    phase2 = phase_mean + 0.5 * phase_diff
    sigma_pump = pump_probe_spot_width
    sigma_probe = pump_probe_spot_width
    bvals = np.array(bfield, copy=False)
    null_result = np.zeros(bvals.shape)  # vector/scalar to match input shape

    def single_pulse_fcn(t_pulse):
        xdiffs1 = pump_probe_dist - drift_velocity1 * t_pulse
        xdiffs2 = pump_probe_dist - drift_velocity2 * t_pulse
        osc_ang_freqs = 2 * np.pi * GFACTORCONSTANT * gfactor * bvals
        drift_signal_factors1 = np.exp(-xdiffs1**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        drift_signal_factors2 = np.exp(-xdiffs2**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        effective_amplitudes1 = pulse_amplitude1 * drift_signal_factors1
        effective_amplitudes2 = pulse_amplitude2 * drift_signal_factors2
        if np.all(drift_signal_factors1 > 1e-6) and \
                np.all(drift_signal_factors2 > 1e-6):  # avoid if unnecessary
            return (
                effective_amplitudes1 * np.exp(-t_pulse / lifetime1)
                    * np.cos(osc_ang_freqs * t_pulse + phase1) +
                effective_amplitudes2 * np.exp(-t_pulse / lifetime2)
                    * np.cos(osc_ang_freqs * t_pulse + phase2))
        else:
            return null_result  # such low amplitude it may as well be zero

    pulsesum = sum([single_pulse_fcn(delay_time + pulsenum * LASER_REPRATE)
                    for pulsenum in range(int(num_pulses))], null_result)

    linear_offset = np.array(offset + slope * bvals)

    result = linear_offset + pulsesum
    return result


def fitfcn_featurevector_two_species_exp_sin_decay(
                                        feature_vector,  # X-COORD
                                        config_dict,
                                        num_pulses,
                                        pulse_amplitude, species_amp_ratio,
                                        lifetime1, lifetime2,
                                        gfactor_mean, gfactor_diff,
                                        phase_mean, phase_diff,
                                        absorption_vs_efield_coeff1,
                                        absorption_vs_efield_coeff2,
                                        lifetime_vs_efield_coeff1,
                                        lifetime_vs_efield_coeff2,
                                        gfactor_vs_efield_coeff,
                                        mobility, pump_probe_dist):
    """
    Important note on function parameters:
    When calling any of the curve fitting procedures in the scandataprocessing
    module, parameters NOT marked free are NOT sent to scipy's curve_fit.
    Rather, a partial function with those parameters 'baked-in' is sent.
    This means non-numerical parameters are allowed as long as they are
    _never_ set as free parameters.
    
    Expected feature vector:
    ([unused], delaytime, efield, bfield,
     pump_probe_dist, wavelength, temperature, set_temperature,
     runID, index_in_run)

    Expected units:
    (times): ps
    (positions): um
    efield: V/cm
    bfield: mT
    wavelength: nm
    temperature: K
    drift_velocity: um/ps
    slope, offset: varies
    """
    pulse_amplitude1 = pulse_amplitude
    pulse_amplitude2 = pulse_amplitude * species_amp_ratio
    gfactor1 = gfactor_mean - 0.5 * gfactor_diff
    gfactor2 = gfactor_mean + 0.5 * gfactor_diff
    phase1 = phase_mean - 0.5 * phase_diff
    phase2 = phase_mean + 0.5 * phase_diff

    feature_vector = np.array(feature_vector, copy=False)
    if len(feature_vector.shape) > 1:
        t = feature_vector[1, :]
        efield = feature_vector[2, :]
        bfield = feature_vector[3, :]
        pump_probe_dist = feature_vector[4, :]
        wavelength = feature_vector[5, :]
        temp = feature_vector[6, :]
        set_temp = feature_vector[7, :]
        runID = feature_vector[8, :]
        index_in_run = feature_vector[9, :]
        vector_length = np.max([len(feature_vector[i, :])
                                   for i in range(1, feature_vector.shape[0])])
        null_result = np.zeros(vector_length)
    else:
        _, t, efield, bfield, \
            pump_probe_dist, wavelength, \
            temp, set_temp, runID, index_in_run = feature_vector
        null_result = 0

#    adjusted_amplitude_1 = pulse_amplitude1 * \
#                           (estimated_temp / set_temp)**absorption_vs_temp_exponent1
    adjusted_amplitude_1 = pulse_amplitude1 * \
                           (1.0 + absorption_vs_efield_coeff1 * np.abs(efield))
    adjusted_amplitude_2 = pulse_amplitude2 * \
                           (1.0 + absorption_vs_efield_coeff2 * np.abs(efield))
    adjusted_lifetime1 = lifetime1 * \
                           (1.0 + lifetime_vs_efield_coeff1 * np.abs(efield))
    adjusted_lifetime2 = lifetime2 * \
                           (1.0 + lifetime_vs_efield_coeff2 * np.abs(efield))
    adjusted_gfactor1 = gfactor1 * \
                        (1.0 + gfactor_vs_efield_coeff * np.abs(efield)**1)
    adjusted_gfactor2 = gfactor2 * \
                        (1.0 + gfactor_vs_efield_coeff * np.abs(efield)**1)
    drift_velocity1 = mobility * efield
    drift_velocity2 = mobility * efield

    # TODO: check bounds, e.g. no negative lifetimes
    # TODO: figure out low-lifetime-overflow problem (exp(-inf))
    osc_ang_freq1 = 2 * np.pi * GFACTORCONSTANT * adjusted_gfactor1 * bfield
    osc_ang_freq2 = 2 * np.pi * GFACTORCONSTANT * adjusted_gfactor2 * bfield
    sigma_probe = 17.5  # probe beam waist in um, +- 0.5um, from Marta's paper
    sigma_pump = sigma_probe  # assuming no diffusion


    def single_pulse_fcn(t_pulse):
        xdiff1 = pump_probe_dist - drift_velocity1 * t_pulse
        xdiff2 = pump_probe_dist - drift_velocity2 * t_pulse
        drift_signal_factor1 = np.exp(-xdiff1**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        drift_signal_factor2 = np.exp(-xdiff2**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        effective_amplitude1 = adjusted_amplitude_1 * drift_signal_factor1
        effective_amplitude2 = adjusted_amplitude_2 * drift_signal_factor2
#       needs to be replaced by an element-by-element filter now that
#       feature vectors are involved
#        if np.any(drift_signal_factor1 > 1e-6) or \
#                np.any(drift_signal_factor2 > 1e-6):  # avoid if unnecessary
        return (
            effective_amplitude1 * np.exp(-t_pulse / adjusted_lifetime1)
                * np.cos(osc_ang_freq1 * t_pulse + phase1) +
            effective_amplitude2 * np.exp(-t_pulse / adjusted_lifetime2)
                * np.cos(osc_ang_freq2 * t_pulse + phase2))
#        else:
#            return null_result

    pulsesum = sum([single_pulse_fcn(t + pulsenum * LASER_REPRATE)
                    for pulsenum in range(int(num_pulses))], null_result)

    # HANDLING OFFSETS
    # TODO: in future versions, should make a factory function that adds
    # an arbitrary number of FIXED parameters for each RunID's slope, offset,
    # and any other run-specific stuff (curve_fit can't handle vector
    # parameters, I think). It can return a function w/ fixed number of 
    # scalar-parameters that is a wrapper to this function, where those
    # parameters are converted to the slope and offset vectors. then we can
    # get the slope for all data pts by 'slope_vector[runID_vector]' (& offset)

#    wrapped_t = t.copy()  # deep copy to avoid changing original, may be slow
#    wrapped_t[t > 1e4] -= LASER_REPRATE
#    linear_offset = offset + slope*wrapped_t

#    linear_offset = offset + slope * t
#    if len(feature_vector.shape) > 1:
#        linear_offset[t > 1e4] -= slope * LASER_REPRATE
#    elif t > 1e4:
#        linear_offset -= slope * LASER_REPRATE

#    print((linear_offset + pulsesum).shape)

#    result = linear_offset + pulsesum
    result = pulsesum
    return result


# %% FIT FUNCTION MODELS
class TwoSpeciesTwoGFactorsTRKRModel(ScanDataModel):
    def __init__(self, **kwargs):
        self.model_name = "Two Species Two GFactors Pump Probe TRKR"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = fitfcn_two_species_twogfactortest_TRKR
        self.model_params = \
            OrderedDict([('config_dict',        {'free': False,  # untouchable
                                                 'initial value': {}}),
                         ('num_pulses',         {'free': False,
                                                 'initial value': 20}),
                         ('pulse_amplitude',    {'free': True,
                                                 'initial value': 0.03,
                                                 'bounds': (0.001, .1)}),
                         ('species_amp_ratio',  {'free': True,
                                                 'initial value': 2.0,
                                                 'bounds': (0, 10.0)}),
                         ('lifetime1',          {'free': False,
                                                 'initial value': 20000,
                                                 'bounds': (15000, np.inf)}),
                         ('lifetime2',          {'free': True,
                                                 'initial value': 3000,
                                                 'bounds': (500, 6000)}),
                         ('osc_period1',        {'free': True,
                                                 'initial value': 560}),
                         ('osc_period2',        {'free': True,
                                                 'initial value': 550}),
                         ('time_offset',        {'free': False,
                                                 'initial value': 0,
                                                 'bounds': (-60, 60)}),
                         ('phase1',             {'free': True,
                                                 'initial value': np.pi}),
                         ('phase2',             {'free': True,
                                                 'initial value': 0}),
                         ('drift_velocity1',    {'free': False,
                                                 'initial value': 0}),
                         ('drift_velocity2',    {'free': True,
                                                 'initial value': 0}),
                         ('pump_probe_spot_width',
                                                {'free': False,
                                                 'initial value': 17.5}),
                         ('pump_probe_dist',    {'free': False,
                                                 'initial value': 0}),
                         ('slope',              {'free': True,
                                                 'initial value': 0}),
                         ('offset',             {'free': True,
                                                 'initial value': 0}),
                        ])
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        self.ignore_weights = False

        # for ScanDataSets: coord spanning the set
        self.fit_result_scan_coord = "Wavelength (nm)" 

        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val


# %% FIT FUNCTION MODELS
class TwoSpeciesPumpProbeTRKRModel(ScanDataModel):
    def __init__(self, **kwargs):
        self.model_name = "Two Species Pump Probe TRKR"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = fitfcn_two_species_exp_sin_decay_TRKR
        self.model_params = \
            OrderedDict([('config_dict',        {'free': False,  # untouchable
                                                 'initial value': {}}),
                         ('num_pulses',         {'free': False,
                                                 'initial value': 20}),
                         ('pulse_amplitude',    {'free': True,
                                                 'initial value': 0.04,
                                                 'bounds': (0, .1)}),
                         ('species_amp_ratio',  {'free': True,
                                                 'initial value': 3.0,
                                                 'bounds': (0, np.inf)}),
                         ('lifetime1',          {'free': False,
                                                 'initial value': 20000,
                                                 'bounds': (15000, np.inf)}),
                         ('lifetime2',          {'free': True,
                                                 'initial value': 2000,
                                                 'bounds': (800, np.inf)}),
                         ('osc_period',         {'free': True,
                                                 'initial value': 800,
                                                 'bounds': (400, 1200)}),
                         ('time_offset',        {'free': False,
                                                 'initial value': 0,
                                                 'bounds': (-60, 60)}),
                         ('phase1',             {'free': False,
                                                 'initial value': np.pi}),
                         ('phase2',             {'free': False,
                                                 'initial value': 0}),
                         ('drift_velocity1',    {'free': False,
                                                 'initial value': 0}),
                         ('drift_velocity2',    {'free': False,
                                                 'initial value': 0}),
                         ('pump_probe_spot_width',
                                                {'free': False,
                                                 'initial value': 17.5}),
                         ('pump_probe_dist',    {'free': False,
                                                 'initial value': 0}),
                         ('slope',              {'free': True,
                                                 'initial value': 0}),
                         ('offset',             {'free': True,
                                                 'initial value': 0}),
                        ])
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        self.ignore_weights = False

        # for ScanDataSets: coord spanning the set
        self.fit_result_scan_coord = "Wavelength (nm)" 

        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val


# %% FIT FUNCTION MODELS
class TwoSpeciesPumpProbeRSAModel(ScanDataModel):
    def __init__(self, **kwargs):
        self.model_name = "Two Species Pump Probe TRKR"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = fitfcn_two_species_exp_sin_decay_TRKR
        self.model_params = \
            OrderedDict([('config_dict',        {'free': False,  # untouchable
                                                 'initial value': {}}),
                         ('num_pulses',         {'free': False,
                                                 'initial value': 20}),
                         ('delay_time',         {'free': False,
                                                 'initial value': -160}),
                         ('pulse_amplitude',    {'free': True,
                                                 'initial value': 0.025,
                                                 'bounds': (0, .05)}),
                         ('species_amp_ratio',  {'free': True,
                                                 'initial value': 0.2,
                                                 'bounds': (0, np.inf)}),
                         ('lifetime1',          {'free': True,
                                                 'initial value': 50000,
                                                 'bounds': (0, np.inf)}),
                         ('lifetime2',          {'free': True,
                                                 'initial value': 10000,
                                                 'bounds': (0, np.inf)}),
                         ('gfactor',            {'free': True,
                                                 'initial value': 0.44}),
                         ('phase_mean',         {'free': True,
                                                 'initial value': np.pi/2,
                                                 'bounds':
                                                     (np.pi/4, 3*np.pi/4)}),
                         ('phase_diff',         {'free': False,
                                                 'initial value': -np.pi}),
                         ('drift_velocity1',    {'free': False,
                                                 'initial value': 0}),
                         ('drift_velocity2',    {'free': False,
                                                 'initial value': 0}),
                         ('pump_probe_spot_width',
                                                {'free': False,
                                                 'initial value': 17.5}),
                         ('pump_probe_dist',    {'free': False,
                                                 'initial value': 0}),
                         ('slope',              {'free': True,
                                                 'initial value': 0}),
                         ('offset',             {'free': True,
                                                 'initial value': 0}),
                        ])
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        self.ignore_weights = False

        # for ScanDataSets: coord spanning the set
        self.fit_result_scan_coord = "Wavelength (nm)" 

        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val


class FeatureVectorTwoSpeciesPumpProbeModel(ScanDataModel):
    def __init__(self, **kwargs):
        self.model_name = "Two Species Pump Probe using Feature Vectors"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = \
            fitfcn_featurevector_two_species_exp_sin_decay
        default_configuration = {'absorption_vs_efield': 'linear',
                                 'lifetime_vs_efield': 'linear',
                                 'gfactor_vs_efield': 'linear'}
        self.model_params = \
            OrderedDict([('config_dict',       {'free': False,  # untouchable
                                                 'initial value':
                                                     default_configuration}),
                         ('num_pulses',         {'free': False,
                                                 'initial value': 20}),
                         ('pulse_amplitude',    {'free': True,
                                                 'initial value': 0.025,
                                                 'bounds': (0, np.inf)}),
                         ('species_amp_ratio',  {'free': False,
                                                 'initial value': 2.0,
                                                 'bounds': (0, np.inf)}),
                         ('lifetime1',          {'free': False,
                                                 'initial value': 20000,
                                                 'bounds': (0, np.inf)}),
                         ('lifetime2',          {'free': False,
                                                 'initial value': 6000,
                                                 'bounds': (0, np.inf)}),
                         ('gfactor_mean',       {'free': False,
                                                 'initial value': 0.4385}),
                         ('gfactor_diff',       {'free': False,
                                                 'initial value': 0.0}),
                         ('phase_mean',         {'free': False,
                                                 'initial value': np.pi/2}),
                         ('phase_diff',         {'free': False,
                                                 'initial value': -np.pi}),
                         ('absorption_vs_efield_coeff1',
                                                {'free': False,
                                                 'initial value': -0e-2}),
                         ('absorption_vs_efield_coeff2',
                                                {'free': False,
                                                 'initial value': 1e-9}),
                         ('lifetime_vs_efield_coeff1',
                                                {'free': False,
                                                 'initial value': 0e-3}),
                         ('lifetime_vs_efield_coeff2',
                                                {'free': False,
                                                 'initial value': 0e-3}),
                         ('gfactor_vs_efield_coeff', 
                                                {'free': False,
                                                 'initial value': -4.0e-4}),
                         ('mobility',           {'free': False,
                                                 'initial value': 1e-4}),
                         ('pump-probe dist',    {'free': False,
                                                 'initial value': 0}),
                        ])
        self.max_fcn_evals = 20000
        self.excluded_intervals = None
        self.ignore_weights = False

        # for ScanDataSets: coord spanning the set
        self.fit_result_scan_coord = "Pump-Probe Distance (um)" 

        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

# 818.0 data:
#initial_params=[20, 0.009, 0.0, 15000, 1.0,
#                0.435, 0.0, -np.pi, 0,
#                0.5, -1.2, -1.2,
#                1.0, 3.2e-4,
#                1e-4, 0, 0],
#param_bounds=[(20,1000), (0, 1), (1.5, 5.0), (1, 1e9), (1, 1e9),
#              (0.3, 0.6), (-0.2, 0.2),
#              (-4*np.pi, 4*np.pi), (-4*np.pi, 4*np.pi),
#              (0.001, 10), (-10.0, 10.0), (-10.0, 10.0),
#              (-10, 10), (0, 1000),
#              (-1e-2, 1e-2), (-1e-6, 1e-6), (-0.01, 0.01)],

# 818.9 data:
#initial_params=[20, 0.025, 2.0, 20260, 8000,
#                0.4385, 0.0, np.pi/2, -np.pi,
#                0.03, -0.8, -1.2,
#                -1.2, 3.5e-4,
#                1e-4, 0, 0],
#param_bounds=[(1,1000), (0, 1), (1.5, 5.0), (1, 1e9), (1e3, 2e4),
#              (0.3, 0.6), (-0.2, 0.2),
#              (-4*np.pi, 4*np.pi), (-4*np.pi, 4*np.pi),
#              (0.01, 0.06), (-10.0, 10.0), (-10.0, 10.0),
#              (-1.8, -1.2), (0, 1e-3),
#              (-1e-2, 1e-2), (-1e-6, 1e-6), (-0.01, 0.01)],

def get_two_species_TRKR_model(**model_kwargs):
    # mostly passthrough? can put options as parameters to override defaults
    model = TwoSpeciesPumpProbeTRKRModel(**model_kwargs)
    return model

def get_two_species_RSA_model(**model_kwargs):
    # mostly passthrough? can put options as parameters to override defaults
    model = TwoSpeciesPumpProbeRSAModel(**model_kwargs)
    return model

def get_one_species_TRKR_model(**model_kwargs):
    # fix a few parameters to make it effectively one species
    model = TwoSpeciesPumpProbeTRKRModel(**model_kwargs)
    model.model_params['species_amp_ratio']['free'] = False
    model.model_params['species_amp_ratio']['initial value'] = 0.0
    model.model_params['lifetime2']['free'] = False
#    model.model_params['phase_mean']['free'] = True
#    model.model_params['phase_mean']['initial_value'] = 0.0
#    model.model_params['phase_mean']['bounds'] = (-2*np.pi, 2*np.pi)
#    model.model_params['phase_diff']['free'] = False
#    model.model_params['phase_diff']['initial_value'] = 0.0
    return model


def get_fvec_two_species_model(**model_kwargs):
    # mostly passthrough? can put options as parameters
    if model_kwargs.get('excluded_intervals') is not None:
        raise ValueError("Warning: feature vector models need to use an " +
                         "alternative implementation of excluded_intervals." +
                         "Use the excluded_intervals parameter in the " +
                         "scandata_list_to_feature_vector_scandata function " +
                         "in the analysis.featurevectors package.")
    model = FeatureVectorTwoSpeciesPumpProbeModel(
                excluded_intervals=None,  # use feature vectors' alternative
                **model_kwargs)
    return model


def get_fvec_one_species_model(**model_kwargs):
    # fix a few parameters to make it effectively one species
    if model_kwargs.get('excluded_intervals') is not None:
        raise ValueError("Warning: feature vector models need to use an " +
                         "alternative implementation of excluded_intervals." +
                         "Use the excluded_intervals parameter in the " +
                         "scandata_list_to_feature_vector_scandata function " +
                         "in the analysis.featurevectors package.")
    model = FeatureVectorTwoSpeciesPumpProbeModel(
                excluded_intervals=None,  # use feature vectors' alternative
                **model_kwargs)
    model.model_params['species_amp_ratio']['free'] = False
    model.model_params['species_amp_ratio']['initial value'] = 0.0
    model.model_params['lifetime1']['initial value'] = 20000
    model.model_params['lifetime2']['initial value'] = 1000
    model.model_params['gfactor_mean']['initial value'] = 0.435
    model.model_params['gfactor_diff']['initial value'] = 0
#    model.model_params['phase_mean']['initial value'] = -np.pi
#    model.model_params['phase_diff']['initial value'] = 0
    return model

#initial_params=[20, 0.009, 0.0, 15000, 1.0,
#                0.435, 0.0, -np.pi, 0,
#                0.5, -1.2, -1.2,
#                1.0, 3.2e-4,
#                1e-4, 0, 0],
