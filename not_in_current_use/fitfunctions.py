# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:32:21 2016

Contains simple functions used by other modules to fit data

Note on mobility dependence of signal:
Assuming gaussian spin packet with center x1 and width s1,
         gaussian probe beam with center x2 and width s2:
Signal expected to follow integral of product of gaussians.
This result can be expressed in the form of a gaussian itself,
varying from 1 at xdiff=x1-x2=0, to 0 as xdiff->infinity. However,
effective sigma^2 = sigma1^2 + sigma2^2. Thus formula used:

drift_signal_factor = np.exp(-xdiff**2/(2*(sigma1**2 + sigma2**2)))

@author: vsih-lab
"""

import numpy as np


#GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_simple_line(t, slope, offset):
    """
    """
    return slope*t + offset


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_1d_gaussian_with_linear_offset(t, amplitude, t0,
                                          sigma, slope, offset):
    """
    """
    return amplitude*np.exp(-0.5*np.power((t - t0)/(sigma), 2.)) + \
                                                            slope*t + offset


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_simple_1d_gaussian(t, amplitude, t0, sigma, offset):
    """
    """
    return amplitude*np.exp(-0.5*np.power((t - t0)/(sigma), 2.)) + offset


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_exp_times_sqrt_decay(t, 
                                pulse_amplitude,
                                lifetime,
                                offset,
                                D,  # diffusion constant
                                D_toffset):  # diffusion_time_offset):
    """
    temp note: this is useful as it is how the gaussian peak
    of a spin distribution should change with time!

    spin propto 1/sqrt(t)
    
    ...or is it ^2/2 in 2d?    
    
    problem: constraint of diffusion time must be small enough
    that (~10000 ps) << D_toffset since spreading is relatively minimal...
    (actually there is a factor exp(-x^2/4D(t+t_off)), but since
    sqrt(4D(t+t_off)) is gaussian spatial width, it must match actual spread
    of packet with time...we can actually fit this with overlap if we need.
    but essentially we need D_toffset to be at least ~20000 or so)

    probably best just to not plot gaussian peak but rather area under
    gaussian peak (integral of exp(-ax^2) is sqrt(pi/a), remember QFT? :D )

    """
    def single_pulse_fcn(t_pulse):
        return pulse_amplitude*np.exp(-t_pulse/lifetime)*\
                    np.power(4*np.pi*D*(t_pulse + D_toffset), -2/2)

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE

    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return offset + thispulse + lastpulse + lastpulse2


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_single_exp_decay(t, pulse_amplitude, lifetime, offset):
    """
    """
    def single_pulse_fcn(t_pulse):
        return pulse_amplitude*np.exp(-t_pulse/lifetime)

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE

    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return offset + thispulse + lastpulse + lastpulse2


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_two_exp_decay(t, pulse_amplitude, lifetime,
                            pulse_amplitude2, lifetime2, offset):
    """
    """
    def single_pulse_fcn(t_pulse):
        return (pulse_amplitude*np.exp(-t_pulse/lifetime) +
                    pulse_amplitude2*np.exp(-t_pulse/lifetime2))

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE

    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return offset + thispulse + lastpulse + lastpulse2


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_two_exp_sin_decay(t, pulse_amplitude, lifetime,
                             pulse_amplitude2, lifetime2,
                             osc_period, phase, slope, offset):
    """
    """
    def single_pulse_fcn(t_pulse):
        return (pulse_amplitude*np.exp(-t_pulse/lifetime) +
                pulse_amplitude2*np.exp(-t_pulse/lifetime2)*np.cos(
                                        2*np.pi*t_pulse/osc_period + phase))

    last_t = t + LASER_REPRATE
    last_t_2 = t + 2*LASER_REPRATE

    thispulse = single_pulse_fcn(t)
    lastpulse = single_pulse_fcn(last_t)
    lastpulse2 = single_pulse_fcn(last_t_2)
    return offset + slope*t + thispulse + lastpulse + lastpulse2


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_two_indep_exp_sin_decay(t, num_pulses,
                                   pulse_amplitude1, pulse_amplitude2,
                                   lifetime1, lifetime2,
                                   osc_period1, osc_period2,
                                   drift_velocity1, drift_velocity2,
                                   phase1, phase2, slope, offset):
    """
    Expected units:
    lifetime: ps
    osc_period: ps
    drift_velocity: um/ps
    """
    def single_pulse_fcn(t_pulse):
        sigma1 = 17.5  # probe beam waist in um, +- 0.5um, from Marta's paper
        sigma2 = sigma1  # assuming no diffusion
        xdiff1 = drift_velocity1*t_pulse
        xdiff2 = drift_velocity2*t_pulse
        drift_signal_factor1 = np.exp(-xdiff1**2/(2*(sigma1**2 + sigma2**2)))
        drift_signal_factor2 = np.exp(-xdiff2**2/(2*(sigma1**2 + sigma2**2)))
        effective_amplitude1 = pulse_amplitude1 * drift_signal_factor1
        effective_amplitude2 = pulse_amplitude2 * drift_signal_factor2
        if np.any(drift_signal_factor1 > 1e-6) or \
                np.any(drift_signal_factor2 > 1e-6):  # avoid if unnecessary
            return (effective_amplitude1*np.exp(-t_pulse/lifetime1)*np.cos(
                                        2*np.pi*t_pulse/osc_period1 + phase1) +
                    effective_amplitude2*np.exp(-t_pulse/lifetime2)*np.cos(
                                        2*np.pi*t_pulse/osc_period2 + phase2))
        else:
            return 0*t_pulse

    pulsesum = sum([single_pulse_fcn(t + pulsenum*LASER_REPRATE)
                    for pulsenum in range(num_pulses)])

    wrapped_t = t.copy()  # deep copy to avoid changing original, may be slow
    wrapped_t[t > 1e4] -= LASER_REPRATE
#    wrapped_t = t[:]  # shallow copy, must avoid changing original!
#    wrapped_t = np.hstack([t[t < 1e4], t[t > 1e4] - 13160])  # assumes ordered
    linear_offset = offset + slope*wrapped_t
    return linear_offset + pulsesum


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_two_opposite_exp_sin_decay(t, num_pulses,
                                      pulse_amplitude1, pulse_amplitude2,
                                      lifetime1, lifetime2,
                                      osc_period,  drift_velocity,
                                      probe_pos, slope, offset):
    """
    Expected units:
    lifetime: ps
    osc_period: ps
    drift_velocity: um/ps
    probe_pos: um
    """
    def single_pulse_fcn(t_pulse):
        sigma1 = 17.5  # probe beam waist in um, +- 0.5um, from Marta's paper
        sigma2 = sigma1  # assuming no diffusion
        xdiff1 = probe_pos - drift_velocity*t_pulse
        xdiff2 = probe_pos - drift_velocity*t_pulse
        drift_signal_factor1 = np.exp(-xdiff1**2/(2*(sigma1**2 + sigma2**2)))
        drift_signal_factor2 = np.exp(-xdiff2**2/(2*(sigma1**2 + sigma2**2)))
        effective_amplitude1 = pulse_amplitude1 * drift_signal_factor1
        effective_amplitude2 = pulse_amplitude2 * drift_signal_factor2
        if np.any(drift_signal_factor1 > 1e-6) or \
                np.any(drift_signal_factor2 > 1e-6):  # avoid if unnecessary
            return (effective_amplitude1*np.exp(-t_pulse/lifetime1)*np.cos(
                                        2*np.pi*t_pulse/osc_period) +
                    effective_amplitude2*np.exp(-t_pulse/lifetime2)*np.cos(
                                        2*np.pi*t_pulse/osc_period - np.pi))
        else:
            return 0*t_pulse

    pulsesum = sum([single_pulse_fcn(t + pulsenum*LASER_REPRATE)
                    for pulsenum in range(num_pulses)])

#    wrapped_t = t.copy()  # deep copy to avoid changing original, may be slow
#    wrapped_t[t > 1e4] -= LASER_REPRATE
#    linear_offset = offset + slope*wrapped_t

    linear_offset = offset + slope * t
    linear_offset[t > 1e4] -= slope * LASER_REPRATE  # in case we forced pos t

    return linear_offset + pulsesum


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_rsa_field_scan(field, num_pulses, delay_time,
                          pulse_amplitude, lifetime, gfactor,
                          field_offset, drift_velocity,
                          phase, slope, offset):
    """
    Expected units:
    delay_time, lifetime: ps
    drift_velocity: um/ps
    field_offset: Tesla
    gfactor: (1/ps) / Tesla
    """
    B = field + field_offset
    def single_pulse_fcn(t_pulse):
        sigma1 = 17.5  # probe beam waist in um, +- 0.5um, from Marta's paper
        sigma2 = sigma1  # assuming no diffusion
        xdiff = drift_velocity*t_pulse
        drift_signal_factor = np.exp(-xdiff**2/(2*(sigma1**2 + sigma2**2)))
        if drift_signal_factor > 1e-6:  # avoid unnecessary calculations
            effective_amplitude = pulse_amplitude * drift_signal_factor
            return effective_amplitude*np.exp(-t_pulse/lifetime)*np.cos(
                        2*np.pi*GFACTORCONSTANT*gfactor*B*t_pulse + phase)
        else:
            return 0*B

    pulsesum = sum([single_pulse_fcn(delay_time + pulsenum*LASER_REPRATE)
                    for pulsenum in range(num_pulses)])
    return offset + slope*field + pulsesum


# %% NEEDS SPHINX DOCUMENTATION
def fitfcn_featurevector_two_heated_species_exp_sin_decay(
                                        feature_vector, num_pulses,
                                        pulse_amplitude, species_amp_ratio,
                                        lifetime1, lifetime2,
                                        gfactor_mean, gfactor_diff,
                                        phase_mean, phase_diff,
                                        efield_heating_coeff,
                                        absorption_vs_temp_exponent1,
                                        absorption_vs_temp_exponent2,
                                        lifetime_vs_temp_exponent,
                                        gfactor_temp_sensitivity,
                                        mobility, slope_vector, offset_vector):
    """
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
        vector_length = np.max([len(feature_vector[i, :])
                                   for i in range(1, 8)])
    else:
        _, t, efield, bfield, \
            pump_probe_dist, wavelength, temp = feature_vector
        vector_length = 1
    # lifetime vs temp model:
    # Ben found that spin lifetime ~ T^-2.2, where T is lattice temp.
    # Also Marta referenced a paper showing roughly linear dependence
    # of electron temp vs field...needs to be updated much more carefully.
    set_temp = 30  # fixed for now, ignoring feature vector. nonlinear weirdness.
    estimated_temp = set_temp + efield_heating_coeff * np.abs(efield)
    adjusted_amplitude_1 = pulse_amplitude1 * \
                           (estimated_temp / set_temp)**absorption_vs_temp_exponent1
    adjusted_amplitude_2 = pulse_amplitude2 * \
                           (estimated_temp / set_temp)**absorption_vs_temp_exponent2
    adjusted_lifetime1 = lifetime1 * \
                         (estimated_temp / set_temp)**lifetime_vs_temp_exponent
    adjusted_lifetime2 = lifetime2 * \
                         (estimated_temp / set_temp)**lifetime_vs_temp_exponent
    adjusted_gfactor1 = gfactor1 - \
                        gfactor_temp_sensitivity * (estimated_temp - set_temp)
    adjusted_gfactor2 = gfactor2 - \
                        gfactor_temp_sensitivity * (estimated_temp - set_temp)
    osc_ang_freq1 = 2 * np.pi * GFACTORCONSTANT * adjusted_gfactor1 * bfield
    osc_ang_freq2 = 2 * np.pi * GFACTORCONSTANT * adjusted_gfactor2 * bfield
    sigma_probe = 17.5  # probe beam waist in um, +- 0.5um, from Marta's paper
    sigma_pump = sigma_probe  # assuming no diffusion

    def single_pulse_fcn(t_pulse):
        xdiff1 = pump_probe_dist - mobility * efield * t_pulse
        xdiff2 = pump_probe_dist - mobility * efield * t_pulse
        drift_signal_factor1 = np.exp(-xdiff1**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        drift_signal_factor2 = np.exp(-xdiff2**2 /
                                      (2 * (sigma_pump**2 + sigma_probe**2)))
        effective_amplitude1 = adjusted_amplitude_1 * drift_signal_factor1
        effective_amplitude2 = adjusted_amplitude_2 * drift_signal_factor2
        if np.any(drift_signal_factor1 > 1e-6) or \
                np.any(drift_signal_factor2 > 1e-6):  # avoid if unnecessary
            return (
                effective_amplitude1 * np.exp(-t_pulse / adjusted_lifetime1)
                    * np.cos(osc_ang_freq1 * t_pulse + phase1) +
                effective_amplitude2 * np.exp(-t_pulse / adjusted_lifetime2)
                    * np.cos(osc_ang_freq2 * t_pulse + phase2))
        else:
            if vector_length > 1:
                return np.zeros(vector_length)
            else:
                return 0

    pulsesum = sum([single_pulse_fcn(t + pulsenum * LASER_REPRATE)
                    for pulsenum in range(int(num_pulses))])

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

