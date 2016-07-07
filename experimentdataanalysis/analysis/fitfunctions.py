# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:32:21 2016

Contains simple functions used by other modules to fit data

@author: vsih-lab
"""

import numpy as np


LASER_REPRATE = 13100  # ps period


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
