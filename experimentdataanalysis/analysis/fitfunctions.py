# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:32:21 2016

Contains simple functions used by other modules to fit data

@author: vsih-lab
"""

import numpy as np


LASER_REPRATE = 13100  # ps period


# %% NEEDS TEST, SPHINX DOCUMENTATION
def fitfcn_simple_1d_gaussian(t, amplitude, t0, sigma, offset):
    """
    """
    return amplitude*np.exp(-np.power((t - t0)/(2.*sigma), 2.)) + offset


# %%
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


# %%
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
