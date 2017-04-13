# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 12:41:18 2017

@author: Michael
"""

from collections import OrderedDict

import numpy as np

from experimentdataanalysis.analysis.scandataprocessing \
    import ScanDataModel, scandata_list_fit


# GLOBAL CONSTANTS
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13158  # ps period


# %%
def get_pulse_sum_vector(spin_lifetime, gfactor, bfield, initial_phase=0):
    """
    Determines effect of summed spins over many pulses and returns the net
    phase and amplitude at zero time delay expected to result. Does not
    take into account any b-field-axis zeeman polarization of spins/nuclei,
    or high-intensity polarization saturation and/or repolarization effects.

    Formulae are from Chris Trowbridge's Ph.D. thesis,
    eqns. 5.12, 5.13, 5.29, and 5.30.

    Does not handle more than one species, so run this function on each
    species individually if possible.

    Expected units:
    polarization: unitless scalar in [0, 1]
    lifetime: ps
    """
    osc_ang_freq = 2 * np.pi * GFACTORCONSTANT * gfactor * bfield
    theta = osc_ang_freq * LASER_REPRATE
    x = LASER_REPRATE / spin_lifetime
    net_polarization = 1. / np.sqrt(1 - 2 * np.exp(-x) * np.cos(theta) +
                                    np.exp(-2 * x))
    net_phase = initial_phase + np.arctan((np.exp(-x) * np.sin(theta)) /
                                          (1 - np.exp(-x) * np.cos(theta)))
    return (net_polarization, net_phase)  # at zero delay, of course


# %%
def fitfcn_trkr_single_species(delay_time,
                               config_dict,
                               pulse_amplitude,
                               gfactor, bfield, spin_lifetime,
                               initial_phase, slope, offset):
    zero_delay_offset = config_dict.get('zero_delay_offset', 0.0)
    trkr_per_unit_polarization = 1.0
    trkr_phase_offset = 0.0
    pos_def_delay = (delay_time + zero_delay_offset) % LASER_REPRATE \
                        - zero_delay_offset
    osc_ang_freq = 2 * np.pi * GFACTORCONSTANT * gfactor * bfield
    net_polarization, net_phase = get_pulse_sum_vector(spin_lifetime,
                                                       gfactor, bfield,
                                                       initial_phase)
    final_phase = (net_phase + pos_def_delay * osc_ang_freq) % (2 * np.pi)
    final_polarization = pulse_amplitude * net_polarization * \
                            np.exp(-pos_def_delay / spin_lifetime)
    signal = trkr_per_unit_polarization * final_polarization * \
                                    np.cos(final_phase + trkr_phase_offset)
    output = signal + delay_time * slope + offset  # NOT pos-definite
    return output


# %%
def fitfcn_rsa_single_species(bfield,
                              config_dict,
                              pulse_amplitude,
                              gfactor, delay_time, spin_lifetime,
                              initial_phase, slope, offset):
    signal = fitfcn_trkr_single_species(delay_time,
                                        config_dict,
                                        pulse_amplitude,
                                        gfactor, bfield, spin_lifetime,
                                        initial_phase, slope=0, offset=0)
    output = signal + bfield * slope + offset
    return output


# %% FIT FUNCTION MODELS
class OneSpecies_NoDNP_TRKR_Model(ScanDataModel):
    def __init__(self, **kwargs):
        self.model_name = "One Species Pump Probe TRKR"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = fitfcn_trkr_single_species
        self.model_params = \
            OrderedDict([('config_dict',        {'free': False,  # untouchable
                                                 'initial value': {}}),
                         ('pulse_amplitude',    {'free': True,
                                                 'initial value': 0.025,
                                                 'bounds': (0, .05)}),
                         ('gfactor',            {'free': True,
                                                 'initial value': 0.44}),
                         ('bfield',             {'free': False,
                                                 'initial value': 0.0}),
                         ('spin_lifetime',      {'free': True,
                                                 'initial value': 20000,
                                                 'bounds': (0, np.inf)}),
                         ('initial_phase',      {'free': False,
                                                 'initial value': 0.0,
                                                 'bounds':
                                                     (-4*np.pi, 4*np.pi)}),
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
class OneSpecies_NoDNP_RSA_Model(ScanDataModel):
    def __init__(self, **kwargs):
        self.model_name = "One Species Pump Probe RSA"
        self.yfield = None  # if not set, scandata.y used as fit 'y' values
        self.fitfunction = fitfcn_rsa_single_species
        self.model_params = \
            OrderedDict([('config_dict',        {'free': False,  # untouchable
                                                 'initial value': {}}),
                         ('pulse_amplitude',    {'free': True,
                                                 'initial value': 0.025,
                                                 'bounds': (0, .05)}),
                         ('gfactor',            {'free': True,
                                                 'initial value': 0.44}),
                         ('delay_time',         {'free': False,
                                                 'initial value': -160.0}),
                         ('spin_lifetime',      {'free': True,
                                                 'initial value': 20000,
                                                 'bounds': (0, np.inf)}),
                         ('initial_phase',      {'free': False,
                                                 'initial value': 0.0,
                                                 'bounds':
                                                     (-4*np.pi, 4*np.pi)}),
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


# %%
# DNP: given (time, field) and (last_time, last_field, last_dnp_field),
#      calculate:
#          - dnp_field resulting from (time, last_field), to account for
#            change since last time point until this point
#          - MAKE SURE IT'S NOT TIME STEP DEPENDENT
#          - NOTE: this is a NON-hysteretical model, in the sense that
#          - only field and d(field) matters
def get_new_dnp_field(elapsed_time, bfield, last_dnp_field,
                      target_dnp_fcn=None, dnp_time_const=300.0):
    if elapsed_time <= 0:
        return last_dnp_field
    if target_dnp_fcn is not None:
        target_dnp_field = target_dnp_fcn(bfield, last_dnp_field)
    else:  # default 
        target_dnp_field = 0  # TEMP
    delta_dnp = target_dnp_field - last_dnp_field
    return target_dnp_field - \
                delta_dnp * np.exp(-elapsed_time / dnp_time_const)


# %%
def extract_dnp_from_field_sweep(bvals, yvals,
                                 ref_field_sweep_period_bvals,
                                 ref_field_sweep_period_yvals,
                                 rsa_period,
                                 y_uncertainty=None,
                                 jump_threshold=0.0,
                                 tracking_threshold_ratio=5.0):
    """
    NOTE: field sweep is assumed periodic, so function will loop around
          and interpolate between provided maximum-b and minimum-b
          (b, y) pairs if they do not span the full rsa_period.
    """
    bvals = 1.0 * np.array(bvals, copy=True)
    yvals = 1.0 * np.array(yvals, copy=True)
    rsa_bvals = 1.0 * np.array(ref_field_sweep_period_bvals, copy=True)
    rsa_yvals = 1.0 * np.array(ref_field_sweep_period_yvals, copy=True)
    # prep: demand monotonic sweep
    b_increasing = (rsa_bvals[1] > rsa_bvals[0])
    last_b = rsa_bvals[0]
    for b in rsa_bvals:
        if b_increasing:
            if b < last_b:
                raise ValueError("ref_field_sweep_period_bvals " +
                                 "must be monotonic!")
        else:
            if b > last_b:
                raise ValueError("ref_field_sweep_period_bvals " +
                                 "must be monotonic!")            
    # prep: if sweep too long, truncate
    if max(rsa_bvals) - min(rsa_bvals) > rsa_period:
        min_bvals_field = max(rsa_bvals) - rsa_period
        allowed_indices = [ind for ind in range(len(rsa_bvals))
                           if rsa_bvals[ind] >= min_bvals_field]
        print("Warning: extract_dnp_from_field_sweep not currently " +
              "equipped to handle field sweep periods greater than " +
              "RSA period, truncating all except last full RSA period." +
              "Kept {} of {} data points.".format(len(allowed_indices),
                                                  len(rsa_bvals)))
        rsa_bvals = rsa_bvals[allowed_indices]
        rsa_yvals = rsa_yvals[allowed_indices]
    # prep: if beginning and end of sweep are exactly rsa_period apart, combine
    if min(rsa_bvals) + rsa_period == max(rsa_bvals):
        min_bvals_index = np.argwhere(rsa_bvals == min(rsa_bvals))[0]
        max_bvals_index = np.argwhere(rsa_bvals == max(rsa_bvals))[-1]
        avg_yval = np.mean([rsa_yvals[min_bvals_index],
                            rsa_yvals[max_bvals_index]])
        rsa_yvals[min_bvals_index] = avg_yval
        new_indices = np.array([ind for ind in range(len(rsa_bvals))
                                if ind != max_bvals_index])
        rsa_bvals = rsa_bvals[new_indices]
        rsa_yvals = rsa_yvals[new_indices]
    rsa_bvals = rsa_bvals % rsa_period
    map_to_sorted = np.argsort(rsa_bvals)
    rsa_bvals = rsa_bvals[map_to_sorted]
    rsa_yvals = rsa_yvals[map_to_sorted]

    def get_rsa_signal(bfield):
        return np.interp(bfield % rsa_period, rsa_bvals, rsa_yvals,
                         period=rsa_period)

    def get_rsa_signal_slope(bfield):
        try:
            len(bfield)
            bfields = bfield
        except TypeError:
            bfields = [bfield]
        results = []
        for bfield in bfields:
            test_spacing = 0.001 * rsa_period
            test_bvals = np.array([bfield - 2 * test_spacing,
                                   bfield - 1 * test_spacing,
                                   bfield,
                                   bfield + 1 * test_spacing,
                                   bfield + 2 * test_spacing])
            test_yvals = get_rsa_signal(test_bvals)
            test_yvals_slope = np.gradient(test_yvals)
            results.append(test_yvals_slope[2] / test_spacing)
        return results

    def find_rsa_match_fields(signal):
        def inverse_interpolate_rsa(match_pair):
            # ensure yvals (the x-coord in interpolation) are _increasing_:
            if rsa_yvals[match_pair[0]] > rsa_yvals[match_pair[1]]:
                match_pair = match_pair[::-1]
            yval_pair = rsa_yvals[match_pair]
            # ensure bvals (the y-coord in interpolation) are periodic:
            bval_pair = rsa_bvals[match_pair]
            if match_pair[0] == max_ind and match_pair[1] == 0:
                bval_pair[1] += rsa_period
            elif match_pair[0] == 0 and match_pair[1] == max_ind:
                bval_pair[1] -= rsa_period
            return np.interp(signal, yval_pair, bval_pair) % rsa_period

        match_field_list = []
        signal_below_threshold = np.argwhere(rsa_yvals < signal)
        signal_above_threshold = np.argwhere(rsa_yvals >= signal)
        for ind in signal_above_threshold:
            max_ind = len(rsa_yvals) - 1
            last_ind = ind - 1 if ind > 0 else np.array([max_ind])
            next_ind = ind + 1 if ind < max_ind else np.array([0])
            if rsa_yvals[ind] == signal:
                match_field_list.append(1.0 * rsa_bvals[int(ind)])
                continue  # false positives will result from following code                
            if last_ind in signal_below_threshold:  # positive slope match
                match_pair = np.hstack([last_ind, ind])  # y must increase
                match_field_list.append(inverse_interpolate_rsa(match_pair))
            if next_ind in signal_below_threshold:  # negative slope match
                match_pair = np.hstack([ind, next_ind])
                match_field_list.append(inverse_interpolate_rsa(match_pair))
        return np.sort(match_field_list)

    fit_dnpcandidates_list = []
    fit_dnpcandidates_error_list = []
    for ind in range(len(bvals)):
        bval, yval = bvals[ind], yvals[ind]

        # try similar yvalues to look for potentially closer candidates
        if y_uncertainty is None:
            y_uncertainty = 0.05 * (max(rsa_yvals) - min(rsa_yvals))
#        yvals_to_test = np.array([yval])

#        # method 1: uncertainties via sampling y-vals
#        mu, sigma = yval, y_uncertainty
#        yvals_to_test = np.random.normal(mu, sigma, 20)
#        yvals_to_test.clip(min=min(rsa_yvals), max=max(rsa_yvals),
#                           out=yvals_to_test)
#        centroid_list = [0.0]
#        for test_yval in yvals_to_test:
#            yval_total_b_candidates = find_rsa_match_fields(test_yval)
#            yval_dnp_candidates = (yval_total_b_candidates - bval) % rsa_period
#            nrows, ncols = len(centroid_list), len(yval_dnp_candidates)
#            distance_matrix = np.zeros((nrows, ncols))
#            for row, centroid in enumerate(centroid_list):
#                for col, candidate in enumerate(yval_dnp_candidates):
#                    distance = \
#                        np.min(np.abs([centroid - candidate - rsa_period,
#                                       centroid - candidate,
#                                       centroid - candidate + rsa_period]))
#                    distance_matrix[row, col] = distance
#            # GOAL: pair up cols & rows to minimize sum of distance of
#            #       the elements (i, j) across each pair i-j. If more columns
#            #       than rows, we then add centroids at each "loser" column
#            #       location (so if recomputed, element at loser column and
#            #       new row would be 0.0). Either way we recompute centroids
#            #       by weighted avg (faster than recomputing from new list):
#            #       (N_vals_in_centroid * centroid val + new value) /
#            #                                          (N_vals_in_centroid + 1)
                

        # method 2: uncertainties via slope:
        yvals_to_test = np.linspace(yval - y_uncertainty / 2,
                                    yval + y_uncertainty / 2,
                                    3)  # must be odd to not lose orig. yval
        yvals_to_test.clip(min=min(rsa_yvals), max=max(rsa_yvals),
                           out=yvals_to_test)
        total_b_candidates = np.hstack([find_rsa_match_fields(test_yval)
                                        for test_yval in yvals_to_test])
        dnp_candidates = (total_b_candidates - bval) % rsa_period
        slope_mags = np.abs(get_rsa_signal_slope(total_b_candidates))
        slope_floor = 2 * y_uncertainty / rsa_period  # caps s_dnp @ period/2
        slope_mags.clip(min=slope_floor, out=slope_mags)
        dnp_uncertainties = y_uncertainty / slope_mags

        # IF NEEDED, ADD K-MEANS CLUSTERING AND COLLAPSE...

        fit_dnpcandidates_list.append(dnp_candidates)
        fit_dnpcandidates_error_list.append(dnp_uncertainties)

    return fit_dnpcandidates_list, fit_dnpcandidates_error_list

#       # ABANDONED TRACKING CODE
#        # Find which new DNP candidates adjacent to last DNP value
#        last_dnp_val = fit_dnpvals[ind - 1]  # looping ok, as = 0 initially
#        current_candidate = min(dnp_candidates)
#        next_candidate = current_candidate
#        # keep incrementing up until we straddle current DNP
#        counter = 0
#        while(next_candidate < last_dnp_val):
#            current_candidate = next_candidate
#            if current_candidate == max(dnp_candidates):  # step: increment
#                dnp_candidates += rsa_period
#                next_candidate = dnp_candidates[0]
#            else:
#                for candidate in dnp_candidates:
#                    if candidate > current_candidate:
#                        next_candidate = candidate
#                        break
#            if counter < 1000:  # DEBUGGING
#                counter += 1
#            else:
#                raise Exception("\n".join(
#                    ["help, stuck!",
#                     "last dnp = {}".format(last_dnp_val),
#                     "current_candidate: {}".format(current_candidate),
#                     "next_candidate: {}".format(next_candidate),
#                     "dnp_candidates: {}".format(dnp_candidates)]))
#        lower_candidate = current_candidate
#        upper_candidate = next_candidate
#        if lower_candidate < 0:  # but can't go below zero!
#            lower_candidate = upper_candidate
#
#        # Decide which DNP candidate is more likely
#        lower_candidate_delta_dnp = np.abs(lower_candidate - last_dnp_val)
#        upper_candidate_delta_dnp = np.abs(upper_candidate - last_dnp_val)
#        if upper_candidate_delta_dnp != 0:
#            delta_dnp_ratio = lower_candidate_delta_dnp / upper_candidate_delta_dnp
#        else:
#            if lower_candidate_delta_dnp != 0:
#                delta_dnp_ratio = np.inf
#            else:
#                delta_dnp_ratio = 1.0
#        log_threshold = np.log(tracking_threshold_ratio)
#        if delta_dnp_ratio == 1.0:  # only one candidate
#            new_dnp_val = lower_candidate
#        elif np.log(delta_dnp_ratio) <= -log_threshold:  # lower >= ~5x closer 
#            new_dnp_val = lower_candidate
#        elif np.log(delta_dnp_ratio) >= +log_threshold:  # upper >= ~5x closer 
#            new_dnp_val = upper_candidate
#        else:
#            if np.abs(yval) - np.abs(last_yval) > 0:  # assume RSA up = DNP up
#                new_dnp_val = upper_candidate
#            else:  # or vice versa
#                new_dnp_val = lower_candidate
#        if np.abs(new_dnp_val - last_dnp_val) < jump_threshold * x_uncertainty:
#            new_dnp_val = last_dnp_val  # RESET
#        fit_dnpvals[ind] = new_dnp_val
#        last_yval = yval

# need to test vs:
#   empty lists
#   repeating bfield values
#   different rsa field ranges, including negative and mixed bounds, weird #s
#   weird asymmetries around the rsa start/stop
#   rsa multiple sweeps to try and pull from


# %%  TESTING SCRIPTS
if False:
# %%
    ref_field_sweep_period_bvals = np.array([110, 115, 120, 125, 130])
    ref_field_sweep_period_yvals = np.array([-1, -1, 10, -1, -1])
    rsa_period = 20

    bvals = [0]
    yvals = [-1]
#    bvals = np.hstack([np.linspace(110, 130, 11),
#                       np.linspace(130, 110, 11),
#                       np.linspace(110, 130, 11),
#                       np.linspace(130, 110, 11)])
#    yvals = np.hstack([00.0 * np.ones(11),
#                       20.0 * np.ones(11),
#                       00.0 * np.ones(11),
#                       20.0 * np.ones(11)])
    extract_dnp_from_field_sweep(bvals, yvals,
                                 ref_field_sweep_period_bvals,
                                 ref_field_sweep_period_yvals,
                                 rsa_period)  # should be +0.5, +19.5


# %%  TESTING SCRIPTS (CONT.)
if False:
# %%
    import matplotlib.pyplot as plt

    time_span = 3000.0  # s
    dnp_time_const = 300.0
    dnp_magnitude = 4.0
    warump_time = 6000.0
    warump_bfield = 110.0

    b_min = 110.0
    b_max = 130.0
    b_offset = 0.0
    b_min += b_offset
    b_max += b_offset

    spin_lifetime = 20000.0
    gfactor = 0.44
    initial_phase = 0.0
    delay_time = -60
    pulse_amplitude = 1.0

    def dnp_field_fcn(bfield, last_dnp_field):  # closure, captures above vars
        effective_bfield = bfield + last_dnp_field
        net_polarization, net_phase = get_pulse_sum_vector(spin_lifetime,
                                                           gfactor,
                                                           effective_bfield,
                                                           initial_phase)

        # DNP style 1 - proportional to spin signal
        equilibrium_dnp_field = dnp_magnitude * net_polarization

#        # DNP style 2 - proportional to capped spin signal
#        equilibrium_dnp_field = \
#            dnp_magnitude * np.min([1.0, 0.7 * net_polarization])

#        # DNP style 3 - proportional to spin signal past threshold
#        equilibrium_dnp_field = \
#            dnp_magnitude * (np.max([1.0, 1.0 * net_polarization]) - 1.0)

#        # DNP constant offset
#        equilibrium_dnp_field += 0.5 * dnp_magnitude

        return equilibrium_dnp_field

    nvals = 2400
    tvals = np.linspace(0, time_span, nvals)  # in SECONDS

    nsteps = 600  # must be factor of nvals
    b_steps = np.hstack([np.linspace(b_min, b_max, int(np.floor(nsteps / 4))),
                         np.linspace(b_max, b_min, int(np.ceil(nsteps / 4))),
                         np.linspace(b_min, b_max, int(np.floor(nsteps / 4))),
                         np.linspace(b_max, b_min, int(np.ceil(nsteps / 4)))])
    bvals = np.hstack([stepval * np.ones(int(nvals / nsteps))
                       for stepval in b_steps])
#    bvals[:int(np.floor(nvals / 2))] = 103

#    # forward/back sweep
#    bvals = np.hstack([np.linspace(100, 130, int(np.floor(nvals / 2))),
#                       np.linspace(130, 100, int(np.ceil(nvals / 2)))])

#    bvals = np.linspace(100, 130, nvals) # linear sweep

#    bvals = 100 * np.ones_like(tvals)  # step fcn
#    bvals[int(nvals/2):] = 120

    initial_dnp = 0.0
    for i in range(10):  # get stable interdependent DNP field & RSA:
        initial_dnp = get_new_dnp_field(warump_time, warump_bfield, initial_dnp,
                                        target_dnp_fcn=dnp_field_fcn,
                                        dnp_time_const=dnp_time_const)
    dnpvals = initial_dnp * np.ones_like(tvals)
    delta_t = tvals[1] - tvals[0]
    for ind, (tval, bval) in enumerate(zip(tvals, bvals)):
        if ind > 0:
            dnpvals[ind] = get_new_dnp_field(delta_t, bval, dnpvals[ind - 1],
                                             target_dnp_fcn=dnp_field_fcn,
                                             dnp_time_const=dnp_time_const)
    signal = fitfcn_trkr_single_species(delay_time, {},
                                        pulse_amplitude, gfactor,
                                        bvals + dnpvals, spin_lifetime,
                                        initial_phase, slope=0, offset=0)
    
    def fancy_plot_vs_time(tvals, bvals, signal, axes):
        segment_xvals = []
        segments = []
        last_xval = -999999
        current_segment_yvals = []
        for xval, yval in zip(bvals[:],
                              signal[:]):
            if abs(xval - last_xval) > 0.0001:  # new segment found
                if last_xval != -999999:  # break off old segment
                    segment_xvals.append(last_xval)
                    segments.append(np.array(current_segment_yvals))
                last_xval = xval
                current_segment_yvals = []
            current_segment_yvals.append(yval)
        current_start_index = 0
        for segment_index, (xval, yvals) in enumerate(zip(segment_xvals, segments)):
            if segment_index % 2 == 0:
                segment_color = 'b'
            else:
                segment_color = 'r'
            indices = current_start_index + np.arange(len(yvals))
            times = tvals[:][indices]
            current_start_index += len(yvals)
            axes.plot(times, yvals, '-', color=segment_color)
        #    axes.plot(xval * np.ones_like(yvals), yvals, 'd', color=segment_color)

    plt.figure()
    plt.subplot(5,1,1)
    plt.plot(tvals, bvals)
    plt.xlabel('Time (s)')
    plt.ylabel('B (mT)')
    plt.subplot(5,1,2)
    plt.plot(tvals, dnpvals)
    plt.xlabel('Time (s)')
    plt.ylabel('DNP (mT)')

    ax = plt.subplot(5,1,3)
    fancy_plot_vs_time(tvals, bvals, signal, ax)
    plt.xlabel('Time (s)')
    plt.ylabel('TRKR (AU)')

    quarterway_index = int(np.floor(nvals / 4))
    halfway_index = int(np.floor(nvals / 2))
    threequarterway_index = int(np.floor(3 * nvals / 4))

    ax = plt.subplot(5,2,7)
    plt.plot(bvals[:quarterway_index],
             signal[:quarterway_index])
    plt.ylabel('TRKR (AU)')
    plt.text(0.5, 0.9, 'Sweep Up #1', transform = ax.transAxes,
             horizontalalignment='center', verticalalignment='top')

    ax = plt.subplot(5,2,9)
    plt.plot(bvals[quarterway_index:halfway_index],
             signal[quarterway_index:halfway_index])
    plt.ylabel('TRKR (AU)')
    plt.xlabel('External Field (mT)')
    plt.text(0.5, 0.9, 'Sweep Down #1', transform = ax.transAxes,
             horizontalalignment='center', verticalalignment='top')

    ax = plt.subplot(5,2,8)
    plt.plot(bvals[halfway_index:threequarterway_index],
             signal[halfway_index:threequarterway_index])
    plt.text(0.5, 0.9, 'Sweep Up #2', transform = ax.transAxes,
             horizontalalignment='center', verticalalignment='top')

    ax = plt.subplot(5,2,10)
    plt.plot(bvals[threequarterway_index:],
             signal[threequarterway_index:])
    plt.xlabel('External Field (mT)')
    plt.text(0.5, 0.9, 'Sweep Down #2', transform = ax.transAxes,
             horizontalalignment='center', verticalalignment='top')

    rsa_period = 1.0 / (LASER_REPRATE * GFACTORCONSTANT * gfactor)
    ref_field_sweep_period_bvals = np.linspace(max(bvals),
                                               max(bvals) + rsa_period,
                                               200)
    ref_field_sweep_period_yvals = \
        fitfcn_trkr_single_species(delay_time, {},
                                   pulse_amplitude, gfactor,
                                   ref_field_sweep_period_bvals, 
                                   spin_lifetime, initial_phase,
                                   slope=0, offset=0)
    reconstructed_dnp, reconstructed_dnp_err = \
        extract_dnp_from_field_sweep(bvals, signal,
                                     ref_field_sweep_period_bvals,
                                     ref_field_sweep_period_yvals,
                                     rsa_period)
    plt.figure()
    ax = plt.subplot(3,1,1)
    plt.plot(tvals, bvals)
    plt.ylabel('B (mT)')
    ax = plt.subplot(3,1,2)
    fancy_plot_vs_time(tvals, bvals, signal, ax)
    plt.ylabel('TRKR (AU)')
    ax = plt.subplot(3,1,3)


    for tval, bval, dnp_list, dnp_err_list in zip(tvals, bvals,
                                                  reconstructed_dnp,
                                                  reconstructed_dnp_err):
        t_array = tval * np.ones_like(dnp_list)
#        plt.plot(t_array, dnp_list, 'd')
        plt.errorbar(t_array, dnp_list, yerr=None,  # dnp_err_list,
                     fmt='bd')
    plt.plot(tvals, dnpvals, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('Est. DNP (mT)')



# %%  TESTING SCRIPTS (CONT.)
if False:
# %%
    bfields = np.arange(110, 150+1, 1)
    amplitudes = []
    phases = []
    for bfield in bfields:
        amplitude, phase = get_pulse_sum_vector(spin_lifetime, gfactor, bfield, initial_phase=0)
        amplitudes.append(amplitude)
        phases.append(phase)
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(bfields, amplitudes)
    plt.subplot(2,1,2)
    plt.plot(bfields, phases)

# %%  TESTING SCRIPTS (CONT.)
if False:
# %%
    pass
    # NEW DATA: since all in -500 to 0 range, ignore exponential decay
    # and possible 2nd species, and just fit to phase of a cosine!


