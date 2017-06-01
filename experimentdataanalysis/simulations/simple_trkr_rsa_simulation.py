# General imports
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# GLOBAL CONSTANTS
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13158  # ps period


# %%
def get_pulse_sum_vector(bfield, spin_lifetime, gfactor, initial_phase=0):
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
def trkr_decaying_cosine(delay_time, total_bfield,
                         pulse_amplitude,
                         gfactor, spin_lifetime,
                         initial_phase, slope, offset):
    trkr_per_unit_polarization = 1.0
    trkr_phase_offset = 0.0
    pos_def_delay = (delay_time + zero_delay_offset) % LASER_REPRATE
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


def generate_TRKR_simulation_params(simulation_constants, seed=None):
    if seed is not None:
        np.random.seed(seed)
    gfactor = simulation_constants['gfactor']
    amplitude_baseline = simulation_constants['amplitude_baseline']
    amplitude_sin_amp = simulation_constants['amplitude_sin_amp']
    amplitude_sin_nperiods = simulation_constants['amplitude_sin_nperiods']
    phase_offset_baseline = simulation_constants['phase_offset_baseline']
    phase_offset_cos_amp = simulation_constants['phase_offset_cos_amp']
    phase_offset_cos_nperiods = simulation_constants['phase_offset_cos_nperiods']
    y_offsets_scale = simulation_constants['y_offsets_scale']
    noise_scale = simulation_constants['noise_scale']
    gfactors = gfactor * np.ones(ndatasets)
    amplitudes = \
        (amplitude_baseline + amplitude_sin_amp *
            np.sin((2 * np.pi * amplitude_sin_nperiods / ndatasets) * np.arange(ndatasets)))
    phase_offsets = \
        (phase_offset_baseline + phase_offset_cos_amp *
            np.sin((2 * np.pi * phase_offset_cos_nperiods / ndatasets) * np.arange(ndatasets)))
    y_offsets = np.random.normal(size=ndatasets, scale=y_offsets_scale)
    noisefcn = lambda: np.random.normal(size=nx, scale=noise_scale)
    noise_layers = [noisefcn() for dataset_index in range(ndatasets)]
    simulation_params = {
        'gfactors': gfactors,
        'amplitudes': amplitudes,
        'phase_offsets': phase_offsets,
        'y_offsets': y_offsets,
        'noise_layers': noise_layers
    }
    return simulation_params


def generate_TRKR_simulation_dataframe(tvals, dataset_bvals, simulation_params,
                                       suppress_plot=False):
    ndatasets = len(dataset_bvals)
    nx = len(tvals)
    gfactors = simulation_params['gfactors']
    amplitudes = simulation_params['amplitudes']
    phase_offsets = simulation_params['phase_offsets']
    y_offsets = simulation_params['y_offsets']
    noise_layers = simulation_params['noise_layers']
    indices_1d = np.arange(nx)
    indices_2d = np.arange(ndatasets)
    scan_1d_results= []
    for index_2d in indices_2d:
        delay_times = tvals
        b_external = dataset_bvals[index_2d]
        amplitude = amplitudes[index_2d]
        gfactor = gfactors[index_2d]
        phase_offset = phase_offsets[index_2d]
        y_offset = y_offsets[index_2d]
        noise = noise_layers[index_2d]
        yvals = fitfcn_cosine(delay_times, b_external,
                              gfactor, amplitude, phase_offset, y_offset)
        noisy_yvals = yvals + noise
        scan_1d_results.append(noisy_yvals)

    X_bvals, X_tvals = np.meshgrid(dataset_bvals, tvals, indexing='ij',
                                   sparse=False, copy=True)  # not sure on ideal settings here
    independent_data_matrices = [X_tvals, X_bvals]
    measured_data = np.array(scan_1d_results)

    plt.imshow(measured_data, interpolation='none', aspect=nx/ndatasets)

    # pandas dataframe conversion
    run_ids = np.zeros(measured_data.size, dtype=np.int)
    indices_2d, indices_1d = np.meshgrid(np.arange(len(dataset_bvals)),
                                         np.arange(len(tvals)),
                                         indexing='ij', sparse=False, copy=True)
    dataframe = pd.DataFrame({'run_id'        : run_ids,
                              'index_2d'      : indices_2d.flatten(),
                              'index_1d'      : indices_1d.flatten(),
                              'b_external'    : X_bvals.flatten(),
                              'probe_delay'   : X_tvals.flatten(),
                              'kerr_rotation' : measured_data.flatten(),
                             })
    dataframe.set_index(['run_id', 'index_2d', 'index_1d'], drop=True, append=False, inplace=True)
    dataframe.sort_index(ascending=True, inplace=True)  # not actually necessary, but nice to be sure
    return dataframe
