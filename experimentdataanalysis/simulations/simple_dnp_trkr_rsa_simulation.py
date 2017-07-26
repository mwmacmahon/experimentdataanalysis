# General imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from experimentdataanalysis.simulations.simple_trkr_rsa_simulation \
    import get_pulse_sum_vector, \
           trkr_decaying_cosine

# GLOBAL CONSTANTS
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13158  # ps period


# function defs for reference
# get_pulse_sum_vector(bfield, spin_lifetime, gfactor, initial_phase=0)
# trkr_decaying_cosine(delay_time, total_bfield,
#                      pulse_amplitude,
#                      gfactor, spin_lifetime,
#                      initial_phase, extra_phase_offset,
#                      slope, offset)

























# %%
def trkr_decaying_cosine(delay_time, total_bfield,
                         pulse_amplitude,
                         gfactor, spin_lifetime,
                         initial_phase, extra_phase_offset,
                         slope, offset):
    zero_delay_offset = 0.0
    pos_def_delay = (delay_time + zero_delay_offset) % LASER_REPRATE
    osc_ang_freq = 2 * np.pi * GFACTORCONSTANT * gfactor * total_bfield
    net_polarization, net_phase = get_pulse_sum_vector(spin_lifetime,
                                                       gfactor, total_bfield,
                                                       initial_phase)
    final_phase = (net_phase + pos_def_delay * osc_ang_freq) % (2 * np.pi)
    final_amplitude = pulse_amplitude * net_polarization * \
                            np.exp(-pos_def_delay / spin_lifetime)
    signal = final_amplitude * np.cos(final_phase + extra_phase_offset)
    output = signal + delay_time * slope + offset  # NOT pos-definite
    return output


# %%
def generate_TRKR_simulation_params(nxvals, ndatasets,
                                    simulation_constants, seed=None):
    if seed is not None:
        np.random.seed(seed)
    # extract the parameters used to determine the parameters
    pulse_amplitude_mean = simulation_constants['pulse_amplitude_mean']
    pulse_amplitude_std_err = simulation_constants['pulse_amplitude_std_err']
    gfactor = simulation_constants['gfactor']
    spin_lifetime_mean = simulation_constants['spin_lifetime_mean']
    spin_lifetime_std_err = simulation_constants['spin_lifetime_std_err']
    initial_phase = simulation_constants['initial_phase']
    extra_phase_offset_baseline = simulation_constants['extra_phase_offset_baseline']
    extra_phase_offset_cos_amp = simulation_constants['extra_phase_offset_cos_amp']
    extra_phase_offset_cos_nperiods = simulation_constants['extra_phase_offset_cos_nperiods']
    slopes_scale = simulation_constants['slopes_scale']
    offsets_scale = simulation_constants['offsets_scale']
    noise_scale = simulation_constants['noise_scale']

    # generate actual parameter values for each simulated dataset
    pulse_amplitudes = pulse_amplitude_mean + \
                        np.random.normal(size=ndatasets, scale=pulse_amplitude_std_err)
    gfactors = gfactor * np.ones(ndatasets)
    spin_lifetimes = spin_lifetime_mean + \
                        np.random.normal(size=ndatasets, scale=spin_lifetime_std_err)
    initial_phases = initial_phase * np.ones(ndatasets)
    extra_phase_offsets = \
        (extra_phase_offset_baseline + extra_phase_offset_cos_amp *
            np.sin((2 * np.pi * extra_phase_offset_cos_nperiods / ndatasets) *
                                                               np.arange(ndatasets)))
    slopes = np.random.normal(size=ndatasets, scale=slopes_scale)
    offsets = np.random.normal(size=ndatasets, scale=offsets_scale)
    noisefcn = lambda nx: np.random.normal(size=nx, scale=noise_scale)
    noise_layers = [noisefcn(nxvals) for dataset_index in range(ndatasets)]
    dataset_model_params_dicts = \
        [{'pulse_amplitude'    : pulse_amplitudes[dataset_index],
          'gfactor'            : gfactors[dataset_index],
          'spin_lifetime'      : spin_lifetimes[dataset_index],
          'initial_phase'      : initial_phases[dataset_index],
          'extra_phase_offset' : extra_phase_offsets[dataset_index],
          'slope'              : slopes[dataset_index],
          'offset'             : offsets[dataset_index]}
         for dataset_index in range(ndatasets)]
    simulation_params = {
        'dataset_model_params_dicts': dataset_model_params_dicts,
        'noise_layers': noise_layers,
        # also save raw parameters for easy lookup
        'pulse_amplitudes': pulse_amplitudes,
        'gfactors': gfactors,
        'spin_lifetimes': spin_lifetimes,
        'initial_phases': initial_phases,
        'extra_phase_offsets': extra_phase_offsets,
        'slopes': slopes,
        'offsets': offsets,
    }
    return simulation_params


# %%
def generate_TRKR_simulation_dataframe(tvals, dataset_bvals, simulation_params,
                                       suppress_plot=False):
    """For function trkr_decaying_cosine() in this module"""
    nx = len(tvals)
    ndatasets = len(dataset_bvals)
    indices_1d = np.arange(nx)
    indices_2d = np.arange(ndatasets)
    dataset_model_params_dicts = simulation_params['dataset_model_params_dicts']
    noise_layers = simulation_params['noise_layers']
    scan_1d_results= []
    for index_2d in indices_2d:
        delay_times = tvals
        dataset_b_external = dataset_bvals[index_2d]
        dataset_params_dict = dataset_model_params_dicts[index_2d]
        yvals = trkr_decaying_cosine(delay_times,
                                     dataset_b_external,
                                     **dataset_params_dict)
        noisy_yvals = yvals + noise_layers[index_2d]
        scan_1d_results.append(noisy_yvals)
    X_bvals, X_tvals = np.meshgrid(dataset_bvals, tvals, indexing='ij',
                                   sparse=False, copy=True)  # not sure on ideal settings here
    independent_data_matrices = [X_tvals, X_bvals]
    measured_data = np.array(scan_1d_results)

    if not suppress_plot:
        plt.figure()
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
