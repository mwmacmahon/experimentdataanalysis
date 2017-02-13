# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:32:32 2016

Designed to run simple toy models of drift-diffusion spin systems
with one or more electron bands.

# TODO: potential speedup from taking same packet for _every_ packet,
# and just translating to the side. This only works without diffusion,
# but perhaps that's not such a bad thing in this case. Consider leaving
# it out.

@author: Michael
"""

import matplotlib.pyplot as plt
import numpy as np


GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period


# GOAL: make this part hyper-efficient by only making a few original
# gaussians and looking up repeats. Perhaps by, for example, making
# any gaussian of common width X by taking the same gaussian and making
# every other gaussian a view of the original with different indexing and
# a different scaling factor. Can 
def gaussian_profile(xvals, A, mu, sigma):
    xvals = np.array(xvals, copy=False)
    A = np.array(A, copy=False)
    mu = np.array(mu, copy=False)
    sigma = np.array(sigma, copy=False)
    gaussian = A * np.exp(-(xvals - mu)**2 / (2 * sigma**2))
    return gaussian


def get_absorption_amp_phase(laser_width, laser_power, laser_wavelength, 
                             species_physics_dict):
    if 'fudge' in species_physics_dict.keys():
        # model in use: Absorption amp, phase identical for species
        area_under_gaussian = laser_power
        amplitude = area_under_gaussian / (laser_width * np.sqrt(2 * np.pi))
        phase = 0.0
#        # this model not in use
#        relative_amp = species_physics_dict['fudge']['relative_amp']
#        relative_phase = species_physics_dict['fudge']['relative_phase']
#        amplitude = laser_power * relative_amp  
#        phase = relative_phase
    else:
        raise NotImplementedError("No actual physics yet, need fudge factors")
    return amplitude, phase


def get_kerr_rotation_amp_phase(laser_width, laser_power, laser_wavelength, 
                                species_physics_dict):
    if 'fudge' in species_physics_dict.keys():
        # model in use: Kerr amp, phase different for species
        relative_amp = species_physics_dict['fudge']['relative_amp']
        relative_phase = species_physics_dict['fudge']['relative_phase']
        area_under_gaussian = laser_power * relative_amp
        amplitude = area_under_gaussian / (laser_width * np.sqrt(2 * np.pi))
        phase = relative_phase
#        # this model not in use
#        amplitude = laser_power
#        phase = 0.0
    else:
        raise NotImplementedError("No actual physics yet, need fudge factors")
    return amplitude, phase


class CarrierPacket():
    def __init__(self, time, laser_position,
                 laser_width, laser_power, laser_wavelength,
                 species_physics_dict):
        # NOTE: all of this COMPLETELY independent of whether we
        # we want to do RSA, TRKR, spatial profiles, etc.
        self.physics_dict = species_physics_dict
        self.generation_time = time
        self.initial_position = laser_position
        self.initial_width = laser_width
        # more like polarization * (# affected carriers), oh well
        self.initial_polarization, self.initial_orientation = \
            get_absorption_amp_phase(laser_width,
                                     laser_power, laser_wavelength, 
                                     self.physics_dict)

    def get_time_elapsed_current_state(self, tvals, evals, bvals):
        """
        dependending on if tvals, evals, or bvals is a vector,
        may return a vector of same size. if none are vectors, returns
        a set of scalar values representing a future state of spin
        packet at the given parameters. Not intended to have more
        than one non-scalar parameters.
        
        multiple tvals: all outputs vectors except packet width (no diffusion)
        multiple evals: all outputs scalars except packet position
        multiple bvals: all outputs scalars except packet orientation
        
        Note: diffusion is NOT implemented in this version.
        """
        tvals = np.array(tvals)
        evals = np.array(evals)
        bvals = np.array(bvals)
        num_vecs = np.count_nonzero([xvals.size != 1
                                     for xvals in [tvals, evals, bvals]])
        
        if num_vecs > 1:
            raise ValueError("no more than one input to " +
                             "get_time_elapsed_current_state can be " +
                             "nonscalar")
        elapsed_tvals = tvals - self.generation_time
        if tvals.size == 1 and elapsed_tvals < 0:
            raise ValueError("attempting to get state of a spin packet " +
                             "before it was created.")
        elif tvals.size > 1 and np.min(elapsed_tvals) < 0:
            raise NotImplementedError("have not yet written code to " +
                                      "handle a time slice spanning both " +
                                      "before AND after pulse creation.")
        electron_temp = self.physics_dict['temperature_vs_efield_fcn'](evals)
        effective_initial_polarization = \
            self.initial_polarization * \
                self.physics_dict['absorption_vs_temp_scaling_fcn'](electron_temp)
        spin_lifetime = \
            self.physics_dict['spin_lifetime'] * \
                self.physics_dict['lifetime_vs_temp_scaling_fcn'](electron_temp)
        gfactor = \
            self.physics_dict['gfactor'] + \
                self.physics_dict['gfactor_vs_temp_offset_fcn'](electron_temp)
        polarization = (effective_initial_polarization *
                        np.exp(-elapsed_tvals / spin_lifetime))
        orientation = (self.initial_orientation +
                       2 * np.pi * elapsed_tvals * bvals * gfactor *
                                                               GFACTORCONSTANT)
        orientation = ((orientation + np.pi/2) % (2*np.pi)) - np.pi/2
        position = (self.initial_position + elapsed_tvals *
                    self.physics_dict['mobility'] * evals)
        width = self.initial_width
        return (elapsed_tvals, polarization, orientation,
                position, width)

    def get_complex_spin_profile(self, tvals, xvals,
                                 applied_E_fields, external_B_fields):
        tvals = np.array(tvals, copy=False)
        xvals = np.array(xvals, copy=False)
        evals = np.array(applied_E_fields, copy=False)
        bvals = np.array(external_B_fields, copy=False)
        elapsed_tvals, packet_polarization, \
            packet_orientation, packet_position, packet_width = \
                self.get_time_elapsed_current_state(tvals, evals, bvals)

        # above variables may be a wide variety of scalars and vectors,
        # but should all work itself out to a single proper vector result
        profile = gaussian_profile(xvals, packet_polarization,
                                   packet_position, packet_width)
        complex_profile = profile * np.exp(1j * packet_orientation)
        return complex_profile

    def get_probe_spin_convolution(self, tvals, xvals,
                                   laser_width, laser_power, laser_wavelength,
                                   applied_E_fields, external_B_fields):
        tvals = np.array(tvals, copy=False)
        xvals = np.array(xvals, copy=False)
        evals = np.array(applied_E_fields, copy=False)
        bvals = np.array(external_B_fields, copy=False)
        elapsed_tvals, packet_polarization, \
            packet_orientation, packet_position, packet_width = \
                self.get_time_elapsed_current_state(tvals, evals, bvals)

        # above variables may be a wide variety of scalars and vectors,
        # but should all work itself out to a single proper vector result.
        # note convolution equation assumes an area under curve of 1.0,
        # so we implicity divided each gaussian by area under gaussian,
        # so the corrected convolution result is the normalized
        # convolution result times the area under each gaussian
        convolution_mu = packet_position
        convolution_sigma = np.sqrt(laser_width**2 + packet_width**2)
        normalized_convolution = gaussian_profile(xvals, 1.0,
                                                  convolution_mu,
                                                  convolution_sigma)
        area_under_laser_gaussian = laser_power
        area_under_spin_gaussian = (packet_polarization *
                                    packet_width * np.sqrt(2 * np.pi))
        probe_spin_convolution = (normalized_convolution *
                                    area_under_laser_gaussian *
                                        area_under_spin_gaussian)
        amp_coeff, phase_offset = \
            get_kerr_rotation_amp_phase(laser_width,
                                        laser_power, laser_wavelength, 
                                        self.physics_dict)
        amplitude = amp_coeff * probe_spin_convolution
        profile = amplitude * np.cos(phase_offset + packet_orientation)
        return profile


class ExperimentChannel():
    def __init__(self, *species_physics_dicts):
        self.species_physics_dict_list = list(species_physics_dicts)
        self.packet_list = []

    def generate_pump_packets(self, time, laser_position,
                              laser_width, laser_power, laser_wavelength):
        for species_dict in self.species_physics_dict_list:
            species_packet = \
                CarrierPacket(time, laser_position,
                              laser_width, laser_power, laser_wavelength,
                              species_dict)
            self.packet_list.append(species_packet)

    def get_valid_packets_for_profile(self, tvals, xvals, evals, bvals):
        """
        rule: only one of the variables should be a vector!
        """
        num_vecs = np.count_nonzero([xvals.size != 1
                                     for xvals in [tvals, xvals,
                                                   evals, bvals]])
        if num_vecs == 0:
            raise ValueError("at least one input to " +
                             "get_valid_packets_for_profile must be " +
                             "nonscalar")
        if num_vecs > 1:
            raise ValueError("no more than one input to " +
                             "get_valid_packets_for_profile can be " +
                             "nonscalar")
        if tvals.size == 1:
            valid_packet_list = [packet
                                 for packet in self.packet_list
                                 if packet.generation_time <= tvals]
        else:
            valid_packet_list = [packet
                                 for packet in self.packet_list
                                 if packet.generation_time <= np.min(tvals)]
            if len(valid_packet_list) < len(self.packet_list):
                print("Warning: current version cannot handle time windows " +
                      "spanning the creation of a pulse packet. Packets " +
                      "within this window are being ignored.")
        return valid_packet_list
        

    def get_complex_spin_profile(self, tvals, xvals,
                                 applied_E_fields, external_B_fields):
        tvals = np.array(tvals, copy=False)
        xvals = np.array(xvals, copy=False)
        evals = np.array(applied_E_fields, copy=False)
        bvals = np.array(external_B_fields, copy=False)
        valid_packet_list = self.get_valid_packets_for_profile(tvals, xvals,
                                                               evals, bvals)
        packet_profiles = [packet.get_complex_spin_profile(tvals, xvals,
                                                           evals, bvals)
                           for packet in valid_packet_list]
        return np.sum(np.vstack(packet_profiles), axis=0)

    def get_probe_profile(self, tvals, xvals,
                          laser_width, laser_power, laser_wavelength,
                          applied_E_fields, external_B_fields):
        tvals = np.array(tvals, copy=False)
        xvals = np.array(xvals, copy=False)
        evals = np.array(applied_E_fields, copy=False)
        bvals = np.array(external_B_fields, copy=False)
        valid_packet_list = self.get_valid_packets_for_profile(tvals, xvals,
                                                               evals, bvals)
        packet_profiles = [packet.get_probe_spin_convolution(tvals, xvals,
                                                             laser_width,
                                                             laser_power,
                                                             laser_wavelength,
                                                             evals, bvals)
                           for packet in valid_packet_list]
        return np.sum(np.vstack(packet_profiles), axis=0)


# %%
if __name__ == "__main__":
# %%
    species1_physics_dict = {'fudge': {'relative_amp': 1.0,
                                       'relative_phase': np.pi},
                             'gfactor': 0.4385,
                             'mobility': 1e-4,  # (um/ps)/(V/cm). 1e-4: 2um/ns/Vapp
                             'spin_lifetime': 20260,
                             'temperature_vs_efield_fcn': 
                                 lambda E: 30 + 0.03 * np.abs(E)**2.0,
                             'absorption_vs_temp_scaling_fcn':
                                 lambda T: (T / 30)**-0.8,
                             'lifetime_vs_temp_scaling_fcn':
                                 lambda T: (T / 30)**-1.2,
                             'gfactor_vs_temp_offset_fcn':
                                 lambda T: (T - 30) * (-3.5e-4)}

    species2_physics_dict = {'fudge': {'relative_amp': 2.0,
                                       'relative_phase': 0},
                             'gfactor': 0.4385,
                             'mobility': 1e-4,  # (um/ps)/(V/cm). 1e-4: 2um/ns/Vapp
                             'spin_lifetime': 8000,
                             'temperature_vs_efield_fcn': 
                                 lambda E: 30 + 0.03 * np.abs(E)**2.0,
                             'absorption_vs_temp_scaling_fcn':
                                 lambda T: (T / 30)**-1.2,
                             'lifetime_vs_temp_scaling_fcn':
                                 lambda T: (T / 30)**-1.2,
                             'gfactor_vs_temp_offset_fcn':
                                 lambda T: (T - 30) * (-3.5e-4)}

    laser_kwargs = {'laser_width': 17.5,  # um, radius
                    'laser_power': 1.0,  # AU
                    'laser_wavelength': 818.9}  # nm
    
    experiment_state_kwargs = {'applied_E_fields': 15,  # V/cm
                               'external_B_fields': 300}  # mT
    
    n_pulses = 20
    laser_pos = 0.0

    vanilla_experiment = ExperimentChannel(species1_physics_dict,
                                           species2_physics_dict)
    current_time = -n_pulses * LASER_REPRATE  # last pulse at t=0
    for i in range(n_pulses):
        current_time += LASER_REPRATE
        vanilla_experiment.generate_pump_packets(current_time, laser_pos,
                                                 **laser_kwargs)


# %% EXPORT TRKR DATA WITH NOISE
if False:
# %%
    tvals = np.linspace(0, 7500, 200)  # ps
    probe_position_list = [0]  # um
#    probe_position_list = np.arange(-20, 100+1, 2)
#    applied_E_fields = [15]  # V/cm
    applied_E_fields = np.arange(-40, 40+1, 4)
    external_B_fields = [300]  # mT
    num_fake_runs = 5
    for run_ind in range(num_fake_runs):
        for probe_position in probe_position_list:
            for current_B_field in external_B_fields:
                for applied_E_field in applied_E_fields:
                    experiment_field_dict = \
                        {'applied_E_fields': applied_E_field,
                         'external_B_fields': current_B_field}
                    probe_profile = vanilla_experiment.get_probe_profile(
                        tvals, probe_position,
                        **laser_kwargs, **experiment_field_dict)
                    probe_noise = 0.0005 * np.random.randn(tvals.size)
                    probe_offset = .005 * (2 * (0.5 - np.random.random()))
                    probe_slope = .005 * np.random.randn()
                    probe_linear_offset = \
                        probe_offset + probe_slope * np.linspace(-0.5, 0.5,
                                                                 tvals.size)
                    probe_profile += probe_noise + probe_linear_offset

#                    plt.plot(tvals, probe_profile, 'bd')
#                    plt.plot(tvals,
#                             probe_profile + probe_noise + probe_linear_offset, 'r')
#                    plt.title('Spin profile w/ probe convolution and Kerr physics')
#                    raise KeyboardInterrupt

                    directory = "C:\\Data\\fake_data\\fake_trkr"
                    filename = "TRKR_MirrorZ_{:.0f}_".format(probe_position) + \
                               "{:.0f}Vcm_".format(float(applied_E_field)) + \
                               "{:.0f}mT_".format(float(current_B_field)) + \
                               "30K_818.9nm_run{:03d}.dat".format(run_ind + 1)
                    filepath = directory + "\\" + filename
        #            headerstr = "\n".join(10*["header header header"])
                    headerstr = '\t'.join(["scancoord","lockin2x"])
                    np.savetxt(filepath,
                               np.vstack([tvals, probe_profile]).T,
                               fmt='%.6e\t%.12e', delimiter='\t', newline='\n',
                               header=headerstr, comments='')


# %% STREAKED POSITION VS TIME OVER MANY, MANY TIME WINDOWS
if False:
# %%
    xvals = np.linspace(-60, 120, 180)
    nrows = 3
    ncols = 8
    n_segments = nrows * ncols
    for seg_ind in range(n_segments):
        n_streaks = 20  
        streak_length = (164 /
                         (nrows * experiment_state_kwargs['external_B_fields']))
        current_time = 0 + streak_length * seg_ind
        greyness = 0.0
        box_min = 0.0
        box_max = 0.0
        for i in range(n_streaks):
            color = np.array([0.0 + 0.8 * (1.0 - i / n_streaks)**1.5,
                              0.1 + 0.6 * (0.0 + i / n_streaks)**1.5,
                              0.2 + 0.8 * (0.0 + i / n_streaks)**1.5])
            color = (1 - greyness) * color + greyness * np.array([1, 1, 1])
            color = tuple(color)
            plt.subplot(nrows, ncols,
                        ncols * (seg_ind % nrows) + (seg_ind // nrows) + 1)
            probe_profile = vanilla_experiment.get_probe_profile(
                current_time + streak_length * i / n_streaks,
                xvals, **laser_kwargs, **experiment_state_kwargs)
            plt.plot(xvals, probe_profile, color=(color))
            plt.ylim([-0.01, 0.01])
            plt.xticks([])
            plt.yticks([])
            title_string = \
                "t: {:.0f} to {:.0f}".format(current_time,
                                             current_time + streak_length)
            plt.title(title_string, {'fontsize': 12})
            if min(probe_profile) < box_min:
                box_min = min(probe_profile)
            if max(probe_profile) > box_max:
                box_max = max(probe_profile)
        plt.plot(xvals, np.full_like(probe_profile, box_min), 'r')
        plt.plot(xvals, np.full_like(probe_profile, box_max), 'r')


# %% COMPARISON OF DIRECT SPATIAL SPIN PROFILE VS MEASURED
if False:
# %%
    xvals = np.linspace(-60, 120, 180)
    n_streaks = 100
    current_time = 400
    streak_length = 300
    for i in range(n_streaks):
        greyness = 0.2
        color = np.array([0.0 + 0.8 * (0.0 + i / n_streaks)**1.5,
                          0.1 + 0.6 * (1.0 - i / n_streaks)**1.5,
                          0.2 + 0.8 * (1.0 - i / n_streaks)**1.5])
        color = (1 - greyness) * color + greyness * np.array([1, 1, 1])
        color = tuple(color)
        plt.subplot(2,1,1)
        complex_spin_profile = vanilla_experiment.get_complex_spin_profile(
            current_time + streak_length * i / n_streaks,
            xvals, **experiment_state_kwargs)
        plt.plot(xvals, np.real(complex_spin_profile), color=(color))
        plt.title('Direct spin profile')
        plt.subplot(2,1,2)
        probe_profile = vanilla_experiment.get_probe_profile(
            current_time + streak_length * i / n_streaks,
            xvals, **laser_kwargs, **experiment_state_kwargs)
        plt.plot(xvals, probe_profile, color=(color))
        plt.title('Spin profile w/ probe convolution and Kerr physics')


# %% COMPARISON OF DIRECT POLAR. VS TIME VS TRKR W/ PROBE CONVOLUTION
if True:
# %%
    tvals = np.linspace(0, 7500, 200)
    probe_position = 0
    plt.subplot(2,1,1)
    complex_spin_profile = vanilla_experiment.get_complex_spin_profile(
        tvals, probe_position, **experiment_state_kwargs)
    plt.plot(tvals, np.real(complex_spin_profile), 'd')
    plt.title('Direct spin profile')
    plt.subplot(2,1,2)
    probe_profile = vanilla_experiment.get_probe_profile(
        tvals, probe_position, **laser_kwargs, **experiment_state_kwargs)
    plt.plot(tvals, probe_profile, 'd')     
    plt.title('Spin profile w/ probe convolution and Kerr physics')


# %% RSA TEST
if False:
# %%
    applied_E_field = 15  # V/cm
    external_B_fields = np.linspace(-0.1, 0.1, 1000)    
    experiment_field_dict = {'applied_E_fields': applied_E_field,
                             'external_B_fields': external_B_fields}

    delay_time = LASER_REPRATE - 40
    probe_position = 0
    plt.subplot(2,1,1)
    complex_spin_profile = vanilla_experiment.get_complex_spin_profile(
        delay_time, probe_position, **experiment_field_dict)
    plt.plot(external_B_fields, np.real(complex_spin_profile), 'd')
    plt.xlim([min(external_B_fields), max(external_B_fields)])
    plt.title('Direct spin profile')
    plt.subplot(2,1,2)
    probe_profile = vanilla_experiment.get_probe_profile(
        delay_time, probe_position, **laser_kwargs, **experiment_field_dict)
    plt.plot(external_B_fields, probe_profile, 'd')
    plt.xlim([min(external_B_fields), max(external_B_fields)])
    plt.title('Spin profile w/ probe convolution and Kerr physics')


# %% EXPORT RSA DATA WITH NOISE
if False:
# %%
    delay = np.array(-100)  # ps
    probe_position = 0  # um
    probe_position_list = np.arange(-20, 101, 2)
    applied_E_fields = [15]  # V/cm
    external_B_fields = np.linspace(-40, 40, 201)  # mT
    num_fake_runs = 2
    for run_ind in range(num_fake_runs):
        for probe_position in probe_position_list:
            for applied_E_field in applied_E_fields:
                experiment_field_dict = {'applied_E_fields': applied_E_field,
                                         'external_B_fields': external_B_fields}
                probe_profile = vanilla_experiment.get_probe_profile(
                    delay, probe_position,
                    **laser_kwargs, **experiment_field_dict)
                probe_noise = 0.0005 * np.random.randn(external_B_fields.size)
                probe_offset = .000 * (2 * (0.5 - np.random.random()))
                probe_slope = .000 * np.random.randn()
                probe_linear_offset = \
                    probe_offset + probe_slope * \
                        np.linspace(-0.5, 0.5, external_B_fields.size)
                probe_profile += probe_noise + probe_linear_offset

#                plt.plot(external_B_fields, probe_profile, '-d')
#                plt.plot(external_B_fields,
#                         probe_profile + probe_noise + probe_linear_offset, 'r')
#                plt.title('Spin profile w/ probe convolution and Kerr physics')
#                raise KeyboardInterrupt

                directory = "C:\\Data\\fake_data\\fake_rsa"
                filename = "RSA_MirrorZ_{:.0f}_".format(probe_position) + \
                           "{:.0f}Vcm_".format(float(applied_E_field)) + \
                           "{:.0f}ps_".format(float(delay)) + \
                           "30K_818.9nm_run{:03d}.dat".format(run_ind + 1)
                filepath = directory + "\\" + filename
    #            headerstr = "\n".join(10*["header header header"])
                headerstr = '\t'.join(["scancoord","lockin2x"])
                np.savetxt(filepath,
                           np.vstack([external_B_fields, probe_profile]).T,
                           fmt='%.6e\t%.12e', delimiter='\t', newline='\n',
                           header=headerstr, comments='')

#    from scipy.interpolate import UnivariateSpline
#    drift_spline_n_pts = 8
#    drift_spline_padding = 8
#    drift_spline_scale = 0.005 * np.random.randn()
#    x = np.linspace(min(tvals), max(tvals), drift_spline_n_pts)
#    y = np.random.randn(drift_spline_n_pts)
#    x_padded = np.linspace(min(tvals), max(tvals),
#                           drift_spline_n_pts + 2 * drift_spline_padding)
#    y_padded = np.hstack([y[0] * np.ones(drift_spline_padding),
#                          y,
#                          y[-1] * np.ones(drift_spline_padding)])
#    s = UnivariateSpline(x, y, s=None)
#    s_padded = UnivariateSpline(x_padded, y_padded, s=None)
##    plt.plot(x, y, 'bd')
##    plt.plot(x_padded, y_padded, 'gd')
#    plt.plot(tvals, s(tvals), 'b')
#    plt.plot(tvals, s_padded(tvals), 'g')
#    probe_drift = drift_spline_scale * s(tvals)



