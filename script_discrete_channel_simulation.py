# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 17:53:06 2016

@author: mwmac
"""

import matplotlib.pyplot as plt
import numpy as np

GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period


def toy_kerr_rotation_curve(x_to_test):
    # define curve
    x = np.linspace(-2, 2, 100)
    x1 = x[np.argwhere(x < -0.2)]
    x2 = x[np.argwhere(np.logical_and(x >= -0.2, x < 0.2))]
    x3 = x[np.argwhere(x >= 0.2)]    
    y1 = np.clip(-1.5/(abs(x1) + 0.12), -3, 0)
    y2 = np.interp(x2, (-0.2, 0.2), (-3, 1))
    y3 = np.clip(1.1/(abs(x3) + 0.8), 0, 1)
    y = np.concatenate([y1, y2, y3])[:,0]
#    plt.plot(x, y, 'd')
    return np.interp(x_to_test, x, y)


def generate_absorption_operator(band_edge, relative_phase):
    """
    Generates a function that takes a spin polarization & orientation profile,
    a laser intensity profile, and that laser's wavelength, then returns
    the new polarization profile using all the fun selection rules kramers
    kronig blah blah rules and stuff ideally. Takes as parameters the
    band edge of the target spin species, and any other parameters needed
    to inform how the absorption should play out.
    """
    def absorption_operator(complex_spin_profile,
                            laser_intensity_profile, laser_wavelength):
        # TODO: insert complexities of adding the new pulse!
        complex_spin_profile += laser_intensity_profile * \
                                    np.exp(1j * relative_phase)
#        if laser_wavelength > band_edge:
#            complex_spin_profile += laser_intensity_profile  # on-axis = real
#        else:
#            complex_spin_profile -= laser_intensity_profile  # on-axis = real

        return complex_spin_profile
    return absorption_operator


def generate_kerr_rotation_operator(band_edge, relative_phase):
    """
    Generates a function that takes a spin polarization & orientation profile,
    a laser intensity profile, and that laser's wavelength, then returns
    the expected kerr rotation signal returned using all the fancy kramers
    kronig beep boop stuff ideally. Takes as parameters the band edge
    of the target spin species, and any other parameters needed
    to inform how the kerr rotation should play out.
    """
    def kerr_rotation_operator(complex_spin_profile,
                               laser_intensity_profile, laser_wavelength):
        # TODO: fancy wavelength and band_edge dependent stuff, even
        # probe-wavelength dependent stuff for a real show off
#        delta_wavelength = laser_wavelength - band_edge
#        complex_spin_profile *= toy_kerr_rotation_curve(delta_wavelength)
#        complex_spin_profile *= np.exp(1j * -1 * relative_phase)
        measured_profile = np.real(complex_spin_profile)
        probe_convolution = laser_intensity_profile * measured_profile
        signal = sum(probe_convolution)
        return signal
    return kerr_rotation_operator


class SpinProfile():
    def __init__(self, xmin, xmax, nboxes,
                 carrier_density, spin_lifetime, gfactor, mobility,
                 transfer_decay_const, diffusion_const, **operator_args):
        self.xvals = np.linspace(xmin, xmax, nboxes)  # pos
        # assuming constant carrier density, variable polarization & angle.
        # but might be different between two spin species
        self.carrier_density = carrier_density  # assume const vs pos, time
        self.polarization = np.zeros_like(self.xvals)  # from 0 to 1
        self.orientation = np.zeros_like(self.xvals)  # angle in radians
        self.spin_lifetime = spin_lifetime
        self.gfactor = gfactor
        self.mobility = mobility
        self.transfer_decay_const = transfer_decay_const
        self.diffusion_const = diffusion_const
        self.absorption_operator = generate_absorption_operator(**operator_args)
        self.kerr_rotation_operator = generate_kerr_rotation_operator(**operator_args)
        # for spatial steps due to drift:
        self.x_grid_spacing = (xmax - xmin) / nboxes
        self.fractional_spatial_step = 0.0

    def get_complex_spin_profile(self):
        return self.carrier_density * self.polarization * \
                                        np.exp(1j * self.orientation)

    def set_complex_spin_profile(self, complex_profile):
        self.polarization = np.abs(complex_profile) / self.carrier_density
        # restrict orientation to phases -pi/2 to 3*pi/2 to avoid cutoff at pi
        self.orientation = np.angle(complex_profile)
        self.orientation = ((self.orientation + np.pi/2) % (2*np.pi)) - np.pi/2

    def generate_packet(self, laser_wavelength,
                        laser_power, laser_width, laser_position):
        laser_intensity_profile = laser_power * \
            np.exp(-(self.xvals - laser_position)**2 / (2*laser_width**2))
        new_complex_spin_profile = \
            self.absorption_operator(self.get_complex_spin_profile(),
                                     laser_intensity_profile, laser_wavelength)
        self.set_complex_spin_profile(new_complex_spin_profile)

    def probe_at_position(self, laser_wavelength,
                          laser_power, laser_width, laser_position):
        laser_intensity_profile = laser_power * \
            np.exp(-(self.xvals - laser_position)**2 / (2*laser_width**2))
        return self.kerr_rotation_operator(self.get_complex_spin_profile(),
                                           laser_intensity_profile,
                                           laser_wavelength)

    def time_evolve(self, time_step, e_field, external_b_field):
        # rotation and decay of polarization
        b_field = external_b_field  # no DNP or SOC atm
        self.polarization *= np.exp(-time_step / self.spin_lifetime)
        self.orientation += (2 * np.pi * self.gfactor * b_field *
                             time_step * GFACTORCONSTANT)

        # restrict oorientation to phases -pi/2 to 3*pi/2 to avoid cutoff at pi
        self.orientation = ((self.orientation + np.pi/2) % (2*np.pi)) - np.pi/2

        # calculate profile of spins transferred into other spin species
        untransferred_ratio = np.exp(-self.transfer_decay_const * time_step)
        transferred_complex_spin_profile = \
            (1.0 - untransferred_ratio) * self.get_complex_spin_profile()
        self.polarization *= untransferred_ratio

        # movement of spins due to electric field if enough time passed
        # remember grid is "boxes" of spins w/ polarization & orientation
        drift_velocity = self.mobility * e_field
        distance_travelled = drift_velocity * time_step
        spatial_step_fraction = distance_travelled / self.x_grid_spacing
        self.fractional_spatial_step += spatial_step_fraction
        if abs(self.fractional_spatial_step) > 1.0:
            new_polarization = np.zeros_like(self.polarization)
            num_steps = int(abs(self.fractional_spatial_step))
            if self.fractional_spatial_step > 0:
                new_polarization[num_steps:] = self.polarization[:-num_steps]
                new_orientation = np.full_like(self.orientation,
                                               self.orientation[0])
                new_orientation[num_steps:] = self.orientation[:-num_steps]
                self.fractional_spatial_step -= num_steps  # positive shift
            else:
                new_polarization[:-num_steps] = self.polarization[num_steps:]
                new_orientation = np.full_like(self.orientation,
                                               self.orientation[-1])
                new_orientation[:-num_steps] = self.orientation[num_steps:]
                self.fractional_spatial_step += num_steps  # negative shift
            self.polarization = new_polarization
            self.orientation = new_orientation

        return transferred_complex_spin_profile



# %%
if __name__ == "__main__":
    # TO KEEP FLOW SIMPLE, FLOW WILL ONLY OCCUR IN INTEGER SPATIAL STEPS
    # THUS HIGH TIME RESOLUTION NECESSITATES HIGH SPATIAL RESOLUTION!

    # general parameters
    xmin, xmax = (-100, 250)  # um
    nboxes = 350
    num_pulses = 10
    time_step = 10.0  # ps

    # field parameters
    electric_field = 15  # V/cm
    external_b_field = 0.3

    # recording parameters
    surpress_plot = False
    if not surpress_plot:
        plt.figure()
        pos_graph_axes = plt.subplot(3,1,1)
        trkr_graph_axes_list = [plt.subplot(3,2,3),
                                plt.subplot(3,2,4),
                                plt.subplot(3,2,5),
                                plt.subplot(3,2,6)]
    probe_position_list = [-10, 0, 20, 30]    
    measurement_window = [-500, 7000]
    time_step = 5.0  # ps
    streak_start = -500
    streak_end = 7000
    streak_duration = streak_end - streak_start
    streak_line_count = 100  # lines in streak, give or take one, limited by time step
    streak_spacing = streak_duration // streak_line_count
    steps_per_pulse = int(round(LASER_REPRATE/time_step))
    lab_time = -(num_pulses - 1)*LASER_REPRATE  # lab_time=0 is last pulse
    num_time_indices = 1 + \
        int(np.ceil((measurement_window[1]-measurement_window[0])/time_step))
    trkr_vs_pos_list = [np.zeros(num_time_indices)
                        for num in range(len(probe_position_list))]
    timeseries_current_index = 0

    # laser parameters
    pump_laser_kwargs = dict(laser_wavelength=818.9,  # nm
                             laser_power=1.0,  # AU
                             laser_width=20.0,  # um
                             laser_position=0.0)  # um
    probe_laser_kwargs = dict(laser_wavelength=818.9,  # nm
                              laser_power=1.0,  # AU
                              laser_width=20.0)  # um
                              # add position manually!

    # create species profiles, will evolve alongside each other
    species1 = SpinProfile(xmin, xmax, nboxes,
                           carrier_density=1.0,
                           spin_lifetime=20000.0,  # ps
                           gfactor=.44,
                           mobility=1e-4,  # (um/ps)/(V/cm). 1e-4 -> 2um/ns/Vapp
                           transfer_decay_const=1e-8,
                           diffusion_const=0,
                           band_edge=818.4,
                           relative_phase=0*np.pi)
    species2 = SpinProfile(xmin, xmax, nboxes,
                           carrier_density=2.0,
                           spin_lifetime=2000.0,  # ps
                           gfactor=.436,
                           mobility=1e-4,  # (um/ps)/(V/cm). 1e-4 -> 2um/ns/Vapp
                           transfer_decay_const=1e-8,
                           diffusion_const=0,
                           band_edge=819.4,
                           relative_phase=-2*np.pi/3)

    def time_evolve_system():
        species1_transferred_complex_spin_profile = \
            species1.time_evolve(time_step, electric_field, external_b_field)
        species2_transferred_complex_spin_profile = \
            species2.time_evolve(time_step, electric_field, external_b_field)
        species1.set_complex_spin_profile(
            species1.get_complex_spin_profile() +
            species2_transferred_complex_spin_profile)
        species2.set_complex_spin_profile(
            species2.get_complex_spin_profile() +
            species1_transferred_complex_spin_profile)
        
    for i in range(num_pulses):
        species1.generate_packet(**pump_laser_kwargs)
        species2.generate_packet(**pump_laser_kwargs)
        for j in range(steps_per_pulse):  # time steps per pulse
            time_evolve_system()
            # compute and record channel state if necessary
            lab_time += time_step
            streak_time = lab_time - streak_start
            record_trkr = all([lab_time >= measurement_window[0],
                               lab_time <= measurement_window[1]])
            record_streak = all([not surpress_plot,
                                 streak_time >= 0,
                                 streak_time < streak_duration,
                                 streak_time % streak_spacing < time_step])
            if record_trkr or record_streak:
                profile1 = species1.polarization * np.cos(species1.orientation)
                profile2 = species2.polarization * np.cos(species2.orientation)
#                profile1 = np.array(
#                    [species1.probe_at_position(**probe_laser_kwargs,
#                                                laser_position=x)
#                     for x in species1.xvals])
#                profile2 = np.array(
#                    [species2.probe_at_position(**probe_laser_kwargs,
#                                                laser_position=x)
#                     for x in species1.xvals])
                sum_profile = profile1 + profile2
            if record_streak and not surpress_plot:
                brightness = 1 - (streak_time / streak_duration)**4.0
                if brightness > 1 or brightness < 0:
                    raise Exception("error: invalid brightness value")
                pos_graph_axes.plot(species1.xvals, sum_profile,
                                    color=(3 * tuple([brightness])))
            if record_trkr:
                for list_ind, pos in enumerate(probe_position_list):
#                    pos_index = np.argwhere(species1.xvals > pos)[0]
#                    series_val = sum_profile[pos_index]
                    series_val = \
                        species1.probe_at_position(**probe_laser_kwargs,
                                                   laser_position=pos) + \
                        species2.probe_at_position(**probe_laser_kwargs,
                                                   laser_position=pos)
                    series = trkr_vs_pos_list[list_ind]
                    if timeseries_current_index < len(series):
                        series[timeseries_current_index] = series_val
                    else:
                        print('oops, t={}, i={}'.format(lab_time,
                                                        timeseries_current_index))
                timeseries_current_index += 1

    # Plot results
    if not surpress_plot:
        times = np.linspace(measurement_window[0],
                            measurement_window[1], num_time_indices)
        posymin, posymax = pos_graph_axes.yaxis.get_view_interval()
        for posnum in range(4):
            axes = trkr_graph_axes_list[posnum]
            axes.plot(times, trkr_vs_pos_list[posnum])
            axes.text(0.5, 0.9, "TRKR @x={}".format(probe_position_list[posnum]),
                      horizontalalignment='center', verticalalignment='center',
                      fontsize=18, transform=axes.transAxes)
            axes.set_yticklabels([])
            if posnum < 2:
                axes.set_xticklabels([])
            pos_graph_axes.plot(2*[probe_position_list[posnum]],
                                [posymin, posymax - .001], 'r')           

#    plt.subplot(2,1,2)
#    plt.plot(1 + species1.polarization, 'r:')
#    plt.plot(1 + species1.orientation / np.pi, 'b:')
#    plt.plot(2 + species2.polarization, 'r:')
#    plt.plot(2 + species2.orientation / np.pi, 'b:')
#    plt.plot(species1.polarization * np.cos(species1.orientation) +
#             species2.polarization * np.cos(species2.orientation), 'gd')


