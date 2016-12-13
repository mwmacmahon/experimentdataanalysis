# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 13:57:37 2016

Designed to run simple toy models of drift-diffusion spin systems
with one or more electron bands.

model: channel is 1D and represented by a gigantic 1D vector. In practice,
we use numpy arrays to access each slice of the giant 1D vector, each slice
corresponding to the population of a given carrier type along each position
x_i along the channel.

@author: mwmac
"""

import matplotlib.pyplot as plt
import numpy as np

GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period

class carrier_packet():
    def __init__(self, identity, laser_wavelength, laser_power,
                 laser_width, initial_orientation, initial_polarization,
                 initial_position, spin_lifetime, gfactor, mobility,
                 diffusion_const, channel_xvals):
        self.identity = identity
        # TODO: actually roll this into polarization, no separate carrier density
        self.carrier_density_gaussian_amplitude = \
            laser_power * self.get_absorption_coeff(laser_wavelength)
        self.carrier_density_gaussian_width = laser_width
        self.orientation = initial_orientation
        self.polarization = initial_polarization  # magnitude, from 0 to 1
        self.position = initial_position
        self.spin_lifetime = spin_lifetime
        self.gfactor = gfactor
        self.mobility = mobility
        self.diffusion_const = diffusion_const
        self.channel_xvals = channel_xvals

        self.update_state(time_step=0,  # just to initialize output attributes
                          electric_field=0,
                          external_b_field=0)

    def get_probe_profile(self, laser_width, laser_power, laser_wavelength):
        probe_convolution = np.exp(-(self.channel_xvals - self.position)**2 /
                                   (2*(laser_width**2 +
                                       self.carrier_density_gaussian_width**2)))
        signal_profile = (self.carrier_density_gaussian_amplitude *
                          self.polarization * np.exp(1j * self.orientation) *
                          probe_convolution)

        # TODO: do this right
        absorption_coeff = self.get_absorption_coeff(laser_wavelength)
        trkr_coeff = absorption_coeff
        return trkr_coeff * signal_profile

    def get_absorption_coeff(self, laser_wavelength):
        if self.identity == "low bandgap spins":
            return 1.0
        elif self.identity == "high bandgap spins":
            return 1.0
        else:
            raise Exception("invalid identity for spin packet '{}'".format(
                                                                self.identity))

    def update_state(self, time_step, **changes):  # time step in ps
        for key, val in changes.items():
            # external parameters, should be refreshed each update
            if key == "electric_field": self.electric_field = val  # V/cm
            if key == "external_b_field": self.external_b_field = val  # Tesla

        # calculated parameters, must be recalculated each update
        self.drift_velocity = self.mobility * self.electric_field
        self.magnetic_field = self.external_b_field  # no DNP or SOC atm

        # time-evolve other parameters according to time-step
        self.polarization = (self.polarization *
                             np.exp(-time_step / self.spin_lifetime))
        self.orientation += (2 * np.pi * GFACTORCONSTANT * time_step *
                              self.magnetic_field * self.gfactor)
        if self.orientation > np.pi:
            self.orientation -= 2 * np.pi
        self.position += self.drift_velocity * time_step
        self.carrier_density_gaussian_width += self.diffusion_const * time_step



#        # handling decays from any source
#        survival_chance_list = np.exp(-time_step / self.spin_lifetime)
#        net_decay_chance = 1.0 - np.product(survival_chance_list)
#        num_decays = min(self.num_carriers,
#                         int(net_decay_chance * self.num_carriers))

#        # scattering 
##        decaying_indices = np.random.permutation(self.num_carriers)[num_decays:]
#        live_indices = np.random.permutation(self.num_carriers)[num_decays:]
#        self.orientation = self.orientation[live_indices]
#        self.position = self.position[live_indices]
#        self.drift_velocity = self.drift_velocity[live_indices]
#        self.num_carriers = len(self.orientation)

#        # finally, recalculate new net polarization
#        self.net_polarization_magnitude = \
#            np.abs(sum(np.exp(1j*self.orientation)) / self.num_carriers)
#        self.net_polarization_angle = \
#            np.angle(sum(np.exp(1j*self.orientation)))

    def __str__(self):
        s = ("packet contains {} {} spins ".format(self.num_carriers,
                                                   self.identity) +
             "of net angle {} and polarization {}".format(
                 self.net_polarization_angle,
                 self.net_polarization_magnitude))
        return s


# %%
# TODO: RSA SIMULATION, make this into a function of B field
# lose ability to use as workspace if not name == main though
def find_TRKR_vs_position(electric_field,  # V/cm
                          external_b_field,  # Tesla
                          probe_laser_wavelength,  # nm
                          pump_laser_wavelength,
                          probe_laser_power,  # AU
                          pump_laser_power,
                          probe_laser_width,  # um
                          pump_laser_width,
                          packet1kwargs,
                          packet2kwargs,
                          suppress_plot=False):
    def get_new_pulse_packets():
        packet1 = carrier_packet(**packet1kwargs)
        packet2 = carrier_packet(**packet2kwargs)
        packet_list = [packet1, packet2]
        return packet_list
    
    plt.figure()
    pos_graph_axes = plt.subplot(3,1,1)
    trkr_graph_axes_list = [plt.subplot(3,2,3),
                            plt.subplot(3,2,4),
                            plt.subplot(3,2,5),
                            plt.subplot(3,2,6)]
    probe_position_list = [-10, 0, 20, 40]

    num_pulses = 15
    measurement_window = [-500, 7000]
    time_step = 5.0  # ps
    streak_start = 6000
    streak_end = 7000
    streak_duration = streak_end - streak_start
    streak_line_count = 100  # lines in streak, give or take one, limited by time step
    streak_spacing = streak_duration // streak_line_count
    steps_per_pulse = int(round(LASER_REPRATE/time_step))
    lab_time = -(num_pulses - 1)*LASER_REPRATE  # lab_time=0 is last pulse
    packet_list = []
    num_time_indices = 1 + \
        int(np.ceil((measurement_window[1]-measurement_window[0])/time_step))
    trkr_vs_pos_list = [np.zeros(num_time_indices)
                        for num in range(len(probe_position_list))]
    timeseries_current_index = 0
    for i in range(num_pulses):  # num laser pulses
        packet_list.extend(get_new_pulse_packets())
        for j in range(steps_per_pulse):  # time steps per pulse
            lab_time += time_step
            for packet in packet_list:
                packet.update_state(time_step,
                                    electric_field=electric_field,
                                    external_b_field=external_b_field)
            # compute and record channel state if necessary
            streak_time = lab_time - streak_start
            record_trkr = all([lab_time >= measurement_window[0],
                               lab_time <= measurement_window[1]])
            record_streak = all([not suppress_plot,
                                 streak_time >= 0,
                                 streak_time < streak_duration,
                                 streak_time % streak_spacing < time_step])
            if record_trkr or record_streak:
                packet_profile_list = \
                    [packet.get_probe_profile(probe_laser_width,
                                              probe_laser_power,
                                              probe_laser_wavelength)
                     for packet in packet_list]
                vector_sum_series = np.real(np.sum(packet_profile_list,
                                                   axis=0))
#                scalar_sum_series = np.sum([np.real(series)
#                                         for series in packet_profile_list],
#                                           axis=0)
            if record_streak and not suppress_plot:
                axes = pos_graph_axes
                brightness = 1 - (streak_time / streak_duration)**4.0
                if brightness > 1 or brightness < 0:
                    raise Exception("error: invalid brightness value")
                axes.plot(channel_xvals, vector_sum_series,
                          color=(3 * tuple([brightness])))
            if record_trkr:
                for list_ind, pos in enumerate(probe_position_list):
                    pos_index = np.argwhere(channel_xvals > pos)[0]
                    series_val = vector_sum_series[pos_index]
                    series = trkr_vs_pos_list[list_ind]
                    if timeseries_current_index < len(series):
                        series[timeseries_current_index] = series_val
                    else:
                        print('oops, t={}, i={}'.format(lab_time,
                                                        timeseries_current_index))
                timeseries_current_index += 1

    # Plot results
    if not suppress_plot:
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

    return packet_list, trkr_vs_pos_list







# %%
if __name__ == "__main__":
    channel_xvals = np.linspace(-50, 250, 200)
    electric_field = 15  # V/cm
    external_b_fields = [0.3]
#    external_b_fields = [0.15, 0.225, 0.3]
#    external_b_fields = [0.1, 0.125, 0.15, 0.175,
#                         0.2, 0.225, 0.25, 0.275, 0.3]
#    external_b_field = 0.3  # Tesla
#    external_b_field = 0.6  # Tesla

    probe_laser_wavelength = 818.9  # nm
    pump_laser_wavelength = 818.9
    probe_laser_power = 1.0  # AU
    pump_laser_power = 1.0
    probe_laser_width = 20.0  # um
    pump_laser_width = 20.0

    packet1kwargs = dict(identity="low bandgap spins",
                         laser_wavelength=pump_laser_wavelength,
                         laser_power=pump_laser_power,
                         laser_width=pump_laser_width,
                         initial_orientation=0*np.pi/3,  # angle, -pi/2 to pi/2
                         initial_polarization=0.25,
                         initial_position=0.0,
                         spin_lifetime=20000.0,
                         gfactor=.44,
                         mobility=1e-4,  # (um/ps)/(V/cm). 1e-4 -> 2um/ns/Vapp
                         diffusion_const=0,
                         channel_xvals=channel_xvals)

    packet2kwargs = dict(identity="high bandgap spins",
                         laser_wavelength=pump_laser_wavelength,
                         laser_power=pump_laser_power,
                         laser_width=pump_laser_width,
                         initial_orientation=-1.5*np.pi/3,  # angle, pi/2 to 3*pi/2
                         initial_polarization=0.5,
                         initial_position=0.0,
                         spin_lifetime=4000.0,
                         gfactor=.43,
                         mobility=1e-4,  # (um/ps)/(V/cm). 1e-4 -> 2um/ns/Vapp
                         diffusion_const=0,
                         channel_xvals=channel_xvals)

    packet_lists = []
    trkr_vs_pos_lists = []
    for external_b_field in external_b_fields:
        new_packet_list, new_trkr_vs_pos_list = \
            find_TRKR_vs_position(electric_field,  # V/cm
                                  external_b_field,  # Tesla
                                  probe_laser_wavelength,  # nm
                                  pump_laser_wavelength,
                                  probe_laser_power,  # AU
                                  pump_laser_power,
                                  probe_laser_width,  # um
                                  pump_laser_width,
                                  packet1kwargs,
                                  packet2kwargs,
                                  suppress_plot=False)


# %%
#if __name__ == "__main__":
#    channel_xvals = np.linspace(-50, 250, 200)
#    electric_field = 20  # V/cm
#    external_b_fields = [0.1, 0.15, 0.2, 0.25, 0.3]
##    external_b_fields = [0.1, 0.125, 0.15, 0.175,
##                         0.2, 0.225, 0.25, 0.275, 0.3]
##    external_b_field = 0.3  # Tesla
##    external_b_field = 0.6  # Tesla
#
#    probe_laser_wavelength = 818.9  # nm
#    pump_laser_wavelength = 818.9
#    probe_laser_power = 1.0  # AU
#    pump_laser_power = 1.0
#    probe_laser_width = 20.0  # um
#    pump_laser_width = 20.0
#
#    packet1kwargs = dict(identity="low bandgap spins",
#                         laser_wavelength=pump_laser_wavelength,
#                         laser_power=pump_laser_power,
#                         laser_width=pump_laser_width,
#                         initial_orientation=0*np.pi/3,  # angle, -pi to pi
#                         initial_polarization=0.5,
#                         initial_position=0.0,
#                         spin_lifetime=20000.0,
#                         gfactor=.44,
#                         mobility=1e-4,  # (um/ps)/(V/cm). 1e-4 -> 2um/ns/Vapp
#                         diffusion_const=0,
#                         channel_xvals=channel_xvals)
#
#    packet2kwargs = dict(identity="high bandgap spins",
#                         laser_wavelength=pump_laser_wavelength,
#                         laser_power=pump_laser_power,
#                         laser_width=pump_laser_width,
#                         initial_orientation=-2*np.pi/3,  # angle, -2pi to pi
#                         initial_polarization=0.5,
#                         initial_position=0.0,
#                         spin_lifetime=20000.0,
#                         gfactor=.43,
#                         mobility=1e-4,  # (um/ps)/(V/cm). 1e-4 -> 2um/ns/Vapp
#                         diffusion_const=0,
#                         channel_xvals=channel_xvals)
#
#    packet_lists = []
#    trkr_vs_pos_lists = []
#    for external_b_field in external_b_fields:
#        new_packet_list, new_trkr_vs_pos_list = \
#            find_TRKR_vs_position(electric_field,  # V/cm
#                                  external_b_field,  # Tesla
#                                  probe_laser_wavelength,  # nm
#                                  pump_laser_wavelength,
#                                  probe_laser_power,  # AU
#                                  pump_laser_power,
#                                  probe_laser_width,  # um
#                                  pump_laser_width,
#                                  packet1kwargs,
#                                  packet2kwargs,
#                                  suppress_plot=False)
#
#
