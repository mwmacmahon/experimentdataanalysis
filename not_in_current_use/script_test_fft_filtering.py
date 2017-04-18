# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 20:09:50 2016

@author: mwmac
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.filters as filters

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata, get_high_pass_filtered_scandata
import experimentdataanalysis.parsing.dataseriesparsing as dsparsing


# %%
# beating TRKR
dirpath="C:\\Data\\august_data\\160902\\GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime"
scandatalist = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(dirpath))
#dataseries = scandatalist[13].dataseries_list[1]  # perfect decaying sine
dataseries = scandatalist[12].dataseries_list[1]
dataseries_min_freq_cutoff = 0.0005

#==============================================================================
# # RSA
# dirpath="C:\\Data\\august_data\\160902\\WavelengthDependence_RSA"
# scandatalist = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(dirpath))
# dataseries = scandatalist[1].dataseries_list[1]
# dataseries_min_freq_cutoff = 20
# 
#==============================================================================


#xvals, oldyvals = dataseries.data(raw=True)
#indices_under_zero = xvals < -10
#xvals[indices_under_zero] += 13160

xvals = np.linspace(-500,7000,225)
#xvals = np.linspace(-1000,14000,550)
#xvals = np.hstack([np.linspace(0,7000,210), np.linspace(12660, 13160, 15)])
#xvals = np.hstack([np.linspace(0,7000,210), np.linspace(12660, 13160, 15),])
#xvals = np.linspace(0,7000,210)
oldyvals = np.exp(-xvals/20000)*np.cos(2*np.pi*xvals/530) + \
           np.exp(-xvals/20000)*np.cos(2*np.pi*xvals/570)

# include several copies of pulse
numpulses = 4
xvals = np.hstack([xvals + 13160*pulsenum
                   for pulsenum in range(numpulses)])
oldyvals = np.hstack([oldyvals
                      for pulsenum in range(numpulses)])


if len(xvals) % 2 > 0:  # ensure even # data pts
    xvals = xvals[1:]
    oldyvals = oldyvals[1:]
inverse_sample_rate = xvals[1] - xvals[0]
f_space_freqs = np.fft.rfftfreq(oldyvals.shape[-1], d=inverse_sample_rate)
f_space_oldyvals = np.fft.rfft(oldyvals)
f_space_newyvals = f_space_oldyvals.copy()
f_space_newyvals[f_space_freqs <= dataseries_min_freq_cutoff] = 0
newyvals = np.fft.irfft(f_space_newyvals)

#plt.figure()
plt.hold(True)

axes1 = plt.subplot(3,1,1)
axes1.plot(xvals, oldyvals, 'b')
axes1.plot(xvals, newyvals, 'r')
axes1.set_xlabel("Delay (ps)")
axes1.set_ylabel("Signal")
axes1.set_autoscalex_on(False)
axes1.set_xlim([min(xvals), max(xvals)])

# FFT plot vs freq: 
axes2 = plt.subplot(3,1,3)
axes2.plot(f_space_freqs[1:], abs(f_space_oldyvals[1:]), 'b')
axes2.plot(f_space_freqs[1:], abs(f_space_newyvals[1:]), 'r')
axes2.set_xlabel("Frequency (1/ps)")
axes2.set_ylabel("abs( FFT(Signal) )")
axes2.set_autoscalex_on(False)
axes2.set_xlim([min(f_space_freqs[1:]), max(f_space_freqs[1:])])

# FFT plot vs period: 
axes3 = plt.subplot(3,1,2)
f_space_periods = 1/f_space_freqs[1:]
f_space_oldys = abs(f_space_oldyvals[1:])
f_space_newys = abs(f_space_newyvals[1:])
axes3.semilogx(f_space_periods, f_space_oldys, 'b')
axes3.semilogx(f_space_periods, f_space_newys, 'r')
axes3.set_xlabel("1/Frequency (ps)")
axes3.set_ylabel("abs( FFT(Signal) )")
axes3.set_autoscalex_on(False)
axes3.set_xlim([min(f_space_periods), max(f_space_periods)])

plt.tight_layout()

# %%

dirpath = ("C:\\Data\\august_data\\160902\\" +
           "WavelengthDependence_RSA")
#dirpath="C:\\Data\\august_data\\160902\\GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime"
scandatalist = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(dirpath))
#dataseries = scandatalist[13].dataseries_list[1]  # perfect decaying sine
dataseries = scandatalist[1].dataseries_list[1]
x_scale_factor = 5e-4
plt.figure()
plt.hold(True)
data = dataseries.yvals()
xvals = [x_scale_factor*index for index in range(len(data))]
plt.plot(xvals, data, 'wd', label='original data')

sigma_list = [0, 0.006]
filtered_data = []
for sigma in sigma_list:
    fdata = filters.gaussian_filter1d(data,
                                      sigma=sigma/x_scale_factor,
                                      mode='reflect')
    plt.plot(xvals, fdata, linewidth=2,
             label='gaussian smoothed, sigma= {}'.format(sigma))
    filtered_data.append(fdata)
plt.legend()

# %%
smoothed_data = filtered_data[0]
trend_data = filtered_data[1]

plt.plot(smoothed_data, 'wd')
plt.plot(smoothed_data - filtered_data[1])

