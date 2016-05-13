# -*- coding: utf-8 -*-
"""
Analyze experiments of form:
1D scan: delay scan
2D scan: mirror (or nothing):
3D scan: voltage

might change to:
1D scan: mirror
2D scan: delay scan:
3D scan: voltage
"""

import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
import experimentdataanalysis.analysis.dataclassfitting as dcfitting
import experimentdataanalysis.analysis.dataclassgraphing as dcgraphing
import experimentdataanalysis.analysis.dataseriesprocessing as dsprocessing
from experimentdataanalysis.analysis.fitfunctions \
    import fitfcn_simple_1d_gaussian, fitfcn_single_exp_decay, \
        fitfcn_two_exp_decay
from experimentdataanalysis.analysis.multidataseriesprocessing \
    import dataseries_iterable_fit, scandata_iterable_fit
import experimentdataanalysis.parsing.dataclassparsing as dcparsing

# %%
if __name__ == "__main__":
# %%
    scandata_list = list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
#        directorypath="C:\\Data\\160506\\Voltage_0_OverlapVsDelay_Channel3_1mT_033XT-B11_819.5nm_30K_2Dscan_DelayTime_MirrorZ"))
        directorypath="C:\\Data\\160506\\Voltage_4_OverlapVsDelay_Channel3_1mT_033XT-B11_819.5nm_30K_2Dscan_DelayTime_MirrorZ"))
#        directorypath="C:\\Data\\160506\\Voltage_8_OverlapVsDelay_Channel3_1mT_033XT-B11_819.5nm_30K_2Dscan_DelayTime_MirrorZ"))
#    scandata_list += list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelCenter_200mT_Channel3_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime"))

#    dataseries_list = [scandata.dataseries[0] for scandata in scandata_list]


# %%
    #fit gaussians to data
    #params = amplitude, t0, sigma, offset
    key_field_ind = 0
    free_params = [True, True, True, False]
    initial_params = [0.001, 0, 20, 0]
    param_bounds = [(-0.1, 0.1), (-100, 100), (10, 60), (0, 0)]
    fit_scandata_list = scandata_iterable_fit(
                            scandata_list, key_field_ind,
                            fitfcn_simple_1d_gaussian,
                            free_params, initial_params, param_bounds,
                            multiprocessing=False)

    #create dataseries of fit parameters & their uncertainties:
    xvalkey = "MiddleScanCoord"
    xvals = [scandata.scaninfo[xvalkey] for scandata in fit_scandata_list]
    fitparams = list(zip(*(scandata.fitdata[key_field_ind].fitparams
                           for scandata in fit_scandata_list)))
    fitparams_dataseries_list = [DataSeries(zip(xvals, fitparam))
                                 for fitparam in fitparams]
    fitparamstds = list(zip(*(scandata.fitdata[key_field_ind].fitparamstds
                              for scandata in fit_scandata_list)))
    fitparamstds_dataseries_list = [DataSeries(zip(xvals, fitparamstd))
                                    for fitparamstd in fitparamstds]
# screwing around with errors for absolute_sigma testing
#    fitparamstds_dataseries_list = [DataSeries(zip(xvals, [val for val in fitparamstd[:-1]] + [fitparamstd[-1]]))
#                                    for fitparamstd in fitparamstds]

    #pick out one particular fit parameter dataseries for further analysis
    key_fit_param_ind = 0
    raw_fit_param_dataseries = \
                        fitparams_dataseries_list[key_fit_param_ind]
    raw_fit_param_uncertainty_dataseries = \
                        fitparamstds_dataseries_list[key_fit_param_ind]

    #get rid of negative xvals (for xvals representing delay time):
    remove_negative_xvals_enabled = True
    zero_delay_offset = 0
    if remove_negative_xvals_enabled:
        fit_param_dataseries = \
                            dsprocessing.get_positive_time_delay_dataseries(
                                    raw_fit_param_dataseries,
                                    zero_delay_offset)
        fit_param_uncertainty_dataseries = \
                            dsprocessing.get_positive_time_delay_dataseries(
                                    raw_fit_param_uncertainty_dataseries,
                                    zero_delay_offset)
    else:
        fit_param_dataseries = raw_fit_param_dataseries
        fit_param_uncertainty_dataseries = raw_fit_param_uncertainty_dataseries

#    #shift xval zero position for fitting convenience:
#    shift_xvals_zero_enabled = True
#    fit_param_dataseries, [fit_param_uncertainty_dataseries] = \
#                            dsprocessing.get_x_offset_dataseries_TRKRstyle(
#                                    fit_param_dataseries,
#                                    [fit_param_uncertainty_dataseries])
    

    #plot parameter dataseries:
    plot_enabled = True
    if plot_enabled:
        plot_xvals = fit_param_dataseries.xvals(unfiltered=True)
        plot_yvals = fit_param_dataseries.yvals(unfiltered=True)
        plot_yerrs = fit_param_uncertainty_dataseries.yvals(unfiltered=True)
        plt.errorbar(plot_xvals, plot_yvals, yerr=plot_yerrs, fmt='.')
#        plt.plot(plot_xvals, plot_yvals, '.')


# %%
    #no loop over dirs implemented yet, so list w/ only one entry:
    fit_param_dataseries_list = [fit_param_dataseries]
    fit_param_weights_list = [fit_param_uncertainty_dataseries]

    #fit lifetimes to each parameter dataseries
    #params = pulse_amplitude, lifetime, offset
    free_params = [True, True, True, True, False]
    initial_params = [0.03, 50, 0.03, 1000, 0]
    param_bounds = [(-1, 1), (1, 100), (-1, 1), (10, 100000), (0, 0)]
    fitdata_list = dataseries_iterable_fit(
                            fit_param_dataseries_list,
                            fitfcn_two_exp_decay,
                            free_params, initial_params, param_bounds,
                            fit_param_weights_list, multiprocessing=False)

    print(fitdata_list[0].fitparams)
    print(fitdata_list[0].fitparamstds)

    #plot parameter dataseries:
    if plot_enabled:
        old_xvals = fitdata_list[0].fitdataseries.xvals(unfiltered=True)
        old_yvals = fitdata_list[0].fitdataseries.yvals(unfiltered=True)
        plot_xvals = np.linspace(min(old_xvals), max(old_xvals), 100)
        plot_yvals = fitfcn_two_exp_decay(plot_xvals,
                                             *fitdata_list[0].fitparams)
        plt.hold(True)
        plt.plot(plot_xvals, plot_yvals, ':')









