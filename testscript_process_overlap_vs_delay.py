# -*- coding: utf-8 -*-
"""
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
        fitfcn_two_exp_decay, fitfcn_exp_times_sqrt_decay, \
        fitfcn_simple_line
from experimentdataanalysis.analysis.multidataseriesprocessing \
    import dataseries_iterable_fit, scandata_iterable_fit, \
        scandata_iterable_sort
import experimentdataanalysis.parsing.dataclassparsing as dcparsing

# %%
class SubDirectoryData:
    def __init__(self, subdirname):
        self.subdirname = subdirname
        self.scandata_list = []
        self.chronological_key = "FastScanIndex"
        self.xval_key = "MiddleScanCoord"
        self.fit_param_dataseries_list = None  # placeholders, not actually needed
        self.fit_param_uncertainty_dataseries_list = None

# %%
if __name__ == "__main__":

# %%
    scandata_list = list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
#        directorypath="C:\\Data\\early_may_good_data"))
        directorypath="C:\\Data\\early_may_good_data"))
#        directorypath="C:\\Data\\160505\\Channel 3\\Voltage_4_Experiment_Channel3_1mT_to4V_033XT-B11_819.5nm_30K_2Dscan_DelayTime_MirrorZ"))
#        directorypath="C:\\Data\\160506\\Voltage_4_OverlapVsDelay_Channel3_1mT_033XT-B11_819.5nm_30K_2Dscan_DelayTime_MirrorZ"))
#        directorypath="C:\\Data\\160506\\Voltage_8_OverlapVsDelay_Channel3_1mT_033XT-B11_819.5nm_30K_2Dscan_DelayTime_MirrorZ"))
#    scandata_list += list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelmobility_200mT_Channel3_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime"))

#    dataseries_list = [scandata.dataseries[0] for scandata in scandata_list]

# %%
# sort scandata by directory into sublists:
subdir_list = []
current_subdirdata = SubDirectoryData("[default]")
last_subdirname = ""
for scandata in scandata_list:
    subdirname = scandata.filepath.split("\\")[-2]
    if subdirname != last_subdirname:
        current_subdirdata = SubDirectoryData(subdirname)
        current_subdirdata.chronological_key = "FastScanIndex"
        current_subdirdata.xval_key = "MiddleScanCoord"
        subdir_list.append(current_subdirdata)
        last_subdirname = subdirname
    current_subdirdata.scandata_list.append(scandata)

# %%
# sort subdirectories by chronological order
# to get rid of misorderings due to file parsing order being alphanumeric
for subdirdata in subdir_list:
    subdir_scandata_list = subdirdata.scandata_list
    primary_key = subdirdata.chronological_key
    secondary_key = subdirdata.chronological_key
    subdirdata.scandata_list, _ = scandata_iterable_sort(subdir_scandata_list,
                                                         primary_key,
                                                         secondary_key)

# %%
# break into smaller pseduo-subdir groups based on chronological order
# and middle coordinate - e.g. x=1,2,3,4,1,2,3,4 -> x=1,2,3,4 ; x=1,2,3,4
break_up_subdirs_w_repeats = True
if break_up_subdirs_w_repeats:
    new_subdir_list = []
    for subdirdata in subdir_list:
        # now, group by x-values
        xvals = [scandata.scaninfo[subdirdata.xval_key]
                 for scandata in subdirdata.scandata_list]
        repeat_counts = [xvals[:ind].count(xvals[ind])
                         for ind in range(len(xvals))]
        sublist_index_lists = []
        # if xvals can be broken into continuous sublists w/ identical sets
        # of values, do so! this test catches even sublists w/ different
        # orderings. Note: this requires a properly sorted scandata list!
        if max(repeat_counts) > 0 and repeat_counts == sorted(repeat_counts):
            for sublist_index in range(max(repeat_counts) + 1):
                index_list = [ind for ind, count in enumerate(repeat_counts)
                              if count == sublist_index]
                sublist_index_lists.append(index_list)
            for new_subdir_indices in sublist_index_lists:
                new_subdirdata = SubDirectoryData(subdirdata.subdirname)
                new_subdirdata.scandata_list = [subdirdata.scandata_list[ind]
                                                for ind in new_subdir_indices]
                new_subdirdata.xval_key = subdirdata.xval_key
                new_subdir_list.append(new_subdirdata)
        else:
            new_subdir_list.append(subdirdata)
    num_old = sum(len(subdir.scandata_list) for subdir in subdir_list)
    num_new = sum(len(subdir.scandata_list) for subdir in new_subdir_list)
    if num_old == num_new:
        subdir_list = new_subdir_list
    else:
        print("Error: lost {} scandata breaking up subdirectories"\
                .format(num_old - num_new))

# %%
for subdirdata in subdir_list:
    # fit gaussians to data
    # params = amplitude, t0, sigma, offset
    key_field_ind = 0
    fastscan_fitfunction = fitfcn_simple_1d_gaussian
    free_params = [True, True, True, True]
    initial_params = [0.001, 0, 20, 0]
    param_bounds = [(0, 0.1), (-200, 200), (5, 200), (-0.01, 0.01)]
    fit_scandata_list = scandata_iterable_fit(
                            subdirdata.scandata_list, key_field_ind,
                            fastscan_fitfunction,
                            free_params, initial_params, param_bounds,
                            multiprocessing=False)

    # create dataseries of fit parameters & their uncertainties vs some xval:
    xval_key = subdirdata.xval_key
    xvals = [scandata.scaninfo[xval_key] for scandata in fit_scandata_list]
    fitparams = list(zip(*(scandata.fitdata[key_field_ind].fitparams
                           for scandata in fit_scandata_list)))
    fitparams_dataseries_list = [DataSeries(zip(xvals, fitparam))
                                 for fitparam in fitparams]
    fitparamstds = list(zip(*(scandata.fitdata[key_field_ind].fitparamstds
                              for scandata in fit_scandata_list)))
    fitparamstds_dataseries_list = [DataSeries(zip(xvals, fitparamstd))
                                    for fitparamstd in fitparamstds]
#    # screwing around with errors for absolute_sigma testing
#    fitparamstds_dataseries_list = [DataSeries(zip(xvals, [val for val in fitparamstd[:-1]] + [fitparamstd[-1]]))
#                                    for fitparamstd in fitparamstds]

    # save this data to subdir data structure:
    subdirdata.scandata_list = fit_scandata_list
    subdirdata.fit_param_dataseries_list = fitparams_dataseries_list
    subdirdata.fit_param_uncertainty_dataseries_list = \
                                                fitparamstds_dataseries_list

# %%
fit_result_info_list = []
fit_result_dataseries_list = []
fit_result_sigmas_list = []
for ind, subdirdata in enumerate(subdir_list):
    # calculate fit result from fit parameters,
    # define function to fit _that_ to

    # fit result: gaussian center location
    result_fitfunction = fitfcn_simple_line
    # params = slope, offset
    free_params = [True, True]
    initial_params = [1, 1]
    param_bounds = [(-1000, 1000), (-1000, 1000)]
    fit_result_dataseries = subdirdata.fit_param_dataseries_list[1]
    fit_result_uncertainty_dataseries = \
                subdirdata.fit_param_uncertainty_dataseries_list[1]
    uncertainty_threshold = 100

#    # fit result: gaussian width 
#    result_fitfunction = fitfcn_simple_line
#    # params = slope, offset
#    free_params = [True, True]
#    initial_params = [0, 0]
#    param_bounds = [(-1000, 1000), (-1000, 1000)]
#    fit_result_dataseries = subdirdata.fit_param_dataseries_list[2]
#    fit_result_uncertainty_dataseries = \
#                subdirdata.fit_param_uncertainty_dataseries_list[2]
#    uncertainty_threshold = 100

#    # fit result: area under gaussian (amplitude * sqrt(width)) (w/o consts)
#    result_fitfunction = fitfcn_two_exp_decay
#    # params = pulse1_amp, lifetime1, pulse2_amp, lifetime2, offset
#    free_params = [True, True, True, True, False]
#    initial_params = [0.05, 100, 0.05, 2000, 0]
#    param_bounds = [(-1, 1), (10, 500), (-1, 1), (10, 1e6), (1, -1)]
#    amplitudes = subdirdata.fit_param_dataseries_list[0]
#    amplitudes_sigma = subdirdata.fit_param_uncertainty_dataseries_list[0]
#    widths = subdirdata.fit_param_dataseries_list[2]
#    widths_sigma = subdirdata.fit_param_uncertainty_dataseries_list[2]
#    fit_result_dataseries = DataSeries(
#        [(x, a*w)
#         for (x, a), (_, w) in zip(amplitudes.datatuples(raw=True),
#                                   widths.datatuples(raw=True))],
#        excluded_intervals=amplitudes.excluded_intervals())
#    # note: uncertainty calc. assumes width, amplitude independent...
#    fit_result_uncertainty_dataseries = DataSeries(
#        [(x, np.sqrt((a*w_sig)**2 + (w*a_sig)**2))
#         for (x, a), (_, w), (_, a_sig), (_, w_sig) in \
#                                 zip(amplitudes.datatuples(raw=True),
#                                     widths.datatuples(raw=True),
#                                     amplitudes_sigma.datatuples(raw=True),
#                                     widths_sigma.datatuples(raw=True))],
#        excluded_intervals=amplitudes.excluded_intervals())
#    uncertainty_threshold = 0.1

    # purge poor quality fits
    badindices = []
    goodindices = []
    for (ind, scandata), (x, result), (_, uncertainty) in zip(
                                        enumerate(subdirdata.scandata_list),
                                        fit_result_dataseries,
                                        fit_result_uncertainty_dataseries):
        fitparams = scandata.fitdata[key_field_ind].fitparams
        fitparamstds = scandata.fitdata[key_field_ind].fitparamstds
        bad_data = False
        if fitparamstds[0] > 0.01:
            bad_data = True
        if uncertainty > uncertainty_threshold:
            bad_data = True
        # other "don't use this data" flags go here...
        if bad_data:
            badindices.append(ind)
        else:
            goodindices.append(ind)
    if len(goodindices) > 0:  # ignore subdirs w/ no data / bad scan params
        fit_result_dataseries = DataSeries(
            [list(fit_result_dataseries.datatuples(raw=True))[i]
             for i in goodindices])
        fit_result_uncertainty_dataseries = DataSeries(
            [list(fit_result_uncertainty_dataseries.datatuples(raw=True))[i]
             for i in goodindices])
    
        # get rid of negative xvals (for xvals representing delay time):
        remove_negative_xvals_enabled = True
        zero_delay_offset = 0
        if remove_negative_xvals_enabled:
            fit_result_dataseries = \
                        dsprocessing.get_positive_time_delay_dataseries(
                                        fit_result_dataseries,
                                        zero_delay_offset)
            fit_result_uncertainty_dataseries = \
                        dsprocessing.get_positive_time_delay_dataseries(
                                        fit_result_uncertainty_dataseries,
                                        zero_delay_offset)
    
    #    # shift xval zero position for fitting convenience:
    #    shift_xvals_zero_enabled = True
    #    fit_result_dataseries, [fit_result_uncertainty_dataseries] = \
    #                            dsprocessing.get_x_offset_dataseries_TRKRstyle(
    #                                    fit_result_dataseries,
    #                                    [fit_result_uncertainty_dataseries])
    
        #add extracted dataseries / uncertainty dataseries to analysis list
        fit_result_info_list.append(subdirdata.scandata_list[0].scaninfo)
        fit_result_dataseries_list.append(fit_result_dataseries)
        fit_result_sigmas_list.append(fit_result_uncertainty_dataseries)

# %%
if True:
    # TEMP
    fit_result_dataseries_list = fit_result_dataseries_list[:150]
    fit_result_sigmas_list = fit_result_sigmas_list[:150]

    # fit lifetimes to each parameter dataseries, altogether for efficiency
    fit_result_fitdata_list = dataseries_iterable_fit(
                                fit_result_dataseries_list,
                                result_fitfunction,
                                free_params, initial_params, param_bounds,
                                fit_result_sigmas_list, max_fcn_evals=20000,
                                multiprocessing=False)

# %%
if True:
    # plot a parameter dataseries logarithmically:
    plot_enabled = True
    semilog_plot = True
    plot_subdir_indices = [8]
    for ind, (info, param_fit, param_series, param_std_series) in \
            enumerate(zip(fit_result_info_list, fit_result_fitdata_list,
                          fit_result_dataseries_list, fit_result_sigmas_list)):
        if plot_enabled and ind in plot_subdir_indices:
            if semilog_plot:
                pltfcn = plt.semilogy
            else:
                pltfcn = plt.plot

            plt.figure()
            plot_xvals = param_series.xvals(unfiltered=True)
            plot_yvals = param_series.yvals(unfiltered=True)
            plot_yerrs = param_std_series.yvals(unfiltered=True)
            pltfcn(plot_xvals, plot_yvals, 'd')

            plt.hold(True)

            old_xvals = param_fit.fitdataseries.xvals(unfiltered=True)
            old_yvals = param_fit.fitdataseries.yvals(unfiltered=True)
            plot_xvals = np.linspace(min(old_xvals), max(old_xvals), 100)
            plot_yvals = \
                result_fitfunction(plot_xvals, *param_fit.fitparams)
            pltfcn(plot_xvals, plot_yvals, ':')

# %%
if True:
    # plot a parameter dataseries:
    plot_enabled = True
    plot_subdir_indices = [8]
    for ind, (info, param_fit, param_series, param_std_series) in \
            enumerate(zip(fit_result_info_list, fit_result_fitdata_list,
                          fit_result_dataseries_list, fit_result_sigmas_list)):
        if plot_enabled and ind in plot_subdir_indices:
            plt.figure()
            plot_xvals = param_series.xvals(unfiltered=True)
            plot_yvals = param_series.yvals(unfiltered=True)
            plot_yerrs = param_std_series.yvals(unfiltered=True)
            plt.errorbar(plot_xvals, plot_yvals, yerr=plot_yerrs, fmt='d')

            plt.hold(True)

            old_xvals = param_fit.fitdataseries.xvals(unfiltered=True)
            old_yvals = param_fit.fitdataseries.yvals(unfiltered=True)
            plot_xvals = np.linspace(min(old_xvals), max(old_xvals), 100)
            plot_yvals = \
                result_fitfunction(plot_xvals, *param_fit.fitparams)
            plt.plot(plot_xvals, plot_yvals, ':')
            plt.title(info["Filepath"].split("\\")[-2])

# %%
if True:
    # sort through gaussian mobilities
    scaninfos = [scaninfo
                 for scaninfo, fitdata in zip(fit_result_info_list,
                                              fit_result_fitdata_list)
                 if fitdata is not None]
    voltages = [scaninfo["Voltage"] for scaninfo in scaninfos]
    mobilities = [fitdata.fitparams[0] for fitdata in fit_result_fitdata_list
               if fitdata is not None]
    mobility_sigmas = [fitdata.fitparamstds[0]
                     for fitdata in fit_result_fitdata_list
                     if fitdata is not None]

    #filter out bad ones:
    badindices = []
    goodindices = []
    for ind, (scaninfo, voltage, mobility, mobility_sigma) in \
                        enumerate(zip(scaninfos, voltages,
                                      mobilities, mobility_sigmas)):
        bad_data = False
        if scaninfo["Channel"] == 3:
            bad_data = True
        if abs(mobility) > 0.1:
            bad_data = True
        if mobility_sigma > 0.1:
            bad_data = True
        # other "don't use this data" flags go here...
        if bad_data:
            badindices.append(ind)
        else:
            goodindices.append(ind)
    voltages = [voltages[i] for i in goodindices]
    mobilities = [mobilities[i] for i in goodindices]
    mobility_sigmas = [mobility_sigmas[i] for i in goodindices]

    plt.errorbar(voltages, mobilities, yerr=mobility_sigmas, fmt='.')
    plt.title("Lifetime vs Electric Field")
    plt.xlabel("Electric field * 500um (V)")
    plt.ylabel("Crudely fitted mobilities")



# %%
if True:
    # sort through lifetimes
    scaninfos = [scaninfo
                 for scaninfo, fitdata in zip(fit_result_info_list,
                                              fit_result_fitdata_list)
                 if fitdata is not None]
    voltages = [scaninfo["Voltage"] for scaninfo in scaninfos]
    lifetimes = [fitdata.fitparams[3] for fitdata in fit_result_fitdata_list
                 if fitdata is not None]
    lifetime_sigmas = [fitdata.fitparamstds[3]
                       for fitdata in fit_result_fitdata_list
                       if fitdata is not None]

    #filter out bad ones:
    badindices = []
    goodindices = []
    for ind, (scaninfo, voltage, lifetime, lifetime_sigma) in \
                                    enumerate(zip(scaninfos, voltages,
                                                  lifetimes, lifetime_sigmas)):
        bad_data = False
        if lifetime > 15000:
            bad_data = True
        if lifetime_sigma > 1000:
            bad_data = True
        # other "don't use this data" flags go here...
        if bad_data:
            badindices.append(ind)
        else:
            goodindices.append(ind)
    voltages = [voltages[i] for i in goodindices]
    lifetimes = [lifetimes[i] for i in goodindices]
    lifetime_sigmas = [lifetime_sigmas[i] for i in goodindices]

    plt.errorbar(voltages, lifetimes, yerr=lifetime_sigmas, fmt='.')
    plt.title("Lifetime vs Electric Field")
    plt.xlabel("Electric field * 500um (V)")
    plt.ylabel("Crudely fitted lifetime")



# %%
    # TODO: plot ln(signal) vs time delay and see if linear?
