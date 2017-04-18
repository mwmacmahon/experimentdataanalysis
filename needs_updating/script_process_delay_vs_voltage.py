# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import SinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer


#==============================================================================
# SCRIPT EXAMPLES - SEE BOTTOM
#==============================================================================

# %% ANALYSIS MODELS
spin_lifetime_model = \
    SinusoidalSpinLifetimeModel(
        max_fcn_evals=200000,
        free_params=[True, True, True, True, True, True, True, True],
        initial_params = [0.002, 50, 0.002, 1000, 810, 0, 0, 0],
        param_bounds = [(-1, 1), (1, 200), (-1, 1), (10, 1e6),
                        (810, 830), (-np.pi, np.pi),
                        (-1e-6, 1e-6), (-0.01, 0.01)],
        error_thresholds=[None, None, 0.2, 2000, None, None, None, None],
        dim2_key="Voltage",
        field_index=0,  # use lockin2x directly instead of, say, area
        excluded_intervals=[[-15, 100], [7000, 15000]])


# %% FILTER FUNCTIONS
def get_filter_fcn_no_super_long_lifetimes(threshold=5000):
    def filter_fcn(scandata):
        if scandata.fitdata_list[0] is None:
            return False
        else:
            return scandata.fitdata_list[0].fitparams[3] < threshold
    return filter_fcn


def get_filter_fcn_no_first_n_scans_in_series(num_ignored=1):
    def filter_fcn(scandata):
        if 'FastScanIndex' in scandata.scaninfo_list[0]:
            return scandata.scaninfo_list[0]['FastScanIndex'] > 1
        else:
            return False
    return filter_fcn


# WARNING: ONLY FOR "C:\Data\160702\delay_scans_-8V_to_8V"
def get_filter_fcn_only_one_channel(target_channel):
    def filter_fcn(scandata):
        try:
            if scandata.scaninfo_list[0]['MiddleScanType'] == 'StageY':
                ypos = scandata.scaninfo_list[0]['MiddleScanCoord']
                if target_channel == 1:
                    return ypos < 3.55
                elif target_channel == 2:
                    return ypos >= 3.55 and ypos < 3.95
                elif target_channel == 3:
                    return ypos >= 3.95
                else:
                    raise Exception('filter_to_one_channel: ' +
                                    'invalid target_channel')
            else:
                print("Warning: filter_to_one_channel: scandata does " +
                      "not conform to expected scan parameters, " +
                      "filtered out.")
                return False
        except KeyError:
            return False
    return filter_fcn


# %% PLOTTING FUNCTION
def plot_scandata(scandata, field_index,
                  label="", fmt="-bd", fit_fmt="xr:"):
    x_vals, y_vals = scandata.dataseries_list[field_index].data()
    error_dataseries = scandata.error_dataseries_list[field_index]
    if error_dataseries is not None:
        _, y_errs = error_dataseries.data()
        plt.errorbar(x_vals, y_vals, yerr=y_errs, label=label, fmt=fmt)
    else:
        plt.plot(x_vals, y_vals, fmt, label=label)
    if scandata.fitdata_list[field_index] is not None:
        x_vals, y_vals = \
            scandata.fitdata_list[field_index].fitdataseries.data()
        plt.plot(x_vals, y_vals, fit_fmt)


# %% OPTIONAL SCRIPTS - TO USE UNCOMMENT if __name__ == "__main__"
#==============================================================================
# ANALYZER FIT FORMULAS:
# (1) delay scans w/ lifetime fits
#         (fields = measurement reading types)
#         (for checkout of fits by eye)
#     scandata_list = \
#         analyzer.collapse_to_scandata_list(filtered=True)
#     field_index = 0: usual lockin2x or w/e (w/ gaussian fit)
#     field_index > 0: other measurements
# 
# 
# (2) lifetime fit data
#         (fields = calculated fit values)
#         (for use of calculated fit values, e.g. lifetimes or gaussian widths)
#     scandata_list = \
#         analyzer.collapse_to_model_fit_scandata_list(new_scan_dimension,
#                                                      filtered=True)
#     field_index = 0: short decay - pulse amplitude
#     field_index = 1: short decay - spin lifetime
#     field_index = 2: long decay - pulse amplitude
#     field_index = 3: long decay - spin lifetime
# 
#==============================================================================
# %%
if __name__ == "__main__":
# %%  ANALYSIS OF DELAY SCANS, LIFETIME VS VOLTAGE
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    analyzer = ScanDataSetsAnalyzer(spin_lifetime_model,
#                                    dirpath="C:\\Data\\160702\\delayscans_-6V_to_6V_noautocenter",
#                                    dirpath="C:\\Data\\160702\\delayscans_-6V_to_6V",
                                    dirpath="C:\\Data\\august_data\\160901\\DelayScansNewWavelength",
                                    uncertainty_value=1e-4,
                                    set_key="Voltage"  # group consecutive runs
                                    )

    # fix improper "StageZ" 2nd coord, change Vapp to V/cm anyway
    for scandataset in analyzer.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.scaninfo_list[0]['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field

#    analyzer.break_up_repeating_scandatasets()  # breaks up sets too much!
    analyzer.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
#    analyzer.apply_transform_to_all_scandata(add_excluded_interval_scandata,
#                                             start=7000, stop=15000)
#    analyzer.apply_transform_to_all_scandata(add_excluded_interval_scandata,
#                                             start=-15, stop=100)
    analyzer.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetimes(threshold=999999))
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    lifetime_scandata_list = \
        analyzer.collapse_to_scandata_list(filtered=False)
    collapsed_scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=False)


# %% OVERVIEW OF FITS
    field_index = 0  # lockin2x
    for scandata in lifetime_scandata_list[:]:
       plot_scandata(scandata, field_index, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    plt.text(8000, 0,
             "Last fit lifetime: {}ns\n     +={}ns".format(
                 scandata.fitdata_list[field_index].fitparams[3]/1000,
                 scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.show()


# %%
    field_index = 3  # dataseries: x:field, y:single long lifetime fit value
    plt.figure()
    plt.hold(True)
    for scandata in collapsed_scandata_list[:]:
        plot_scandata(scandata, field_index, fmt="bd")

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Applied Voltage (V)  |  Electric Field (V/500um)")
    plt.ylabel("Spin lifetime")  
    plt.show()

    
    scandata_list = \
        analyzer.collapse_to_scandata_list(filtered=True)

# %%  MULTICHANNEL ANALYSIS OF "C:\\Data\\160702\\delay_scans_-8V_to_8V"
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    analyzer = ScanDataSetsAnalyzer(spin_lifetime_model,
                                    dirpath="C:\\Data\\160702\\delay_scans_-8V_to_8V",
                                    uncertainty_value=1e-4,
                                    set_key="MiddleScanCoord"  # group by channel
                                    )
#    analyzer.break_up_repeating_scandatasets()  # breaks up sets too much!
    analyzer.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
#    analyzer.apply_transform_to_all_scandata(add_excluded_interval_scandata,
#                                             start=7000, stop=15000)
#    analyzer.apply_transform_to_all_scandata(add_excluded_interval_scandata,
#                                             start=-15, stop=100)
    analyzer.fit_all_scandata_to_model(multiprocessing=True)

    # DEFINE SHARED FILTERS TO USE EVERY TIME FOR CONVENIENCE
    # note: important to _pass_ analyzer, don't use snapshot at this location!
    #       probably okay due to shallow copying, but don't risk it
    def clear_and_apply_shared_filters_to_analyzer(analyzer):  
        analyzer.clear_filters_from_each_scandataset()
        analyzer.add_filter_to_each_scandataset(
            get_filter_fcn_no_super_long_lifetimes(threshold=5000))
        analyzer.add_filter_to_each_scandataset(
            get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))

    # APPLY FILTERS, GET FIRST CHANNEL LIFETIME FIT DATA, PLOT CHANNEL 1
    clear_and_apply_shared_filters_to_analyzer(analyzer)
    analyzer.add_filter_to_each_scandataset(
        get_filter_fcn_only_one_channel(target_channel=1))
    scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=True)
    field_index = 3  # dataseries: x:voltage, y:single long lifetime fit value
    plt.figure()
    plt.hold(True)
    for scandata in scandata_list:
        plot_scandata(scandata, field_index, label="Channel 1", fmt=":rd")

    # REPEAT KEY STEPS FOR CHANNEL 2
    clear_and_apply_shared_filters_to_analyzer(analyzer)
    analyzer.add_filter_to_each_scandataset(
        get_filter_fcn_only_one_channel(target_channel=2))
    scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=True)
    for scandata in scandata_list:
        plot_scandata(scandata, field_index, label="Channel 2 (wired)", fmt=":bd")

    # REPEAT KEY STEPS FOR CHANNEL 3
    clear_and_apply_shared_filters_to_analyzer(analyzer)
    analyzer.add_filter_to_each_scandataset(
        get_filter_fcn_only_one_channel(target_channel=3))
    scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=True)
    for scandata in scandata_list:
        plot_scandata(scandata, field_index, label="Channel 3", fmt=":gd")

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Applied Voltage (V)  |  Electric Field (V/500um)")
    plt.ylabel("Spin polarization")
    plt.legend(loc='center right')
    plt.show()

    # SAVE SCANDATA
