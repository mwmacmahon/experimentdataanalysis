# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata, add_excluded_interval_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import GaussianModel, SpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer



#==============================================================================
# SCRIPT EXAMPLES - SEE BOTTOM
#==============================================================================

# %% ANALYSIS MODELS
gaussian_model = \
    GaussianModel(max_fcn_evals=80000,
                  error_thresholds=[0.1, 30, 30, None, None])
lifetime_model = \
    SpinLifetimeModel(max_fcn_evals=80000,
                      free_params=[True, True, True, True, True],
                      error_thresholds=[0.2, 200, 0.2, 2000, None],
                      dim2_key="Voltage",
                      field_index=0)  # use amplitude instead of area

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
                  label="", fmt="-bd", fit_fmt="xr:", unfiltered=True):
    x_vals, y_vals = \
        scandata.dataseries_list[field_index].\
                                        datalists(unfiltered=unfiltered)
    error_dataseries = scandata.error_dataseries_list[field_index]
    if error_dataseries is not None:
        _, y_errs = error_dataseries.datalists(unfiltered=unfiltered)
        plt.errorbar(x_vals, y_vals, yerr=y_errs, label=label, fmt=fmt)
    else:
        plt.plot(x_vals, y_vals, fmt, label=label)
    if scandata.fitdata_list[field_index] is not None:
        x_vals, y_vals = \
            scandata.fitdata_list[field_index].fitdataseries.datalists()
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
    analyzer = ScanDataSetsAnalyzer(gaussian_model,
                                    dirpath="C:\\Data\\160702\\delayscans_-6V_to_6V_reference_overlapscans",
                                    uncertainty_value=1e-4)
#    analyzer.break_up_repeating_scandatasets()  # breaks up sets too much!
    analyzer.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS TO GAUSSIAN FIT DATA
    analyzer.add_filter_to_each_scandataset(
        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    # 3rd coord is "Voltage" after MirrorZ and DelayTime
    gaussian_fit_scandata_list = \
        analyzer.collapse_to_scandata_list(filtered=True)
    collapsed_gaussian_fit_scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="Voltage",
                                                     filtered=True)

    # LOAD AND FIT GAUSSIAN FIT DATA TO LIFETIMES
    analyzer2 = ScanDataSetsAnalyzer(
        lifetime_model,
        scandata_list=collapsed_gaussian_fit_scandata_list,
        set_key=None)  # one scandataset per voltage?
    analyzer2.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
    analyzer2.apply_transform_to_all_scandata(add_excluded_interval_scandata,
                                             start=7000, stop=15000)
    analyzer2.apply_transform_to_all_scandata(add_excluded_interval_scandata,
                                             start=-15, stop=100)
    analyzer2.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS, PLOT CHANNEL LIFETIME FIT DATA
    analyzer2.add_filter_to_each_scandataset(
        get_filter_fcn_no_super_long_lifetimes(threshold=5000))
    analyzer2.add_filter_to_each_scandataset(
        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    # no 4th coord after voltage needed, not 4D!
    lifetime_scandata_list = \
        analyzer2.collapse_to_scandata_list(filtered=True)
    collapsed_lifetime_scandata_list = \
        analyzer2.collapse_to_model_fit_scandata_list(
                                                new_scan_type="[Unknown]",
                                                filtered=True)
    plt.figure()
    plt.hold(True)

#==============================================================================
#     # OVERLAP SCANS: PLOT, LABEL, DISPLAY
#     field_index = 0  # lockin output
#     for scandata in gaussian_fit_scandata_list:
#         plot_scandata(scandata, field_index, fmt="-bd")
#     plt.xlabel("Pump-Probe Distance (um)")
#     plt.ylabel("Spin Lifetime (ps)")  # assumes field_index = 0
#     plt.show()
#==============================================================================

#==============================================================================
#     # GAUSSIAN FIT PROPERTIES VS TIME DELAY: PLOT, LABEL, DISPLAY
#     field_index = 0  # gaussian amplitude
# #    field_index = 1  # gaussian width
# #    field_index = 2  # gaussian center pos
# #    field_index = 3  # gaussian area (w/ lifetime fit)
#     for scandata in lifetime_scandata_list:
#         plot_scandata(scandata, field_index, fmt="-bd")
#     plt.xlabel("Pump-Probe Time Delay (ps)")
#     plt.ylabel("Gaussian Fit Amplitude (AU)")  # assumes field_index = 0
#     plt.show()
#==============================================================================

    # SPIN LIFETIME FIT PROPERTIES VS VOLTAGE: PLOT, LABEL, DISPLAY
#    field_index = 0  # short decay - pulse amplitude
#    field_index = 1  # short decay - spin lifetime
#    field_index = 2  # long decay - pulse amplitude
    field_index = 3  # long decay - spin lifetime
    for scandata in collapsed_lifetime_scandata_list:
       plot_scandata(scandata, field_index, fmt="-bd")
    plt.xlabel("Applied Voltage (V)  |  Electric Field (V/500um)")
    plt.ylabel("Spin Lifetime (ps)")  # assumes field_index = 3
    plt.show()


    # PICK: WHAT TO DROP INTO WORKSPACE    
    scandata_list = gaussian_fit_scandata_list
#    scandata_list = lifetime_scandata_list
#    scandata_list = collapsed_lifetime_scandata_list
