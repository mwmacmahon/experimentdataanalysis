# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata, get_high_pass_filtered_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import IndependentSinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer, sort_scandata_into_sets
import experimentdataanalysis.parsing.dataseriesparsing as dsparsing


#==============================================================================
# SCRIPT EXAMPLES - SEE BOTTOM
#==============================================================================

# ANALYSIS MODELS
# MODEL 1: Two decaying cosines w/ relative phase. Long lifetime fixed.
#          Start beyond pump transient (<200ps lifetime) signal - hidden 3rd.
# params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
#          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
#          phase1, phase2, slope, offset
spin_lifetime_model = \
    IndependentSinusoidalSpinLifetimeModel(
        field_index=0,  # lockin2x
        max_fcn_evals=20000,
        free_params=[False, True, True, True, False,
                     True, True, False, False,
                     True, True, True, True],
        initial_params = [40, 0.01, 0.01, 10000, 20260,
                          550, 554.62, 0, 0,
                          np.pi/3, -np.pi/2, 0, 0],
        param_bounds = [(1, 1000), (0, 1), (0, 1),
                        (1, 1e9), (1e3, 1e9),
                        (500, 600), (500, 600),
                        (0, 0), (0, 0),
                        (-np.pi, np.pi), (-2*np.pi, np.pi),
                        (-1e-4, 1e-4), (-0.01, 0.01)],
        error_thresholds=[None, None, None, None, None,
                          None, None, None, None,
                          None, None, None, None],
        dim2_key="Electric Field (V/cm)",  # look at results of fit vs field
        excluded_intervals=[[-15, 400]])
#        excluded_intervals=[[-15, 400], [7000, 15000]])


# IMPORTANT RESULTS FOR MODEL 1 CONSIDERATION:
# one good RSA fit to long species, @818.0nm, no voltage applied:
# lifetime: 20.26ns +- 0.08ns,
# freq_per_T ~= 5.41e8 Ts (have to grab again) -> g=.00615? uhh
# \-> osc_period @300mT: 554.620ps +- .045ps
# field offset: -0.7565mT +- 0.0056mT

# ANALYSIS MODELS
# MODEL 2: Two decaying cosines w/ relative phase, as before. This time
#          no fixing of long lifetime or oscillation factor, and short
#          lifetime signal is assumed to be very low by starting params
# params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
#          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
#          phase1, phase2, slope, offset
spin_lifetime_model_2 = \
    IndependentSinusoidalSpinLifetimeModel(
        field_index=0,  # lockin2x
        max_fcn_evals=2000,
        free_params=[False, True, False, False, False,
                     True, False, False, False,
                     True, False, True, True],
        initial_params = [40, 0.01, 0, 20260, 20260,
                          808, 808, 0, 0,
                          -np.pi, np.pi/2, 0, 0],
        param_bounds = [(1, 1000), (0, 1), (0, 1),
                        (1, 1e6), (1, 1e6),
                        (750, 850), (750, 850),
                        (0, 0), (0, 0),
                        (-2*np.pi, np.pi), (-2*np.pi, np.pi),
                        (-1e-4, 1e-4), (-0.01, 0.01)],
        error_thresholds=[None, None, None, None, None,
                          None, None, None, None,
                          None, None, None, None],
        dim2_key="Electric Field (V/cm)",  # look at results of fit vs field
#        excluded_intervals=[[-15, 400]])
        excluded_intervals=[[-15, 400], [7000, 15000]])



# FILTER FUNCTIONS
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


# PLOTTING FUNCTION
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


# %%
if __name__ == "__main__":
# %%  ANALYSIS OF 818.9nm DELAY SCANS
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER

#    dirpath = ("C:\\Data\\160702\\" +
#               "delayscans_-6V_to_6V_noautocenter")
#               "delayscans_-6V_to_6V")
#    dirpath = ("C:\\Data\\august_data\\160901\\" +
#               "DelayScansNewWavelength")  # requires other assumptions in model! do seperately below...
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "WavelengthDependence_TRKR_300mT")
#               "WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "LowV_818.0nm_WavelengthDependence_TRKRvsV_200mT")
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    fixed_uncertainty = 1e-4  # manually set uncertainty of data set
    model = spin_lifetime_model
    sort_key = None  # just one ScanDataSet

    scandata_list = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))
    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    analyzer = ScanDataSetsAnalyzer(scandataset_list)


    # fix improper "StageZ" 2nd coord, change Vapp to V/cm anyway
    for scandataset in analyzer.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.scaninfo_list[0]['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field
                scaninfo['Electric Field (V/cm)'] = field

#    analyzer.break_up_repeating_scandatasets()  # breaks up sets too much!
    analyzer.apply_transform_to_all_scandata(get_high_pass_filtered_scandata,
                                              min_freq_cutoff=0)
    analyzer.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
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
    for scandata in lifetime_scandata_list[0:10]:
       plot_scandata(scandata, field_index, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    if scandata.fitdata_list[field_index] is not None:
        plt.text(8000, 0,
                 "Last fit lifetime: {}ns\n     +={}ns".format(
                     scandata.fitdata_list[field_index].fitparams[3]/1000,
                     scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.show()


# %%
    field_indices = [8,9]  # dataseries: x:field, y:short/long phase
    plt.figure()
    plt.hold(True)
    for field_index in field_indices:
        for scandata in collapsed_scandata_list[:]:
            plot_scandata(scandata, field_index, fmt="d",
                          label=scandata.fields[field_index])

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Electric Field (V/cm)")
    plt.ylabel("")  
    plt.legend()
    plt.show()

# %%  ANALYSIS OF 818.0nm DELAY SCANS
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER

#    dirpath = ("C:\\Data\\160702\\" +
#               "delayscans_-6V_to_6V_noautocenter")
#               "delayscans_-6V_to_6V")
#    dirpath = ("C:\\Data\\august_data\\160901\\" +
#               "DelayScansNewWavelength")  # requires other assumptions in model! do seperately below...
    dirpath = ("C:\\Data\\august_data\\160902\\" +
#               "GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "WavelengthDependence_TRKR_300mT")
#               "WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
               "LowV_818.0nm_WavelengthDependence_TRKRvsV_200mT")
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    fixed_uncertainty = 1e-4  # manually set uncertainty of data set
    model = spin_lifetime_model_2
    sort_key = None  # just one ScanDataSet

    scandata_list = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))
    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    analyzer2 = ScanDataSetsAnalyzer(scandataset_list)

    # fix improper "StageZ" 2nd coord, change Vapp to V/cm anyway
    for scandataset in analyzer2.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.scaninfo_list[0]['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field
                scaninfo['Electric Field (V/cm)'] = field

#    analyzer.break_up_repeating_scandatasets()  # breaks up sets too much!
    analyzer2.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
    analyzer2.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetime1s(threshold=999999))
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    lifetime_scandata_list = \
        analyzer2.collapse_to_scandata_list(filtered=False)
    collapsed_scandata_list = \
        analyzer2.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                      filtered=False)


# %% OVERVIEW OF FITS
    field_index = 1  # lockin2x
    for scandata in lifetime_scandata_list[:]:
       plot_scandata(scandata, field_index, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    if scandata.fitdata_list[field_index] is not None:
        plt.text(8000, 0,
                 "Last fit lifetime: {}ns\n     +={}ns".format(
                     scandata.fitdata_list[field_index].fitparams[3]/1000,
                     scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.show()

    print(scandata.fitdata_list[field_index].fitparams)

    scandata_list = lifetime_scandata_list


# %%
    field_indices = [0]  # dataseries: x:field, y:short/long phase
    plt.figure()
    plt.hold(True)
    for field_index in field_indices:
        for scandata in collapsed_scandata_list[:]:
            plot_scandata(scandata, field_index, fmt="d",
                          label=scandata.fields[field_index])

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Applied Voltage (V)  |  Electric Field (V/500um)")
    plt.ylabel("")  
    plt.legend()
    plt.show()

    if scandata.fitdata_list[field_index] is not None:
        print(scandata.fitdata_list[field_index].fitparams)

    scandata_list = collapsed_scandata_list

