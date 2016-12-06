# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata, get_high_pass_filtered_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import RSAFieldScanModel, IndependentSinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer


#==============================================================================
# SCRIPT EXAMPLES - SEE BOTTOM
#==============================================================================

# ANALYSIS MODELS
# MODEL 1: RSA for only longest lifetime species. 
# params = num_pulses, delay_time,
#          pulse_amplitude, lifetime, freq_per_T,
#          field_offset, drift_velocity, phase, slope, offset
rsa_field_scan_model = \
    RSAFieldScanModel(
        field_index=0,  # lockin2x
        max_fcn_evals=40000,
        free_params=[False, False,
                     True, True, True,
                     False, False, True, True, True],
        initial_params=[40, -160,
                        0.005, 2e4, 0.006,
                        -0.0012, 0, -np.pi, 0, 0.001],
        param_bounds=[(1, 1000), (-1000, 10000),
                      (0, 1), (1, 1e9), (1e-4, 0.1),
                      (-0.1, 0.1), (-0.1, 0.1),
                      (-2*np.pi, 2*np.pi), (-0.5, 0.5), (-0.1, 0.1)],
        error_thresholds=[None, None,
                          None, None, None,
                          None, None, None, None, None],
        dim2_key="Wavelength")

# MODEL 2: Two decaying cosines w/ relative phase. Long lifetime fixed.
#          Start beyond pump transient (<200ps lifetime) signal - hidden 3rd.
# params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
#          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
#          phase1, phase2, slope, offset
spin_lifetime_model = \
    IndependentSinusoidalSpinLifetimeModel(
        field_index=0,  # lockin2x
        max_fcn_evals=20000,
        free_params=[False, True, True, True, False,
                     True, False, False, False,
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
        dim2_key="Wavelength",
        excluded_intervals=[[-15, 400]])
#        excluded_intervals=[[-15, 400], [7000, 15000]])

# IMPORTANT RESULTS FOR MODEL 2 CONSIDERATION:
# one good RSA fit to long species, @818.0nm, no voltage applied:
# lifetime: 20.26ns +- 0.08ns,
# freq_per_T ~= 5.41e8 Ts (have to grab again) -> g=.00615? uhh
# \-> osc_period @300mT: 554.620ps +- .045ps
# field offset: -0.7565mT +- 0.0056mT
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


def get_filter_fcn_this_voltage_only(voltage=0):
    def filter_fcn(scandata):
        if 'Voltage' in scandata.scaninfo_list[0]:
            return scandata.scaninfo_list[0]['Voltage'] == voltage
        else:
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
# %%  ANALYSIS OF DELAY SCANS, LIFETIME VS VOLTAGE
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    analyzer = ScanDataSetsAnalyzer(spin_lifetime_model,
#                                    dirpath="C:\\Data\\160702\\delayscans_-6V_to_6V_noautocenter",
#                                    dirpath="C:\\Data\\160702\\delayscans_-6V_to_6V",
#                                    dirpath="C:\\Data\\august_data\\160901\\DelayScansNewWavelength",
#                                    dirpath="C:\\Data\\august_data\\160902\\WavelengthDependence_TRKR_300mT\\WavelengthDependence_TRKR_300mT_033XT-A5_817.7nm_30K_2Dscan_Voltage_DelayTime",
                                    dirpath="C:\\Data\\august_data\\160902\\WavelengthDependence_TRKR_300mT",
#                                    dirpath="C:\\Data\\august_data\\160902\\WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime",
#                                    dirpath="C:\\Data\\august_data\\160902\\GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime",
#                                    dirpath="C:\\Data\\august_data\\160902\\BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime",
                                    uncertainty_value=1e-4,
                                    set_key="Voltage"
                                    )

    # fix improper "StageZ" 2nd coord, change Vapp to V/cm anyway
    for scandataset in analyzer.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.scaninfo_list[0]['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field

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
    field_index = 0  # x:delay, y:lockin2x
    for scandata in lifetime_scandata_list[:]:
        plot_scandata(scandata, field_index, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    plt.text(7100, 0.01,
             "Last fit short lifetime:\n{:4.1f}ns +={:4.1f}ns".format(
                 scandata.fitdata_list[field_index].fitparams[3]/1000,
                 scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.text(7100, -0.01,
             "Last fit long lifetime:\n{:4.1f}ns +={:4.1f}ns".format(
                 scandata.fitdata_list[field_index].fitparams[4]/1000,
                 scandata.fitdata_list[field_index].fitparamstds[4]/1000))
    plt.show()
    print(scandata.fitdata_list[field_index].fitparams)

    scandata_list = lifetime_scandata_list


# %%
    field_index = 2  # x:wavelength, y:short lifetime
    plt.figure()
    plt.hold(True)
    for scandata in collapsed_scandata_list[:]:
        if scandata.scaninfo_list[0]['Voltage'] == 0:
            plot_scandata(scandata, field_index, fmt=":bd",
                          label="0 V/cm")
        elif scandata.scaninfo_list[0]['Voltage'] == 0.75:
            plot_scandata(scandata, field_index, fmt=":rd",
                          label="15 V/cm")

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Short Spin Lifetime (ps)")
    plt.legend()
    plt.title("Results of fit to TRKR @300mT,\nlong lifetime & g-factor held constant from RSA")
    plt.show()

    scandata_list = collapsed_scandata_list


# %%  ANALYSIS OF RSA FIELD SCAN DATA
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    analyzer2 = ScanDataSetsAnalyzer(rsa_field_scan_model,
                                    dirpath="C:\\Data\\august_data\\160902\\WavelengthDependence_RSA",
                                    uncertainty_value=1e-4,
                                    set_key=None  # all runs grouped
                                    )

#    analyzer2.apply_transform_to_all_scandata(get_high_pass_filtered_scandata,
#                                              min_freq_cutoff=0.01)
    analyzer2.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer2.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetimes(threshold=999999))
#    analyzer2.add_filter_to_each_scandataset(
#        get_filter_fcn_this_voltage_only(voltage=0))

    lifetime_scandata_list = \
        analyzer2.collapse_to_scandata_list(filtered=False)
    collapsed_scandata_list = \
        analyzer2.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=False)


# %% OVERVIEW OF RSA FITS
    field_index = 0  # x:field, y:lockin2x
    for scandata in lifetime_scandata_list[0:9]:
       plot_scandata(scandata, field_index, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    plt.text(0.065, 0.002,
             "Last fit lifetime:\n{:4.1f}ns +={:4.1f}ns".format(
                 scandata.fitdata_list[field_index].fitparams[3]/1000,
                 scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.text(0.065, -0.002,
             "Last fit freq_per_T:\n{:6.3f}/T +={:6.3f}/T".format(
                 scandata.fitdata_list[field_index].fitparams[4],
                 scandata.fitdata_list[field_index].fitparamstds[4]))
    plt.show()
    print(scandata.fitdata_list[field_index].fitparams)


# %%
    field_index = 2  # x:wavelength, y:lifetime
    plt.figure()
    plt.hold(True)
    for scandata in collapsed_scandata_list[:]:
        plot_scandata(scandata, field_index, fmt="d",
                      label=scandata.fields[field_index])

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("")
    plt.legend()
    plt.title("RSA fit results")
    plt.show()

    scandata_list = collapsed_scandata_list

