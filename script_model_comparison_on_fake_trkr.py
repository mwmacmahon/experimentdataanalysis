# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 13:34:40 2016

@author: Michael
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.filters as filters
import scipy.stats as stats

from experimentdataanalysis.analysis.scandataprocessing \
    import make_scandata_time_delay_positive, \
           make_scandata_phase_continuous, \
           gaussian_smooth_scandata, \
           process_scandata_fields
from experimentdataanalysis.analysis.scandatamodels \
    import RSAFieldScanModel, IndependentSinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer, sort_scandata_into_sets
from experimentdataanalysis.parsing.scandataparsing import \
        fetch_dir_as_unfit_scandata_iterator


GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
#GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar


# FIT RESULT PROCESSING

def add_mean_and_std_info_to_scandata(scandata, indices_to_use):
    def add_mean_and_std_info(xvals, yvals, yerrvals,
                              indices_to_use, scandatainfo,
                              field_name, *args, **kwargs):
        yvals_to_use = yvals[indices_to_use]
        mean_and_std_results = {field_name + '_mean': np.mean(yvals_to_use),
                                field_name + '_std': np.std(yvals_to_use)}
        scandatainfo.update(mean_and_std_results)
        return None, None, None  # only change yvals

    xyyerr_fcn = add_mean_and_std_info
    xyyerr_fcn_args = [indices_to_use, scandata.info]
    xyyerr_fcn_kwargs = {}
    process_scandata_fields(scandata, xyyerr_fcn, xyyerr_fcn_args,
                            xyyerr_fcn_kwargs)
                                        

def add_linear_fit_results_info_to_scandata(scandata, indices_to_use):
    def add_linear_fit_results_info(xvals, yvals, yerrvals,
                                    indices_to_use, scandatainfo,
                                    field_name, *args, **kwargs):
        xvals_to_use = xvals[indices_to_use]
        yvals_to_use = yvals[indices_to_use]
        slope, intercept, _, _, _ = stats.linregress(xvals_to_use,
                                                     yvals_to_use)
        key_xvals = np.array([818.0, 818.9])
        fit_yvals = key_xvals*slope + intercept
        linear_fit_results = {field_name + '_linear_fit_slope': slope,
                              field_name + '_linear_fit_intercept': intercept,
                              field_name + '_linear_fit_at_818.0nm': fit_yvals[0],
                              field_name + '_linear_fit_at_818.9nm': fit_yvals[1]}
        scandatainfo.update(linear_fit_results)
        return None, None, None  # only change yvals

    xyyerr_fcn = add_linear_fit_results_info
    xyyerr_fcn_args = [indices_to_use, scandata.info]
    xyyerr_fcn_kwargs = {}
    process_scandata_fields(scandata, xyyerr_fcn, xyyerr_fcn_args,
                            xyyerr_fcn_kwargs)


def add_mean_and_std_info(scandata, field_name):
    mean_and_std_results = {'mean': np.mean(yvals),
                            'std': np.std(yvals)}
    scandata.info.update(mean_and_std_results)
    return scandata


def add_linear_fit_results_info(scandata, field_name):
    # warning: error not used! must prune bad fits elsewhere
    xvals, yvals = scandata.get_field_xy(field_name)
    slope, intercept, _, _, _ = stats.linregress(xvals, yvals)
    key_xvals = np.array([818.0, 818.9])
    fit_yvals = key_xvals*slope + intercept
    linear_fit_results = {field_name + '_linear_fit_slope': slope,
                          field_name + '_linear_fit_intercept': intercept,
                          field_name + '_linear_fit_at_818.0nm': fit_yvals[0],
                          field_name + '_linear_fit_at_818.9nm': fit_yvals[1]}
    scandata.info.update(linear_fit_results)
    return scandata


# FILTER FUNCTIONS
def get_filter_fcn_no_super_long_lifetime(threshold=5000):
    def filter_fcn(scandata):
        if scandata.fitdata_list[0] is None:
            return False
        else:
            return scandata.fitdata_list[0].fitparams[3] < threshold
    return filter_fcn


def get_filter_fcn_no_first_n_scans_in_series(num_ignored=1):
    def filter_fcn(scandata):
        if 'FastScanIndex' in scandata.info:
            return scandata.info['FastScanIndex'] > 1
        else:
            return False
    return filter_fcn


def get_filter_fcn_this_voltage_only(voltage=0):
    def filter_fcn(scandata):
        if 'Voltage' in scandata.info:
            return scandata.info['Voltage'] == voltage
        else:
            return False
    return filter_fcn


# PLOTTING FUNCTION
def plot_scandata(scandata, field_name, model=None,
                  label="", fmt="-bd", fit_fmt="xr:"):
    x_vals, y_vals, y_errs = scandata.get_field_xyyerr(field_name)
    if y_errs is not None:
        plt.errorbar(x_vals, y_vals, yerr=y_errs, label=label, fmt=fmt)
    else:
        plt.plot(x_vals, y_vals, fmt, label=label)
    if scandata.get_field_fitdata(field_name) is not None:
        y_vals = scandata.get_field_fitdata(field_name).fityvals
        if model is not None:
            x_vals = np.linspace(min(x_vals), max(x_vals), 1000)
            params = scandata.get_field_fitdata(field_name).fitparams
            y_vals = model.fitfunction(x_vals, *params)
        plt.plot(x_vals, y_vals, fit_fmt)


def plot_rsa_fit_scandata(field_name, fit_rsa_scandata_list):
    plt.figure()
    plt.hold(True)
    try:
        for scandata in fit_rsa_scandata_list[:]:
            plot_scandata(scandata, field_name, fmt="bd")
    except TypeError as err:
        raise err
        scandata = fit_rsa_scandata_list
        plot_scandata(scandata, field_name, fmt="bd")
    plt.xlabel("Field (T)")
    plt.ylabel("Kerr Rotation (AU)")
    plt.text(0.065, 0.002,
             "Last fit lifetime:\n{:4.1f}ns +={:4.1f}ns".format(
                 scandata.get_field_fitdata(field_name).fitparams[3]/1000,
                 scandata.get_field_fitdata(field_name).fitparamstds[3]/1000))
    plt.text(0.065, -0.002,
             "Last fit gfactor:\n{:6.3f}/T +={:6.3f}/T".format(
                 scandata.get_field_fitdata(field_name).fitparams[4],
                 scandata.get_field_fitdata(field_name).fitparamstds[4]))
    plt.show()
    print(scandata.get_field_fitdata(field_name).fitparams)


def plot_trkr_fit_scandata(field_name, fit_trkr_scandata_list):
    plt.figure()
    plt.hold(True)
    try:
        for scandata in fit_trkr_scandata_list[:]:
            plot_scandata(scandata, field_name, fmt="bd")
    except TypeError:
        scandata = fit_trkr_scandata_list
        plot_scandata(scandata, field_name, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    plt.text(7100, 0.01,
             "Last fit short lifetime:\n{:4.1f}ns +={:4.1f}ns".format(
                 scandata.get_field_fitdata(field_name).fitparams[3]/1000,
                 scandata.get_field_fitdata(field_name).fitparamstds[3]/1000))
    plt.text(7100, -0.01,
             "Last fit long lifetime:\n{:4.1f}ns +={:4.1f}ns".format(
                 scandata.get_field_fitdata(field_name).fitparams[4]/1000,
                 scandata.get_field_fitdata(field_name).fitparamstds[4]/1000))
    plt.show()
    print(scandata.get_field_fitdata(field_name).fitparams)


def plot_fit_param_scandata(field_names, fit_results_scandata_list):
    plt.figure()
    plt.hold(True)
    for field_name in field_names:
        try:
            for scandata in fit_results_scandata_list[:]:
                plot_scandata(scandata, field_name, fmt="bd", label=field_name)
        except TypeError:
            scandata = fit_results_scandata_list
            plot_scandata(scandata, field_name, fmt="bd", label=field_name)
    # LABEL AND DISPLAY GRAPH
    plt.legend(loc='best')
    plt.show()


# %%  ANALYSIS OF FAKE DELAY SCANS VS B FIELD

    # MODEL 1: Fixed long lifetime
    # params = num_pulses, pulse_amplitude1, pulse_amplitude2,
    #          lifetime1, lifetime2, osc_period,
    #          drift_velocity, slope, offset
    simple_trkr_model = \
        IndependentSinusoidalSpinLifetimeModel(
            field_name="lockin2x",
            max_fcn_evals=20000,
            free_params=[False, True, True,
                         True, True, True,
                         False, True, True],
            initial_params=[40, 0.04, 0.04,
                            20000, 2000, 600,
                            0, 0, 0],
            param_bounds=[(1,1000), (0, 1), (0, 1),
                          (1, 1e9), (1, 1e9), (200, 1200),
                          (-1e3, 1e3), (-1e-6, 1e-6), (-0.01, 0.01)],
            free_params=[False, True, True,
                         False, True, True, True,
                         False, False,
                         False, False, True, True],
            fit_result_scan_coord="Pump-Probe Distance (um)",
            excluded_intervals=[[-15, 400]],
#            excluded_intervals=[[-15, 400], [7000, 15000]],
            BField=300)

    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    dirpath = ("C:\\Data\\fake_data\\fake_trkr")
    fixed_uncertainty = 1e-3  # manually set uncertainty of data set
    model = trkr_model_0V
    sort_key = "Magnetic Field (mT)"  # one scandataset for each field value

    # filename parsing pattern: if string in element[0] is found in filepath
    # separated from other characters by '_'s, will record adjacent number
    # and store in scandata's info dict under the key element[1], as a float
    # if possible.
    # e.g. TRKR_15Vcm_230mT.dat -> {"Electric Field (V/cm)": 15.0,
    #                               "Magnetic Field (mT)": 230.0}
    in_filepath_element_keyword_list = [("Vcm", "Electric Field (V/cm)"),
                                        ("mT", "Magnetic Field (mT)"),
                                        ("K", "SetTemperature (K)"),
                                        ("nm", "Wavelength (nm)"),
                                        ("run", "RunIndex")]
    # for this one, if element[0] found, next element stored w/ key element[1]
    in_filepath_next_element_keyword_list = [("MirrorZ",
                                              "Pump-Probe Distance (um)")]
    filepath_parsing_keyword_lists = [[],
                                      in_filepath_next_element_keyword_list,
                                      in_filepath_element_keyword_list]
    scandata_list = \
        list(fetch_dir_as_unfit_scandata_iterator(
                    directorypath=dirpath,
                    key_field="lockin2x",
                    key_field_error_val=fixed_uncertainty,
                    parsing_keywordlists=filepath_parsing_keyword_lists))

    # ---------------
    # TEMPORARY, FOR SPEED:
#    scandata_list = scandata_list[0:4]
    # ---------------

    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    field_name = model.field_name  # for future use

    # Change drift velocity based on voltage. Assumes each set has same
    # voltage for every scan!
    for scandataset in scandataset_list:
        e_field_list = [scandata.info['Electric Field (V/cm)']
                        for scandata in scandataset.scandata_list]
        b_field_list = [scandata.info['Magnetic Field (mT)']
                        for scandata in scandataset.scandata_list]
        set_voltage = voltage_list[0]
        b_field = b_field_list[0]
        if any(voltage != set_voltage for voltage in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            # Drift velocity per voltage
            drift_velocity_per_volt = 2e-3  # in um/(ps*V)
            drift_velocity = drift_velocity_per_volt*set_voltage  # in um/ps
            scandataset.model.free_params[6] = False  # drift_velocity1
            scandataset.model.initial_params[6] = drift_velocity

        if any(bval != b_field for bval in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            scandata.model.Bfield = b_field

            scandataset.model.initial_params[7] = drift_velocity
            scandataset.model.initial_params[8] = drift_velocity

            guess_period_long = long_period_init + 2.7 * abs(set_voltage)
            scandataset.model.free_params[6] = True  # let it fit from here
            scandataset.model.initial_params[6] = guess_period_long

    analyzer2 = ScanDataSetsAnalyzer(scandataset_list)

#    # smooth over data with a 40ps wide gaussian convolution filter
#    analyzer2.apply_transform_to_all_scandata(
#                                    gaussian_smooth_scandata,
#                                    gaussian_width=40)
    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
    analyzer2.apply_transform_to_all_scandata(
                                    gaussian_smooth_scandata,
                                    fields_to_process=[field_name],
                                    gaussian_width=600,
                                    edge_handling='reflect',
                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer2.apply_transform_to_all_scandata(make_scandata_time_delay_positive,
                                              zero_delay_offset=-15,
                                              neg_time_weight_multiplier=5.0)

    # scandatasets don't share models, can't multiprocess in this version:
    analyzer2.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer2.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetime(threshold=999999))
#    analyzer2.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    fit_trkr_scandata_list = \
        analyzer2.collapse_to_scandata_list(filtered=False)
    trkr_fit_results_scandata_list = \
        analyzer2.collapse_to_model_fit_scandata_list(filtered=False)


# %% OVERVIEW OF FITS
    field_name = "lockin2x"  # x:delay, y:lockin2x
    plot_trkr_fit_scandata(field_name, fit_trkr_scandata_list[3:4])
    scandata_list = fit_trkr_scandata_list


# %%
    field_name = "lockin2x"  # x:wavelength, y:short lifetime
    plt.figure()
    plt.hold(True)
    for scandata in trkr_fit_results_scandata_list[:]:
        if scandata.info['Voltage'] == 0:
            plot_scandata(scandata, field_name, fmt=":bd",
                          label="0 V/cm")
        elif scandata.info['Voltage'] == 0.75:
            plot_scandata(scandata, field_name, fmt=":rd",
                          label="15 V/cm")
    if field_name == 1:  # long pulse amplitude guess from avg. RSA data
        plt.plot(rsa_vs_wavelength_fit_results['wavelengths'],
                 rsa_vs_wavelength_fit_results['pulse_amplitude'],
                 'gd:', label="RSA fit in ROI")
    elif field_name == 3:  # long lifetime
        plt.plot(rsa_vs_wavelength_fit_results['wavelengths'],
                 rsa_vs_wavelength_fit_results['lifetime'],
                 'gd:', label="RSA fit in ROI")
    elif field_name == 4:  # long oscillation period
        plt.plot(rsa_vs_wavelength_fit_results['wavelengths'],
                 rsa_vs_wavelength_fit_results['osc_period_300mT'],
                 'gd:', label="RSA fit in ROI")

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Wavelength (nm)")
    plt.ylabel(scandata.fields[field_name])
    plt.legend(loc='best')
    plt.title("Results of fit to TRKR @300mT,\nlong lifetime & g-factor held constant from RSA")
    plt.show()

    scandata_list = trkr_fit_results_scandata_list
