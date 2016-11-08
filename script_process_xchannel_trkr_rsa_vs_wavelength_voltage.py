# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.filters as filters
import scipy.stats as stats

from experimentdataanalysis.analysis.scandataprocessing \
    import get_positive_time_delay_scandata, \
           get_continuous_phase_scandata, \
           get_gaussian_smoothed_scandata, \
           aggregate_and_process_scandata_into_dict
from experimentdataanalysis.analysis.scandatamodels \
    import RSAFieldScanModel, IndependentSinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer, sort_scandata_into_sets
from experimentdataanalysis.parsing.scandataparsing import \
        fetch_dir_as_unfit_scandata_iterator


GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar


# FIT RESULT PROCESSING

def get_mean_and_std_dict(xvals, yvals, yerrorvals):
    mean_and_std_results = {'mean': np.mean(yvals),
                            'std': np.std(yvals)}
    return mean_and_std_results


def get_linear_fit_results_dict(xvals, yvals, yerrorvals):
    # warning: error not used! must prune bad fits elsewhere
    slope, intercept, _, _, _ = stats.linregress(xvals, yvals)
    key_xvals = np.array([818.0, 818.9])
    fit_yvals = key_xvals*slope + intercept
    linear_fit_results = {'linear_fit_slope': slope,
                          'linear_fit_intercept': intercept,
                          'linear_fit_at_818.0nm': fit_yvals[0],
                          'linear_fit_at_818.9nm': fit_yvals[1]}
    return linear_fit_results


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


# %%
if __name__ == "__main__":
# %%  ANALYSIS OF RSA FIELD SCAN DATA
    # MODEL 1: RSA for only longest lifetime species. 
    # params = num_pulses, delay_time,
    #          pulse_amplitude, lifetime, gfactor,
    #          field_offset, drift_velocity, phase, slope, offset
    rsa_field_scan_model = \
        RSAFieldScanModel(
            field_name="lockin2x",
            max_fcn_evals=40000,
            free_params=[False, False,
                         True, True, True,
                         True, False, True, True, True],
            initial_params=[40, -160,
                            0.007, 2e4, 0.431,
                            -0.0012, 0, -np.pi, 0, 0],
            param_bounds=[(1, 1000), (-1000, 10000),
                          (0, 1), (1, 1e9), (0.3, 0.6),
                          (-0.1, 0.1), (-0.1, 0.1),
                          (-2*np.pi, 2*np.pi), (-0.5, 0.5), (-0.1, 0.1)],
            error_thresholds=[None, None,
                              None, None, None,
                              None, None, None, None, None],
            fit_result_scan_coord="Wavelength")


    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "WavelengthDependence_RSA")
    fixed_uncertainty = 1e-3  # manually set uncertainty of data set
    model = rsa_field_scan_model
    sort_key = None  # all scans grouped together
    scandata_list = list(fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))
    original_scandata_list_rsa = scandata_list  # don't lose this...
    # ---------------
    # TEMPORARY, FOR SPEED:
#    scandata_list = scandata_list[5:6]
    # ---------------
    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    field_name = model.field_name  # for future use
    analyzer = ScanDataSetsAnalyzer(scandataset_list)
#    # smooth over data with a 40ps wide gaussian convolution filter
#    analyzer.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    gaussian_width=40)
    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
#    analyzer.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    field_names_to_process=[field_name],
#                                    gaussian_width=0.01,
#                                    edge_handling='reflect',
#                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
    analyzer.fit_all_scandata_to_model(multiprocessing=True)

    # use smoothed-over difference of raw data and fits for estimate of drift,
    # since can't use smoothed-over raw data with RSA as it doesn't average to
    # zero in a smooth fashion. Modify data in-place and fit again.
    for scandataset in analyzer.scandataset_list:
        for scandata in scandataset.scandata_list:
            yvals = scandata.y
            fityvals = scandata.fitdata.fityvals
            smoothedoffset = filters.gaussian_filter1d(yvals - fityvals,
                                                       sigma=10,
                                                       mode='reflect')    
            scandata.y = yvals - smoothedoffset
#            plt.figure()
#            plt.subplot(2,1,1)
#            plt.plot(xvals, yvals, label='before drift correction')
#            plt.subplot(2,1,2)
#            plt.plot(xvals, scandata.dataseries_list[field_name].yvals(),
#                     label='after drift correction')
    analyzer.fit_all_scandata_to_model(multiprocessing=True)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetime(threshold=999999))
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_this_voltage_only(voltage=0))

    fit_rsa_scandata_list = \
        analyzer.collapse_to_scandata_list(filtered=False)
    rsa_fit_results_scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=False)


# %% OVERVIEW OF RSA FITS
    field_name = "lockin2x"  # x:field, y:lockin2x
    plot_rsa_fit_scandata(field_name, fit_rsa_scandata_list[:])
    scandata_list = fit_rsa_scandata_list

# %%
    field_names = ["gfactor"]
    plot_fit_param_scandata(field_names, rsa_fit_results_scandata_list)
    plt.xlabel("Wavelength (nm)")
    plt.title("RSA fit results")
    plt.show()

    scandata_list = rsa_fit_results_scandata_list


# %% MODEL 1 EXTRACTED DATA REVIEW
    # First, extract all data from fit results into a comprehensive dict
    # This gives us an estimate of several factors easily looked up in dict
    # as a function of wavelength, and leaving out "bad fits":
    # - "wavelengths": the wavelengths of which the other fields are a function
    # - "pulse_amplitude": spin amplitude generated per pulse, long species(?)
    # - "lifetime": spin lifetime, long species (?)
    # - "gfactor": precession factor, long species (?). wavelength-dependent!?
    # - "field_offset": how off-center the RSA was from B=0, mostly external B
    # Furthermore, these flags can be added:
    # "_mean" brings a scalar: the mean of the values across wavelengths
    # "_std" brings a scalar: the std around said mean value across wavelengths
    # "_error" brings up the error from the diagonal elements of 
    #           the fit covariance matrix for each wavelength, and
    # "_linear_fit_at_818.0nm" [or 818.9nm] brings the linear interpolation
    #           estimate for the value at those two wavelengths.

    # First, trim indices to those in the region of interest, namely wavelength
    # less than 818.9, as above that plot results differ significantly, and
    # later datasets will be at 818.0 and 818.9, and finally the 818.9 data
    # set is too low quality to use its results
    indices_to_use = []
    full_wavelength_list = rsa_fit_results_scandata_list[0].x
    for index, wavelength in enumerate(full_wavelength_list):
        if wavelength >= 818.0 and wavelength < 818.9:
            indices_to_use.append(index)
    fit_results_scandata = rsa_fit_results_scandata_list[0]  # just one entry
    rsa_vs_wavelength_fit_results = aggregate_and_process_scandata_into_dict(
                                fit_results_scandata,
                                process_fcn_list=[get_mean_and_std_dict,
                                                  get_linear_fit_results_dict],
                                xcoord_indices_to_use=indices_to_use,
                                xcoord_name="wavelengths")

     # plot w/ superimposed linear fit interpolated values, mean, std,
    field_name = "lifetime"
    plt.errorbar(rsa_vs_wavelength_fit_results['wavelengths'],
                 rsa_vs_wavelength_fit_results[field_name],
                 rsa_vs_wavelength_fit_results[field_name + '_error'],
                 label=field_name, fmt='bd')
    mean = rsa_vs_wavelength_fit_results[field_name + '_mean']
    std = rsa_vs_wavelength_fit_results[field_name + '_std']
    plt.plot([817.9, 819.0], np.ones((2))*mean, 'b-', label='mean value')
    plt.plot([817.9, 819.0], np.ones((2))*(mean + std), 'b:')
    plt.plot([817.9, 819.0], np.ones((2))*(mean - std), 'b:', label='+- one std')
    plt.plot([818.0], rsa_vs_wavelength_fit_results[field_name +
                                                    '_linear_fit_at_818.0nm'],
             'rd', label='linear fit interpolation')
    plt.plot([818.9], rsa_vs_wavelength_fit_results[field_name +
                                                    '_linear_fit_at_818.9nm'],
             'rd')
    plt.title(field_name)
    plt.xlim(817.9, 819.0)
    plt.legend(loc='best')

    # extract parameter averages/stds, since generally same in region:
    print('Pulse amplitude (AU): {:4g} +- {:4g}'.format(
            rsa_vs_wavelength_fit_results["pulse_amplitude_mean"],
            rsa_vs_wavelength_fit_results["pulse_amplitude_std"]))
    print('Spin lifetime (ps): {:4g} +- {:4g}'.format(
            rsa_vs_wavelength_fit_results["lifetime_mean"],
            rsa_vs_wavelength_fit_results["lifetime_std"]))
    print('g-factor: {:4g} +- {:4g}'.format(
            rsa_vs_wavelength_fit_results["gfactor_mean"],
            rsa_vs_wavelength_fit_results["gfactor_std"]))
    print('Oscillation period @200mT (ps): {:4g} +- {:4g}'.format(
            rsa_vs_wavelength_fit_results["osc_period_200mT_mean"],
            rsa_vs_wavelength_fit_results["osc_period_200mT_std"]))
    print('Oscillation period @300mT (ps): {:4g} +- {:4g}'.format(
            rsa_vs_wavelength_fit_results["osc_period_300mT_mean"],
            rsa_vs_wavelength_fit_results["osc_period_300mT_std"]))

    # OUTPUT:
    #Pulse amplitude (AU): 0.00774479 +- 0.00110934
    #Spin lifetime (ps): 20012.9 +- 380.617
    #g-factor: 0.431126 +- 0.00133554
    #Oscillation period @200mT (ps): 828.641 +- 2.56808
    #Oscillation period @300mT (ps): 552.427 +- 1.71205
    # plug pulse amplitude, lifetime, 300mT period into next model,
    # using param bounds as +- 4 standard deviations
    long_pulseamp_init = rsa_vs_wavelength_fit_results["pulse_amplitude_mean"]
#    long_pulseamp_bounds = (long_pulseamp_init - 4*np.std(rsa_pulse_amplitude),
#                            long_pulseamp_init + 4*np.std(rsa_pulse_amplitude))
    long_pulseamp_bounds = (0, 1)
    long_lifetime_init = rsa_vs_wavelength_fit_results["lifetime_mean"]
    rsa_lifetime_std = rsa_vs_wavelength_fit_results["lifetime_std"]
#    long_lifetime_bounds = (long_lifetime_init - 4*np.std(rsa_lifetime_std),
#                            long_lifetime_init + 4*np.std(rsa_lifetime_std))
    long_lifetime_bounds = (0.5*long_lifetime_init, 1.5*long_lifetime_init)
    long_period_init = rsa_vs_wavelength_fit_results["osc_period_300mT_mean"]
#    long_period_bounds = (long_period_init - 4*np.std(rsa_period_300mT),
#                          long_period_init + 4*np.std(rsa_period_300mT))
    long_period_bounds = (500, 600)

# %%  ANALYSIS OF DELAY SCANS VS WAVELENGTH

    # MODEL 2: Two decaying cosines w/ relative phase. Long lifetime fixed.
    #          Start beyond pump transient (<200ps lifetime) signal
    # params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
    #          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
    #          phase1, phase2, slope, offset
    trkr_model_0V = \
        IndependentSinusoidalSpinLifetimeModel(
            field_name="lockin2x",
            max_fcn_evals=20000,
            free_params=[False, True, True, True, False,
                         True, False, False, False,
                         True, True, True, True],
            initial_params = [40, 0.01, long_pulseamp_init,
                              10000, long_lifetime_init,
                              550, long_period_init, 0, 0,
                              2*np.pi/3, -2*np.pi/3, 0, 0],
            param_bounds = [(1, 1000), (0, 1), long_pulseamp_bounds,
                            (1, 1e6), long_lifetime_bounds,
                            (500, 600), long_period_bounds,
                            (0, 0), (0, 0),
                            (-np.pi, np.pi), (-2*np.pi, np.pi),
                            (-1e-4, 1e-4), (-0.01, 0.01)],
            error_thresholds=[None, None, None, None, None,
                              None, None, None, None,
                              None, None, None, None],
            fit_result_scan_coord="Wavelength",
            excluded_intervals=[[-15, 400]],
#            excluded_intervals=[[-15, 400], [7000, 15000]],
            BField=0.3)
    
    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "WavelengthDependence_TRKR_300mT")
    fixed_uncertainty = 1e-3  # manually set uncertainty of data set
    model = trkr_model_0V
    sort_key = "Voltage"  # one scandataset for each voltage

    scandata_list = list(fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))

    # ---------------
    # TEMPORARY, FOR SPEED:
#    scandata_list = scandata_list[0:4]
    # ---------------

    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    field_name = model.field_name  # for future use

    # Change drift velocity based on voltage. Assumes each set has same
    # voltage for every scan!
    for scandataset in scandataset_list:
        voltage_list = [scandata.info['Voltage']
                        for scandata in scandataset.scandata_list]
        set_voltage = voltage_list[0]
        if any(voltage != set_voltage for voltage in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            # Drift velocity per voltage
            drift_velocity_per_volt = 2e-3  # in um/(ps*V)
            drift_velocity = drift_velocity_per_volt*set_voltage  # in um/ps
            scandataset.model.free_params[7] = False  # drift_velocity1
            scandataset.model.free_params[8] = False  # drift_velocity2
            scandataset.model.initial_params[7] = drift_velocity
            scandataset.model.initial_params[8] = drift_velocity

            guess_period_long = long_period_init + 2.7 * abs(set_voltage)
            scandataset.model.free_params[6] = True  # let it fit from here
            scandataset.model.initial_params[6] = guess_period_long

    analyzer2 = ScanDataSetsAnalyzer(scandataset_list)

    # since now sorted by "Voltage", can change to electric field as well as
    # fix improper "StageZ" 2nd coord, change units to V/cm anyway
    for scandataset in analyzer2.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.info['Voltage']*20
            scandata.info['MiddleScanType'] = 'Electric Field (V/cm)'
            scandata.info['MiddleScanCoord'] = field

#    # smooth over data with a 40ps wide gaussian convolution filter
#    analyzer2.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    gaussian_width=40)
    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
    analyzer2.apply_transform_to_all_scandata(
                                    get_gaussian_smoothed_scandata,
                                    field_names_to_process=[field_name],
                                    gaussian_width=600,
                                    edge_handling='reflect',
                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer2.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
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
        analyzer2.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=False)


# %% OVERVIEW OF FITS
    field_name = 0  # x:delay, y:lockin2x
    plot_trkr_fit_scandata(field_name, fit_trkr_scandata_list[4:8])
    scandata_list = fit_trkr_scandata_list


# %%
    field_name = 1  # x:wavelength, y:short lifetime
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


# %% MODEL 2 EXTRACTED DATA REVIEW
    # First extract all data from fit results into a comprehensive dict
    # This gives us an estimate of several factors easily looked up in dict
    # as a function of wavelength, and leaving out "bad fits":
    # - "wavelengths": the wavelengths of which the other fields are a function
    # - "long/short_amplitude": spin amplitude generated per pulse
    # - "long/short_lifetime": spin lifetime
    # - "long/short_osc_period": precession period, wavelength-dependent!?
    # - "long/short_drift_velocity": generally fixed based on known mobility
    # - "long/short_phase": relative phase between long & short pulses
    # Furthermore, these flags can be added:
    # "_mean" brings a scalar: the mean of the values across wavelengths
    # "_std" brings a scalar: the std around said mean value across wavelengths
    # "_error" brings up the error from the diagonal elements of 
    #           the fit covariance matrix for each wavelength, and
    # "_linear_fit_at_818.0nm" [or 818.9nm] brings the linear interpolation
    #           estimate for the value at those two wavelengths.

    # First, trim indices to those in the region of interest, namely wavelength
    # less than 818.9, as above that plot results differ significantly, and
    # later datasets will be at 818.0 and 818.9, and finally the 818.9 data
    # set is too low quality to use its results
    indices_to_use = []
    full_wavelength_list = trkr_fit_results_scandata_list[0].x
    for index, wavelength in enumerate(full_wavelength_list):
        if wavelength >= 818.0 and wavelength <= 818.9:
            indices_to_use.append(index)
    trkr_vs_wavelength_fit_results_0Vcm = {}
    trkr_vs_wavelength_fit_results_15Vcm = {}
    for i in range(2):  # for 0V/cm, 15V/cm data:
        if i == 0:
            fit_results_dict = trkr_vs_wavelength_fit_results_0Vcm
        else:
            fit_results_dict = trkr_vs_wavelength_fit_results_15Vcm
        fit_results_scandata = trkr_fit_results_scandata_list[i]
        fit_results_dict.update(
            aggregate_and_process_scandata_into_dict(
                                fit_results_scandata,
                                process_fcn_list=[get_mean_and_std_dict,
                                                  get_linear_fit_results_dict],
                                xcoord_indices_to_use=indices_to_use,
                                xcoord_name="wavelengths"))

         # plot w/ superimposed linear fit interpolated values, mean, std,
        plt.figure()
        field_name = "short_amplitude"
        if field_name + '_error' in fit_results_dict.keys():
            plt.errorbar(fit_results_dict['wavelengths'],
                         fit_results_dict[field_name],
                         fit_results_dict[field_name + '_error'],
                         label=field_name, fmt='bd')
        else:
            plt.plot(fit_results_dict['wavelengths'],
                         fit_results_dict[field_name],
                         'bd', label=field_name)
        mean = fit_results_dict[field_name + '_mean']
        std = fit_results_dict[field_name + '_std']
        plt.plot([817.9, 819.0], np.ones((2))*mean, 'b-', label='mean value')
        plt.plot([817.9, 819.0], np.ones((2))*(mean + std), 'b:')
        plt.plot([817.9, 819.0], np.ones((2))*(mean - std), 'b:', label='+- one std')
        plt.plot([818.0], fit_results_dict[field_name +
                                                        '_linear_fit_at_818.0nm'],
                 'rd', label='linear fit interpolation')
        plt.plot([818.9], fit_results_dict[field_name +
                                                        '_linear_fit_at_818.9nm'],
                 'rd')
        plt.xlim(817.9, 819.0)
        plt.legend(loc='best')
        if i == 0:
            plt.title('TRKR fit results, 0 V/cm')
        else:
            plt.title('TRKR fit results, 15 V/cm')

        # extract parameter averages/stds, since generally same in region:
        if i == 0:
            print("--- TRKR FIT RESULTS AT 0 V/cm ---")
        else:
            print("--- TRKR FIT RESULTS AT 15 V/cm ---")
        print('Short pulse amplitude (AU): {:4g} +- {:4g}'.format(
                fit_results_dict["short_amplitude_mean"],
                fit_results_dict["short_amplitude_std"]))
        print('Long pulse amplitude (AU): {:4g} +- {:4g}'.format(
                fit_results_dict["long_amplitude_mean"],
                fit_results_dict["long_amplitude_std"]))
        print('Short spin lifetime (ps): {:4g} +- {:4g}'.format(
                fit_results_dict["short_lifetime_mean"],
                fit_results_dict["short_lifetime_std"]))
        print('Long spin lifetime (ps): {:4g} +- {:4g}'.format(
                fit_results_dict["long_lifetime_mean"],
                fit_results_dict["long_lifetime_std"]))
        print('Short g-factor: {:4g} +- {:4g}'.format(
                fit_results_dict["short_gfactor_mean"],
                fit_results_dict["short_gfactor_std"]))
        print('Long g-factor: {:4g} +- {:4g}'.format(
                fit_results_dict["long_gfactor_mean"],
                fit_results_dict["long_gfactor_std"]))
        print('Short drift velocity: {:4g} +- {:4g}'.format(
                fit_results_dict["short_drift_velocity_mean"],
                fit_results_dict["short_drift_velocity_std"]))
        print('Long drift velocity: {:4g} +- {:4g}'.format(
                fit_results_dict["long_drift_velocity_mean"],
                fit_results_dict["long_drift_velocity_std"]))

    # OUTPUT:
    # short pulse amplitude: assume wavelength-dependent, by some scaling factor
    short_pulseamp_init_818_0_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "short_amplitude_linear_fit_at_818.0nm"]
    short_pulseamp_init_818_9_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "short_amplitude_linear_fit_at_818.9nm"]
    short_pulseamp_init_818_0_15Vcm = \
        trkr_vs_wavelength_fit_results_15Vcm[
                                    "short_amplitude_linear_fit_at_818.0nm"]
    short_pulseamp_init_818_9_15Vcm = \
        trkr_vs_wavelength_fit_results_15Vcm[
                                    "short_amplitude_linear_fit_at_818.9nm"]
    short_pulseamp_bounds_818_0 = (0, 1)
    short_pulseamp_bounds_818_9 = (0, 1)

    long_pulseamp_init_818_0_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "long_amplitude_linear_fit_at_818.0nm"]
    long_pulseamp_init_818_9_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "long_amplitude_linear_fit_at_818.9nm"]
    long_pulseamp_init_818_0_15Vcm = \
        trkr_vs_wavelength_fit_results_15Vcm[
                                    "long_amplitude_linear_fit_at_818.0nm"]
    long_pulseamp_init_818_9_15Vcm = \
        trkr_vs_wavelength_fit_results_15Vcm[
                                    "long_amplitude_linear_fit_at_818.9nm"]
    long_pulseamp_bounds_818_0 = (0, 1)
    long_pulseamp_bounds_818_9 = (0, 1)

    # short lifetime: 0V/cm data is linear, but 15V/cm data is unreliable due
    #                 to first data point, and rest are ~flat. Use average
    #                 of those 3 points (indices 1,2,3) for both wavelengths -
    #                 just above the low end of 0V/cm curve vs wavelength.
    #                 Of course, we don't expect lifetime to change vs
    #                 wavelength anyway...so let's use the 818.9nm data since
    #                 since that's where this species is seen strongest.
    short_lifetime_init_818_0_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "short_lifetime_linear_fit_at_818.9nm"]
    short_lifetime_init_818_9_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "short_lifetime_linear_fit_at_818.9nm"]
    short_lifetime_init_818_0_15Vcm = \
        np.mean(trkr_vs_wavelength_fit_results_15Vcm["short_lifetime"][1:3])
    short_lifetime_init_818_9_15Vcm = \
        np.mean(trkr_vs_wavelength_fit_results_15Vcm["short_lifetime"][1:3])
    short_lifetime_bounds_818_0 = (1, 1e9)
    short_lifetime_bounds_818_9 = (1, 1e9)

    # long lifetime: The starting values from RSA of just under 20ns never
    #                changed during TRKR iteration when not fixed, despite the
    #                RSA values definitely changing vs. wavelength. Presumably,
    #                RSA is more sensitive to changes in lifetime, so the TRKR
    #                data was unhelpful. Of course, we don't expect lifetime
    #                to change vs wavelength anyway...let's use the 818.0nm
    #                data since that's seemingly not obscured by other species.
    #                For 15V/cm data...just use the same, I suppose, since the
    #                TRKR data isn't distinguishable any more than for 0V/cm.
    long_lifetime_init_818_0_0Vcm = \
        rsa_vs_wavelength_fit_results["lifetime_linear_fit_at_818.0nm"]
#        trkr_vs_wavelength_fit_results_0Vcm["long_lifetime_linear_fit_at_818.0nm"]
    long_lifetime_init_818_9_0Vcm = \
        rsa_vs_wavelength_fit_results["lifetime_linear_fit_at_818.0nm"]
#        trkr_vs_wavelength_fit_results_0Vcm["long_lifetime_linear_fit_at_818.9nm"]
    long_lifetime_init_818_0_15Vcm = \
        rsa_vs_wavelength_fit_results["lifetime_mean"]
#        trkr_vs_wavelength_fit_results_15Vcm["long_lifetime_linear_fit_at_818.0nm"]
    long_lifetime_init_818_9_15Vcm = \
        rsa_vs_wavelength_fit_results["lifetime_mean"]
#        trkr_vs_wavelength_fit_results_15Vcm["long_lifetime_linear_fit_at_818.9nm"]
    long_lifetime_bounds_818_0 = (1, 1e9)
    long_lifetime_bounds_818_9 = (1, 1e9)

    # short species g-factor: Use 818.9nm data - where that species has
    #                         the most influence on the measured results
    #                        Note we must convert to 200mT for 818.0nm scan.
    short_period_init_818_0_0Vcm = \
        1.5*trkr_vs_wavelength_fit_results_0Vcm[
                                    "short_osc_period_linear_fit_at_818.9nm"]
    short_period_init_818_9_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "short_osc_period_linear_fit_at_818.9nm"]
    short_period_init_818_0_15Vcm = \
        1.5*trkr_vs_wavelength_fit_results_15Vcm[
                                    "short_osc_period_linear_fit_at_818.9nm"]
    short_period_init_818_9_15Vcm = \
        trkr_vs_wavelength_fit_results_15Vcm[
                                    "short_osc_period_linear_fit_at_818.9nm"]
    short_period_bounds_818_0 = (400, 1000)
    short_period_bounds_818_9 = (400, 1000)

    # long species g-factor: Use 818.0nm data - where that species has
    #                        the most influence on the measured results.
    #                        Note we must convert to 200mT for 818.0nm scan.
    long_period_init_818_0_0Vcm = \
        1.5*trkr_vs_wavelength_fit_results_0Vcm[
                                    "long_osc_period_linear_fit_at_818.0nm"]
    long_period_init_818_9_0Vcm = \
        trkr_vs_wavelength_fit_results_0Vcm[
                                    "long_osc_period_linear_fit_at_818.0nm"]
    long_period_init_818_0_15Vcm = \
        1.5*trkr_vs_wavelength_fit_results_15Vcm[
                                    "long_osc_period_linear_fit_at_818.0nm"]
    long_period_init_818_9_15Vcm = \
        trkr_vs_wavelength_fit_results_15Vcm[
                                    "long_osc_period_linear_fit_at_818.0nm"]
    long_period_bounds_818_0 = (400, 1000)
    long_period_bounds_818_9 = (400, 1000)


# %%  ANALYSIS OF 818.0nm DELAY SCANS
    
    # MODEL 2: Two decaying cosines w/ relative phase, as before. This time
    #          no fixing of long lifetime or oscillation factor, and short
    #          lifetime signal is assumed to be very low by starting params
    # params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
    #          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
    #          phase1, phase2, slope, offset
    trkr_model_3 = \
        IndependentSinusoidalSpinLifetimeModel(
            field_name="lockin2x",
            max_fcn_evals=2000,
            free_params=[False, False, True, False, True,
                         False, True, False, False,
                         False, True, True, True],
            initial_params = [40,
                              short_pulseamp_init_818_0_0Vcm,
                              long_pulseamp_init_818_0_0Vcm,
                              short_lifetime_init_818_0_0Vcm,
                              long_lifetime_init_818_0_0Vcm,
                              short_period_init_818_0_0Vcm,
                              long_period_init_818_0_0Vcm,
#                              0, 0.01, 20260, 20260,
#                              808, 808,
                              0, 0,
                              -np.pi, np.pi/2, 0, 0],
            param_bounds = [(1, 1000),
                            short_pulseamp_bounds_818_0,
                            long_pulseamp_bounds_818_0,
                            short_lifetime_bounds_818_0,
                            long_lifetime_bounds_818_0,
                            short_period_bounds_818_0,
                            long_period_bounds_818_0,
#                            (0, 1), (0, 1), (1, 1e6), (1, 1e6),
#                            (750, 850), (750, 850),
                            (0, 0), (0, 0),
                            (-2*np.pi, np.pi), (-2*np.pi, np.pi),
                            (-1e-4, 1e-4), (-0.01, 0.01)],
            error_thresholds=[None, None, None, None, None,
                              None, None, None, None,
                              None, None, None, None],
            fit_result_scan_coord="Electric Field (V/cm)",  # look at results of fit vs field
            excluded_intervals=[[-15, 400]],
#            excluded_intervals=[[-15, 400], [7000, 15000]],
            BField=0.3)


    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
#    dirpath = ("C:\\Data\\160702\\" +
#               "delayscans_-6V_to_6V_noautocenter")
#               "delayscans_-6V_to_6V")
#    dirpath = ("C:\\Data\\august_data\\160901\\" +
#               "DelayScansNewWavelength")  # requires other assumptions in model! do seperately below...
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "LowV_818.0nm_WavelengthDependence_TRKRvsV_200mT")
#               "GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "WavelengthDependence_TRKR_300mT")
#               "WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    fixed_uncertainty = 1e-3  # manually set uncertainty of data set
    model = trkr_model_3
    sort_key = "Voltage"  # just one ScanDataSet

    scandata_list = list(fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))
    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)

    # Change drift velocity based on voltage. Assumes each set has same
    # voltage for every scan!
    field_list = []
    model_initial_params_list = []
    for scandataset in scandataset_list:
        voltage_list = [scandata.info['Voltage']
                        for scandata in scandataset.scandata_list]
        set_voltage = voltage_list[0]
        if any(voltage != set_voltage for voltage in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            # Drift velocity per voltage
            drift_velocity_per_volt = 2e-3  # in um/(ps*V)
            drift_velocity = drift_velocity_per_volt*set_voltage  # in um/ps
            scandataset.model.free_params[7] = False  # drift_velocity1
            scandataset.model.free_params[8] = False  # drift_velocity2
            scandataset.model.initial_params[7] = drift_velocity
            scandataset.model.initial_params[8] = drift_velocity

            # Other parameters: set starting guess as linear
            # interpolation between 0V/cm and 15V/cm value
            current_field = set_voltage * 20
            abs_field = abs(current_field)
            if abs_field > 15:
                abs_field = 15
            for index in range(1, 7):
                if index == 1:
                    value_at_0Vcm = short_pulseamp_init_818_0_0Vcm
                    value_at_15Vcm = short_pulseamp_init_818_0_15Vcm
                elif index == 2:
                    value_at_0Vcm = long_pulseamp_init_818_0_0Vcm
                    value_at_15Vcm = long_pulseamp_init_818_0_15Vcm
                elif index == 3:
                    value_at_0Vcm = short_lifetime_init_818_0_0Vcm
                    value_at_15Vcm = short_lifetime_init_818_0_15Vcm
                elif index == 4:
                    value_at_0Vcm = long_lifetime_init_818_0_0Vcm
                    value_at_15Vcm = long_lifetime_init_818_0_15Vcm
                elif index == 5:
                    value_at_0Vcm = short_period_init_818_0_0Vcm
                    value_at_15Vcm = short_period_init_818_0_15Vcm
                elif index == 6:
                    value_at_0Vcm = long_period_init_818_0_0Vcm
                    value_at_15Vcm = long_period_init_818_0_15Vcm
                value_per_field = (value_at_15Vcm - value_at_0Vcm)/15
                new_value = value_at_0Vcm + value_per_field * abs_field
                scandataset.model.initial_params[index] = new_value
#                print('---')
#                print('0V/cm value: {}'.format(value_at_0Vcm))
#                print('15V/cm value: {}'.format(value_at_15Vcm))
#                print('@{}V/cm, new value: {}'.format(current_field, new_value))
            if current_field not in field_list:
                field_list.append(current_field)
                current_initial_params = \
                    dict(zip(scandataset.model.model_params,
                             scandataset.model.initial_params))
                model_initial_params_list.append(current_initial_params)

    analyzer4 = ScanDataSetsAnalyzer(scandataset_list)

    # since now sorted by "Voltage", can change to electric field as well as
    # fix improper "StageZ" 2nd coord, change units to V/cm anyway
    for scandataset in analyzer4.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.info['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field
                scaninfo['Electric Field (V/cm)'] = field

    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
#    analyzer4.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    field_names_to_process=[field_name],
#                                    gaussian_width=600,
#                                    edge_handling='reflect',
#                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer4.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15,
                                             neg_time_weight_multiplier=1.0)
                        #  low multiplier, not enough data to avoid bad fits

    # scandatasets don't share models, can't multiprocess in this version:
    analyzer4.fit_all_scandata_to_model(multiprocessing=False)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer4.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetime(threshold=999999))
#    analyzer4.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    fit_trkr_scandata_list_3 = \
        analyzer4.collapse_to_scandata_list(filtered=False)

    # regroup scandatasets into one big set to get single collapsed scandata:
    analyzer4.regroup_scandatasets(new_model=model,
                                   sort_key=None)
    trkr_fit_results_scandata_list_3 = \
        analyzer4.collapse_to_model_fit_scandata_list(
                                                    new_scan_type="[Unknown]",
                                                    filtered=False)

    # collapsed list: jiggle phase to keep phase continuous via np.unwrap()
    trkr_fit_results_scandata_list_3 = [get_continuous_phase_scandata(scandata)
                                for scandata in trkr_fit_results_scandata_list_3]


# %%
    field_name = 0  # lockin2x
    for scandata in fit_trkr_scandata_list_3[6:7]:
       plot_scandata(scandata, field_name, model=model, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    if scandata.get_field_fitdata(field_name) is not None:
        plt.text(8000, 0,
                 "Last fit lifetime: {:.3f}ns\n                         +-{:.3f}ns".format(
                     scandata.get_field_fitdata(field_name).fitparams[3]/1000,
                     scandata.get_field_fitdata(field_name).fitparamstds[3]/1000))
    plt.show()

    print('Voltage: {}'.format(scandata.info['Voltage']))
    if scandata.get_field_fitdata(field_name) is not None:
        print(scandata.get_field_fitdata(field_name).fitparams)

    scandata_list = fit_trkr_scandata_list_3


# %%
    field_names = [5]  # dataseries: x:field, y:short/long phase
    plt.figure()
    plt.hold(True)
    for field_name in field_names:
        for scandata in trkr_fit_results_scandata_list_3[:]:
            fit_param_name = scandata.fields[field_name]
            plot_scandata(scandata, field_name,
                          fmt="d", label=fit_param_name)
        try:
            initial_vals = [initial_params[fit_param_name]
                            for initial_params in model_initial_params_list]
            plt.plot(field_list, initial_vals, 'd:',
                     label="{}: Initial Guesses".format(fit_param_name))
        except KeyError:
            print("No fit parameter by name " + fit_param_name +
                  " found in model. Initial guess fit omitted")

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Electric Field (V/cm)")
    plt.ylabel("")  
    plt.legend(loc='best')
    plt.show()


## %% OVERVIEW OF FITS
#    field_name = 0  # lockin2x
#    for scandata in fit_trkr_scandata_list_3[:]:
#       plot_scandata(scandata, field_name, fmt="bd")
#    plt.xlabel("Delay (ps)")
#    plt.ylabel("Kerr Rotation (AU)")
#    if scandata.get_field_fitdata(field_name) is not None:
#        plt.text(8000, 0,
#                 "Last fit lifetime: {}ns\n     +={}ns".format(
#                     scandata.get_field_fitdata(field_name).fitparams[3]/1000,
#                     scandata.get_field_fitdata(field_name).fitparamstds[3]/1000))
#    plt.show()
#
#    print('Voltage: {}'.format(scandata.info['Voltage']))
#    if scandata.get_field_fitdata(field_name) is not None:
#        print(scandata.get_field_fitdata(field_name).fitparams)
#
#    scandata_list = fit_trkr_scandata_list_2
#
#
## %%
#    field_names = [5]  # dataseries: x:field, y:short/long phase
#    plt.figure()
#    plt.hold(True)
#    for field_name in field_names:
#        for scandata in trkr_fit_results_scandata_list_3[:]:
#            plot_scandata(scandata, field_name, fmt="d",
#                          label=scandata.fields[field_name])
#
#    # LABEL AND DISPLAY GRAPH
#    plt.xlabel("Applied Voltage (V)  |  Electric Field (V/500um)")
#    plt.ylabel("")  
#    plt.legend(loc='best')
#    plt.show()
#
#    if scandata.get_field_fitdata(field_name) is not None:
#        print(scandata.get_field_fitdata(field_name).fitparams)
#
#    scandata_list = collapsed_scandata_list


# %%  ANALYSIS OF 818.9nm DELAY SCANS
    # IMPORTANT RESULTS FOR MODEL CONSIDERATION:
    # one good RSA fit to long species, @818.0nm, no voltage applied:
    # lifetime: 20.26ns +- 0.08ns,
    # gfactor ~= 5.41e8 Ts (have to grab again) -> g=.00615? uhh
    # \-> osc_period @300mT: 554.620ps +- .045ps
    # field offset: -0.7565mT +- 0.0056mT
    # note fitting model 2 shows ~1% increase in osc_period moving to +-2V!

    # MODEL 3: Two decaying cosines w/ relative phase. Long lifetime fixed.
    #          Start beyond pump transient (<200ps lifetime) signal
    # params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
    #          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
    #          phase1, phase2, slope, offset
    trkr_model_2_0V = \
        IndependentSinusoidalSpinLifetimeModel(
            field_name="lockin2x",
            max_fcn_evals=20000,
            free_params=[False, True, True, True, False,
                         True, True, False, False,
                         True, True, True, True],
            initial_params = [40,
                              short_pulseamp_init_818_9_0Vcm,
                              long_pulseamp_init_818_9_0Vcm,
                              short_lifetime_init_818_9_0Vcm,
                              long_lifetime_init_818_9_0Vcm,
                              short_period_init_818_9_0Vcm,
                              long_period_init_818_9_0Vcm,
#                              0.015, 0.008, 4000, 20260,
#                              555.5, 554.62,
                              0, 0,
                              0, np.pi, 0, 0],
            param_bounds = [(1, 1000),
                              short_pulseamp_bounds_818_9,
                              long_pulseamp_bounds_818_9,
                              short_lifetime_bounds_818_9,
                              long_lifetime_bounds_818_9,
                              short_period_bounds_818_9,
                              long_period_bounds_818_9,
#                            (0, 0.1), (0.003, 0.1), (1, 1e4), (1e3, 1e7),
#                            (500, 600), (500, 600),
                            (0, 0), (0, 0),
                            (-4*np.pi, 4*np.pi), (-4*np.pi, 4*np.pi),
                            (-1e-4, 1e-4), (-0.01, 0.01)],
            error_thresholds=[None, None, None, None, None,
                              None, None, None, None,
                              None, None, None, None],
            fit_result_scan_coord="Electric Field (V/cm)",  # look at results of fit vs field
            excluded_intervals=[[-15, 400]],
#            excluded_intervals=[[-15, 400], [7000, 15000]],
            BField=0.3)


    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
#    dirpath = ("C:\\Data\\160702\\" +
#               "delayscans_-6V_to_6V_noautocenter")
#               "delayscans_-6V_to_6V")
#    dirpath = ("C:\\Data\\august_data\\160901\\" +
#               "DelayScansNewWavelength")  # requires other assumptions in model! do seperately below...
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "LowV_818.0nm_WavelengthDependence_TRKRvsV_200mT")
#               "WavelengthDependence_TRKR_300mT")
#               "WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    fixed_uncertainty = 1e-3  # manually set uncertainty of data set
    model = trkr_model_2_0V
    sort_key = "Voltage"  # one scandataset for each voltage

    scandata_list = list(fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))

    # ---------------
    # TEMPORARY, FOR SPEED:
#    scandata_list = scandata_list[5:6]
    # ---------------

    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    field_name = model.field_name  # for future use

    # Change drift velocity based on voltage. Assumes each set has same
    # voltage for every scan!
    field_list = []
    model_initial_params_list = []
    for scandataset in scandataset_list:
        voltage_list = [scandata.info['Voltage']
                        for scandata in scandataset.scandata_list]
        set_voltage = voltage_list[0]
        if any(voltage != set_voltage for voltage in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            # Drift velocity per voltage
            drift_velocity_per_volt = 2e-3  # in um/(ps*V)
            drift_velocity = drift_velocity_per_volt*set_voltage  # in um/ps
            scandataset.model.free_params[7] = False  # drift_velocity1
            scandataset.model.free_params[8] = False  # drift_velocity2
            scandataset.model.initial_params[7] = drift_velocity
            scandataset.model.initial_params[8] = drift_velocity

            # Other parameters: set starting guess as linear
            # interpolation between 0V/cm and 15V/cm value
            for index in range(1, 7):
                if index == 1:
                    value_at_0Vcm = short_pulseamp_init_818_9_0Vcm
                    value_at_15Vcm = short_pulseamp_init_818_9_15Vcm
                elif index == 2:
                    value_at_0Vcm = long_pulseamp_init_818_9_0Vcm
                    value_at_15Vcm = long_pulseamp_init_818_9_15Vcm
                elif index == 3:
                    value_at_0Vcm = short_lifetime_init_818_9_0Vcm
                    value_at_15Vcm = short_lifetime_init_818_9_15Vcm
                elif index == 4:
                    value_at_0Vcm = long_lifetime_init_818_9_0Vcm
                    value_at_15Vcm = long_lifetime_init_818_9_15Vcm
                elif index == 5:
                    value_at_0Vcm = short_period_init_818_9_0Vcm
                    value_at_15Vcm = short_period_init_818_9_15Vcm
                elif index == 6:
                    value_at_0Vcm = long_period_init_818_9_0Vcm
                    value_at_15Vcm = long_period_init_818_9_15Vcm
                value_per_field = (value_at_15Vcm - value_at_0Vcm)/15
                current_field = set_voltage * 20
                abs_field = abs(current_field)
                if abs_field > 15:
                    abs_field = 15
                new_value = value_at_0Vcm + value_per_field * abs_field
                scandataset.model.initial_params[index] = new_value
#                print('---')
#                print('0V/cm value: {}'.format(value_at_0Vcm))
#                print('15V/cm value: {}'.format(value_at_15Vcm))
#                print('@{}V/cm, new value: {}'.format(current_field, new_value))
                if current_field not in field_list:
                    field_list.append(current_field)
                    current_initial_params = scandataset.model.initial_params
                    model_initial_params_list.append(current_initial_params)

    analyzer3 = ScanDataSetsAnalyzer(scandataset_list)


    # since now sorted by "Voltage", can change to electric field as well as
    # fix improper "StageZ" 2nd coord, change units to V/cm anyway
    for scandataset in analyzer3.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.info['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field
                scaninfo['Electric Field (V/cm)'] = field

    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
#    analyzer3.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    field_names_to_process=[field_name],
#                                    gaussian_width=600,
#                                    edge_handling='reflect',
#                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer3.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15,
                                             neg_time_weight_multiplier=5.0)
    

    # scandatasets don't share models, can't multiprocess in this version:
    analyzer3.fit_all_scandata_to_model(multiprocessing=False)


    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer3.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetime(threshold=999999))
#    analyzer3.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    fit_trkr_scandata_list_2 = \
        analyzer3.collapse_to_scandata_list(filtered=False)

    # regroup scandatasets into one big set to get single collapsed scandata:
    analyzer3.regroup_scandatasets(new_model=model,
                                  sort_key=None)
    trkr_fit_results_scandata_list_2 = \
        analyzer3.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=False)

    # collapsed list: jiggle phase to keep phase continuous via np.unwrap()
    trkr_fit_results_scandata_list_2 = [get_continuous_phase_scandata(scandata)
                               for scandata in trkr_fit_results_scandata_list_2]




