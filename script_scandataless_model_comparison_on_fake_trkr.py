# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 13:34:40 2016

@author: Michael
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.filters as filters
import scipy.stats as stats

from experimentdataanalysis.analysis.dataclasses import ScanData
from experimentdataanalysis.analysis.scandataprocessing \
    import make_scandata_time_delay_positive, \
           make_scandata_phase_continuous, \
           gaussian_smooth_scandata, \
           process_scandata_fields, \
           generic_curve_fit, \
           scandata_model_fit
from experimentdataanalysis.analysis.scandatamodels \
    import FeatureVectorsTwoLifetimesOppositePhaseTRKRModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import sort_scandata_into_sets, \
           fit_scandataset_list_to_model, \
           collapse_scandataset_to_model_fit_scandata_list
from experimentdataanalysis.parsing.scandataparsing import \
        fetch_dir_as_unfit_scandata_iterator


#GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar


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
                  label="", fmt="-bd", fit_fmt="xr-"):
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
#    plt.figure()
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


def plot_trkr_fit_scandata(field_name, fit_trkr_scandata_list,
                           fmt="bd", fit_fmt="xr-"):
#    plt.figure()
    plt.hold(True)
    try:
        for scandata in fit_trkr_scandata_list[:]:
            plot_scandata(scandata, field_name,
                          fmt=fmt, fit_fmt=fit_fmt)
    except TypeError:
        scandata = fit_trkr_scandata_list
        plot_scandata(scandata, field_name,
                      fmt=fmt, fit_fmt=fit_fmt)
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    plt.show()
    print(scandata.get_field_fitdata(field_name).fitparams)


def plot_fit_param_scandata(field_names, fit_results_scandata_list,
                            fmt="bd", fit_fmt="xr-"):
    plt.figure()
    plt.hold(True)
    for field_name in field_names:
        try:
            for scandata in fit_results_scandata_list[:]:
                plot_scandata(scandata, field_name, label=field_name,
                              fmt=fmt, fit_fmt=fit_fmt)
        except TypeError:
            scandata = fit_results_scandata_list
            plot_scandata(scandata, field_name, label=field_name,
                          fmt=fmt, fit_fmt=fit_fmt)
    # LABEL AND DISPLAY GRAPH
    plt.legend(loc='best')
    plt.show()


# %%  ANALYSIS OF FAKE DELAY SCANS VS B FIELD

    # MODEL
    # params = num_pulses, pulse_amplitude, species_amp_ratio,
    #          lifetime1, lifetime2, gfactor, mobility,
    #          slope, offset
    feature_vector_model = \
        FeatureVectorsTwoLifetimesOppositePhaseTRKRModel(
            field_name="measurement",
            max_fcn_evals=10000,
            free_params=[False, True, True,
                         False, True, True, False,
                         False, False],
            initial_params=[40, 0.04, 2.0,
                            20000, 2000, 0.44, 1e-4,
                            0, 0],
            param_bounds=[(1,1000), (0, 1), (-100, 100),
                          (1, 1e9), (1, 1e9), (0.3, 0.6), (1e-6, 1),
                          (-1e-6, 1e-6), (-0.01, 0.01)],
            fit_result_scan_coord="Pump-Probe Distance (um)",
            excluded_intervals=None,  # below options not usable on f-vectors
#            excluded_intervals=[[-15, 400]],
#            excluded_intervals=[[-15, 400], [7000, 15000]],
            b_field=300,
            ignore_weights=True)
    excluded_time_intervals = [[-15, 400]]  # put here instead

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

    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
    dirpath = ("C:\\Data\\fake_data\\fake_trkr")
    fixed_uncertainty = 1e-3  # manually set uncertainty of data set
    model = feature_vector_model
    # feature_vectors: (measurement, delaytime, efield, bfield,
    #                   pump_probe_dist, wavelength, temperature)
    fvec_fields = [('measurement', "lockin2x"),
                   ('time', "scancoord"),
                   ('efield', "Electric Field (V/cm)"),
                   ('bfield', "Magnetic Field (mT)"),
                   ('pump_probe_dist', "Pump-Probe Distance (um)"),
                   ('wavelength', "Wavelength (nm)"),
                   ('temperature', "SetTemperature (K)")]
    timefield = 'scancoord'
    yfield = 'lockin2x'
    field_name = model.field_name  # for future use

    # Fetch scandata, start with one big scandataset
    scandata_list = \
        list(fetch_dir_as_unfit_scandata_iterator(
                    directorypath=dirpath,
                    yfield=yfield,
                    yfield_error_val=fixed_uncertainty,
                    parsing_keywordlists=filepath_parsing_keyword_lists))
    
    for scandata in scandata_list:
        scandata.info['Wavelength (nm)'] = 818.9
        scandata.info['SetTemperature (K)'] = 30.0


# %%
    # feature_vectors: (measurement, delaytime, efield, bfield,
    #                   pump_probe_dist, wavelength, temperature)
    feature_vector_array = None
    for scandata in scandata_list:
        indices_to_use_mask = np.logical_and.reduce(
            np.vstack([np.logical_or(getattr(scandata, timefield) < t_min,
                                     getattr(scandata, timefield) > t_max)
                       for t_min, t_max in excluded_time_intervals]))
        measurement = getattr(scandata, yfield)[indices_to_use_mask]
        scandata_feature_vector_array = np.zeros((len(measurement),
                                                  len(fvec_fields)))
        scandata_feature_vector_array[:, 0] = measurement
        for fvec_ind, (element_name, field_name) in enumerate(fvec_fields):
            try:
                scandata_feature_vector_array[:, fvec_ind] = \
                    getattr(scandata, field_name)[indices_to_use_mask]
            except AttributeError:
                try:
                    scandata_feature_vector_array[:, fvec_ind] = \
                        scandata.info[field_name] * np.ones(len(measurement))
                except KeyError:
                    raise KeyError("unable to find " + element_name +
                                   " in fields or info dict of scandata")
        if feature_vector_array is not None:
            feature_vector_array = np.vstack([feature_vector_array,
                                              scandata_feature_vector_array])
        else:
            feature_vector_array = scandata_feature_vector_array
    feature_vector_indices = np.arange(len(feature_vector_array)) + 1
    feature_vector_scandata = ScanData(['feature_vector',
                                        'measurement', 'index'],
                                       [feature_vector_array,
                                        feature_vector_array[:, 0],
                                        feature_vector_indices])


# %%
#    fitdata = generic_curve_fit(feature_vector_array,
#                                feature_vector_array[:, 0], None,
#                                model.fitfunction,
#                                model.free_params, model.initial_params, 
#                                model.param_bounds, model.max_fcn_evals)
    fitdata = scandata_model_fit(feature_vector_scandata, model)

    for param_label, param in zip(fitdata.fitparamlabels, fitdata.fitparams):
        print("{}: {}".format(param_label, param))





# %%
#    # FOR TESTING: CUT SIZE OF DATA DOWN
#    original_scandata_list = scandata_list
#    scandata_list = scandata_list[::10]

    # one scandataset for each (b-field value, pump-probe-pos)
    sort_keys = ["Magnetic Field (mT)", "Pump-Probe Distance (um)"]
    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_keys)

    for scandataset in scandataset_list:
        # SET OSCILLATION PERIOD MODEL PARAMETER INITIAL GUESS
        b_field_list = [scandata.info['Magnetic Field (mT)']
                        for scandata in scandataset.scandata_list]
        set_b_field = b_field_list[0]
        if any(b_field != set_b_field for b_field in b_field_list):
            print("No common magnetic field in ScanDataSet, " +
                  "cannot set exact oscillation period in fit model.")
        else:
            scandataset.model.b_field = set_b_field
            assumed_gfactor = 0.44
            scandataset.model.initial_params[5] = \
                (GFACTORCONSTANT * assumed_gfactor * set_b_field)**(-1)
        # SET DRIFT VELOCITY MODEL PARAMETER FIXED VALUE
        e_field_list = [scandata.info['Electric Field (V/cm)']
                        for scandata in scandataset.scandata_list]
        set_e_field = e_field_list[0]
        if any(e_field != set_e_field for e_field in e_field_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set exact drift velocity in fit model.")
        else:
            # Drift velocity per voltage
            mobility_coeff = 1e-4  # in um/(ps*V/cm)
            drift_velocity = mobility_coeff * set_e_field  # in um/ps
            scandataset.model.free_params[6] = False  # drift_velocity
            scandataset.model.initial_params[6] = drift_velocity
        # SET PUMP PROBE DISTANCE MODEL PARAMETER FIXED VALUE
        distance_list = [scandata.info['Pump-Probe Distance (um)']
                         for scandata in scandataset.scandata_list]
        set_distance = distance_list[0]
        if any(distance != set_distance for distance in distance_list):
            print("No common pump probe distance in ScanDataSet, " +
                  "cannot set exact drift velocity in fit model.")
        else:
            # Drift velocity per voltage
            scandataset.model.free_params[7] = False  # probe_pos
            scandataset.model.initial_params[7] = set_distance

#    # smooth over data with a 40ps wide gaussian convolution filter
#    for scandataset in scandataset_list:
#        scandataset.apply_transform_to_scandata(gaussian_smooth_scandata,
#                                                gaussian_width=40)

    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
    for scandataset in scandataset_list:
        scandataset.apply_transform_to_scandata(
                                    gaussian_smooth_scandata,
                                    fields_to_process=[field_name],
                                    gaussian_width=600,
                                    edge_handling='reflect',
                                    subtract_smoothed_data_from_original=True)

    # add 13160ps to all negative delay times
    for scandataset in scandataset_list:
        scandataset.apply_transform_to_scandata(
                                    make_scandata_time_delay_positive,
                                    zero_delay_offset=-15,
                                    neg_time_weight_multiplier=5.0)

    # scandatasets don't share models, can't multiprocess in this version:
    fit_scandataset_list_to_model(scandataset_list, multiprocessing=False)
#    for scandataset in scandataset_list:
#        scandataset.purge_failed_fit_scandata()
    fit_trkr_scandata_list = [scandata
                              for scandataset in scandataset_list
                              for scandata in scandataset.scandata_list
                              if scandata.fitdata is not None]
    trkr_fit_results_scandata_list = \
        collapse_scandataset_to_model_fit_scandata_list(scandataset_list)


# %% OVERVIEW OF FITS
    plot_trkr_fit_scandata(field_name, fit_trkr_scandata_list[25:30])
#    scandata_list = fit_trkr_scandata_list

# %%
    param_name = "amplitude1"
    plt.figure()
    plt.hold(True)
    for scandata in trkr_fit_results_scandata_list[:]:
        plot_scandata(scandata, param_name, fmt=":bd",
                      label="")

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("")
    plt.ylabel(param_name)
#    plt.legend(loc='best')
    plt.title("")
    plt.show()

#    scandata_list = trkr_fit_results_scandata_list

# %% GRAND MODEL FIT COMPARISON STUFF START

# %% nicely shows fits are great within 1 std dev
    for scandata in fit_trkr_scandata_list[::10]:
        fitdata = scandata.fitdata_lockin2x
        freeindices = np.array(fitdata.freeparamindices)
        fitparams = np.array(fitdata.fitparams)[freeindices]
        fitparamstds = np.array(fitdata.fitparamstds)[freeindices]
        paramlabels = np.array(fitdata.fitparamlabels)[freeindices]
        covmat = scandata.fitdata_lockin2x.covariancematrix
        if not np.all(np.linalg.eigvals(covmat) >= 0):
            print("error, covariance matrix " +
                  "not positive semidefinite, skipping scandata")
            continue

        distribution = np.random.multivariate_normal(fitparams, covmat, 100)
        plot2dindices = [0, 2]

        axes = plt.subplot(1,2,1)
#        num_steps = 5
#        xmin, ymin = (fitparams - 3 * fitparamstds)[plot2dindices]
#        xmax, ymax = (fitparams + 3 * fitparamstds)[plot2dindices]
#        xstep = (xmax - xmin) / (num_steps - 1)
#        ystep = (ymax - ymin) / (num_steps - 1)
#        x, y = np.mgrid[xmin:xmax+xstep:xstep, ymin:ymax+ystep:ystep]
#        param_mesh = np.empty(x.shape + (len(fitparams),))
#        param_mesh[:, :] = fitparams
#        param_mesh[:, :, plot2dindices[0]] = x
#        param_mesh[:, :, plot2dindices[1]] = y
#        rv = stats.multivariate_normal(fitparams, covmat, allow_singular=True)
#        plt.contourf(x, y, rv.pdf(param_mesh), 10)

        axes.add_patch(plt.Rectangle((fitparams - fitparamstds)[plot2dindices],
                                     *(2 * fitparamstds)[plot2dindices],
                                     fill=False))
        plt.plot(*fitparams[plot2dindices], 'ro', markersize=20)
        axes.plot(*distribution[:, plot2dindices].T, 'bd')
#        plt.xlim([xmin, xmax])
#        plt.ylim([ymin, ymax])
        plt.xlabel(paramlabels[plot2dindices[0]])
        plt.ylabel(paramlabels[plot2dindices[1]])
        plt.subplot(1,2,2)
        plt.plot(*scandata.xy, 'bo')
        for fit_param_set in distribution[:, :]:
            yvals = fitdata.partialfcn(scandata.x, *fit_param_set)
            if np.max(np.abs(yvals)) < 2 * np.max(np.abs(scandata.y)):
                plt.plot(scandata.x, yvals, 'g:')

# %%
    for scandata in fit_trkr_scandata_list[:80:10]:
        fitdata = scandata.fitdata_lockin2x
        freeindices = np.array(fitdata.freeparamindices)
        fitparams = np.array(fitdata.fitparams)[freeindices]
        fitparamstds = np.array(fitdata.fitparamstds)[freeindices]
        paramlabels = np.array(fitdata.fitparamlabels)[freeindices]
        covmat = scandata.fitdata_lockin2x.covariancematrix
        if not np.all(np.linalg.eigvals(covmat) >= 0):
            print("error, covariance matrix " +
                  "not positive semidefinite, skipping scandata")
            continue

        plot2dindices = [0, 2]
        distribution = np.random.multivariate_normal(fitparams, covmat, 20)
        axes = plt.subplot(1,2,1)
        axes.add_patch(plt.Rectangle((fitparams - fitparamstds)[plot2dindices],
                                     *(2 * fitparamstds)[plot2dindices],
                                     fill=False))
#        plt.plot(*fitparams[plot2dindices], 'ro', markersize=20)
        axes.plot(*distribution[:, plot2dindices].T, 'bd')
        plt.xlabel(paramlabels[plot2dindices[0]])
        plt.ylabel(paramlabels[plot2dindices[1]])
        plt.subplot(1,2,2)
        plt.plot(*scandata.xy, 'bo')
        for fit_param_set in distribution[:, :]:
            yvals = fitdata.partialfcn(scandata.x, *fit_param_set)
            if np.max(np.abs(yvals)) < 2 * np.max(np.abs(scandata.y)):
                plt.plot(scandata.x, yvals, 'g:')

