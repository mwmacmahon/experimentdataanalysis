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
LASER_REPRATE = 13160  # ps period



# PLOTTING FUNCTIONS
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


# %%
if __name__ == '__main__':
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
            ignore_weights=False)
    excluded_time_intervals = [[0, 400],
                               [LASER_REPRATE - 15, 15000]]  # put here instead

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
    #                   pump_probe_dist, wavelength, temperature,
    #                   runID, index_in_run)
    timefield = 'scancoord'
    yfield = 'lockin2x'  # or '= model.field_name'
    fvec_fields = [('measurement', yfield),
                   ('time', timefield),
                   ('efield', "Electric Field (V/cm)"),
                   ('bfield', "Magnetic Field (mT)"),
                   ('pump_probe_dist', "Pump-Probe Distance (um)"),
                   ('wavelength', "Wavelength (nm)"),
                   ('temperature', "SetTemperature (K)"),
                   ('runID', None),
                   ('index_in_run', None)]

    # Fetch scandata, start with one big scandataset
    scandata_list = \
        list(fetch_dir_as_unfit_scandata_iterator(
                    directorypath=dirpath,
                    yfield=yfield,
                    yfield_error_val=fixed_uncertainty,
                    parsing_keywordlists=filepath_parsing_keyword_lists))
    
    for scandata in scandata_list:
        scandata.info['Wavelength (nm)'] = 818.9  # since filename reading iffy
        scandata.info['SetTemperature (K)'] = 30.0  # UNLESS HAVE ACTUAL TEMP
        gaussian_smooth_scandata(scandata,
                                 fields_to_process=[yfield],
                                 gaussian_width=600,
                                 edge_handling='reflect',
                                 subtract_smoothed_data_from_original=True)
        make_scandata_time_delay_positive(scandata,
                                          zero_delay_offset=-15,
                                          neg_time_weight_multiplier=5.0)


# %%
    # feature_vectors: (measurement, delaytime, efield, bfield,
    #                   pump_probe_dist, wavelength, temperature,
    #                   runID, index_in_run)
    feature_vector_array = None
    measurement_error_array = None
    for scandata_index, scandata in enumerate(scandata_list):
        indices_to_use_mask = np.logical_and.reduce(
            np.vstack([np.logical_or(getattr(scandata, timefield) < t_min,
                                     getattr(scandata, timefield) > t_max)
                       for t_min, t_max in excluded_time_intervals]))
        scandata_nfvecs = np.count_nonzero(indices_to_use_mask)
        scandata_feature_vector_array = np.zeros((scandata_nfvecs,
                                                  len(fvec_fields)))
        for fvec_ind, (element_name, field_name) in enumerate(fvec_fields):
            if field_name is None:  # meaning field value not found any dict
                continue
            try:
                scandata_feature_vector_array[:, fvec_ind] = \
                    getattr(scandata, field_name)[indices_to_use_mask]
            except AttributeError:
                try:
                    scandata_feature_vector_array[:, fvec_ind] = \
                        scandata.info[field_name] * np.ones(scandata_nfvecs)
                except KeyError:
                    raise KeyError("unable to find " + element_name +
                                   " in fields or info dict of scandata")
        # last two feature_vector indices: runID, index_in_run
        scandata_feature_vector_array[:, -2] = \
            scandata_index * np.ones(scandata_nfvecs)
        scandata_feature_vector_array[:, -1] = np.arange(scandata_nfvecs)
        scandata_measurement_error = \
            getattr(scandata, fvec_fields[0][1] + '_error')[indices_to_use_mask]
        if feature_vector_array is not None:
            feature_vector_array = np.vstack([feature_vector_array,
                                              scandata_feature_vector_array])
            if measurement_error_array is not None:  # only append if not None
                measurement_error_array = np.hstack([measurement_error_array,
                                                     scandata_measurement_error])
                
        else:
            feature_vector_array = scandata_feature_vector_array
            measurement_error_array = scandata_measurement_error
    fvec_scandata = ScanData(['feature_vector', 'measurement',
                              'runID', 'index_in_run'],
                             [feature_vector_array,
                              feature_vector_array[:, 0],
                              feature_vector_array[:, -2],
                              feature_vector_array[:, -1]])
    fvec_scandata.yerr = measurement_error_array


# %%
#    fitdata = generic_curve_fit(feature_vector_array,
#                                feature_vector_array[:, 0], None,
#                                model.fitfunction,
#                                model.free_params, model.initial_params, 
#                                model.param_bounds, model.max_fcn_evals)
    fitdata = scandata_model_fit(fvec_scandata, model)

    for param_label, param, param_std in zip(fitdata.fitparamlabels,
                                             fitdata.fitparams,
                                             fitdata.fitparamstds):
        print("{}: {:.6g} +- {:.6g}".format(param_label, param, param_std))




# %%


# trkr data:

# WITHOUT SMOOTHING
#num_pulses: 40 +- 0
#pulse_amplitude: 0.0228067 +- 2.78938e-05
#species_amp_ratio: 1.78913 +- 0.00354041
#lifetime1: 20000 +- 0
#lifetime2: 2011.62 +- 4.93839
#gfactor: 0.439998 +- 4.4904e-06
#mobility: 0.0001 +- 0
#slope: 0 +- 0
#offset: 0 +- 0
#
#fitdata.meansquarederror
#Out[25]: 7.2138053492442867e-08


# WITH SMOOTHING
#num_pulses: 40 +- 0
#pulse_amplitude: 0.0228049 +- 2.79123e-05
#species_amp_ratio: 1.79655 +- 0.00355929
#lifetime1: 20000 +- 0
#lifetime2: 2002.31 +- 4.90805
#gfactor: 0.440006 +- 4.50148e-06
#mobility: 0.0001 +- 0
#slope: 0 +- 0
#offset: 0 +- 0
#
#fitdata.meansquarederror
#Out[29]: 6.0788889090491028e-08































#
#
#
## %%
##    # FOR TESTING: CUT SIZE OF DATA DOWN
##    original_scandata_list = scandata_list
##    scandata_list = scandata_list[::10]
#
#    # one scandataset for each (b-field value, pump-probe-pos)
#    sort_keys = ["Magnetic Field (mT)", "Pump-Probe Distance (um)"]
#    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_keys)
#
#    for scandataset in scandataset_list:
#        # SET OSCILLATION PERIOD MODEL PARAMETER INITIAL GUESS
#        b_field_list = [scandata.info['Magnetic Field (mT)']
#                        for scandata in scandataset.scandata_list]
#        set_b_field = b_field_list[0]
#        if any(b_field != set_b_field for b_field in b_field_list):
#            print("No common magnetic field in ScanDataSet, " +
#                  "cannot set exact oscillation period in fit model.")
#        else:
#            scandataset.model.b_field = set_b_field
#            assumed_gfactor = 0.44
#            scandataset.model.initial_params[5] = \
#                (GFACTORCONSTANT * assumed_gfactor * set_b_field)**(-1)
#        # SET DRIFT VELOCITY MODEL PARAMETER FIXED VALUE
#        e_field_list = [scandata.info['Electric Field (V/cm)']
#                        for scandata in scandataset.scandata_list]
#        set_e_field = e_field_list[0]
#        if any(e_field != set_e_field for e_field in e_field_list):
#            print("No common voltage in ScanDataSet, " +
#                  "cannot set exact drift velocity in fit model.")
#        else:
#            # Drift velocity per voltage
#            mobility_coeff = 1e-4  # in um/(ps*V/cm)
#            drift_velocity = mobility_coeff * set_e_field  # in um/ps
#            scandataset.model.free_params[6] = False  # drift_velocity
#            scandataset.model.initial_params[6] = drift_velocity
#        # SET PUMP PROBE DISTANCE MODEL PARAMETER FIXED VALUE
#        distance_list = [scandata.info['Pump-Probe Distance (um)']
#                         for scandata in scandataset.scandata_list]
#        set_distance = distance_list[0]
#        if any(distance != set_distance for distance in distance_list):
#            print("No common pump probe distance in ScanDataSet, " +
#                  "cannot set exact drift velocity in fit model.")
#        else:
#            # Drift velocity per voltage
#            scandataset.model.free_params[7] = False  # probe_pos
#            scandataset.model.initial_params[7] = set_distance
#
##    # smooth over data with a 40ps wide gaussian convolution filter
##    for scandataset in scandataset_list:
##        scandataset.apply_transform_to_scandata(gaussian_smooth_scandata,
##                                                gaussian_width=40)
#
#    # drift subtraction:
#    # subtract from data: data times a 400ps wide gaussian convolution filter                        
#    for scandataset in scandataset_list:
#        scandataset.apply_transform_to_scandata(
#                                    gaussian_smooth_scandata,
#                                    fields_to_process=[yfield],
#                                    gaussian_width=600,
#                                    edge_handling='reflect',
#                                    subtract_smoothed_data_from_original=True)
#
#    # add 13160ps to all negative delay times
#    for scandataset in scandataset_list:
#        scandataset.apply_transform_to_scandata(
#                                    make_scandata_time_delay_positive,
#                                    zero_delay_offset=-15,
#                                    neg_time_weight_multiplier=5.0)
#
#    # scandatasets don't share models, can't multiprocess in this version:
#    fit_scandataset_list_to_model(scandataset_list, multiprocessing=False)
##    for scandataset in scandataset_list:
##        scandataset.purge_failed_fit_scandata()
#    fit_trkr_scandata_list = [scandata
#                              for scandataset in scandataset_list
#                              for scandata in scandataset.scandata_list
#                              if scandata.fitdata is not None]
#    trkr_fit_results_scandata_list = \
#        collapse_scandataset_to_model_fit_scandata_list(scandataset_list)
#
#
## %% OVERVIEW OF FITS
#    plot_trkr_fit_scandata(yfield, fit_trkr_scandata_list[25:30])
##    scandata_list = fit_trkr_scandata_list
#
## %%
#    param_name = "amplitude1"
#    plt.figure()
#    plt.hold(True)
#    for scandata in trkr_fit_results_scandata_list[:]:
#        plot_scandata(scandata, param_name, fmt=":bd",
#                      label="")
#
#    # LABEL AND DISPLAY GRAPH
#    plt.xlabel("")
#    plt.ylabel(param_name)
##    plt.legend(loc='best')
#    plt.title("")
#    plt.show()
#
##    scandata_list = trkr_fit_results_scandata_list
#
## %% GRAND MODEL FIT COMPARISON STUFF START
#
## %% nicely shows fits are great within 1 std dev
#    for scandata in fit_trkr_scandata_list[::10]:
#        fitdata = scandata.fitdata_lockin2x
#        freeindices = np.array(fitdata.freeparamindices)
#        fitparams = np.array(fitdata.fitparams)[freeindices]
#        fitparamstds = np.array(fitdata.fitparamstds)[freeindices]
#        paramlabels = np.array(fitdata.fitparamlabels)[freeindices]
#        covmat = scandata.fitdata_lockin2x.covariancematrix
#        if not np.all(np.linalg.eigvals(covmat) >= 0):
#            print("error, covariance matrix " +
#                  "not positive semidefinite, skipping scandata")
#            continue
#
#        distribution = np.random.multivariate_normal(fitparams, covmat, 100)
#        plot2dindices = [0, 2]
#
#        axes = plt.subplot(1,2,1)
##        num_steps = 5
##        xmin, ymin = (fitparams - 3 * fitparamstds)[plot2dindices]
##        xmax, ymax = (fitparams + 3 * fitparamstds)[plot2dindices]
##        xstep = (xmax - xmin) / (num_steps - 1)
##        ystep = (ymax - ymin) / (num_steps - 1)
##        x, y = np.mgrid[xmin:xmax+xstep:xstep, ymin:ymax+ystep:ystep]
##        param_mesh = np.empty(x.shape + (len(fitparams),))
##        param_mesh[:, :] = fitparams
##        param_mesh[:, :, plot2dindices[0]] = x
##        param_mesh[:, :, plot2dindices[1]] = y
##        rv = stats.multivariate_normal(fitparams, covmat, allow_singular=True)
##        plt.contourf(x, y, rv.pdf(param_mesh), 10)
#
#        axes.add_patch(plt.Rectangle((fitparams - fitparamstds)[plot2dindices],
#                                     *(2 * fitparamstds)[plot2dindices],
#                                     fill=False))
#        plt.plot(*fitparams[plot2dindices], 'ro', markersize=20)
#        axes.plot(*distribution[:, plot2dindices].T, 'bd')
##        plt.xlim([xmin, xmax])
##        plt.ylim([ymin, ymax])
#        plt.xlabel(paramlabels[plot2dindices[0]])
#        plt.ylabel(paramlabels[plot2dindices[1]])
#        plt.subplot(1,2,2)
#        plt.plot(*scandata.xy, 'bo')
#        for fit_param_set in distribution[:, :]:
#            yvals = fitdata.partialfcn(scandata.x, *fit_param_set)
#            if np.max(np.abs(yvals)) < 2 * np.max(np.abs(scandata.y)):
#                plt.plot(scandata.x, yvals, 'g:')
#
## %%
#    for scandata in fit_trkr_scandata_list[:80:10]:
#        fitdata = scandata.fitdata_lockin2x
#        freeindices = np.array(fitdata.freeparamindices)
#        fitparams = np.array(fitdata.fitparams)[freeindices]
#        fitparamstds = np.array(fitdata.fitparamstds)[freeindices]
#        paramlabels = np.array(fitdata.fitparamlabels)[freeindices]
#        covmat = scandata.fitdata_lockin2x.covariancematrix
#        if not np.all(np.linalg.eigvals(covmat) >= 0):
#            print("error, covariance matrix " +
#                  "not positive semidefinite, skipping scandata")
#            continue
#
#        plot2dindices = [0, 2]
#        distribution = np.random.multivariate_normal(fitparams, covmat, 20)
#        axes = plt.subplot(1,2,1)
#        axes.add_patch(plt.Rectangle((fitparams - fitparamstds)[plot2dindices],
#                                     *(2 * fitparamstds)[plot2dindices],
#                                     fill=False))
##        plt.plot(*fitparams[plot2dindices], 'ro', markersize=20)
#        axes.plot(*distribution[:, plot2dindices].T, 'bd')
#        plt.xlabel(paramlabels[plot2dindices[0]])
#        plt.ylabel(paramlabels[plot2dindices[1]])
#        plt.subplot(1,2,2)
#        plt.plot(*scandata.xy, 'bo')
#        for fit_param_set in distribution[:, :]:
#            yvals = fitdata.partialfcn(scandata.x, *fit_param_set)
#            if np.max(np.abs(yvals)) < 2 * np.max(np.abs(scandata.y)):
#                plt.plot(scandata.x, yvals, 'g:')

