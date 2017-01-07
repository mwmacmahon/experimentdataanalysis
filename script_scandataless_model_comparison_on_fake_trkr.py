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
    import FeatureVectorsTwoLifetimesOppositePhaseTRKRModel, \
           FeatureVectorsIndependentSinusoidalTRKRModel
from experimentdataanalysis.parsing.scandataparsing import \
        fetch_dir_as_unfit_scandata_iterator


#GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period

# filename parsing pattern: if string in element[0] is found in filepath
# separated from other characters by '_'s, will record adjacent number
# and store in scandata's info dict under the key element[1], as a float
# if possible.
# e.g. TRKR_15Vcm_230mT.dat -> {"Electric Field (V/cm)": 15.0,
#                               "Magnetic Field (mT)": 230.0}
in_filepath_element_keyword_list = [("Vcm", "Electric Field (V/cm)"),
                                    ("mT", "Magnetic Field (mT)"),
                                    ("K", "Set Temperature (K)"),
                                    ("nm", "Wavelength (nm)"),
                                    ("ps", "Delay Time (ps)"),
                                    ("run", "RunIndex"),
                                    ("x", "Voltage (V)")]
# for this one, if element[0] found, next element stored w/ key element[1]
in_filepath_next_element_keyword_list = [("MirrorZ",
                                          "Pump-Probe Distance (um)")]
FILEPATH_PARSING_KEYWORD_LISTS = [[],
                                  in_filepath_next_element_keyword_list,
                                  in_filepath_element_keyword_list]


# HELPER FUNCTIONS
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


def scandata_list_to_feature_vector_scandata(scandata_list, fvec_fields,
                                             excluded_x_intervals=None):
    """
    [description goes here]
    Note: excluded_intervals are checked versus the scandata's xfield
    """
    feature_vector_array = None
    measurement_error_array = None
    for scandata_index, scandata in enumerate(scandata_list):
        if excluded_x_intervals is not None and len(excluded_x_intervals) > 0:
            indices_to_use_mask = np.logical_and.reduce(
                np.vstack([np.logical_or(scandata.x < x_min,
                                         scandata.x > x_max)
                           for x_min, x_max in excluded_x_intervals]))
        else:
            indices_to_use_mask = np.ones(len(scandata), dtype=np.bool)
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
        scandata_measurement_error = getattr(scandata,
                                             fvec_fields[0][1] + '_error',
                                             None)  # default to unknown error
        if scandata_measurement_error is not None:
            scandata_measurement_error = \
                scandata_measurement_error[indices_to_use_mask]
        if feature_vector_array is None:  # if examining first scandata
            feature_vector_array = scandata_feature_vector_array
            measurement_error_array = scandata_measurement_error
        else:
            feature_vector_array = np.vstack([feature_vector_array,
                                              scandata_feature_vector_array])
            if measurement_error_array is not None:  # only append if not None
                measurement_error_array = np.hstack([measurement_error_array,
                                                     scandata_measurement_error])

    fvec_scandata = ScanData(['feature_vector', 'measurement',
                              'runID', 'index_in_run'],
                             [feature_vector_array[:, :],
                              feature_vector_array[:, 0],
                              feature_vector_array[:, -2],
                              feature_vector_array[:, -1]])
    if measurement_error_array is not None:
        fvec_scandata.measurement_error = measurement_error_array
    return fvec_scandata


def split_scandata_into_training_and_test_sets(fvec_scandata, test_fraction,
                                               test_by_run=True):
    if test_by_run:
        num_all_targets = len(scandata_list)
    else:
        num_all_targets = len(fvec_scandata)
    num_test_targets = np.int(np.ceil(test_fraction * num_all_targets))
    test_targets = np.array([True] * num_test_targets +
                            [False] * (num_all_targets - num_test_targets))
    np.random.shuffle(test_targets)
    if test_by_run:
        test_indices = test_targets[fvec_scandata.runID.astype(np.int)]
    else:
        test_indices = test_targets
    training_indices = np.logical_not(test_indices)
    test_set_fvec_scandata = \
        ScanData(fvec_scandata.fields,
                 [getattr(fvec_scandata, field)[test_indices]
                  for field in fvec_scandata.fields])
    training_set_fvec_scandata = \
        ScanData(fvec_scandata.fields,
                 [getattr(fvec_scandata, field)[training_indices]
                  for field in fvec_scandata.fields])
    if fvec_scandata.yerr is not None:
        test_set_fvec_scandata.yerr = fvec_scandata.yerr[test_indices]
        training_set_fvec_scandata.yerr = fvec_scandata.yerr[training_indices]
    return training_set_fvec_scandata, test_set_fvec_scandata


# %%
if __name__ == '__main__':
# %%
    # GENERAL OPTIONS
    test_fraction = 0.1
    test_by_run = True  # TEST BY RUN VS INDIVIDUAL POINTS

# %%
    # DATA TO FIT
#    dirpath = ("C:\\Data\\fake_data\\fake_trkr")
#    dirpath = ("C:\\Data\\fake_data\\fake_rsa")
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "BestDataFromWavelengthDependence_TRKRvsV_300mT_" +
               "033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    # DATA OPTIONS:
    # manipulations on each scandata, e.g. filling missing info:
    def update_scandata_info(scandata):
        scandata.info['Wavelength (nm)'] = 818.9  # since filename reading iffy
        scandata.info['Set Temperature (K)'] = 30.0  # UNLESS HAVE ACTUAL TEMP
        if 'Voltage (V)' in scandata.info:
            efield = 20 * scandata.info['Voltage (V)']
            scandata.info['Electric Field (V/cm)'] = efield
        if 'Pump-Probe Distance (um)' not in scandata.info.keys():
            scandata.info['Pump-Probe Distance (um)'] = 0.0

    # tweaks to each model, e.g. if dataset has super high amplitude:
    def update_model_params(models_to_compare):
        feature_vector_model_1 = models_to_compare[0]
        # params = num_pulses, pulse_amplitude, species_amp_ratio,
        #          lifetime1, lifetime2, gfactor, phase,
        #          species1_efield_heating_coeff, species2_efield_heating_coeff
        #          mobility, slope, offset
        feature_vector_model_1.free_params = \
                                    [False, True, False,
                                     False, True, True, False,
                                     False, True,
                                     False, False, False]
        feature_vector_model_1.initial_params = \
                                    [20, 0.004, 2.0,
                                     20000, 10000, 0.44, 0,
                                     0.2, 0.5,
                                     1e-4, 0, 0]
        feature_vector_model_1.param_bounds[1] = (0, 1)  # pulse amp
        feature_vector_model_1.param_bounds[2] = (1.5, 3.0)  # species ratio
        feature_vector_model_1.param_bounds[4] = (1e3, 2e4)  # short lifetime
        feature_vector_model_2 = models_to_compare[1]
        # params = num_pulses, pulse_amplitude, species_amp_ratio,
        #          lifetime1, lifetime2, gfactor1, gfactor2,
        #          phase1, phase2, mobility, slope, offset
        feature_vector_model_2.free_params = \
                                      [False, True, True,
                                       False, False, False, True,
                                       False, True, True, False, False]
        feature_vector_model_2.initial_params = \
                                      [20, 0.0045, 1.0,
                                       20260, 20260, 0.44, 0.42,
                                       0, -1*np.pi/3, 1e-4, 0, 0]
        feature_vector_model_2.param_bounds[1] = (1e-3, 1)  # pulse amp
        feature_vector_model_2.param_bounds[6] = (0.42, 0.44)  # gfactor2

    # manually set uncertainty of data
    fixed_uncertainty = 1e-3
#    fixed_uncertainty = None  # ...or don't

    # excluded intervals of data's xfield (assumed scan coordinate)
    # (replaces model-specific excluded intervals)
    # (remember, applied post-data-filters, e.g. forcing positive delay times)
    excluded_intervals = [[-100, 400],  # for simulated TRKR
                          [LASER_REPRATE - 15, 15000]]
#    excluded_intervals = [[-100, 400],  # for actual TRKR @818.9, no neg data yet
#                          [7000, 15000]]


    # DEFINE MAP OF DATA TO FEATURE VECTORS
    # feature_vectors: (measurement, delaytime, efield, bfield,
    #                   pump_probe_dist, wavelength, temperature,
    #                   runID, index_in_run)
#    timefield = 'Delay Time (ps)'  # RSA DATA
#    bfieldfield = 'scancoord'  # RSA DATA
    timefield = 'scancoord'  # TRKR DATA
    bfieldfield = 'Magnetic Field (mT)'  # TRKR DATA
    yfield = 'lockin2x'
    fvec_fields = [('measurement', yfield),
                   ('time', timefield),
                   ('efield', "Electric Field (V/cm)"),
                   ('bfield', bfieldfield),
                   ('pump_probe_dist', "Pump-Probe Distance (um)"),
                   ('wavelength', "Wavelength (nm)"),
#                   ('temperature', "temperature"),  # right now ignores measured temp...
                   ('temperature', "Set Temperature (K)"),  # make set_temp?
                   ('runID', None),
                   ('index_in_run', None)]

    # Model 1: simplified two-opposite-pulses model
    # params = num_pulses, pulse_amplitude, species_amp_ratio,
    #          lifetime1, lifetime2, gfactor, phase,
    #          species1_efield_heating_coeff, species2_efield_heating_coeff
    #          mobility, slope, offset
    feature_vector_model_1 = \
        FeatureVectorsTwoLifetimesOppositePhaseTRKRModel(
            field_name="measurement",
            max_fcn_evals=20000,
            free_params=[False, True, True,
                         True, True, True, True,
                         True, True,
                         False, False, False],
#            initial_params=([40] + \
#                            list((1.0 + np.random.randn(7) *
#                                   np.array([0.2, 0.2,  # std per param
#                                             0.2, 0.2, 0.01, 0.0,
#                                             0.0, 0.0])) *
#                                 np.array([0.022, 1.8,  # correct param
#                                           20000, 2000, 0.44, 0.0,
#                                           0.0, 0.0])) + \
#                            [1e-4, 0, 0]),
            initial_params=[20, 0.022, 1.8,  # correct params
                            20000, 6000, 0.44, 0,
                            0.5, 0.5,
                            1e-4, 0, 0],
            param_bounds=[(1,1000), (0, 1), (0, 100),
                          (1, 1e9), (1, 1e9), (0.3, 0.6), (-2*np.pi, 2*np.pi),
                          (0, 4.0), (0, 4.0),
                          (0, 1), (-1e-6, 1e-6), (-0.01, 0.01)],
            fit_result_scan_coord="Pump-Probe Distance (um)",
            excluded_intervals=None,  # use feature vectors' alternative
            b_field=300,
            ignore_weights=False)

    # DEFINE MODELS TO COMPARE
    # Model 2: two totally-adjustable sines
    # params = num_pulses, pulse_amplitude, species_amp_ratio,
    #          lifetime1, lifetime2, gfactor1, gfactor2,
    #          phase1, phase2, mobility, slope, offset
    feature_vector_model_2 = \
        FeatureVectorsIndependentSinusoidalTRKRModel(
            field_name="measurement",
            max_fcn_evals=20000,
            free_params=[False, True, True,
                         False, True, False, True,
                         False, True, False, False, False],
            initial_params=[20, 0.04, 1.0,
                            20000, 18000, 0.44, 0.43,
                            0, -1*np.pi/3, 1e-4, 0, 0],
            param_bounds=[(1,1000), (0, 1), (0, 100),
                          (1, 1e9), (1, 1e9),
                          (0.3, 0.6), (0.3, 0.6),
                          (-np.pi, np.pi), (-4*np.pi, 4*np.pi),
                          (0, 1), (-1e-6, 1e-6), (-0.01, 0.01)],
            fit_result_scan_coord="Pump-Probe Distance (um)",
            excluded_intervals=None,  # use feature vectors' alternative
            b_field=300,
            ignore_weights=False)

    models_to_compare = [feature_vector_model_1, feature_vector_model_2]
    update_model_params(models_to_compare)
 
    # FETCH DATA AND APPLY CORRECTIONS/FILTERS
    scandata_list = \
        list(fetch_dir_as_unfit_scandata_iterator(
                    directorypath=dirpath,
                    yfield=yfield,
                    yfield_error_val=fixed_uncertainty,
                    parsing_keywordlists=FILEPATH_PARSING_KEYWORD_LISTS))
    for scandata in scandata_list:
        update_scandata_info(scandata)
        gaussian_smooth_scandata(scandata,
                                 fields_to_process=[yfield],
                                 gaussian_width=600,
                                 edge_handling='reflect',
                                 subtract_smoothed_data_from_original=True)
        make_scandata_time_delay_positive(scandata,
                                          zero_delay_offset=-15,
                                          neg_time_weight_multiplier=5.0)

    # CONVERT TO A FEATURE VECTOR SCANDATA
    fvec_scandata = \
        scandata_list_to_feature_vector_scandata(scandata_list, fvec_fields,
                                                 excluded_intervals)

    # SPLIT INTO TRAINING-SET AND TESTING-SET SCANDATA
    training_set_fvec_scandata, test_set_fvec_scandata = \
        split_scandata_into_training_and_test_sets(fvec_scandata,
                                                   test_fraction, test_by_run)

    for model in models_to_compare:
        fitdata = scandata_model_fit(training_set_fvec_scandata, model)
        training_set_fvec_scandata.fit_result = \
            model.fitfunction(training_set_fvec_scandata.x.T, *fitdata.fitparams)
        test_set_fvec_scandata.fit_result = \
            model.fitfunction(test_set_fvec_scandata.x.T, *fitdata.fitparams)
    
        print('---')
        print("model: {}".format(model.model_type))
        print("Fitted parameters from training set:")
        print("({} parameters from {} data points)".format(
                len(fitdata.freeparamindices), len(training_set_fvec_scandata)))
        for param_label, param, param_std in zip(fitdata.fitparamlabels,
                                                 fitdata.fitparams,
                                                 fitdata.fitparamstds):
            if param_std != 0.0:
                print(" -{}: {:.4g} +- {:.4g}".format(param_label,
                                                      param, param_std))
        training_error = np.mean((training_set_fvec_scandata.y -
                                  training_set_fvec_scandata.fit_result)**2)
        test_error = np.mean((test_set_fvec_scandata.y -
                              test_set_fvec_scandata.fit_result)**2)
        print('Goodness of Fit (mean-ssd):')
        print(' training error: {:.4g}'.format(training_error))
        print('  testing error: {:.4g}'.format(test_error))

        # add some useful attributes to model for cross-analysis
        model.fvec_scandata = fvec_scandata
        model.training_set_fvec_scandata = training_set_fvec_scandata
        model.test_set_fvec_scandata = test_set_fvec_scandata
        model.fitdata = fitdata


# %% plot results on one scandata
    scandata_index_to_plot = 20
    plt.figure()
    for model_ind, model in enumerate(models_to_compare):
        fitparams = model.fitdata.fitparams
        scandata_mask = (model.fvec_scandata.runID == scandata_index_to_plot)
        fvecs = model.fvec_scandata.x[scandata_mask]
        yvals = model.fvec_scandata.y[scandata_mask]

        # find xvals, a little roundabout
        original_scandata = scandata_list[scandata_index_to_plot]
        for fvec_elem_index, (_, original_field) in enumerate(fvec_fields):
            if original_field == original_scandata.xfield:
                xvals = model.fvec_scandata.x[scandata_mask, fvec_elem_index]
                break
        
        axes = plt.subplot(len(models_to_compare), 1, model_ind + 1)
        axes.set_title(model.model_type)
        axes.plot(xvals, yvals, 'bd')
        axes.plot(xvals, model.fitfunction(fvecs.T, *fitparams), 'r-')


#        test_scandata = model.test_set_fvec_scandata
#        training_scandata = model.training_set_fvec_scandata
#        test_scandata = \
#            test_scandata[test_scandata.runID == scandata_index_to_plot]
#        training_scandata = \
#            training_scandata[training_scandata.runID == scandata_index_to_plot]
#        xvals = scandata_list[scandata_index_to_plot].x  # assumes none dropped


