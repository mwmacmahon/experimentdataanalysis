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
from experimentdataanalysis.parsing.scandataparsing import \
        fetch_dir_as_unfit_scandata_iterator

import fit_models_two_species_pump_probe as fit_models


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
                                    ("V", "Voltage (V)"),
                                    ("x", "SecondScanCoord")]

# for this one, if element[0] found, next element stored w/ key element[1]
in_filepath_next_element_keyword_list = [("MirrorZ",
                                          "Pump-Probe Distance (um)"),
                                          ("Voltage", "Voltage (V)")]
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

#        scandata_measurement_error = getattr(scandata,
#                                             fvec_fields[0][1] + '_error',
#                                             None)  # default to unknown error
#        if scandata_measurement_error is not None:
#            scandata_measurement_error = \
#                scandata_measurement_error[indices_to_use_mask]

        # TEMPORARY: MESSING WITH WEIGHTS
        # note: for some bizarre reason making it the mean of the ABS of the
        # data ruins the effect - making the data 'weight' independent of
        # magnitude.
#        scandata_measurement_error = \
#           np.ones(scandata_nfvecs) * (1.0 / np.mean(scandata.y))
        scandata_measurement_error = \
           np.ones(scandata_nfvecs) * (1.0 / np.mean(np.abs(scandata.y)))**2

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


def plot_scandata_fvec_fit(scandata_index, scandata_list,
                           axes, model, forced_fitparams=None):
    if forced_fitparams is not None:
        fitparams = forced_fitparams
    else:
        fitparams = model.fitdata.fitparams
    scandata_mask = (model.fvec_scandata.runID == scandata_index)
    fvecs = model.fvec_scandata.x[scandata_mask]
    yvals = model.fvec_scandata.y[scandata_mask]
#    yerrvals = model.fvec_scandata.measurement_error[scandata_mask]

    # find xvals, a little roundabout
    xfield = scandata_list[scandata_index].xfield
    for fvec_elem_index, (_, original_field) in enumerate(fvec_fields):
        if original_field == xfield:
            xvals = model.fvec_scandata.x[scandata_mask, fvec_elem_index]
            break

    axes.set_title(model.model_name)
    axes.plot(xvals, yvals, 'bd', markersize=5)
    axes.plot(xvals, model.fitfunction(fvecs.T, *fitparams), 'r-')
    plt.xticks([])
    plt.yticks([])


# %%
if __name__ == '__main__':
# %%
    # GENERAL OPTIONS
    test_fraction = 0.3
    test_by_run = True  # TEST BY RUN VS INDIVIDUAL POINTS

    # DATA TO FIT
    dirpath = ("C:\\Data\\fake_data\\fake_trkr")
#    dirpath = ("C:\\Data\\fake_data\\fake_trkr_huge_slopes")
#    dirpath = ("C:\\Data\\fake_data\\fake_rsa")
#    dirpath = ("C:\\Data\\august_data\\160902\\" +
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_" +
#               "033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    # DATA OPTIONS:
    # manipulations on each scandata, e.g. filling missing info:
    def update_scandata_info(scandata):
        if 'TRKRvsV' in scandata.info['Filepath']:
            scandata.info['Voltage (V)'] = scandata.info['SecondScanCoord']
        if 'Set Temperature (K)' not in scandata.info.keys():
            scandata.info['Set Temperature (K)'] = 30.0  # ASSUME DEFAULT
        if 'temperature' not in scandata.fields:
            scandata.info['temperature'] = scandata.info['Set Temperature (K)']
        if 'Electric Field (V/cm)' not in scandata.info:
            if 'Voltage (V)' in scandata.info:
                efield = 20 * scandata.info['Voltage (V)']
                scandata.info['Electric Field (V/cm)'] = efield
            else:
                scandata.info['Electric Field (V/cm)'] = 0.0
        if 'Pump-Probe Distance (um)' not in scandata.info.keys():
            scandata.info['Pump-Probe Distance (um)'] = 0.0

    # manually set uncertainty of data
    fixed_uncertainty = 2e-4
#    fixed_uncertainty = None  # ...or don't

    # excluded intervals of data's xfield (assumed scan coordinate)
    # (replaces model-specific excluded intervals)
    # (remember, applied post-data-filters, e.g. forcing positive delay times)
    excluded_intervals = [[-100, 100],
#                          [LASER_REPRATE - 15, 15000]]  # neg dalyta on to 15ps
                          [7000, 15000]]  # no negative data at all


    # DEFINE MAP OF DATA TO FEATURE VECTORS
    # feature_vectors: (measurement, delaytime, efield, bfield,
    #                   pump_probe_dist, wavelength, temperature,
    #                   set_temperature, runID, index_in_run)
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
                   ('temperature', "temperature"),  # right now ignores measured temp...
                   ('set_temperature', "Set Temperature (K)"),  # make set_temp?
                   ('runID', None),
                   ('index_in_run', None)]

    models_to_compare = [fit_models.get_fvec_two_species_model()]
 
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
    
        # PRINT FIT RESULTS
        fitted_coord_strings = []
        fixed_coord_strings = []
        zipdata = zip(fitdata.fitparamlabels,
                      fitdata.fitparams,
                      fitdata.fitparamstds)
        for param_ind, (param_label, param, param_std) in enumerate(zipdata):
            if param_ind not in fitdata.freeparamindices:
                try:
                    fixed_coord_strings.append(
                        " -{}: {:.4g}".format(param_label, param))
                except (TypeError, ValueError):
                    fixed_coord_strings.append(
                        " -{}: {}".format(param_label, param))
            else:
                fitted_coord_strings.append(
                    " -{}: {:.4g} +- {:.4g}".format(param_label,
                                                    param, param_std))
        print('---')
        print("model: {}".format(model.model_name))
        print("Fixed parameters:")
        for param_string in fixed_coord_strings:
            print(param_string)
        print("Fitted parameters from training set:")
        print("({} parameters from {} data points)".format(
                len(fitdata.freeparamindices), len(training_set_fvec_scandata)))
        for param_string in fitted_coord_strings:
            print(param_string)
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


# %%
if False:
# %%
    model = models_to_compare[0]
    fitparams = None  # grab from model
    # OPTIONAL: force fit params
#    fitparams = \
#        [
#            20.0, 0.025, 2.0, 20260, 8000,  # for best 818.9 data 09/02/2016
#            0.4385, 0.0, np.pi/2, -np.pi,
#            0.03, -0.8, -1.2, -1.2, 0.00035,
#            0.0001, 0.0, 0.0
##            20.0, 0.024, 2.0, 20260, 8000,  # for fake_trkr 1/9/2017
##            0.44, 0.0, np.pi/2, -np.pi,
##            0.03, 1.0, -1.5, 0.0004,
##            0.0001, 0.0, 0.0
#        ]
    n_rows = 3
    n_cols = 4
    for scandata_ind in range(len(scandata_list)):
        subplot_ind = scandata_ind % (n_rows * n_cols) + 1
        if subplot_ind == 1:
            plt.figure()
        axes = plt.subplot(n_rows, n_cols, subplot_ind)
        plot_scandata_fvec_fit(scandata_ind, scandata_list,
                               axes, model, fitparams)
        


# %%
if False:
# %% plot results on one scandata
    scandata_index_to_plot = 7
#    for item in scandata_list[scandata_index_to_plot].info.items():
#        print(item)
    current_scandata = scandata_list[scandata_index_to_plot]
    if 'temperature' in scandata_list[scandata_index_to_plot].fields:
        current_temp = np.mean(current_scandata.temperature)
    else:
        current_temp = current_scandata.info['temperature']
    current_field = current_scandata.info['Electric Field (V/cm)']
    print('temp: {} K'.format(current_temp))
    print("field: {} V/cm".format(current_field))
    plt.figure()
    for model_ind, model in enumerate(models_to_compare):
        axes = plt.subplot(len(models_to_compare), 1, model_ind + 1)
        fitfunction = model.fitfunction
        fitparams = model.fitdata.fitparams
        # OPTIONAL: force fit params
#        fitparams = \
#            [
#                20.0, 0.025, 2.0, 20260, 8000,  # for best 818.9 data 09/02/2016
#                0.4385, 0.0, np.pi/2, -np.pi,
#                0.03, -0.8, -1.2, -1.2, 0.00035,
#                0.0001, 0.0, 0.0
#                20.0, 0.024, 2.0, 20260, 8000,  # for fake_trkr 1/9/2017
#                0.44, 0.0, np.pi/2, -np.pi,
#                0.03, 1.0, -1.5, 0.0004,
#                0.0001, 0.0, 0.0
#            ]
        scandata_mask = (model.fvec_scandata.runID == scandata_index_to_plot)
        fvecs = model.fvec_scandata.x[scandata_mask]
        yvals = model.fvec_scandata.y[scandata_mask]
        yerrvals = model.fvec_scandata.measurement_error[scandata_mask]

        # find xvals, a little roundabout
        original_scandata = scandata_list[scandata_index_to_plot]
        for fvec_elem_index, (_, original_field) in enumerate(fvec_fields):
            if original_field == original_scandata.xfield:
                xvals = model.fvec_scandata.x[scandata_mask, fvec_elem_index]
                break

        axes.set_title(model.model_name)
        axes.plot(xvals, yvals, 'bd')
        axes.plot(xvals, model.fitfunction(fvecs.T, *fitparams), 'r-')
        plt.title('field = {} V/cm, temp = {} K'.format(
                    current_field, current_temp))


    plt.title('field = {} V/cm, temp = {} K'.format(
                current_field, current_temp))

#        test_scandata = model.test_set_fvec_scandata
#        training_scandata = model.training_set_fvec_scandata
#        test_scandata = \
#            test_scandata[test_scandata.runID == scandata_index_to_plot]
#        training_scandata = \
#            training_scandata[training_scandata.runID == scandata_index_to_plot]
#        xvals = scandata_list[scandata_index_to_plot].x  # assumes none dropped


