# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 12:49:27 2017

@author: Michael
"""

import numpy as np

from experimentdataanalysis.analysis.dataclasses import ScanData


def scandata_list_to_fvec_scandata(scandata_list, fvec_fields,
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


def split_fvec_scandata_by_training_and_test(fvec_scandata, test_fraction,
                                             test_by_runID=True):
    if test_by_runID:
        num_all_targets = int(max(fvec_scandata.runID)) + 1
    else:
        num_all_targets = len(fvec_scandata)
    num_test_targets = int(np.ceil(np.floor(test_fraction * num_all_targets)))
    test_targets = np.array([True] * num_test_targets +
                            [False] * (num_all_targets - num_test_targets))
    np.random.shuffle(test_targets)
    if test_by_runID:
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
