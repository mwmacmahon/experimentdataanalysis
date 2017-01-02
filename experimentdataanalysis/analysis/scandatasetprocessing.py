# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:47:08 2016

@author: Michael
"""

import numpy as np

from experimentdataanalysis.analysis.dataclasses import FitData, ScanData
from experimentdataanalysis.analysis.scandataprocessing \
    import scandata_fit, scandata_list_fit, scandata_iterable_sort
from experimentdataanalysis.analysis.scandatamodels \
    import ScanDataModel

from copy import deepcopy


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class ScanDataSet:
    """
    Stores a list of ScanData, a model used to fit them, and the results of
    said fits once they are performed. Note this is a mutable data structure,
    and ScanData themselves are mutable (they will, for example, gain a
    FitData key_field after fitting).
    """
    def __init__(self, setname, model):
        self.setname = setname
        self.model = model
        self.chronological_key = "FastScanIndex"
        self.fit_result_scan_coord = model.fit_result_scan_coord
        self.model_param_array_list = None
        self.model_param_uncertainty_array_list = None
        self.model_param_arrays_scan_coord = None
        self.scandata_list = []

    def __setattr__(self, name, value):
        super().__setattr__(name, value)  # passthrough attribute set, but...
        if name == "scandata_list":  # if already fit, reparse fit params
            if self.model_param_arrays_scan_coord is not None:
                self.parse_fit_scandata(self.model_param_arrays_scan_coord)

    def sort_scandata_list(self):
        if len(self.scandata_list) > 1:
            # fix scandataset ordering errors arising from alphanumeric sorting
            primary_key = self.chronological_key
            secondary_key = self.fit_result_scan_coord
            scandata_list_tuple, _ = \
                scandata_iterable_sort(self.scandata_list,
                                       self.model.field_name,
                                       primary_key, secondary_key)
            self.scandata_list = list(scandata_list_tuple)

    def fit_scandata_to_model(self, multiprocessing=False,
                              purge_failed_fits=True):
        # must filter out any pts with error = 0, will kill weighted fit
        self.purge_zero_error_scandata(self.model.field_name)
        scandata_list_fit(self.scandata_list,
                          self.model.field_name,
                          self.model.fitfunction,
                          self.model.free_params,
                          self.model.initial_params,
                          self.model.param_bounds,
                          self.model.max_fcn_evals,
                          self.model.excluded_intervals,
                          self.model.ignore_weights,
                          multiprocessing)
        for scandata in self.scandata_list:
            fitdata = scandata.get_field_fitdata(self.model.field_name)
            if fitdata is not None:
                setattr(scandata, 'fitdata_' + self.model.field_name,
                        FitData(fitdata.fitfunction, fitdata.partialfcn,
                                fitdata.fitparams, fitdata.fitparamstds,
                                self.model.model_params,
                                fitdata.fityvals,
                                fitdata.freeparamindices,
                                fitdata.covariancematrix,
                                fitdata.meansquarederror))
        self.parse_fit_scandata()  # will purge scandata

    def purge_zero_error_scandata(self, field_name):
        def has_zero_error(scandata):
            yerrvals = scandata.get_field_yerr(field_name)
            if yerrvals is None:
                return True
            elif 0 in yerrvals:
                return True
            else:
                return None

        fail_message = ("Warning: ScanData value found with 0 uncertainty " +
                        "in DataSeries, so ScanData was thrown out as " +
                        "incompatible with weighted fit")
        self.purge_scandata_list(field_name, has_zero_error, fail_message)

    def purge_failed_fit_scandata(self, field_name):
        def has_no_fitdata(scandata):
            return scandata.get_field_fitdata(field_name) is None

        fail_message = ("Warning: Fit failed in ScanDataSet, erasing " +
                        "associated ScanData so all ScanData have " +
                        "associated fits during further analysis.")
        self.purge_scandata_list(field_name, has_no_fitdata, fail_message)

    def purge_scandata_list(self, test_field_name,
                            scandata_failed_test_fcn, fail_message):
        filtered_scandata_list = []
        for scandata in self.scandata_list:
            if not scandata_failed_test_fcn(scandata):
                filtered_scandata_list.append(scandata)
            else:
                print(fail_message)
                if 'Filepath' in scandata.info.keys():
                    print("Associated filepath:")
                    print(scandata.info['Filepath'])
        self.scandata_list = filtered_scandata_list

    def parse_fit_scandata(self, fit_result_scan_coord_override=None):
        # TODO: just make this a scandata
        field_name = self.model.field_name
        # if no fitdata for this field parameter, cancel parse
        fitted_scandata_list = \
            [scandata
             for scandata in self.scandata_list
             if scandata.get_field_fitdata(field_name) is not None]
        if len(fitted_scandata_list) == 0:
            return
#        # purge failed fits, as they will cause issues further down.
#        self.purge_failed_fit_scandata(field_name)
        # all remaining scandata should have FitData and be processable:
        # create dataseries of fit parameters & their uncertainties vs xval:
        if fit_result_scan_coord_override is None:
            xfield = self.fit_result_scan_coord
        else:
            xfield = fit_result_scan_coord_override
        x_array = np.array([scandata.info[xfield]
                            for scandata in fitted_scandata_list])
        fitparams = []
        fitparamstds = []
        fitparams = list(zip(*(scandata.get_field_fitdata(field_name).fitparams
                               for scandata in fitted_scandata_list)))
        fitparam_array_list = [np.array(fitparam_series)
                               for fitparam_series in fitparams]
        fitparamstds = \
            list(zip(*(scandata.get_field_fitdata(field_name).fitparamstds
                       for scandata in fitted_scandata_list)))
        fitparamstd_array_list = [np.array(fitparamstd_series)
                                  for fitparamstd_series in fitparamstds]

        # save these results to scandataset:
        self.model_param_array_list = [x_array] + fitparam_array_list
        self.model_param_uncertainty_array_list = fitparamstd_array_list
        self.model_param_arrays_scan_coord = xfield

    def apply_transform_to_scandata(self, transform_fcn, *args, **kwargs):
        new_scandata_list = []
        for scandata in self.scandata_list:
            return_val = transform_fcn(scandata,*args, **kwargs)
            if return_val is None:  # if transform modifies scandata in-place
                new_scandata_list.append(scandata)
            else:  # if transform returns new scandata
                new_scandata_list.append(return_val)
        self.scandata_list = new_scandata_list

    def get_scandata_list(self):
        return self.scandata_list

    def get_model_param_array_list(self):
        if self.model_param_array_list is None:
            return []
        else:
            return self.model_param_array_list.copy()

    def get_model_param_uncertainty_array_list(self):
        if self.model_param_uncertainty_array_list is None:
            return []
        else:
            return self.model_param_uncertainty_array_list.copy()

    def extract_scandata(self, fit_result_scan_coord_override=None):
        if fit_result_scan_coord_override is not None:
            self.parse_fit_scandata(fit_result_scan_coord_override)
        fit_result_scan_coord = self.model_param_arrays_scan_coord
        # extract actual data from model, fields are model output params
        params = self.get_model_param_array_list()  # param 1 = x
        param_sigmas = \
            self.get_model_param_uncertainty_array_list()
        if len(params) == 1 or len(param_sigmas) == 0:
            print('extract_scandata: no successful fits to extract!')
            return None
        field_names, field_arrays, field_error_arrays = \
            self.model.all_model_fields(params, param_sigmas)
        field_names = ([fit_result_scan_coord] + field_names +  # x field name
                       [field_name + '_error' for field_name in field_names])
        field_arrays += field_error_arrays  # already starts with x-array!

        # get new scaninfo to reflect new coord types:
        new_scaninfo = {'fit_scandata_info_dicts': []}
        for scandata in self.scandata_list:
            if scandata.fitdata is not None:
                info_copy = deepcopy(scandata.info)  # take a snapshot
                new_scaninfo['fit_scandata_info_dicts'].append(info_copy)
        return ScanData(field_names,
                        field_arrays,
                        new_scaninfo,
                        xfield=None,  # 1st field is already x-coord!
                        yfield=None)  # defaults to first model param


# %%
def collapse_scandataset_to_model_fit_scandata_list(
                                        scandataset_list,
                                        fit_result_scan_coord_override=None,
                                        ignore_empty_scandatasets=True):
    scandata_list = [scandataset.extract_scandata(fit_result_scan_coord_override)
                     for scandataset in scandataset_list]
    if ignore_empty_scandatasets:
        return [scandata for scandata in scandata_list
                if scandata is not None
                if len(scandata.x) > 0]
    else:
        return scandata_list


# %%
def fit_scandataset_list_to_model(scandataset_list, multiprocessing=False):
    if not multiprocessing:
        for scandataset in scandataset_list:
            scandataset.fit_scandata_to_model()
    else:
        # cobble together all scandata across all and process in parallel
        # if using weights, must kill error = 0 pts, will ruin weighted fit
        scandata_fit_list = []
        merged_scandata_list = []
        for scandataset in scandataset_list:
            model = scandataset.model
            if not model.ignore_weights:
                scandataset.purge_zero_error_scandata(model.field_name)
            for scandata in scandataset.scandata_list:
                scandata_fit_list.append([scandata,
                                          model.field_name,
                                          model.fitfunction,
                                          model.free_params,
                                          model.initial_params,
                                          model.param_bounds,
                                          model.max_fcn_evals,
                                          model.excluded_intervals,
                                          model.ignore_weights])
                merged_scandata_list.append(scandata)
        if len(scandata_fit_list) > 0:
            scandata_list_fit(*zip(*scandata_fit_list),
                              multiprocessing=multiprocessing)
        else:
            print("fit_all_scandata_to_model: No ScanData left after " +
                  "purging ScanData with no weights in model's fit field " +
                  "and ignore_weights flag set to False.")
        # tell scandatasets to examine their newly fit scandata!
        for scandataset in scandataset_list:
            scandataset.parse_fit_scandata()


# %%
def sort_scandata_into_sets(scandata_list, model, sort_keys=["Filepath"]):
    """
    Takes a list of ScanData and sorts them into sets based on one or more
    sort keys given. For example, if sort_key is "Voltage", will sort each
    ScanData into sets with shared scaninfo entry {..."Voltage": X, ...}
    where X is different for each set.

    Returns a list of ScanDataSets, each with the analysis model provided.
    """
    if sort_keys == None:  # base case: no sort keys
        return sort_scandata_into_sets_single_key(scandata_list, model, None)
    elif isinstance(sort_keys, str):  # base case: sort_keys is a string
        return sort_scandata_into_sets_single_key(scandata_list,
                                                  model, sort_keys)
    elif len(sort_keys) == 1:  # base case, only one sort key
        return sort_scandata_into_sets_single_key(scandata_list,
                                                  model, sort_keys[0])
    else:  # just grab set where all but last key already used, recursively
        first_scandataset_list = sort_scandata_into_sets(scandata_list,
                                                         model, sort_keys[:-1])
        scandataset_list = []
        for scandataset in first_scandataset_list:
            scandataset_list.extend(  # subdivide _each_ set by _last_ key
                sort_scandata_into_sets_single_key(scandataset.scandata_list,
                                                   model, sort_keys[-1]))
        return scandataset_list


def sort_scandata_into_sets_single_key(scandata_list, model, sort_key):
    if len(scandata_list) == 0:
        raise TypeError("sort_scandata_into_sets: " +
                        "empty ScanData list provided")
    try:  # test key existence & if values associated with sort_key numerical
        for scandata in scandata_list:
            float(scandata.info[sort_key])
    except KeyError:
        if sort_key is None:
            pass  # all scandata put into one set, ignore sorting
        else:
            raise KeyError("sort_scandata_into_sets: invalid sort_key")
    except ValueError:  # sort key non-numeric, so pre-sort alphanumerically
        scandata_list, _ = scandata_iterable_sort(scandata_list, 0,
                                                  sort_key, sort_key,
                                                  numeric_sort=False)
    else:  # sort key numeric, so pre-sort numerically
        scandata_list, _ = scandata_iterable_sort(scandata_list, 0,
                                                  sort_key, sort_key,
                                                  numeric_sort=True)

    # Divide ScanData into ScanDataSets. Since pre-sorted, create a new
    # set whenever a new value is found.
    scandataset_list = []
    scandata_list = list(scandata_list)
    current_scandataset = ScanDataSet("[default]", ScanDataModel())
    last_setname = ""
    for scandata in scandata_list:
        try:  # use sort key's value to group ScanData into ScanDataSets
            setname = scandata.info[sort_key]
            if sort_key == "Filepath":  # special case
                setname = setname.split("\\")[-2]
        except KeyError:
            if sort_key == None:
                setname = "All ScanData"                    
            else:
                print("Warning: ScanDataSet sort key not found in " +
                      "ScanData! Sets may not be properly grouped")
                setname = "key_error_set"
        if setname != last_setname:
            current_scandataset.sort_scandata_list()
            current_scandataset = ScanDataSet(setname, model.copy())
            scandataset_list.append(current_scandataset)
            last_setname = setname
        current_scandataset.scandata_list.append(scandata)
    current_scandataset.sort_scandata_list()  # sort last scandataset   

    return scandataset_list     


def regroup_scandatasets(scandataset_list, new_model, sort_key="Filepath"):
    all_scandata_list = sum((scandataset.scandata_list
                             for scandataset in scandataset_list), [])
    new_scandataset_list = sort_scandata_into_sets(all_scandata_list,
                                                   model=new_model,
                                                   sort_key=sort_key)
    for scandataset in new_scandataset_list:  # reparse fit results
        scandataset.parse_fit_scandata()
    return new_scandataset_list

