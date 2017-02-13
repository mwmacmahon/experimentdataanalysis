# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:47:08 2016

@author: Michael
"""

import numpy as np

from experimentdataanalysis.analysis.dataclasses import FitData, ScanData
from experimentdataanalysis.analysis.scandataprocessing \
    import ScanDataModel, \
           scandata_list_fit, scandata_list_model_fit, \
           scandata_iterable_sort

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
        self.model_param_label_list = None
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
                                       self.model.yfield,
                                       primary_key, secondary_key)
            self.scandata_list = list(scandata_list_tuple)

    def fit_scandata_to_model(self, multiprocessing=False):
        fit_field = self.model.yfield
        if fit_field is None:
            fit_field = self.scandata_list[0].yfield
        # must filter out any pts with error = 0, will kill weighted fit
        if not self.model.ignore_weights:
            usable_scandata_list = self.get_nonzero_error_scandata(fit_field)
        else:
            usable_scandata_list = self.scandata_list
        scandata_list_model_fit(usable_scandata_list, self.model,
                                multiprocessing=multiprocessing)
#        scandata_list_fit(self.scandata_list,
#                          fit_field,
#                          self.model.fitfunction,
#                          self.model.get_fit_params_is_free_param(),
#                          self.model.get_fit_params_initial_values(),
#                          self.model.get_fit_params_bounds(),
#                          self.model.max_fcn_evals,
#                          self.model.excluded_intervals,
#                          self.model.ignore_weights,
#                          multiprocessing)
        for scandata in usable_scandata_list:
            fitdata = scandata.get_field_fitdata(fit_field)
            if fitdata is not None:
                setattr(scandata, 'fitdata_' + fit_field,
                        FitData(fitdata.fitfunction, fitdata.partialfcn,
                                fitdata.fitparams, fitdata.fitparamstds,
                                self.model.get_fit_params_labels(),
                                fitdata.fityvals,
                                fitdata.freeparamindices,
                                fitdata.covariancematrix,
                                fitdata.meansquarederror))
        self.parse_fit_scandata()  # will purge scandata

    def get_nonzero_error_scandata(self, field_name):
        def has_zero_error(scandata):
            yerrvals = scandata.get_field_yerr(field_name)
            if yerrvals is None:
                return True
            elif 0 in yerrvals:
                return True
            else:
                return None

        nonzero_error_scandata_list = [scandata
                                       for scandata in self.scandata_list
                                       if not has_zero_error(scandata)]
        if len(nonzero_error_scandata_list) < len(self.scandata_list):
            print("Warning: ScanData value found with 0 uncertainty " +
                  "in DataSeries, so ScanData was thrown out as " +
                  "incompatible with weighted fit")
        return nonzero_error_scandata_list

#    def purge_zero_error_scandata(self, field_name):
#        def has_zero_error(scandata):
#            yerrvals = scandata.get_field_yerr(field_name)
#            if yerrvals is None:
#                return True
#            elif 0 in yerrvals:
#                return True
#            else:
#                return None
#
#        fail_message = ("Warning: ScanData value found with 0 uncertainty " +
#                        "in DataSeries, so ScanData was thrown out as " +
#                        "incompatible with weighted fit")
#        self.purge_scandata_list(field_name, has_zero_error, fail_message)
        
#    def purge_failed_fit_scandata(self, field_name):
#        def has_no_fitdata(scandata):
#            return scandata.get_field_fitdata(field_name) is None
#
#        fail_message = ("Warning: Fit failed in ScanDataSet, erasing " +
#                        "associated ScanData so all ScanData have " +
#                        "associated fits during further analysis.")
#        self.purge_scandata_list(field_name, has_no_fitdata, fail_message)

#    def purge_scandata_list(self, test_field_name,
#                            scandata_failed_test_fcn, fail_message):
#        filtered_scandata_list = []
#        for scandata in self.scandata_list:
#            if not scandata_failed_test_fcn(scandata):
#                filtered_scandata_list.append(scandata)
#            else:
#                print(fail_message)
#                if 'Filepath' in scandata.info.keys():
#                    print("Associated filepath:")
#                    print(scandata.info['Filepath'])
#        self.scandata_list = filtered_scandata_list

    def parse_fit_scandata(self, fit_result_scan_coord_override=None):
        # TODO: just make this a scandata
        yfield = self.model.yfield
        if yfield is None:
            yfield = self.scandata_list[0].yfield
        # if no fitdata for this field parameter, cancel parse
        fitted_scandata_list = \
            [scandata
             for scandata in self.scandata_list
             if scandata.get_field_fitdata(yfield) is not None]
        if len(fitted_scandata_list) == 0:
            return
        # all scandata in list should have FitData and be processable:
        # create array of fit parameters & their uncertainties vs xval:
        if fit_result_scan_coord_override is None:
            xfield = self.fit_result_scan_coord
        else:
            xfield = fit_result_scan_coord_override
        x_array = np.array([scandata.info[xfield]
                            for scandata in fitted_scandata_list])
        first_fitdata = fitted_scandata_list[0].get_field_fitdata(yfield)
        fitparamlabels = list(first_fitdata.fitparamlabels)
        fitparamlabel_array_list = [np.array(fitparamlabel_series)
                                    for fitparamlabel_series in fitparamlabels]
        fitparams = list(zip(*(scandata.get_field_fitdata(yfield).fitparams
                               for scandata in fitted_scandata_list)))
        fitparam_array_list = [np.array(fitparam_series)
                               for fitparam_series in fitparams]
        fitparamstds = \
            list(zip(*(scandata.get_field_fitdata(yfield).fitparamstds
                       for scandata in fitted_scandata_list)))
        fitparamstd_array_list = [np.array(fitparamstd_series)
                                  for fitparamstd_series in fitparamstds]

        # save these results to scandataset:
        self.model_param_label_list = [xfield] + fitparamlabel_array_list
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

    def get_model_param_label_list(self):
        if self.model_param_label_list is None:
            return []
        else:
            return self.model_param_label_list.copy()

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
        # extract actual data from model, fields are model output params
        # note "param 1" = xcoord of new scandata, fit_result_scan_coord
        param_labels = self.get_model_param_label_list()  
        params = self.get_model_param_array_list()
        param_sigmas = \
            self.get_model_param_uncertainty_array_list()
        if len(params) == 1 or len(param_sigmas) == 0:
            print('extract_scandata: no successful fits to extract!')
            return None

        field_names = param_labels + [param_name + '_error'  # no x-error
                                      for param_name in param_labels[1:]]
        field_arrays = params + param_sigmas  # already excludes x-error

        # get new scaninfo to reflect new coord types:
        new_scaninfo = {'fit_scandata_info_dicts': []}
        fit_field = self.model.yfield
        if fit_field is None:
            fit_field = self.scandata_list[0].yfield
        fitted_scandata_list = \
            [scandata
             for scandata in self.scandata_list
             if scandata.get_field_fitdata(fit_field) is not None]
        for scandata in fitted_scandata_list:
            # alternative: just keep all scandata? could bloat though
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
                usable_scandata_list = \
                    scandataset.get_nonzero_error_scandata(model.yfield)
            else:
                usable_scandata_list = scandataset.scandata_list
            for scandata in usable_scandata_list:
                scandata_fit_list.append([scandata,
                                          model.yfield,
                                          model.fitfunction,
                                          model.get_fit_params_is_free_param(),
                                          model.get_fit_params_initial_values(),
                                          model.get_fit_params_bounds(),
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
        new_scandataset_list = []
        for scandataset in first_scandataset_list:
            # subdivide _each_ set by _last_ key
            subdivided_scandataset_list = \
                sort_scandata_into_sets_single_key(scandataset.scandata_list,
                                                   model, sort_keys[-1])
            old_setname = scandataset.setname
            for new_scandataset in subdivided_scandataset_list:
                new_scandataset.setname = \
                    str(old_setname) + ", " + str(new_scandataset.setname)
            new_scandataset_list.extend(subdivided_scandataset_list)
        return new_scandataset_list


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
#            print("Warning: sort key not found in ScanData, skipping...")
#            raise KeyError("sort_scandata_into_sets: invalid sort_key")
            pass
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

