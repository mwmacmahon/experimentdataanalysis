# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:47:08 2016

@author: Michael
"""

import numpy as np

from experimentdataanalysis.analysis.dataclasses import ScanData
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
    
    add_filter / add_filters add filters that limit which scandatasets are
    returned by get_
    """
    def __init__(self, setname, model):
        self.setname = setname
        self.model = model
        self.scandata_list = []
        self.scandata_filter_fcn_list = []
        self.chronological_key = "FastScanIndex"
        self.fit_result_scan_coord = model.fit_result_scan_coord
        self.model_param_array_list = None
        self.model_param_uncertainty_array_list = None
        if model is not None:
            self.scandata_filter_fcn_list += model.get_model_filter_fcns()

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

    def fit_scandata_to_model(self, multiprocessing=False):
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
        self.parse_fit_scandata()

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
        self.model_param_array_list = None  # now invalid!
        self.model_param_uncertainty_array_list = None

    def parse_fit_scandata(self, fit_result_scan_coord_override=None):
        # purge failed fits, as they will cause issues further down.
        field_name = self.model.field_name
        self.purge_failed_fit_scandata(field_name)
        # all remaining scandata should have FitData and be processable:
        # create dataseries of fit parameters & their uncertainties vs xval:
        if fit_result_scan_coord_override is None:
            xfield = self.fit_result_scan_coord
        else:
            xfield = fit_result_scan_coord_override
        x_array = np.array([scandata.info[xfield]
                            for scandata in self.scandata_list])
        fitparams = []
        fitparamstds = []
        fitparams = list(zip(*(scandata.get_field_fitdata(field_name).fitparams
                               for scandata in self.scandata_list)))
        fitparam_array_list = [np.array(fitparam_series)
                               for fitparam_series in fitparams]
        fitparamstds = \
            list(zip(*(scandata.get_field_fitdata(field_name).fitparamstds
                       for scandata in self.scandata_list)))
        fitparamstd_array_list = [np.array(fitparamstd_series)
                                  for fitparamstd_series in fitparamstds]

        # save these results to scandataset:
        self.model_param_array_list = [x_array] + fitparam_array_list
        self.model_param_uncertainty_array_list = fitparamstd_array_list

    def apply_transform_to_scandata(self, transform_fcn, *args, **kwargs):
        new_scandata_list = []
        for scandata in self.scandata_list:
            return_val = transform_fcn(scandata,*args, **kwargs)
            if return_val is None:  # if transform modifies scandata in-place
                new_scandata_list.append(scandata)
            else:  # if transform returns new scandata
                new_scandata_list.append(return_val)
        self.scandata_list = new_scandata_list

    def add_scandata_filters(self, *scandata_filter_fcn_list):
        for scandata_filter_fcn in scandata_filter_fcn_list:
            self.add_scandata_filter(scandata_filter_fcn)

    def add_scandata_filter(self, scandata_filter_fcn):
        """
        Adds a filter to an internally kept filter list. Filters are functions
        that accept a scandata and return a boolean - whether scandata passes
        filter or not. E.g. a simple filter might be:
        
        def only_fitted_scandata_filter(scandata):
            return all(fitdata != None for fitdata in scandata.fitdata_list)
        
        Note this is usually done AFTER fitting all scandata, so fit results
        can be used in filtering.
        """
        # insert checks for valid filter fcns?
        self.scandata_filter_fcn_list.append(scandata_filter_fcn)

    def clear_scandata_filters(self):
        self.scandata_filter_fcn_list = []
        if self.model is not None:
            self.scandata_filter_fcn_list += self.model.get_model_filter_fcns()

    def get_filtered_scandata_indices(self):
        """
        Returns a list containing all indices in scandata_list whose
        corresponding scandata meet the requirements off all filters
        in self.scandata_filter_fcn_list
        """
        filtered_indices = []
        for ind, scandata in enumerate(self.scandata_list):
            ok_flag = True
            for scandata_filter_fcn in self.scandata_filter_fcn_list:
                ok_flag = ok_flag and scandata_filter_fcn(scandata)
            if ok_flag:
                filtered_indices.append(ind)
        return filtered_indices

    def get_scandata_list(self, filtered=True):
        if filtered:
            return [scandata
                    for ind, scandata in enumerate(self.scandata_list)
                    if ind in self.get_filtered_scandata_indices()]
        else:
            return self.scandata_list.copy()

    def get_model_param_array_list(self, filtered=True):
        if filtered:
            filtered_indices = self.get_filtered_scandata_indices()
            if len(filtered_indices) > 0:
                return [series[filtered_indices]
                        for series in self.model_param_array_list]
            else:
                print('Warning: all scandata have been filtered out! ' +
                      'get_model_param_array_list filters ignored.')
                return self.model_param_array_list.copy()
        else:
            return self.model_param_array_list.copy()

    def get_model_param_uncertainty_array_list(self, filtered=True):
        if filtered:
            filtered_indices = self.get_filtered_scandata_indices()
            if len(filtered_indices) > 0:
                return [series[filtered_indices]
                        for series in self.model_param_uncertainty_array_list]
            else:
                print('Warning: all scandata have been filtered out! ' +
                      'get_model_param_uncertainty_array_list filters ' +
                      'ignored.')
                return self.model_param_uncertainty_array_list.copy()
        else:
            return self.model_param_uncertainty_array_list.copy()

    def extract_scandata(self, filtered=True,
                         fit_result_scan_coord_override=None):
        if fit_result_scan_coord_override is not None:
            self.parse_fit_scandata(fit_result_scan_coord_override)
            fit_result_scan_coord = fit_result_scan_coord_override
        else:
            fit_result_scan_coord = self.model.fit_result_scan_coord
        # extract actual data from model, fields are model output params
        params = self.get_model_param_array_list(filtered)  # param 1 = x
        param_sigmas = \
            self.get_model_param_uncertainty_array_list(filtered)
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
        for scandata in self.get_scandata_list(filtered):
            info_copy = deepcopy(scandata.info)  # take a snapshot
            new_scaninfo['fit_scandata_info_dicts'].append(info_copy)
        return ScanData(field_names,
                        field_arrays,
                        new_scaninfo,
                        x_field_name=None,  # 1st field is already x-coord!
                        y_field_name=None)  # defaults to first model param


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class ScanDataSetsAnalyzer:
    """
    Object that takes ownership of a given list of ScanDataSets, analyzing
    them in parallel and potentially condensing them into ScanData representing
    fit results of entire ScanDataSets.
    """
    def __init__(self, scandataset_list):
        if any(not isinstance(item, ScanDataSet) for item in scandataset_list):
            raise ValueError("ScanDataSetsAnalyzer: requires " +
                             "ScanDataSet list as input")
        else:
            self.scandataset_list = list(scandataset_list)

    def fit_all_scandata_to_model(self, multiprocessing=False):
        if not multiprocessing:
            for scandataset in self.scandataset_list:
                scandataset.fit_scandata_to_model()
        else:
            # cobble together all scandata across all and process in parallel
            # if using weights, must kill error = 0 pts, will ruin weighted fit
            scandata_fit_list = []
            merged_scandata_list = []
            for scandataset in self.scandataset_list:
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
            for scandataset in self.scandataset_list:
                scandataset.parse_fit_scandata()

    def regroup_scandatasets(self, new_model, sort_key="Filepath"):
        """
        Resort ScanDataSets into different ScanDataSets using the new
        model and sort key. If FitData is still valid for new model
        (i.e. same parameters in same order, same output fields),
        can still collapse new sets into scandata.
        
        Note all ScanDataSets must share the same new model!
        """
        self.scandataset_list = regroup_scandatasets(self.scandataset_list,
                                                     new_model, sort_key)

    def break_up_repeating_scandatasets(self):
        """
        break into smaller scandatasets based on repeated xval patterns
        e.g. {x=1,2,3,4,1,2,3,4} -> {x=1,2,3,4} ; {x=1,2,3,4}
        WARNING: may break up more complex patterns too much.
        be careful of use.
        """
        new_scandataset_list = []
        for scandataset in self.scandataset_list:
            # now, group by x-values
            xvals = [scandata.info[scandataset.fit_result_scan_coord]
                     for scandata in scandataset.scandata_list]
            repeat_counts = [xvals[:ind].count(xvals[ind])
                             for ind in range(len(xvals))]
            sublist_index_lists = []
            # if xvals can be broken into continuous sublists w/ identical sets
            # of values, do so! this test catches even sublists w/ different
            # orderings. Note: this requires a properly sorted scandata list!
            if max(repeat_counts) > 0 and \
                                        repeat_counts == sorted(repeat_counts):
                for sublist_index in range(max(repeat_counts) + 1):
                    index_list = [ind
                                  for ind, count in enumerate(repeat_counts)
                                  if count == sublist_index]
                    sublist_index_lists.append(index_list)
                for new_subdir_indices in sublist_index_lists:
                    new_scandataset = ScanDataSet(scandataset.setname,
                                                  scandataset.model)
                    new_scandataset.scandata_list = \
                        [scandataset.scandata_list[ind]
                         for ind in new_subdir_indices]
                    new_scandataset.fit_result_scan_coord = \
                                            scandataset.fit_result_scan_coord
                    new_scandataset_list.append(new_scandataset)
            else:
                new_scandataset_list.append(scandataset)
        num_old = sum(len(subdir.scandata_list)
                      for subdir in self.scandataset_list)
        num_new = sum(len(subdir.scandata_list)
                      for subdir in new_scandataset_list)
        if num_old == num_new:
            self.scandataset_list = new_scandataset_list
        else:
            print("Error: lost {} scandata breaking up subdirectories"\
                  .format(num_old - num_new))

    def apply_transform_to_all_scandata(self, transform_fcn, *args, **kwargs):
        for scandataset in self.scandataset_list:
            scandataset.apply_transform_to_scandata(transform_fcn,
                                                    *args, **kwargs)

    def clear_filters_from_each_scandataset(self):
        for scandataset in self.scandataset_list:
            scandataset.clear_scandata_filters()

    def add_filter_to_each_scandataset(self, filter_fcn):
        for scandataset in self.scandataset_list:
            scandataset.add_scandata_filter(filter_fcn)

    def collapse_to_scandata_list(self, filtered=True):
        return sum((scandataset.get_scandata_list(filtered)
                    for scandataset in self.scandataset_list),
                   [])

    def collapse_to_model_fit_scandata_list(self, filtered=True,
                                            fit_result_scan_coord_override=None,
                                            ignore_empty_scandatasets=True):
        scandata_list = [scandataset.extract_scandata(filtered,
                                                      fit_result_scan_coord_override)
                         for scandataset in self.scandataset_list]
        if ignore_empty_scandatasets:
            return [scandata for scandata in scandata_list
                    if scandata is not None
                    if len(scandata.x) > 0]
        else:
            return scandata_list

#==============================================================================
#     def extract_model_key_field(self, key_field_fcn_name, filtered=True):
#         """
#         Pulls out the desired model key_field from each ScanDataSet,
#         returning the following lists:
#         1. a list containing the key_field dataseries from each ScanDataSet
#         2. a list containing the corresponding uncertainty dataseries from
#            each ScanDataSet
#         3. a list containing [references to] each ScanDataSet. Note these
#            are mutable, and liable to changes from elsewhere!
#         """
#         dataseries_list = []
#         dataseries_uncertainty_list = []
#         scandataset_list = []
#         for scandataset in self.scandataset_list:
#             key_field_fcn = \
#                 scandataset.model.__getkey_field__(key_field_fcn_name)
#             dataseries, uncertainty_dataseries = key_field_fcn(
#                 scandataset.model_param_array_list,
#                 scandataset.model_param_uncertainty_array_list)
#             dataseries_list.append(dataseries)
#             dataseries_uncertainty_list.append(uncertainty_dataseries)
#             scandataset_list.append(scandataset)
#         return dataseries_list, dataseries_uncertainty_list, scandataset_list
#==============================================================================

#==============================================================================
#     def fit_model_key_field(self, key_field_fcn_name,
#                             scandata_filtered=True,
#                             scandataset_filtered=True):
#         """
#         [like extract_model_key_field, but returns a single dataseries
#          corresponding to the fit of model key_fields instead of the
#          key_fields themselves. to get those, just use
#          extract_model_key_field]
#         """
#         raise NotImplementedError("fit_model_key_field not implemented yet")
# 
#==============================================================================


# %%
def sort_scandata_into_sets(scandata_list, model, sort_key="Filepath"):
    """
    Takes a list of ScanData and sorts them into sets based on the sort
    key given. For example, if sort_key is "Voltage", will sort each
    ScanData into sets with shared scaninfo entry {..."Voltage": X, ...}
    where X is different for each set.

    Returns a list of ScanDataSets, each with the analysis model provided.
    """
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