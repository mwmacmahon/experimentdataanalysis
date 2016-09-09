# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:47:08 2016

@author: Michael
"""

from experimentdataanalysis.analysis.dataclasses \
    import ScanData, DataSeries
from experimentdataanalysis.analysis.multidataseriesprocessing \
    import scandata_iterable_fit, scandata_iterable_sort
from experimentdataanalysis.analysis.scandatamodels \
    import ScanDataModel
import experimentdataanalysis.parsing.dataseriesparsing as dsparsing


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class ScanDataSet:
    """
    Stores a list of ScanData, a model used to fit them, and the results of
    said fits once they are performed. Note this is a mutable data structure,
    and ScanData themselves are mutable (they will, for example, gain a
    FitData attribute after fitting).
    
    add_filter / add_filters add filters that limit which scandatasets are
    returned by get_
    """
    def __init__(self, setname, model):
        self.setname = setname
        self.model = model
        self.scandata_list = []
        self.scandata_filter_fcn_list = []
        self.chronological_key = "FastScanIndex"
        self.dim2_key = model.dim2_key
        self.model_param_dataseries_list = None
        self.model_param_uncertainty_dataseries_list = None
        if model is not None:
            self.scandata_filter_fcn_list += model.get_model_filter_fcns()

    def sort_scandata_list(self):
        if len(self.scandata_list) > 1:
            # fix scandataset ordering errors arising from alphanumeric sorting
            primary_key = self.chronological_key
            secondary_key = self.dim2_key
            scandata_list_tuple, _ = \
                scandata_iterable_sort(self.scandata_list,
                                       self.model.field_index,
                                       primary_key, secondary_key)
            self.scandata_list = list(scandata_list_tuple)

    def fit_scandata_to_model(self):
        # must filter out any pts with error = 0, will kill weighted fit
        field_index = self.model.field_index
        self.purge_zero_error_scandata(field_index)
        fit_scandata_list = \
                scandata_iterable_fit(self.scandata_list,
                                      field_index,
                                      self.model.fitfunction,
                                      self.model.free_params,
                                      self.model.initial_params,
                                      self.model.param_bounds,
                                      self.model.max_fcn_evals,
                                      multiprocessing=False)
        self.parse_fit_scandata_list(self.model.field_index,
                                     fit_scandata_list)

    def purge_zero_error_scandata(self, field_index):
        def has_zero_error(scandata):
            return 0 in scandata.error_dataseries_list[field_index].yvals()

        fail_message = ("Warning: ScanData value found with 0 uncertainty " +
                        "in  DataSeries, so ScanData was thrown out as " +
                        "incompatible with weighted fit")
        self.purge_scandata_list(field_index, has_zero_error, fail_message)

    def purge_failed_fit_scandata(self, field_index):
        def has_no_fitdata(scandata):
            return scandata.fitdata_list[field_index] is None

        fail_message = ("Warning: Fit failed in ScanDataSet, erasing " +
                        "associated ScanData so all ScanData have " +
                        "associated fits during further analysis.")
        self.purge_scandata_list(field_index, has_no_fitdata, fail_message)

    def purge_scandata_list(self, test_field_index,
                            scandata_failed_test_fcn, fail_message):
        filtered_scandata_list = []
        for scandata in self.scandata_list:
            if not scandata_failed_test_fcn(scandata):
                filtered_scandata_list.append(scandata)
            else:
                print(fail_message)
                if 'Filepath' in scandata.scaninfo_list[test_field_index]:
                    print("Associated filepath:")
                    print(scandata.scaninfo_list[test_field_index]['Filepath'])
        self.scandata_list = filtered_scandata_list
        self.model_param_dataseries_list = None  # now invalid!
        self.model_param_uncertainty_dataseries_list = None

    def parse_fit_scandata_list(self, field_index, fit_scandata_list):
        # purge failed fits, as they will cause issues further down.
        self.scandata_list = fit_scandata_list
        self.purge_failed_fit_scandata(field_index)
        # all remaining scandata should have FitData and be processable:
        # create dataseries of fit parameters & their uncertainties vs xval:
        dim2_key = self.dim2_key
        xvals = [scandata.scaninfo_list[0][dim2_key]
                 for scandata in fit_scandata_list]
        fitparams = []
        fitparamstds = []
        fitparams = list(zip(*(scandata.fitdata_list[field_index].fitparams
                               for scandata in self.scandata_list)))
        fitparams_dataseries_list = [DataSeries(zip(xvals, fitparam))
                                     for fitparam in fitparams]
        fitparamstds = list(zip(*(
            scandata.fitdata_list[field_index].fitparamstds
            for scandata in self.scandata_list
                                  )))
        fitparamstds_dataseries_list = [DataSeries(zip(xvals, fitparamstd))
                                        for fitparamstd in fitparamstds]

        # save these results to scandataset:
        self.model_param_dataseries_list = fitparams_dataseries_list
        self.model_param_uncertainty_dataseries_list = \
                                                fitparamstds_dataseries_list

    def apply_transform_to_scandata(self, transform_fcn, *args, **kwargs):
        new_scandata_list = []
        for scandata in self.scandata_list:
            new_scandata_list.append(transform_fcn(scandata,
                                                   *args, **kwargs))
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

    def get_model_param_dataseries_list(self, filtered=True):
        if filtered:
            filtered_indices = self.get_filtered_scandata_indices()
            if len(filtered_indices) > 0:
                return [series.copy_subset(filtered_indices)
                        for series in self.model_param_dataseries_list]
            else:
                print('Warning: all scandata have been filtered out! ' +
                      'get_model_param_dataseries_list filters ignored.')
                return self.model_param_dataseries_list.copy()
        else:
            return self.model_param_dataseries_list.copy()

    def get_model_param_uncertainty_dataseries_list(self, filtered=True):
        if filtered:
            filtered_indices = self.get_filtered_scandata_indices()
            if len(filtered_indices) > 0:
                return [series.copy_subset(filtered_indices)
                        for series in \
                            self.model_param_uncertainty_dataseries_list]
            else:
                print('Warning: all scandata have been filtered out! ' +
                      'get_model_param_uncertainty_dataseries_list filters ' +
                      'ignored.')
                return self.model_param_uncertainty_dataseries_list.copy()
        else:
            return self.model_param_uncertainty_dataseries_list.copy()

    def extract_scandata(self, new_scan_type, filtered=True):
        params = self.get_model_param_dataseries_list(filtered)
        param_sigmas = \
            self.get_model_param_uncertainty_dataseries_list(filtered)
        fields, field_dataseries_list, field_uncertainty_dataseries_list = \
            self.model.all_model_fields(params, param_sigmas)
        if len(self.get_scandata_list(filtered)) > 0:
            scandata_for_scaninfo = self.get_scandata_list(filtered)[0]
        else:  # even if no scandata survive filters, still grab a scaninfo
            scandata_for_scaninfo = self.scandata_list[0]
        scaninfo_list = [scandata_for_scaninfo.scaninfo_list[0].copy()
                         for field in fields]
        # update scaninfo to reflect new coord types:
        updated_scaninfo_list = []
        for scaninfo in scaninfo_list:
            new_scaninfo = scaninfo.copy()  # rule #1: don't modify in place!
            new_scaninfo['FastScanType'] = scaninfo['MiddleScanType']
            new_scaninfo['MiddleScanType'] = new_scan_type
            try:
                new_scaninfo['MiddleScanCoord'] = scaninfo[new_scan_type]
            except KeyError:
                new_scaninfo['MiddleScanCoord'] = "[Unknown]"
            updated_scaninfo_list.append(new_scaninfo)
        return ScanData(fields,
                        updated_scaninfo_list,
                        field_dataseries_list,
                        field_uncertainty_dataseries_list,
                        [None for field in fields])

# %% NEEDS TESTS, SPHINX DOCUMENTATION
class ScanDataSetsAnalyzer:
    """
    Load a directory into a list of ScanDataSets. Each csv/tsv file will be
    added to a ScanDataSet corresponding to its parent directory. Can also
    simply load in a series of scandata, but they will be sorted by

    Note this means in the case of nested directories, a parent directory's
    related ScanDataSet will only contain the free csv/tsv files directly
    within that directory, not those belonging to its subdirectories.
    """
    def __init__(self, model, dirpath="", *, scandata_list=[],
                 uncertainty_value=None, set_key="Filepath"):
        self.set_model = model
        self.scandataset_list = []
        if scandata_list:
            scandata_list = list(scandata_list)
        elif dirpath:
            scandata_list = \
                list(dsparsing.fetch_dir_as_unfit_scandata_iterator(
                     directorypath=dirpath))
        else:
            raise TypeError("No data source provided to ScanDataSetsAnalyzer")
        if uncertainty_value:
            scandata_list = [dsparsing.set_scandata_error(scandata,
                                                          model.field_index,
                                                          uncertainty_value)
                             for scandata in scandata_list]
        self.load_scandata_list(scandata_list, model, set_key)

    def load_scandata_list(self, scandata_list, model, set_key="Filepath"):
        if len(scandata_list) == 0:
            raise TypeError("Empty scandata list provided to " +
                            "ScanDataSetsAnalyzer.")
        try:  # if numerical values for this sort key, pre-sort scandata list
            for scandata in scandata_list:
                float(scandata.scaninfo_list[0][set_key])
        except KeyError:
            if set_key is None:
                pass  # all scandata put into one set, ignore sorting
            else:
                raise KeyError("Invalid set_key for ScanDataSetsAnalyzer")
        except ValueError:
            scandata_list, _ = scandata_iterable_sort(scandata_list, 0,
                                                      set_key, set_key,
                                                      numeric_sort=False)
        else:
            scandata_list, _ = scandata_iterable_sort(scandata_list, 0,
                                                      set_key, set_key,
                                                      numeric_sort=True)
        scandata_list = list(scandata_list)
        current_scandataset = ScanDataSet("[default]", ScanDataModel())
        last_setname = ""
        for scandata in scandata_list:
            try:  # use sort key's value to group ScanData into ScanDataSets
                setname = scandata.scaninfo_list[0][set_key]
                if set_key == "Filepath":  # special case
                    setname = setname.split("\\")[-2]
            except KeyError:
                if set_key == None:
                    setname = "All ScanData"                    
                else:
                    print("Warning: ScanDataSet sort key not found in " +
                          "ScanData! Sets may not be properly grouped")
                    setname = "key_error_set"
            if setname != last_setname:
                current_scandataset.sort_scandata_list()
                current_scandataset = ScanDataSet(setname, model)
                self.scandataset_list.append(current_scandataset)
                last_setname = setname
            current_scandataset.scandata_list.append(scandata)
        current_scandataset.sort_scandata_list()  # sort last scandataset        

    def fit_all_scandata_to_model(self, multiprocessing=False):
        if not multiprocessing:
            for scandataset in self.scandataset_list:
                scandataset.fit_scandata_to_model()
        else:
            model = self.scandataset_list[0].model
            field_index = model.field_index
            model_list = [scandataset.model
                          for scandataset in self.scandataset_list]
            if any(ref_model != model for ref_model in model_list):
                raise TypeError("Error: multiprocessing in a " +
                                "ScanDataSetAnalyzer requires " +
                                "all ScanDataSets to share an " +
                                "identical model")
            # cobble together a giant list of scandata, process
            # together, then redistribute to scandatasets...
            # must filter out any pts with error = 0, will kill weighted fit
            merged_scandata_list = []
            for scandataset in self.scandataset_list:
                scandataset.purge_zero_error_scandata(field_index)
                merged_scandata_list += list(scandataset.scandata_list)
            merged_fit_scandata_list = \
                scandata_iterable_fit(merged_scandata_list,
                                      model.field_index,
                                      model.fitfunction,
                                      model.free_params,
                                      model.initial_params,
                                      model.param_bounds,
                                      model.max_fcn_evals,
                                      multiprocessing=True)
            start_ind = 0
            for scandataset in self.scandataset_list:
                end_ind = start_ind + len(scandataset.scandata_list)
                fit_scandata_list = \
                    merged_fit_scandata_list[start_ind:end_ind]
                scandataset.parse_fit_scandata_list(field_index,
                                                    fit_scandata_list)
                start_ind = end_ind

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
            xvals = [scandata.scaninfo_list[0][scandataset.dim2_key]
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
                                                  self.set_model)
                    new_scandataset.scandata_list = \
                        [scandataset.scandata_list[ind]
                         for ind in new_subdir_indices]
                    new_scandataset.dim2_key = scandataset.dim2_key
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

    def collapse_to_model_fit_scandata_list(self, new_scan_type,
                                            filtered=True,
                                            ignore_empty_scandatasets=True):
        scandata_list = [scandataset.extract_scandata(new_scan_type,
                                                      filtered)
                         for scandataset in self.scandataset_list]
        if ignore_empty_scandatasets:
            return [scandata for scandata in scandata_list
                    if len(scandata.dataseries_list[0]) > 0]
        else:
            return scandata_list

#==============================================================================
#     def extract_model_attribute(self, attribute_fcn_name, filtered=True):
#         """
#         Pulls out the desired model attribute from each ScanDataSet,
#         returning the following lists:
#         1. a list containing the attribute dataseries from each ScanDataSet
#         2. a list containing the corresponding uncertainty dataseries from
#            each ScanDataSet
#         3. a list containing [references to] each ScanDataSet. Note these
#            are mutable, and liable to changes from elsewhere!
#         """
#         dataseries_list = []
#         dataseries_uncertainty_list = []
#         scandataset_list = []
#         for scandataset in self.scandataset_list:
#             attribute_fcn = \
#                 scandataset.model.__getattribute__(attribute_fcn_name)
#             dataseries, uncertainty_dataseries = attribute_fcn(
#                 scandataset.model_param_dataseries_list,
#                 scandataset.model_param_uncertainty_dataseries_list)
#             dataseries_list.append(dataseries)
#             dataseries_uncertainty_list.append(uncertainty_dataseries)
#             scandataset_list.append(scandataset)
#         return dataseries_list, dataseries_uncertainty_list, scandataset_list
#==============================================================================

#==============================================================================
#     def fit_model_attribute(self, attribute_fcn_name,
#                             scandata_filtered=True,
#                             scandataset_filtered=True):
#         """
#         [like extract_model_attribute, but returns a single dataseries
#          corresponding to the fit of model attributes instead of the
#          attributes themselves. to get those, just use
#          extract_model_attribute]
#         """
#         raise NotImplementedError("fit_model_attribute not implemented yet")
# 
#==============================================================================
