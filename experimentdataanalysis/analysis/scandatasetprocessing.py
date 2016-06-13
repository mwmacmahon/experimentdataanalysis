# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:47:08 2016

@author: Michael
"""

from itertools import repeat
import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
import experimentdataanalysis.analysis.dataclassfitting as dcfitting
import experimentdataanalysis.analysis.dataclassgraphing as dcgraphing
import experimentdataanalysis.analysis.dataseriesprocessing as dsprocessing
from experimentdataanalysis.analysis.fitfunctions \
    import fitfcn_simple_1d_gaussian, fitfcn_single_exp_decay, \
        fitfcn_two_exp_decay, fitfcn_exp_times_sqrt_decay, \
        fitfcn_simple_line, fitfcn_1d_gaussian_with_linear_offset
from experimentdataanalysis.analysis.multidataseriesprocessing \
    import dataseries_iterable_fit, scandata_iterable_fit, \
        scandata_iterable_sort
import experimentdataanalysis.parsing.dataclassparsing as dcparsing


# %% NEEDS TESTS, SPHINX DOCUMENTATION
class GaussianModel:
    """
    Stores the parameters and formulas needed to fit scandata to a
    1D model and extract attributes. Should not store any data, although
    it can be modified after creation (e.g. to increase max_fcn_evals
    attribute or adjust the experimental uncertainty assumed).

    Can be inherited from to change the fit model, may only need to change
    __init__ method if first 3 fit parameters are left the same.
    """
    def __init__(self, **kwargs):
        self.field_index = 0
        self.fitfunction = fitfcn_1d_gaussian_with_linear_offset
        # params = amplitude, x0, sigma, slope, offset
        self.free_params = [True, True, True, True, True]
        self.initial_params = [0.001, 0, 20, 0, 0]
        self.param_bounds = [(0, 1), (-200, 200),
                             (5, 200), (-0.01, 0.01), (-0.01, 0.01)]
        self.uncertainty_level = 0.01
        self.max_fcn_evals = 20000
        # keyword args override defaults
        for key, val in kwargs.items():
            self.__dict__[key] = val

    def gaussian_amplitudes(self, model_param_dataseries_list,
                            model_param_uncertainty_dataseries_list):
        amplitudes = model_param_dataseries_list[0]
        amplitudes_sigma = model_param_uncertainty_dataseries_list[0]
        return amplitudes, amplitudes_sigma

    def gaussian_widths(self, model_param_dataseries_list,
                        model_param_uncertainty_dataseries_list):
        widths = model_param_dataseries_list[2]
        widths_sigma = model_param_uncertainty_dataseries_list[2]
        return widths, widths_sigma

    def gaussian_centers(self, model_param_dataseries_list,
                         model_param_uncertainty_dataseries_list):
        centers = model_param_dataseries_list[1]
        centers_sigma = model_param_uncertainty_dataseries_list[1]
        return centers, centers_sigma

    def gaussian_areas(self, model_param_dataseries_list,
                       model_param_uncertainty_dataseries_list):
        amplitudes = model_param_dataseries_list[0]
        amplitudes_sigma = model_param_uncertainty_dataseries_list[0]
        widths = model_param_dataseries_list[2]
        widths_sigma = model_param_uncertainty_dataseries_list[2]
        # note: uncertainty calc. assumes width, amplitude independent...
        areas = amplitudes*widths
        areas_sigma = \
            ((amplitudes*widths_sigma)**2 + (widths*amplitudes_sigma)**2)**0.5
        return areas, areas_sigma


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
        self.xval_key = "MiddleScanCoord"
        self.model_param_dataseries_list = None
        self.model_param_uncertainty_dataseries_list = None

    def sort_scandata_list(self):
        if len(self.scandata_list) > 1:
            # fix scandataset ordering errors arising from alphanumeric sorting
            primary_key = self.chronological_key
            secondary_key = self.xval_key
            self.scandata_list, _ = \
                scandata_iterable_sort(self.scandata_list,
                                       primary_key, secondary_key)

    def fit_scandata_to_model(self):
        self._fit_scandata_in_set(self.model.field_index,
                                 self.model.fitfunction,
                                 self.model.free_params,
                                 self.model.initial_params,
                                 self.model.param_bounds,
                                 self.model.uncertainty_level,
                                 self.model.max_fcn_evals)

    def _fit_scandata_in_set(self, field_index, scandata_fitfunction,
                             free_params, initial_params, param_bounds,
                             uncertainty_level=0.01, max_fcn_evals=20000):
        scan_uncertainties_list = \
            [DataSeries(zip(scandata.dataseries[0].xvals(raw=True),
                            repeat(uncertainty_level)),
                        excluded_intervals = \
                            scandata.dataseries[0].excluded_intervals())
             for scandata in self.scandata_list]
        fit_scandata_list = scandata_iterable_fit(
                                self.scandata_list,
                                field_index, scandata_fitfunction,
                                free_params, initial_params, param_bounds,
                                scan_uncertainties_list,
                                max_fcn_evals, multiprocessing=False)

        # create dataseries of fit parameters & their uncertainties vs xval:
        xval_key = self.xval_key
        xvals = [scandata.scaninfo[xval_key] for scandata in fit_scandata_list]
        fitparams = list(zip(*(scandata.fitdata[field_index].fitparams
                               for scandata in fit_scandata_list)))
        fitparams_dataseries_list = [DataSeries(zip(xvals, fitparam))
                                     for fitparam in fitparams]
        fitparamstds = list(zip(*(scandata.fitdata[field_index].fitparamstds
                                  for scandata in fit_scandata_list)))
        fitparamstds_dataseries_list = [DataSeries(zip(xvals, fitparamstd))
                                        for fitparamstd in fitparamstds]

        # save these results to scandataset:
        self.scandata_list = fit_scandata_list
        self.model_param_dataseries_list = fitparams_dataseries_list
        self.model_param_uncertainty_dataseries_list = \
                                                fitparamstds_dataseries_list

    def add_scandata_filters(self, *scandata_filter_fcn_list):
        for scandata_filter_fcn in scandata_filter_fcn_list:
            self.add_scandata_filter(scandata_filter_fcn)

    def add_scandata_filter(self, scandata_filter_fcn):
        """
        Adds a filter to an internally kept filter list. Filters are functions
        that accept a scandata and return a boolean - whether scandata passes
        filter or not. E.g. a simple filter might be:
        
        def only_fitted_scandata_filter(scandata):
            return scandata.fitdata != None
        
        Note this is usually done AFTER fitting all scandata, so fit results
        can be used in filtering.
        """
        # insert checks for valid filter fcns?
        self.scandata_filter_fcn_list.append(scandata_filter_fcn)

    def clear_scandata_filters(self):
        self.scandata_filter_fcn_list = []

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
        return [scandata
                for ind, scandata in enumerate(self.scandata_list)
                if ind in self.get_filtered_scandata_indices()]

    def get_model_param_dataseries_list(self, filtered=True):
        return [dataseries
                for ind, dataseries in enumerate(
                    self.model_param_dataseries_list)
                if ind in self.get_filtered_scandata_indices()]

    def get_model_param_uncertainty_dataseries_list(self, filtered=True):
        return [dataseries
                for ind, dataseries in enumerate(
                    self.model_param_uncertainty_dataseries_list)
                if ind in self.get_filtered_scandata_indices()]

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
    def __init__(self, model, dirpath="", *, scandata_list=[]):
        self.set_model = model
        self.scandataset_list = []
        self.scandataset_filter_fcn_list = []
        if scandata_list:
            self.load_scandata_list(list(scandata_list), self.set_model)
        elif dirpath:
            scandata_list = \
                list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
                     directorypath=dirpath))
            self.load_scandata_list(scandata_list, self.set_model)

    def load_scandata_list(self, scandata_list, model):
        current_scandataset = ScanDataSet("[default]", None)
        last_setname = ""
        for scandata in scandata_list:
            setname = scandata.filepath.split("\\")[-2]
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
            # speed up by sending unprocessed scandatasets and returning
            # processed scandatasets - not super efficient, so may have to
            # instead cobble together a giant list of scandata, process
            # together, then redistribute to scandatasets...a huge ordeal
            raise NotImplementedError("Error: multiprocessing not " +
                                      "implemented for ScanDataSets yet.")

    def break_up_repeating_scandatasets(self):
        """
        break into smaller scandatasets based on repeated xval patterns
        e.g. {x=1,2,3,4,1,2,3,4} -> {x=1,2,3,4} ; {x=1,2,3,4}
        """
        new_scandataset_list = []
        for scandataset in self.scandataset_list:
            # now, group by x-values
            xvals = [scandata.scaninfo[scandataset.xval_key]
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
                    new_scandataset.xval_key = scandataset.xval_key
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

    def add_scandataset_filters(self, *scandataset_filter_fcn_list):
        for scandataset_filter_fcn in scandataset_filter_fcn_list:
            self.add_scandataset_filter(scandataset_filter_fcn)

    def add_scandataset_filter(self, scandataset_filter_fcn):
        """
        Adds a filter to an internally kept filter list. Filters are functions
        that accept a scandataset and return a boolean - whether scandataset
        passes filter or not. E.g. a simple filter might be:
        
        def only_fitted_scandataset_filter(scandata):
            return len(scandataset.get_filtered_scandata_indices) > 0
        """
        # insert checks for valid filter fcns?
        self.scandataset_filter_fcn_list.append(scandataset_filter_fcn)

    def clear_scandataset_filters(self):
        self.scandataset_filter_fcn_list = []

    def get_filtered_indices(self):
        """
        Returns a list containing all indices in scandata_list whose
        corresponding scandata meet the requirements off all filters
        in self.scandata_filter_fcn_list
        """
        filtered_indices = []
        for ind, scandataset in enumerate(self.scandataset_list):
            ok_flag = True
            for scandataset_filter_fcn in self.scandataset_filter_fcn_list:
                ok_flag = ok_flag and scandataset_filter_fcn(scandataset)
            if ok_flag:
                filtered_indices.append(ind)
        return filtered_indices

    def get_scandataset_list(self, filtered=True):
        return [scandataset
                for ind, scandataset in enumerate(self.scandataset_list)
                if ind in self.get_filtered_indices()]

    def extract_model_attribute(self, attribute_fcn_name,
                                scandata_filtered=True,
                                scandataset_filtered=True):
        """
        Pulls out the desired model attribute from each ScanDataSet,
        returning the following lists:
        1. a list containing the attribute dataseries from each ScanDataSet
        2. a list containing the corresponding uncertainty dataseries from
           each ScanDataSet
        3. a list containing [references to] each ScanDataSet. Note these
           are mutable, and liable to changes from elsewhere!
        
        Note scandataset_filtered only works if fits have been performed
        on each scandataset! If no fit performed, set will be filtered out.
        """
        dataseries_list = []
        dataseries_uncertainty_list = []
        scandataset_list = []
        for scandataset in self.scandataset_list:
            attribute_fcn = \
                scandataset.model.__getattribute__(attribute_fcn_name)
            dataseries, uncertainty_dataseries = attribute_fcn(
                scandataset.model_param_dataseries_list,
                scandataset.model_param_uncertainty_dataseries_list)
            dataseries_list.append(dataseries)
            dataseries_uncertainty_list.append(uncertainty_dataseries)
            scandataset_list.append(scandataset)
        return dataseries_list, dataseries_uncertainty_list, scandataset_list

    def fit_model_attribute(self, attribute_fcn_name,
                            scandata_filtered=True,
                            scandataset_filtered=True):
        """
        [like extract_model_attribute, but returns a single dataseries
         corresponding to the fit of model attributes instead of the
         attributes themselves. to get those, just use
         extract_model_attribute]
        """
        raise NotImplementedError("fit_model_attribute not implemented yet")
