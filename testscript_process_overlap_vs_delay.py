# -*- coding: utf-8 -*-
"""
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
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer, ScanDataSet, GaussianModel


# %% ALTERNATE SCRIPT
if __name__ == "__main__":

# %% CRUNCH DATA
    analyzer = ScanDataSetsAnalyzer(GaussianModel(uncertainty_level=0.01,
                                                  max_fcn_evals=20000),
                                    "C:\\Data\\early_may_good_data")
    analyzer.break_up_repeating_scandatasets()
    analyzer.fit_all_scandata_to_model(multiprocessing=False)
    centers_list, centers_sigma_list, scandataset_list = \
        analyzer.extract_model_attribute("gaussian_centers")


# %%


# %% FILTER AND PLOT
    x_list = []
    y_list = []
    yerr_list = []
    plt.figure()
    plt.hold(True)
    use_flag = True
    for centers, centers_sigma, scandataset \
                    in zip(centers_list, centers_sigma_list, scandataset_list):
        x_vals, y_vals = centers.datalists()
        _, y_errs = centers_sigma.datalists()
        for x, y, yerr in zip(x_vals, y_vals, y_errs):
            use_flag = True
            if yerr > 5:
                use_flag = False
            if scandataset.scandata_list[0].scaninfo['Voltage'] != 0:
                use_flag = False
            if use_flag:
                x_list.append(x)
                y_list.append(y)
                yerr_list.append(yerr)
                plt.errorbar(x_list, y_list, yerr=yerr_list, fmt='.')
        
#    plt.errorbar(x_list, y_list, yerr=yerr_list, fmt='.')
