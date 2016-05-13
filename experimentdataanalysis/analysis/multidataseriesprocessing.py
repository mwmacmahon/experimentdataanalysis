# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:35:06 2016

@author: vsih-lab
"""

from experimentdataanalysis.analysis.dataseriesprocessing \
    import dataseries_fit
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
from experimentdataanalysis.analysis.generalutilities \
    import multiprocessable_map


# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_iterable_fit(scandata_iterable, dataseries_index, fitfunction,
                          free_params, initial_params, param_bounds,
                          max_fcn_evals=20000, multiprocessing=False):
    """
    """
    scandata_list = list(scandata_iterable)
    new_fitdata_list = dataseries_iterable_fit(
                                [scandata.dataseries[dataseries_index]
                                 for scandata in scandata_list],
                                fitfunction, free_params,
                                initial_params, param_bounds, None,
                                max_fcn_evals, multiprocessing)
    new_scandata_list = []
    for ind, scandata in enumerate(scandata_list):
        new_scandata_fitdatalist = list(scandata.fitdata)
        new_scandata_fitdatalist[dataseries_index] = new_fitdata_list[ind]
        new_scandata_list.append(ScanData(scandata.filepath,
                                          scandata.scaninfo.copy(),
                                          scandata.fields,
                                          scandata.dataseries,
                                          tuple(new_scandata_fitdatalist)))
                                  
    return new_scandata_list


# %% NEEDS TEST, SPHINX DOCUMENTATION
def dataseries_iterable_fit(dataseries_iterable, fitfunction,
                            free_params, initial_params, param_bounds,
                            weights_dataseries_iterable=None,
                            max_fcn_evals=20000, multiprocessing=False):
    """
    """
    #package for processing
    dataseries_list = list(dataseries_iterable)
    if weights_dataseries_iterable is None:
        weights_dataseries_list = [None]*len(dataseries_list)
    else:
        weights_dataseries_list = list(weights_dataseries_iterable)
    input_args_list = [[dataseries, fitfunction,
                        free_params, initial_params,
                        param_bounds, weights_dataseries, max_fcn_evals]
                       for dataseries, weights_dataseries in
                       zip(dataseries_list, weights_dataseries_list)]
    output_iter = multiprocessable_map(dataseries_fit,
                                       input_args_list, multiprocessing)

    #unpackage:
    fitdata_list = list(output_iter)
    return fitdata_list
