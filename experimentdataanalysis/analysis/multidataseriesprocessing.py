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
def scandata_iterable_fit(scandata_iterable, field_index, fitfunction,
                          free_params, initial_params, param_bounds,
                          max_fcn_evals=20000, multiprocessing=False):
    """
    """
    scandata_list = list(scandata_iterable)
    new_fitdata_list = dataseries_iterable_fit(
                                [scandata.dataseries_list[field_index]
                                 for scandata in scandata_list],
                                fitfunction, free_params,
                                initial_params, param_bounds,
                                [scandata.error_dataseries_list[field_index]
                                 for scandata in scandata_list],
                                max_fcn_evals, multiprocessing)
    new_scandata_list = []
    for ind, scandata in enumerate(scandata_list):
        new_scandata_fitdata_list = list(scandata.fitdata_list)
        new_scandata_fitdata_list[field_index] = new_fitdata_list[ind]
        new_scandata_list.append(ScanData(scandata.fields,
                                          [scaninfo.copy() for scaninfo in \
                                                      scandata.scaninfo_list],
                                          scandata.dataseries_list,
                                          scandata.error_dataseries_list,
                                          new_scandata_fitdata_list))
                                  
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

# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_iterable_sort(scandata_iterable, field_index,
                           primary_key, secondary_key):
    """
    """
    # Subfunction to use as sort key
    def scandatasortfcn(scandata_2_tuple):
        scandata, _ = scandata_2_tuple
        scaninfo = scandata.scaninfo_list[field_index]
        try:
            return (float(scaninfo[primary_key]),
                    float(scaninfo[secondary_key]))
        except KeyError:
            try:
                return (float(scaninfo[primary_key]),
                        99999999)
            except KeyError:
                try:
                    return (99999999,
                            float(scaninfo[secondary_key]))
                except KeyError:
                    return (99999999, 99999999)
                except ValueError:
                    print("scandata_iterable_sort: Numerical sort keys only!")
                    return (99999999, 99999999)
            except ValueError:
                print("scandata_iterable_sort: Numerical sort keys only!")
                return (99999999, 99999999)
        except ValueError:
            print("scandata_iterable_sort: Numerical sort keys only!")
            return (99999999, 99999999)
        except AttributeError:
            print("scandata_iterable_sort: ScanData expected as list element.")
            return (99999999, 99999999)

    scandata_list = list(scandata_iterable)
    index_ordering = range(len(scandata_list))
    scandata_list, index_ordering = zip(*sorted(
                                        zip(scandata_list, index_ordering),
                                        key=scandatasortfcn))
    return scandata_list, index_ordering
