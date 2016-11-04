# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:35:06 2016

Handles fitting lists or iterables of scandata/dataseries, including
multiprocessing support. Calls dataseriesprocessing.py for actual crunching
of data, even when farming it out to multiprocessing_map.

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
                          max_fcn_evals=20000, multiprocessing=False,
                          excluded_intervals=None):
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
                                max_fcn_evals, multiprocessing,
                                excluded_intervals)
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
                            max_fcn_evals=20000, multiprocessing=False,
                            excluded_intervals=None):
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
                        param_bounds, weights_dataseries,
                        max_fcn_evals, excluded_intervals]
                       for dataseries, weights_dataseries in
                       zip(dataseries_list, weights_dataseries_list)]
    output_list = multiprocessable_map(dataseries_fit,
                                       input_args_list, multiprocessing)
    return output_list

# %% NEEDS TEST, SPHINX DOCUMENTATION
def scandata_iterable_sort(scandata_iterable, field_index,
                           primary_key, secondary_key, numeric_sort=True):
    """
    right now, no mixed numeric/non-numeric keys, because effort
    """
    # Subfunction to use as sort key
    def scandatasortfcn_strings(scandata_2_tuple):
        scandata, _ = scandata_2_tuple
        scaninfo = scandata.scaninfo_list[field_index]
        try:
            return (str(scaninfo[primary_key]),
                    str(scaninfo[secondary_key]))
        except KeyError:
            try:
                return (str(scaninfo[primary_key]),
                        "")
            except KeyError:
                try:
                    return ("",
                            str(scaninfo[secondary_key]))
                except KeyError:
                    return ("", "")
        except AttributeError:
            print("scandata_iterable_sort: ScanData expected as list element.")
            return ("", "")

    # Subfunction to use as sort key
    def scandatasortfcn_numerical(scandata_2_tuple):
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
                    print("scandata_iterable_sort: numerical_sort flag on, " +
                          "numerical sort keys only!")
                    return (99999999, 99999999)
            except ValueError:
                print("scandata_iterable_sort: numerical_sort flag on, " +
                      "numerical sort keys only!")
                return (99999999, 99999999)
        except ValueError:
            print("scandata_iterable_sort: numerical_sort flag on, " +
                  "numerical sort keys only!")
            return (99999999, 99999999)
        except AttributeError:
            print("scandata_iterable_sort: ScanData expected as list element.")
            return (99999999, 99999999)

    scandata_list = list(scandata_iterable) 
    index_ordering = range(len(scandata_list))
    if numeric_sort:
        key_fcn = scandatasortfcn_numerical
    else:
        key_fcn = scandatasortfcn_strings

    scandata_list, index_ordering = zip(*sorted(
                                        zip(scandata_list, index_ordering),
                                        key=key_fcn))
    return scandata_list, index_ordering


# %%
def aggregate_and_process_scandata_into_dict(scandata, process_fcn_list=None,
                                             xcoord_indices_to_use=None,
                                             xcoord_name=None):
    """
    Accepts a scandata, and converts each dataseries into a dict for easy
    lookup of results by field name. Optionally takes a list of process
    functions - where each takes the xvals, yvals & error of a given field
    and outputs a dict of analysis results - and appends the results of those
    functions to the aggregated data dict, with the field name added to the
    front of each process function key, and repeats for all fields in scandata
        (e.g. "linear_fit_slope" -> "[field name 1]_linear_fit_slope",
                                    "[field name 2]_linear_fit_slope", ...)
    """
    if xcoord_indices_to_use is None:
        xcoord_indices_to_use = list(range(len(scandata.fields)))
    if xcoord_name is None:
        xcoord_name = "xvalues"

    xvalues = scandata.dataseries_list[0].xvals()[xcoord_indices_to_use]

    # construct dictionary containing 1d arrays of field values (vs xvalues):
    aggregated_data = {'xcoord_name': xcoord_name,
                       xcoord_name: xvalues,
                       'field_names': list(scandata.fields).copy()}
    aggregated_data.update(
        zip(scandata.fields,
            [series.yvals()[xcoord_indices_to_use]
             for series in scandata.dataseries_list]))
    # append fit parameter uncertainties:
    field_param_error_names = \
        [param_name + "_error" for param_name in scandata.fields]
    aggregated_data.update(
        zip(field_param_error_names,
            [series.yvals()[xcoord_indices_to_use]
             for series in scandata.error_dataseries_list]))
    # perform process fit to attributes
    if process_fcn_list is not None:
        for field_name in scandata.fields:
            for process_fcn in process_fcn_list:
                process_result = \
                    process_fcn(xvalues,
                                aggregated_data[field_name],
                                aggregated_data[field_name + '_error'])
                for key, value in process_result.items():
                    new_key = field_name + "_" + key
                    aggregated_data[new_key] = value

    return aggregated_data

