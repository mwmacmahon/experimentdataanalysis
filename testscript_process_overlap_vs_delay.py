# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import GaussianModel, SpinLifetimeModel, SimplerSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer


# %% ALTERNATE SCRIPT
if __name__ == "__main__":

# %% LOAD AND FIT DATA
    gaussian_model = \
        GaussianModel(max_fcn_evals=80000,
                      error_thresholds=[.01, 20, 20, None, None])
    analyzer = ScanDataSetsAnalyzer(gaussian_model,
#                                    dirpath="C:\\Data\\early_may_good_data",
                                    dirpath="C:\\Data\\early_may_good_data\\160506",
                                    uncertainty_value=5e-4)
    analyzer.break_up_repeating_scandatasets()
    analyzer.fit_all_scandata_to_model(multiprocessing=True)

# %% FILTER DATA
    # Start by clearing existing filters, only necessary on subsequent runs
    analyzer.clear_filters_from_each_scandataset()

    def filter_out_all_but_channel_2_data(scandata):
        return all(scaninfo['Channel'] == 2
                   for scaninfo in scandata.scaninfo_list)

    def filter_out_all_but_channel_3_data(scandata):
        return all(scaninfo['Channel'] == 3
                   for scaninfo in scandata.scaninfo_list)

    def filter_out_nonzero_voltage(scandata):
        return all(scaninfo['Voltage'] == 0
                   for scaninfo in scandata.scaninfo_list)

    # Add all desired filters, not including built-in model fit thresholds
#    analyzer.add_filter_to_each_scandataset(filter_out_nonzero_voltage)
#    analyzer.add_filter_to_each_scandataset(filter_out_all_but_channel_2_data)
#    analyzer.add_filter_to_each_scandataset(filter_out_all_but_channel_3_data)

# %% EXTRACT FILTERED FIT DATA
    gaussian_fits_scandata_list = \
                            analyzer.collapse_to_scandata_list(filtered=True)
    collapsed_gaussian_fit_scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(filtered=True)

# %% LOAD AND FIT GAUSSIAN_FIT DATA TO LIFETIMES
#    lifetime_model = \
#        SimplerSpinLifetimeModel(max_fcn_evals=80000,
#                          error_thresholds=[0.5, 50000, None])
    lifetime_model = \
        SpinLifetimeModel(max_fcn_evals=80000,
                          error_thresholds=[0.2, 200, 0.2, 2000, None],
                          xval_key="Voltage")
    analyzer2 = ScanDataSetsAnalyzer(
        lifetime_model,
        scandata_list=collapsed_gaussian_fit_scandata_list,
        sort_key="Voltage")
    analyzer2.apply_transform_to_all_scandata(get_positive_time_delay_scandata)
    analyzer2.fit_all_scandata_to_model(multiprocessing=True)

# %% FILTER AGAIN?
    # Start by clearing existing filters, only necessary on subsequent runs
    analyzer2.clear_filters_from_each_scandataset()

    def filter_out_under_5_data_points(scandata):
        return all(len(dataseries) >= 5
                   for dataseries in scandata.dataseries_list)

    # Add all desired filters, not including built-in model fit thresholds
    analyzer2.add_filter_to_each_scandataset(filter_out_under_5_data_points)

# %% EXTRACT FILTERED FIT DATA
    lifetime_fit_scandata_list = \
                            analyzer2.collapse_to_scandata_list(filtered=True)
    collapsed_lifetime_fit_scandata_list = \
        analyzer2.collapse_to_model_fit_scandata_list(filtered=True)

# %% TESTING...
#    scandata_list = gaussian_fits_scandata_list
#    scandata_list = collapsed_gaussian_fit_scandata_list
#    scandata_list = lifetime_fit_scandata_list
    scandata_list = collapsed_lifetime_fit_scandata_list

# %% PLOT OVERLAP DATA
#    field_index = 0  # gaussian amplitude
#    field_index = 1  # gaussian width
#    field_index = 2  # gaussian center pos
#    field_index = 3  # gaussian area

#    field_index = 0  # pulse amplitude 1
#    field_index = 1  # spin lifetime 1
#    field_index = 2  # pulse amplitude 2
    field_index = 3  # spin lifetime 2

#    x_list = []
#    y_list = []
#    yerr_list = []
    plt.figure()
    plt.hold(True)
    for scandata in scandata_list:
        x_vals, y_vals = scandata.dataseries_list[field_index].datalists()
        _, y_errs = scandata.error_dataseries_list[field_index].datalists()
        plt.errorbar(x_vals, y_vals, yerr=y_errs, fmt='-d')
        if scandata.fitdata_list[field_index] is not None:
            x_vals, y_vals = \
                scandata.fitdata_list[field_index].fitdataseries.datalists()
            plt.plot(x_vals, y_vals, ':')

#    plt.errorbar(x_list, y_list, yerr=yerr_list, fmt='.')
