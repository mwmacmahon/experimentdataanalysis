# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import GaussianModel, SpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer


# %% ALTERNATE SCRIPT
if __name__ == "__main__":

# %% LOAD AND FIT DATA
    gaussian_model = \
        GaussianModel(max_fcn_evals=80000,
                      error_thresholds=[0.1, 30, 30, None, None])
    analyzer = ScanDataSetsAnalyzer(gaussian_model,
#                                    dirpath="C:\\Data\\early_may_good_data",
#                                    dirpath="C:\\Data\\160702\\delay_scans_vs_overlap-delay_scans\\1V_overlap",
                                    dirpath="C:\\Data\\160702\\delayscans_-6V_to_6V_reference_overlapscans",
                                    uncertainty_value=1e-4)
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
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="Voltage",
                                                     filtered=True)

# %% LOAD AND FIT GAUSSIAN_FIT DATA TO LIFETIMES
    lifetime_model = \
        SpinLifetimeModel(max_fcn_evals=80000,
                          free_params=[True, True, True, True, True],
                          error_thresholds=[0.2, 200, 0.2, 2000, None],
                          dim2_key="Voltage",
                          field_index=0)  # use amplitude instead of area
    analyzer2 = ScanDataSetsAnalyzer(
        lifetime_model,
        scandata_list=collapsed_gaussian_fit_scandata_list,
        set_key=None)  # single unified set, no 4th coordinate in particular
    analyzer2.apply_transform_to_all_scandata(get_positive_time_delay_scandata)
    analyzer2.fit_all_scandata_to_model(multiprocessing=True)

# %% FILTER AGAIN?
    # Start by clearing existing filters, only necessary on subsequent runs
    analyzer2.clear_filters_from_each_scandataset()

    def filter_out_under_5_data_points(scandata):
        return all(len(dataseries) >= 5
                   for dataseries in scandata.dataseries_list)

    # Add all desired filters, not including built-in model fit thresholds
#    analyzer2.add_filter_to_each_scandataset(filter_out_under_5_data_points)

# %% EXTRACT FILTERED FIT DATA
    lifetime_fit_scandata_list = \
        analyzer2.collapse_to_scandata_list(filtered=True)
    collapsed_lifetime_fit_scandata_list = \
        analyzer2.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=True)

# %% scandata_list: gaussian scans w/ gaussian fits
#                   fields = measurement reading types
    scandata_list = gaussian_fits_scandata_list
    field_index = 0  # usual lockin2x or w/e (w/ gaussian fit)


# %% scandata_list: gaussian fit results w/ lifetime fits
#                   fields = calculated fit values
    scandata_list = lifetime_fit_scandata_list
    field_index = 0  # gaussian amplitude
#    field_index = 1  # gaussian width
#    field_index = 2  # gaussian center pos
#    field_index = 3  # gaussian area (w/ lifetime fit)


# %% scandata_list: lifetime fit results
#                   fields = calculated fit values
    scandata_list = collapsed_lifetime_fit_scandata_list
#    field_index = 0  # short decay - pulse amplitude
#    field_index = 1  # short decay - spin lifetime
#    field_index = 2  # long decay - pulse amplitude
    field_index = 3  # long decay - spin lifetime

#    scandata_list = lifetime_fit_scandata_list


# %% PLOT OVERLAP DATA
    plt.figure()
    plt.hold(True)
    for scandata in scandata_list[0:10]:
        x_vals, y_vals = scandata.dataseries_list[field_index].datalists()
        error_dataseries = scandata.error_dataseries_list[field_index]
        if error_dataseries is not None:
            _, y_errs = error_dataseries.datalists()
            plt.errorbar(x_vals, y_vals, yerr=y_errs, fmt='-bd')
        else:
            plt.plot(x_vals, y_vals, '-bd')
        if scandata.fitdata_list[field_index] is not None:
            x_vals, y_vals = \
                scandata.fitdata_list[field_index].fitdataseries.datalists()
            plt.plot(x_vals, y_vals, 'xr:')

#    plt.errorbar(x_list, y_list, yerr=yerr_list, fmt='.')
