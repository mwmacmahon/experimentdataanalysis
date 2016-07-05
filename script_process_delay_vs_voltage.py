# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt

from experimentdataanalysis.analysis.dataseriesprocessing \
    import get_positive_time_delay_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import GaussianModel, SpinLifetimeModel, SinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer


# %% ALTERNATE SCRIPT
if __name__ == "__main__":

# %% LOAD AND FIT DATA
    spin_lifetime_model = \
        SpinLifetimeModel(max_fcn_evals=80000,
                          free_params=[True, True, True, True, True],
                          initial_params = [0.05, 50, 0.05, 2000, 0],
                          param_bounds = [(0, 1), (1, 200),
                                          (0, 1), (10, 1e6),
                                          (-0.1, 0.1)],
                          error_thresholds=[0.2, 200, 0.2, 2000, 0.1],
                          xval_key="Voltage",
                          field_index=0)
#    spin_lifetime_model = \
#        SinusoidalSpinLifetimeModel(
#            max_fcn_evals=80000,
#            error_thresholds=[0.2, 200, 0.2, 2000, 400, None, None],
#            xval_key="Voltage")
    analyzer = ScanDataSetsAnalyzer(spin_lifetime_model,
#                                    dirpath="C:\\Data\\160702\\delay_scans_vs_overlap-delay_scans\\1V_delay",
                                    dirpath="C:\\Data\\160702\\delay_scans_vs_overlap-delay_scans\\4V_delay",
                                    uncertainty_value=1e-4)
    analyzer.break_up_repeating_scandatasets()
    analyzer.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
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
    trkr_fits_scandata_list = \
        analyzer.collapse_to_scandata_list(filtered=False)
    collapsed_trkr_fit_scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(filtered=False)

# %% scandata_list: delay scans w/ lifetime fits
#                   fields = measurement reading types
    scandata_list = trkr_fits_scandata_list
    field_index = 0  # usual lockin2x or w/e (w/ gaussian fit)


# %% scandata_list: lifetime fit data
#                   fields = calculated fit values
    scandata_list = collapsed_trkr_fit_scandata_list
#    field_index = 0  # short decay - pulse amplitude
#    field_index = 1  # short decay - spin lifetime
#    field_index = 2  # long decay - pulse amplitude
    field_index = 3  # long decay - spin lifetime


# %% PLOT OVERLAP DATA
    plt.figure()
    plt.hold(True)
    for scandata in scandata_list:
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

