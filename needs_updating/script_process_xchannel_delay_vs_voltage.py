# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import numpy as np

from experimentdataanalysis.analysis.dataseriesprocessing \
    import process_dataseries_and_error_in_scandata, \
        get_positive_time_delay_scandata, \
        get_continuous_phase_scandata, \
        get_gaussian_smoothed_scandata
from experimentdataanalysis.analysis.scandatamodels \
    import IndependentSinusoidalSpinLifetimeModel
from experimentdataanalysis.analysis.scandatasetprocessing \
    import ScanDataSetsAnalyzer, sort_scandata_into_sets
import experimentdataanalysis.parsing.dataseriesparsing as dsparsing


#==============================================================================
# SCRIPT EXAMPLES - SEE BOTTOM
#==============================================================================





# FILTER FUNCTIONS
def get_filter_fcn_no_super_long_lifetimes(threshold=5000):
    def filter_fcn(scandata):
        if scandata.fitdata_list[0] is None:
            return False
        else:
            return scandata.fitdata_list[0].fitparams[3] < threshold
    return filter_fcn


def get_filter_fcn_no_first_n_scans_in_series(num_ignored=1):
    def filter_fcn(scandata):
        if 'FastScanIndex' in scandata.scaninfo_list[0]:
            return scandata.scaninfo_list[0]['FastScanIndex'] > 1
        else:
            return False
    return filter_fcn


# WARNING: ONLY FOR "C:\Data\160702\delay_scans_-8V_to_8V"
def get_filter_fcn_only_one_channel(target_channel):
    def filter_fcn(scandata):
        try:
            if scandata.scaninfo_list[0]['MiddleScanType'] == 'StageY':
                ypos = scandata.scaninfo_list[0]['MiddleScanCoord']
                if target_channel == 1:
                    return ypos < 3.55
                elif target_channel == 2:
                    return ypos >= 3.55 and ypos < 3.95
                elif target_channel == 3:
                    return ypos >= 3.95
                else:
                    raise Exception('filter_to_one_channel: ' +
                                    'invalid target_channel')
            else:
                print("Warning: filter_to_one_channel: scandata does " +
                      "not conform to expected scan parameters, " +
                      "filtered out.")
                return False
        except KeyError:
            return False
    return filter_fcn


# PLOTTING FUNCTION
def plot_scandata(scandata, field_index, model=None,
                  label="", fmt="-bd", fit_fmt="xr:"):
    x_vals, y_vals = scandata.dataseries_list[field_index].data()
    error_dataseries = scandata.error_dataseries_list[field_index]
    if error_dataseries is not None:
        _, y_errs = error_dataseries.data()
        plt.errorbar(x_vals, y_vals, yerr=y_errs, label=label, fmt=fmt)
    else:
        plt.plot(x_vals, y_vals, fmt, label=label)
    if scandata.fitdata_list[field_index] is not None:
        x_vals, y_vals = \
            scandata.fitdata_list[field_index].fitdataseries.data()
        if model is not None:
            x_vals = np.linspace(min(x_vals), max(x_vals), 1000)
            params = scandata.fitdata_list[field_index].fitparams
            y_vals = model.fitfunction(x_vals, *params)
        plt.plot(x_vals, y_vals, fit_fmt)


# %%
if __name__ == "__main__":
# %%  ANALYSIS OF 818.9nm DELAY SCANS

    # ANALYSIS MODELS
    # IMPORTANT RESULTS FOR MODEL 1 CONSIDERATION:
    # one good RSA fit to long species, @818.0nm, no voltage applied:
    # lifetime: 20.26ns +- 0.08ns,
    # freq_per_T ~= 5.41e8 Ts (have to grab again) -> g=.00615? uhh
    # \-> osc_period @300mT: 554.620ps +- .045ps
    # field offset: -0.7565mT +- 0.0056mT
    # note fitting model 2 shows ~1% increase in osc_period moving to +-2V!
    # MODEL 1: Two decaying cosines w/ relative phase. Long lifetime fixed.
    #          Start beyond pump transient (<200ps lifetime) signal
    # params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
    #          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
    #          phase1, phase2, slope, offset
    spin_lifetime_model = \
        IndependentSinusoidalSpinLifetimeModel(
            field_index=0,  # lockin2x
            max_fcn_evals=20000,
            free_params=[False, True, True, True, False,
                         True, True, False, False,
                         True, True, True, True],
            initial_params = [40, 0.015, 0.008, 4000, 20260,
                              555.5, 554.62, 0, 0,
                              0, np.pi, 0, 0],
            param_bounds = [(1, 1000), (0, 0.1), (0.003, 0.1),
                            (1, 1e4), (1e3, 1e7),
                            (500, 600), (500, 600),
                            (0, 0), (0, 0),
                            (-4*np.pi, 4*np.pi), (-4*np.pi, 4*np.pi),
                            (-1e-4, 1e-4), (-0.01, 0.01)],
            error_thresholds=[None, None, None, None, None,
                              None, None, None, None,
                              None, None, None, None],
            dim2_key="Electric Field (V/cm)",  # look at results of fit vs field
            excluded_intervals=[[-15, 400]])
#            excluded_intervals=[[-15, 400], [7000, 15000]])


    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
#    dirpath = ("C:\\Data\\160702\\" +
#               "delayscans_-6V_to_6V_noautocenter")
#               "delayscans_-6V_to_6V")
#    dirpath = ("C:\\Data\\august_data\\160901\\" +
#               "DelayScansNewWavelength")  # requires other assumptions in model! do seperately below...
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "LowV_818.0nm_WavelengthDependence_TRKRvsV_200mT")
#               "WavelengthDependence_TRKR_300mT")
#               "WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    fixed_uncertainty = 1e-4  # manually set uncertainty of data set
    model = spin_lifetime_model
    sort_key = "Voltage"  # one scandataset for each voltage

    scandata_list = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))

    # ---------------
    # TEMPORARY, FOR SPEED:
#    scandata_list = scandata_list[5:6]
    # ---------------

    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)
    field_index = model.field_index  # for future use

    # Change drift velocity based on voltage. Assumes each set has same
    # voltage for every scan!
    for scandataset in scandataset_list:
        voltage_list = [scandata.scaninfo_list[0]['Voltage']
                        for scandata in scandataset.scandata_list]
        set_voltage = voltage_list[0]
        if any(voltage != set_voltage for voltage in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            # Drift velocity per voltage
            drift_velocity_per_volt = 2e-3  # in um/(ps*V)
            drift_velocity = drift_velocity_per_volt*set_voltage  # in um/ps
            scandataset.model.free_params[7] = False  # drift_velocity1
            scandataset.model.free_params[8] = False  # drift_velocity2
            scandataset.model.initial_params[7] = drift_velocity
            scandataset.model.initial_params[8] = drift_velocity

            # Initial fit guess: oscillation period per voltage
#            guess_period_short = 550 + 7.5 * abs(set_voltage)
#            scandataset.model.free_params[5] = True  # let it fit from here
#            scandataset.model.initial_params[5] = guess_period_short

            guess_period_long = 554.62 + 2.7 * abs(set_voltage)
            scandataset.model.free_params[6] = True  # let it fit from here
            scandataset.model.initial_params[6] = guess_period_long

    analyzer = ScanDataSetsAnalyzer(scandataset_list)


    # since now sorted by "Voltage", can change to electric field as well as
    # fix improper "StageZ" 2nd coord, change units to V/cm anyway
    for scandataset in analyzer.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.scaninfo_list[0]['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field
                scaninfo['Electric Field (V/cm)'] = field

#    # smooth over data with a 40ps wide gaussian convolution filter
#    analyzer.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    gaussian_width=40)
    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
    analyzer.apply_transform_to_all_scandata(
                                    get_gaussian_smoothed_scandata,
                                    field_indices_to_process=[field_index],
                                    gaussian_width=600,
                                    edge_handling='reflect',
                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)
    

    # scandatasets don't share models, can't multiprocess in this version:
    analyzer.fit_all_scandata_to_model(multiprocessing=False)


    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetimes(threshold=999999))
#    analyzer.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    lifetime_scandata_list = \
        analyzer.collapse_to_scandata_list(filtered=False)

    # regroup scandatasets into one big set to get single collapsed scandata:
    analyzer.regroup_scandatasets(new_model=model,
                                  sort_key=None)
    collapsed_scandata_list = \
        analyzer.collapse_to_model_fit_scandata_list(new_scan_type="[Unknown]",
                                                     filtered=False)

    # collapsed list: jiggle phase to keep phase continuous via np.unwrap()
    collapsed_scandata_list = [get_continuous_phase_scandata(scandata)
                               for scandata in collapsed_scandata_list]


# %%
    field_index = 0  # lockin2x
    for scandata in lifetime_scandata_list[5:6]:
       plot_scandata(scandata, field_index, model=model, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    if scandata.fitdata_list[field_index] is not None:
        plt.text(8000, 0,
                 "Last fit lifetime: {:.3f}ns\n                         +-{:.3f}ns".format(
                     scandata.fitdata_list[field_index].fitparams[3]/1000,
                     scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.show()

    print('Voltage: {}'.format(scandata.scaninfo_list[0]['Voltage']))
    if scandata.fitdata_list[field_index] is not None:
        print(scandata.fitdata_list[field_index].fitparams)

    scandata_list = lifetime_scandata_list


# %%
    field_indices = [0,1]  # dataseries: x:field, y:short/long phase
    plt.figure()
    plt.hold(True)
    for field_index in field_indices:
        for scandata in collapsed_scandata_list[:]:
            plot_scandata(scandata, field_index,
                          fmt="d", label=scandata.fields[field_index])

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Electric Field (V/cm)")
    plt.ylabel("")  
    plt.legend()
    plt.show()


# %%  ANALYSIS OF 818.0nm DELAY SCANS
    
    # MODEL 2: Two decaying cosines w/ relative phase, as before. This time
    #          no fixing of long lifetime or oscillation factor, and short
    #          lifetime signal is assumed to be very low by starting params
    # params = num_pulses, pulse_amp1, pulse_amp2, lifetime1, lifetime2,
    #          osc_period1, osc_period2, drift_voltage1, drift_voltage2,
    #          phase1, phase2, slope, offset
    spin_lifetime_model_2 = \
        IndependentSinusoidalSpinLifetimeModel(
            field_index=0,  # lockin2x
            max_fcn_evals=2000,
            free_params=[False, False, True, False, True,
                         False, True, False, False,
                         False, True, True, True],
            initial_params = [40, 0, 0.01, 20260, 20260,
                              808, 808, 0, 0,
                              -np.pi, np.pi/2, 0, 0],
            param_bounds = [(1, 1000), (0, 1), (0, 1),
                            (1, 1e6), (1, 1e6),
                            (750, 850), (750, 850),
                            (0, 0), (0, 0),
                            (-2*np.pi, np.pi), (-2*np.pi, np.pi),
                            (-1e-4, 1e-4), (-0.01, 0.01)],
            error_thresholds=[None, None, None, None, None,
                              None, None, None, None,
                              None, None, None, None],
            dim2_key="Electric Field (V/cm)",  # look at results of fit vs field
            excluded_intervals=[[-15, 400]])
#            excluded_intervals=[[-15, 400], [7000, 15000]])


    # LOAD DATA, ORGANIZE, AND FIT IN ANALYZER
#    dirpath = ("C:\\Data\\160702\\" +
#               "delayscans_-6V_to_6V_noautocenter")
#               "delayscans_-6V_to_6V")
#    dirpath = ("C:\\Data\\august_data\\160901\\" +
#               "DelayScansNewWavelength")  # requires other assumptions in model! do seperately below...
    dirpath = ("C:\\Data\\august_data\\160902\\" +
               "LowV_818.0nm_WavelengthDependence_TRKRvsV_200mT")
#               "GoodDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "WavelengthDependence_TRKR_300mT")
#               "WavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")
#               "BestDataFromWavelengthDependence_TRKRvsV_300mT_033XT-A5_818.9nm_30K_2Dscan_Voltage_DelayTime")

    fixed_uncertainty = 1e-4  # manually set uncertainty of data set
    model = spin_lifetime_model_2
    sort_key = "Voltage"  # just one ScanDataSet

    scandata_list = list(dsparsing.fetch_dir_as_unfit_scandata_iterator(
                                     directorypath=dirpath,
                                     key_field="lockin2x",
                                     key_field_error_val=fixed_uncertainty))
    scandataset_list = sort_scandata_into_sets(scandata_list, model, sort_key)

    # Change drift velocity based on voltage. Assumes each set has same
    # voltage for every scan!
    for scandataset in scandataset_list:
        voltage_list = [scandata.scaninfo_list[0]['Voltage']
                        for scandata in scandataset.scandata_list]
        set_voltage = voltage_list[0]
        if any(voltage != set_voltage for voltage in voltage_list):
            print("No common voltage in ScanDataSet, " +
                  "cannot set drift velocity.")
        else:
            # Drift velocity per voltage
            drift_velocity_per_volt = 2e-3  # in um/(ps*V)
            drift_velocity = drift_velocity_per_volt*set_voltage  # in um/ps
            scandataset.model.free_params[7] = False  # drift_velocity1
            scandataset.model.free_params[8] = False  # drift_velocity2
            scandataset.model.initial_params[7] = drift_velocity
            scandataset.model.initial_params[8] = drift_velocity


    analyzer2 = ScanDataSetsAnalyzer(scandataset_list)

    # since now sorted by "Voltage", can change to electric field as well as
    # fix improper "StageZ" 2nd coord, change units to V/cm anyway
    for scandataset in analyzer2.scandataset_list:
        for scandata in scandataset.scandata_list:
            field = scandata.scaninfo_list[0]['Voltage']*20
            for scaninfo in scandata.scaninfo_list:
                scaninfo['MiddleScanType'] = 'Electric Field (V/cm)'
                scaninfo['MiddleScanCoord'] = field
                scaninfo['Electric Field (V/cm)'] = field

    # drift subtraction:
    # subtract from data: data times a 400ps wide gaussian convolution filter                        
#    analyzer2.apply_transform_to_all_scandata(
#                                    get_gaussian_smoothed_scandata,
#                                    field_indices_to_process=[field_index],
#                                    gaussian_width=600,
#                                    edge_handling='reflect',
#                                    subtract_smoothed_data_from_original=True)
    # add 13160ps to all negative delay times
    analyzer2.apply_transform_to_all_scandata(get_positive_time_delay_scandata,
                                             zero_delay_offset=-15)

    # scandatasets don't share models, can't multiprocess in this version:
    analyzer2.fit_all_scandata_to_model(multiprocessing=False)

    # APPLY FILTERS AND EXTRACT FITTED SCANDATA AND SCANDATA OF FITS
#    analyzer2.add_filter_to_each_scandataset(
#        get_filter_fcn_no_super_long_lifetimes(threshold=999999))
#    analyzer2.add_filter_to_each_scandataset(
#        get_filter_fcn_no_first_n_scans_in_series(num_ignored=1))
    lifetime_scandata_list2 = \
        analyzer2.collapse_to_scandata_list(filtered=False)

    # regroup scandatasets into one big set to get single collapsed scandata:
    analyzer2.regroup_scandatasets(new_model=model,
                                   sort_key=None)
    collapsed_scandata_list2 = \
        analyzer2.collapse_to_model_fit_scandata_list(
                                                    new_scan_type="[Unknown]",
                                                    filtered=False)

    # collapsed list: jiggle phase to keep phase continuous via np.unwrap()
    collapsed_scandata_list2 = [get_continuous_phase_scandata(scandata)
                                for scandata in collapsed_scandata_list2]


# %% OVERVIEW OF FITS
    field_index = 0  # lockin2x
    for scandata in lifetime_scandata_list2[:]:
       plot_scandata(scandata, field_index, fmt="bd")
    plt.xlabel("Delay (ps)")
    plt.ylabel("Kerr Rotation (AU)")
    if scandata.fitdata_list[field_index] is not None:
        plt.text(8000, 0,
                 "Last fit lifetime: {}ns\n     +={}ns".format(
                     scandata.fitdata_list[field_index].fitparams[3]/1000,
                     scandata.fitdata_list[field_index].fitparamstds[3]/1000))
    plt.show()

    print('Voltage: {}'.format(scandata.scaninfo_list[0]['Voltage']))
    if scandata.fitdata_list[field_index] is not None:
        print(scandata.fitdata_list[field_index].fitparams)

    scandata_list = lifetime_scandata_list


# %%
    field_indices = [5]  # dataseries: x:field, y:short/long phase
    plt.figure()
    plt.hold(True)
    for field_index in field_indices:
        for scandata in collapsed_scandata_list2[:]:
            plot_scandata(scandata, field_index, fmt="d",
                          label=scandata.fields[field_index])

    # LABEL AND DISPLAY GRAPH
    plt.xlabel("Applied Voltage (V)  |  Electric Field (V/500um)")
    plt.ylabel("")  
    plt.legend()
    plt.show()

    if scandata.fitdata_list[field_index] is not None:
        print(scandata.fitdata_list[field_index].fitparams)

    scandata_list = collapsed_scandata_list

