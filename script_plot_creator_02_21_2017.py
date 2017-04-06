
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.filters as filters
import scipy.stats as stats

from experimentdataanalysis.analysis.dataclasses import ScanData
from experimentdataanalysis.analysis.scandataprocessing \
    import make_scandata_time_delay_positive, \
           make_scandata_phase_continuous, \
           gaussian_smooth_scandata, \
           process_scandata_fields, \
           scandata_model_fit, \
           scandata_iterable_sort
from experimentdataanalysis.analysis.featurevectors import \
    scandata_list_to_fvec_scandata, \
    split_fvec_scandata_by_training_and_test
from experimentdataanalysis.analysis.scandatasetprocessing \
    import sort_scandata_into_sets, fit_scandataset_list_to_model
from experimentdataanalysis.parsing.scandataparsing import \
        fetch_dir_as_unfit_scandata_iterator

import fit_models_two_species_pump_probe as two_species_model


#GFACTORCONSTANT = 0.013996  # 1/(ps*Tesla), = bohr magneton/2*pi*hbar
GFACTORCONSTANT = 1.3996e-5  # 1/(ps*mTesla), = bohr magneton/2*pi*hbar
LASER_REPRATE = 13160  # ps period

# FILENAME-TO-INFO-DICT PARSING RULES
# 1. If first string found, register second string as
#    tag containing third string/value
#        e.g. if keyword_list contains ("warmup", "Warmup?", True):
#             "...warmup..." -> {"Warmup?": True}
this_element_keyword_list = [("TRKR", "IsTRKR?", True),
                             ("RSA", "IsRSA?", True)]

# 2. If string in element[0] is found in filepath separated from other
#    characters by '_'s, will record adjacent number and store in
#    scandata's info dict under the key element[1], as a float
#    if possible.
#    e.g. TRKR_15Vcm_230mT.dat -> {"Electric Field (V/cm)": 15.0,
#                                  "Magnetic Field (mT)": 230.0}
filepath_element_keyword_list = [("Vcm", "Electric Field (V/cm)"),
                                 ("mT", "Magnetic Field (mT)"),
                                 ("K", "Set Temperature (K)"),
                                 ("nm", "Wavelength (nm)"),
                                 ("ps", "Delay Time (ps)"),
                                 ("run", "RunIndex"),
                                 ("sAt", "Time At Field + Junk"),
#                                 ("V", "Voltage (V)"),
                                 ("x", "SecondScanCoord")]

# for this one, if element [0] is found,
# next element stored w/ key given by elements [1][0], [1][1], [1][2], etc.
filepath_next_element_keyword_list = [("Ind", "FirstScanIndex"),
                                      ("2Dscan", ["SecondScanType",
                                                  "FirstScanType"]),
                                      ("Voltage", "Voltage (V)"),
                                      ("MirrorZ", "Pump-Probe Distance (um)")]
FILEPATH_PARSING_KEYWORD_LISTS = [this_element_keyword_list,
                                  filepath_next_element_keyword_list,
                                  filepath_element_keyword_list]




yfield = 'lockin1x'
fixed_uncertainty = None


# DIFFERENT PLOTS:
# 1. TRKR vs field, use:
#    - TRKR dir, 
#    - 'SecondScanCoord' to sort,
#    - 0::3 OR 1::3 OR 2::3 indices from scandata_list (for run #1/2/3)
# 2. RSA scans forward/backward/both @ one delay, use:
#    - RSA dir,
#    - 'Delay Time (ps)' to sort
#    - 0:3 OR 3:6 OR 6:9 for -500, -280, -60 respectively
# 3. RSA forward OR backward OR both @ all delays
#    - RSA dir,
#    - 'Delay Time (ps)' to sort
#    - 0::3 OR 1::3 OR 2::3 for forwards, backwards, or both respectively


dirpath = ("C:\\Data\\feb2017_data_part2\\" +
               "170217\\" +
                   "RSA_FieldScanTesting_From-300mT_Increasing_Narrow_-160ps___818.60nm_10K_2Dscan_Voltage_BExternal")
#                   "TRKR_FieldScanTesting_From-300mT_SlowSweepUpFrom242mT___818.58nm_10K_2Dscan_BExternal_DelayTime")
#                   "FieldScanTesting_EarlyScans")
#               "170219\\phaselocktesting\\" +
#                   "RSA_PhaseTracking")
#                   "TRKR_PhaseTracking_from110mT___816.48nm_10K_2Dscan_BExternal_DelayTime")

scandata_list = \
    list(fetch_dir_as_unfit_scandata_iterator(
                directorypath=dirpath,
                yfield=yfield,
                yfield_error_val=fixed_uncertainty,
                parsing_keywordlists=FILEPATH_PARSING_KEYWORD_LISTS))
print("Extracted {} ScanData from dir path".format(len(scandata_list)))

scandata_list, temp = scandata_iterable_sort(scandata_list,
#                                             'Delay Time (ps)',
                                            'SecondScanCoord',
                                             'FirstScanIndex',
#                                             'Filepath',
#                                             'Filepath',
                                             numeric_sort=True)
#scandata_list = scandata_list[0:3]
#scandata_list = scandata_list[3:6]
#scandata_list = scandata_list[6:9]
#scandata_list = scandata_list[0::3]
#scandata_list = scandata_list[7:8]

# %% PLOT DATA

#nrows = 3  # SET NUM ROWS, GET # COLS
#ncols = int(np.ceil(len(scandata_list) / nrows))
ncols = 1  # SET NUM COLS, GET # ROWS
nrows = int(np.ceil(len(scandata_list) / ncols))

xlims = [min(scandata_list[0].x), max(scandata_list[0].x)]
ylims = [min(min(scandata.y) for scandata in scandata_list),
         max(max(scandata.y) for scandata in scandata_list)]
plt.figure()
for ind, scandata in enumerate(scandata_list):
    col = int(np.ceil((ind + 1) / nrows))  # (row, col) starting at (1,1)
    row = ind + 1 - nrows * (col - 1)
    subplot_index = col + ncols * (row - 1)
    plt.subplot(nrows, ncols, subplot_index)

    xvals, yvals = scandata.xy
    plt.plot(xvals, yvals, 'bd')
#    for xy_ind, (xval, yval) in enumerate(zip(xvals, yvals)):
#        scale = 1.0 * xy_ind / len(xvals)
#        color = (1.0 * scale, 0, 1 - 1.0 * scale)
#        plt.plot(xval, yval, 'o',
#                 markerfacecolor=color, markeredgecolor=color)


    plt.plot(xlims, [0, 0], 'k:')
    plt.xlim(xlims)
    plt.ylim(ylims)

    if col == ncols: # on last column, last row is unique to that column
        col_nrows = len(scandata_list) - nrows * (ncols - 1)
    else:
        col_nrows = nrows

    if row < col_nrows:
        plt.xticks([])
    else:
        FCTR = 1
#        plt.xticks([min(scandata.x), max(scandata.x)])
#        plt.xticks([round(min(scandata.x)*FCTR)/FCTR,
#                    round(max(scandata.x*FCTR))/FCTR])

    if col > 1:
        plt.yticks([])
    else:
        plt.yticks([0])


    # TITLING
    if row == 1 and col == 1:
        plt.title("Lessening buildup time, but slowing sweep")

#    plt.title("B = {}mT".format(scandata.info['SecondScanCoord']))

#    plt.title("B = {}mT (jumped from 0)".format(scandata.info['Magnetic Field (mT)']))

#    plt.title("Scanning {}mT to {}mT".format(  # label every plot by min/max x
#                scandata.x[0], scandata.x[-1]))

#    if row == 1: # label RSA individually @ different delays
#        plt.title('Scanning 110mT to 130mT\n@{}ps'.format(
#                    scandata.info['Delay Time (ps)']))
#    elif row == 2:
#        plt.title('@{}ps'.format(scandata.info['Delay Time (ps)']))
#    elif row == 3:
#        plt.title('@{}ps'.format(scandata.info['Delay Time (ps)']))

#    if row == 1: # label RSA individually, different scan types
#        plt.title('@{}ps\nScanning 110mT to 130mT'.format(
#                    scandata.info['Delay Time (ps)']))
#    elif row == 2:
#        plt.title('Scanning 130mT to 110mT')
#    elif row == 3:
#        plt.title('Scanning 110mT to 130mT to 110mT')

#    if row == 1:  # label 2D scan vs field
#        row_end_scandata = \
#            scandata_list[min(len(scandata_list) - 1, ind + nrows - 1)]
#        plt.title("B: {}mT to {}mT".format(
#            scandata.info['SecondScanCoord'],
#            row_end_scandata.info['SecondScanCoord']))


# %% OPTIONAL: OVERLAY RSA DATA ONTO TRKR DATA
assert(rsa_scandata_list is not None)  # create elsewhere by running above code

for ind, scandata in enumerate(scandata_list):  # go back through TRKRs
    col = int(np.ceil((ind + 1) / nrows))  # (row, col) starting at (1,1)
    row = ind + 1 - nrows * (col - 1)
    subplot_index = col + ncols * (row - 1)
    plt.subplot(nrows, ncols, subplot_index)

    trkr_bfield = scandata.info['SecondScanCoord']
    for rsa_scandata in rsa_scandata_list:
        rsa_delay_time = rsa_scandata.info['Delay Time (ps)']
        field_match_indices = np.argwhere(rsa_scandata.x == trkr_bfield)
        yvals = rsa_scandata.y[field_match_indices]
        if len(yvals) > 0:
            xvals = rsa_delay_time * np.ones_like(yvals)
            plt.plot(xvals, yvals, 'rd')
    plt.xlim([-510, 210])
    plt.ylim([-0.07, 0.14])

# %% ALTERNATIVE DATA PLOTTING

plt.figure()
for ind, scandata in enumerate(scandata_list):
    plt.subplot(len(scandata_list), 1, ind + 1)
    plt.plot(*scandata.xy, 'bd')
