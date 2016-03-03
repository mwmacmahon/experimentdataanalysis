# -*- coding: utf-8 -*-
"""
Analyze experiments of form:
1D scan: delay scan
2D scan: mirror (or nothing):
3D scan: voltage

might change to:
1D scan: mirror
2D scan: delay scan:
3D scan: voltage
"""

import matplotlib.pyplot as plt

import experimentdataanalysis.analysis.dataclassfitting as dcfitting
import experimentdataanalysis.analysis.dataclassgraphing as dcgraphing
import experimentdataanalysis.parsing.dataclassparsing as dcparsing

# %%
if __name__ == "__main__":
# %%
    rawscandata_iter = dcparsing.fetch_dir_as_unfit_scandata_iterator(
        directorypath="C:\\Data\\febdata\\Experiment - Channel 2")
#        directorypath="C:\Data\decdata\Channel 3 Run 1")
#        directorypath="C:\\Data\\decdata\\representative")
    scandata_list = list(dcfitting.fit_scandata_iterable(
        rawscandata_iter,
        timeseriesfitfunction=dcfitting.fit_timeseries_with_one_decaying_cos,
        fit_drift=True, multiprocessing=True))

# %%
    dcgraphing.graph_scandata_iterable(scandata_list)

# %%
    fig, ax = plt.subplots()
    voltages = []
    voltages2 = []
    lifetimes = []
    lifetimes2 = []
    lifetimestds = []
    lifetimestds2 = []
    for scandata in scandata_list:
        if scandata.fitdata is not None:
            field = scandata.scaninfo['MiddleScanCoord']
            if scandata.fitdata.fitparamstds[4] < 1000:
                voltage = scandata.scaninfo['Voltage']
                lifetime = scandata.fitdata.fitparams[4]
                lifetimestd = scandata.fitdata.fitparamstds[4]
                if field == 0.08:
                    voltages.append(voltage)
                    lifetimes.append(lifetime)
                    lifetimestds.append(lifetimestd)
                else:
                    voltages2.append(voltage)
                    lifetimes2.append(lifetime)
                    lifetimestds2.append(lifetimestd)
    ax.errorbar(voltages, lifetimes, yerr=lifetimestds, fmt='ro',
                label="B=80mT")
    ax.errorbar(voltages2, lifetimes2, yerr=lifetimestds2, fmt='bo',
                label="B=160mT")
    ax.legend(loc='best', fancybox=True, framealpha=0.5)
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Net Spin Lifetime (ps)")
    ax.set_title("Spin Lifetime vs Voltage - Active Channel")

# %%
    fig, ax = plt.subplots()
    voltages = []
    voltages2 = []
    gfactors = []
    gfactors2 = []
    gfactorstds = []
    gfactorstds2 = []
    scalefactor = 11.3713  # =hbar/bohr magneton in ps*T
    for scandata in scandata_list:
        if scandata.fitdata is not None:
            field = scandata.scaninfo['MiddleScanCoord']
            voltage = scandata.scaninfo['Voltage']
            freq = abs(scandata.fitdata.fitparams[5])
            gfactor = scalefactor*freq/field
            frequencystd = scandata.fitdata.fitparamstds[5]
            gfactorstd = scalefactor*frequencystd/field
            if (abs(gfactor - 0.44) < 0.1 and gfactorstd < 0.1):
                if field == 0.08:
                    voltages.append(voltage)
                    gfactors.append(gfactor)
                    gfactorstds.append(gfactorstd)
                else:
                    voltages2.append(voltage)
                    gfactors2.append(gfactor)
                    gfactorstds2.append(gfactorstd)
    ax.errorbar(voltages, gfactors, yerr=gfactorstds, fmt='ro',
                label="B=80mT")
    ax.errorbar(voltages2, gfactors2, yerr=gfactorstds2, fmt='bo',
                label="B=160mT")
    ax.legend(loc='best', fancybox=True, framealpha=0.5)
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Measured g-factor")
    ax.set_title("Spin g-factor vs Voltage - Active Channel")

# %%
#    rawscandata = dcparsing.fetch_csv_as_unfit_scandata(
#        filepath="C:\\Data\\decdata\\representative\\ch3run1_7V.dat")
#    scandata = dcfitting.fit_scandata(
#        rawscandata,
#        timeseriesfitfunction=dcfitting.fit_timeseries_with_one_decaying_cos)
#    dcgraphing.plot_scandata(scandata)
