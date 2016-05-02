# -*- coding: utf-8 -*-
"""
Analyze experiments of form:
1D scan: delay scan
(lots of runs)
"""

import matplotlib.pyplot as plt

import experimentdataanalysis.analysis.dataclassfitting as dcfitting
import experimentdataanalysis.analysis.dataclassgraphing as dcgraphing
import experimentdataanalysis.parsing.dataclassparsing as dcparsing

# %%
if __name__ == "__main__":
# %%
    scandata_list = list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
        directorypath="C:\\Data\\aprilfoolsxinlin"))
#    scandata_list += list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelCenter_200mT_Channel3_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime"))
    scandata_list = list(dcfitting.fit_scandata_iterable(
        scandata_list,
#        dataseriesfitfunction=None,
        dataseriesfitfunction=dcfitting.fit_dataseries_to_two_decays,
        fit_drift=False, multiprocessing=False))

# %%
#    _ = dcgraphing.graph_scandata_iterable(
#            scandata_list, plot_options={'suppress_legend': True})

# %%
    primarykey = 'MiddleScanCoord'
    secondarykey = 'Voltage'

    def scandatasortfcn(scandata):
        try:
            return (scandata.scaninfo[primarykey],
                    scandata.scaninfo[secondarykey])
        except KeyError:
            try:
                return (scandata.scaninfo[primarykey], 0)
            except KeyError:
                try:
                    return (0, scandata.scaninfo[secondarykey])
                except KeyError:
                    print("No valid sort keys found in scandata tags!")
                    return (100000, 100000)
        except AttributeError:
            print("Scandata expected as list element.")
            return (100000, 1000000)

    scandata_list = sorted(scandata_list, key=scandatasortfcn)

# %%
    scandata_list = [scandata for scandata in scandata_list
                     if scandata.scaninfo['MiddleScanCoord'] == 0.16]

# %%
    plotfield = "lockin2x"
    ytype = 'Voltage'
    data2d = []
    data2d_scandata = []
    data2d_scandata_indices = []
    for scandata in list(scandata_list):
        for ind, field in enumerate(scandata.fields):
            if field == plotfield:
                data2d.append(scandata.dataseries[ind].yvals(unfiltered=True))
                data2d_scandata.append(scandata)
                data2d_scandata_indices.append(ind)

    # scale by laserpower?
#    for ind, datarow in enumerate(data2d):
#        rowlaserpower = None
#        for ind2, field in enumerate(data2d_scandata[ind]):
#            if field == "laserpower":
#                rowlaserpower = data2d_scandata[ind].dataseries[ind2]
#        if rowlaserpower is not None:
#            for val, laserpower in zip(datarow, rowlaserpower):
#                val = val / laserpower

    ind = data2d_scandata_indices[0]
    xmin = min(data2d_scandata[0].dataseries[ind].xvals(unfiltered=True))
    xmax = max(data2d_scandata[0].dataseries[ind].xvals(unfiltered=True))
    ymin = min(scandata.scaninfo[ytype]
               for scandata in data2d_scandata)
    ymax = max(scandata.scaninfo[ytype]
               for scandata in data2d_scandata)

    yvals = []
    omitted_indices = []
    for ind, scandata in enumerate(data2d_scandata):
        if scandata.scaninfo[ytype] not in yvals:
            yvals.append(scandata.scaninfo[ytype])
        else:
            omitted_indices.append(ind)
    data2d = [element for ind, element in enumerate(data2d)
              if ind not in omitted_indices]
    data2d_scandata = [element for ind, element in enumerate(data2d_scandata)
                       if ind not in omitted_indices]
    data2d_scandata_indices = \
        [element for ind, element in enumerate(data2d_scandata_indices)
         if ind not in omitted_indices]

    # TODO: MAKE NOT HARD-SET FOR APS MEETING VALUES:
    ymin = ymin*20  # temp:E-Field
    ymax = ymax*20  # temp:E-Field
    imageplot = plt.imshow(data2d,
                           interpolation="nearest",
                           extent=[xmin, xmax, ymax, ymin])
    plt.xlabel('Pump-Probe Delay (ps)')
    plt.ylabel('Electric Field (V/cm)')  # instead of ytype
    imageplot.axes.set_aspect('auto')
    plt.colorbar(imageplot, label="Kerr Rotation (AU)")

# %%
    voltages = []
    voltages2 = []
    voltages3 = []
    voltages4 = []
    lifetimes = []
    lifetimes2 = []
    lifetimes3 = []
    lifetimes4 = []
    lifetimestds = []
    lifetimestds2 = []
    lifetimestds3 = []
    lifetimestds4 = []
    for scandata in scandata_list:
        if scandata.fitdata[0] is not None:
            if scandata.scaninfo['MiddleScanType'] == 'BExternal':
                field = scandata.scaninfo['MiddleScanCoord']
            else:
                field = 0.2
            if scandata.fitdata[0].fitparamstds[4] < 1000:
                voltage = scandata.scaninfo['Voltage']
                lifetime = scandata.fitdata[0].fitparams[4]
                lifetimestd = scandata.fitdata[0].fitparamstds[4]
                if field == 0:
                    voltages4.append(voltage)
                    lifetimes4.append(lifetime)
                    lifetimestds4.append(lifetimestd)
                elif field == 0.08:
                    voltages.append(voltage)
                    lifetimes.append(lifetime)
                    lifetimestds.append(lifetimestd)
                elif field == 0.16:
                    voltages2.append(voltage)
                    lifetimes2.append(lifetime)
                    lifetimestds2.append(lifetimestd)
                elif field == 0.2:
                    voltages3.append(voltage)
                    lifetimes3.append(lifetime)
                    lifetimestds3.append(lifetimestd)
    if voltages or voltages2 or voltages3 or voltages4:
        fig, ax = plt.subplots()
    if voltages:
        ax.errorbar(voltages, lifetimes, yerr=lifetimestds, fmt='ro',
                    label="B=80mT")
    if voltages2:
        ax.errorbar(voltages2, lifetimes2, yerr=lifetimestds2, fmt='bo',
                    label="B=160mT")
    if voltages3:
        ax.errorbar(voltages3, lifetimes3, yerr=lifetimestds3, fmt='bo',
                    label="B=200mT")
    if voltages4:
        ax.errorbar(voltages4, lifetimes4, yerr=lifetimestds4, fmt='bo',
                    label="B=?mT")
    if voltages or voltages2 or voltages3 or voltages4:
        ax.legend(loc='best', fancybox=True, framealpha=0.5)
        ax.set_xlabel("Electric Field (V/cm)")
        ax.set_ylabel("Net Spin Lifetime (ps)")
        ax.set_title("Spin Lifetime vs Voltage - Adjacent Channel")

# %%
#    fig, ax = plt.subplots()
#    voltages = []
#    voltages2 = []
#    gfactors = []
#    gfactors2 = []
#    gfactorstds = []
#    gfactorstds2 = []
#    scalefactor = 11.3713  # =hbar/bohr magneton in ps*T
#    for scandata in scandata_list:
#        if scandata.fitdata is not None:
#            field = scandata.scaninfo['MiddleScanCoord']
#            voltage = scandata.scaninfo['Voltage']
#            freq = abs(scandata.fitdata[0].fitparams[5])
#            gfactor = scalefactor*freq/field
#            frequencystd = scandata.fitdata[0].fitparamstds[5]
#            gfactorstd = scalefactor*frequencystd/field
#            if (abs(gfactor - 0.44) < 0.1 and gfactorstd < 0.1):
#                if field == 0.08:
#                    voltages.append(voltage)
#                    gfactors.append(gfactor)
#                    gfactorstds.append(gfactorstd)
#                else:
#                    voltages2.append(voltage)
#                    gfactors2.append(gfactor)
#                    gfactorstds2.append(gfactorstd)
#    ax.errorbar(voltages, gfactors, yerr=gfactorstds, fmt='ro',
#                label="B=80mT")
#    ax.errorbar(voltages2, gfactors2, yerr=gfactorstds2, fmt='bo',
#                label="B=160mT")
#    ax.legend(loc='best', fancybox=True, framealpha=0.5)
#    ax.set_xlabel("Voltage (V)")
#    ax.set_ylabel("Measured g-factor")
#    ax.set_title("Spin g-factor vs Voltage - Active Channel")

# %%
#    rawscandata = dcparsing.fetch_csv_as_unfit_scandata(
#        filepath="C:\\Data\\decdata\\representative\\ch3run1_7V.dat")
#    scandata = dcfitting.fit_scandata(
#        rawscandata,
#        dataseriesfitfunction=dcfitting.fit_dataseries_with_one_decaying_cos)
#    dcgraphing.plot_scandata(scandata)
