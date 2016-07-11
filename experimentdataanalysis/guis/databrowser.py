# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:27:08 2016

@author: Michael
"""

import os, pkg_resources

from matplotlib.backends.backend_qt4agg \
    import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from PyQt4 import QtCore, QtGui, uic

from experimentdataanalysis.analysis.dataclasses \
    import FitData, FitFunc, ScanData, DataSeries
import experimentdataanalysis.analysis.dataclassgraphing as dcgraphing
import experimentdataanalysis.parsing.dataclassparsing as dcparsing
from experimentdataanalysis.guis.guistarter import QApplicationStarter


# %%
class DataBrowserWindow(QtGui.QMainWindow):
    """
    Qt application that lets user browse and fit saved experiment data.
    """
    def __init__(self, appoutput, app_saved_state=None):
        """
        Initializes the DataBrowserWindow. Takes an AppOutput object
        or other QObject with the function "receive_output([1 pos. param])"
        used to export any "return values", and optionally takes a
        app_saved_state object emitted from a previous instance.
        """
        super(DataBrowserWindow, self).__init__()
        resource_package = __name__
        gui_filepath = pkg_resources.resource_filename(resource_package,
                                                       "ui_databrowser.ui")
        uic.loadUi(gui_filepath, self)
        self.ignore_input = True
        self.canvas = MPLCanvas(dpi=75)
        self.figure_container.addWidget(self.canvas)
        # Set up combo boxes, load fit functions
        for cmbbox in [self.cmb_primarysort, self.cmb_secondarysort]:
            cmbbox.addItem('Channel')
            cmbbox.addItem('FastScanIndex')
            cmbbox.addItem('MiddleScanCoord')
            cmbbox.addItem('SetTemperature')
            cmbbox.addItem('Voltage')
            cmbbox.addItem('Wavelength')
        self.cmb_primarysort.setCurrentIndex(0)
        self.cmb_secondarysort.setCurrentIndex(1)
        # Set up application signals
        self.btn_loadfile.pressed.connect(
                                        self.callback_add_file)
        self.btn_loaddir.pressed.connect(
                                        self.callback_load_dir)
        self.btn_excludescandata.pressed.connect(
                                        self.callback_exclude_scandata)
        self.btn_clearallscandata.pressed.connect(
                                        self.callback_clear_all)
        self.list_scandata.itemSelectionChanged.connect(
                                        self.callback_update_selection)
        self.btn_sortdata.pressed.connect(
                                        self.callback_sort_scan_list)
        self.btn_setprimary.pressed.connect(
                                        self.callback_set_primary_datatype)
        self.btn_plot1d.toggled.connect(
                                        self.callback_update_plots)
        self.btn_plot2d.toggled.connect(
                                        self.callback_update_plots)
        self.btn_plotfitparam.toggled.connect(
                                        self.callback_update_plots)
        self.cmb_datatype.activated.connect(
                                        self.callback_update_datatype)
        self.btn_fitstart.pressed.connect(
                                        self.callback_perform_fit)
        self.btn_fitstore.pressed.connect(
                                        self.callback_store_fit)
        self.btn_fiterasetemp.pressed.connect(
                                        self.callback_erase_temp_fit)
        self.btn_fitdelete.pressed.connect(
                                        self.callback_delete_fit)
# FROM SERIESEDITOR - FOR REFERENCE ONLY
#        self.btn_seriesedit.pressed.connect(self.edit_scandata)
#        self.btn_seriesdelete.pressed.connect(self.delete_scandata)
#        self.btn_seriesadd.pressed.connect(self.add_scandata)
#        self.radio_linear.toggled.connect(self.update_state)
#        self.radio_log.toggled.connect(self.update_state)
#        self.edit_min.editingFinished.connect(self.update_state)
#        self.edit_max.editingFinished.connect(self.update_state)
#        self.edit_numvalues.editingFinished.connect(self.update_state)
#        self.check_combinescandata.toggled.connect(self.update_state)
#        self.check_breaksubsections.toggled.connect(self.update_state)
#        self.slider_numsubsections.valueChanged.connect(self.update_state)
#        self.check_shuffle.toggled.connect(self.update_state)
#        self.btn_finish.pressed.connect(self.close)
        self.connect(self,  # emitter
                     QtCore.SIGNAL('output_signal(PyQt_PyObject)'),  # signal
                     appoutput.receive_output)  # receiver
        self.current_scan_list = []
        self.temporary_fitdata_list = []
        self.last_scaninfo = None
        if app_saved_state is not None:
            for scandata in app_saved_state['current_scan_list']:
                self.add_scandata_to_list(scandata)
            if len(self.current_scan_list) > 0:
                self.list_scandata.setCurrentRow(0)
        self.current_scan_list_changed_since_last_update = True
        self.update_state()
        self.ignore_input = False
        self.statusBar.showMessage("Ready")
        self.show()
        # plot once more if data pre-loaded since plot can be wonky un-shown:
        if len(self.current_scan_list) > 0:
            self.statusBar.showMessage("Plotting...")
            self.current_scan_list_changed_since_last_update = True
            self.update_state()
            self.statusBar.showMessage("Ready")

    def closeEvent(self, event):
        """
        Close window and send the output receiver a copy of the output
        data and saved state.
        """
        app_saved_state = {'current_scan_list': self.current_scan_list}
        output_data = tuple(self.current_scan_list)
        self.emit(QtCore.SIGNAL('output_signal(PyQt_PyObject)'),
                  tuple([app_saved_state, output_data]))
        super(DataBrowserWindow, self).closeEvent(event)

    def callback_load_dir(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Loading...")
        try:
            scandata_list = \
                list(dcparsing.fetch_dir_as_unfit_scandata_iterator())
        # TODO: catch exceptions!
        except:
            self.ignore_input = False
            raise
        self.clear_all_scandata(suppress_update=True)
        for scandata in scandata_list:
            self.add_scandata_to_list(scandata)
        self.list_scandata.setCurrentRow(0)
        self.update_state()
        self.callback_sort_scan_list()
        self.update_state()
        self.statusBar.showMessage("Ready")
        self.ignore_input = False

    def callback_add_file(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Loading...")
        try:
            scandata = dcparsing.fetch_csv_as_unfit_scandata()
        # TODO: catch exceptions!
        except:
            self.ignore_input = False
            raise
        self.add_scandata_to_list(scandata)
        if len(self.current_scan_list) == 1:
            self.list_scandata.setCurrentRow(0)
        self.update_state()
        self.statusBar.showMessage("Ready")
        self.ignore_input = False

    def callback_exclude_scandata(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                self.list_scandata.takeItem(rownum)
                self.current_scan_list.pop(rownum)
                self.temporary_fitdata_list.pop(rownum)
                self.current_scan_list_changed_since_last_update = True
                self.update_state()
        self.ignore_input = False

    def callback_clear_all(self, suppress_update=False):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.clear_all_scandata(suppress_update)
        self.ignore_input = False

    def callback_update_selection(self):
        if self.ignore_input:
            return
        # Force plot update if plotting 1D or 2D
        if self.btn_plot1d.isChecked() or self.btn_plot2d.isChecked():
            self.current_scan_list_changed_since_last_update = True
        self.update_state(new_current_selection=True)

    def callback_update_plots(self, toggled_on):
        # ignore "toggled off" signals in addition to input during processing
        if self.ignore_input or not toggled_on:
            return
        # Force plot update
        self.current_scan_list_changed_since_last_update = True
        self.update_state()

    def callback_update_datatype(self):
        if self.ignore_input:
            return
        # Force plot update
        self.current_scan_list_changed_since_last_update = True
        self.update_state()

    def callback_set_primary_datatype(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        if self.get_active_scandata() is not None:
            self.statusBar.showMessage("Processing...")
            newfield = self.cmb_datatype.currentText()
            old_scandata_list = self.current_scan_list[:]  # DEEP COPY
            new_scandata_list = []
            for scandata in old_scandata_list:
                if newfield in scandata.fields:
                    new_scandata_list.append(scandata)
            if old_scandata_list != new_scandata_list:
                yndialog = QtGui.QMessageBox()
                yndialog.setText("Warning - the following scans lack this " +
                                 "field and will be removed:")
                yndialog.setInformativeText(
                    "\n".join(self.list_scandata.item(ind).text()
                              for ind, scandata in enumerate(old_scandata_list)
                              if scandata not in new_scandata_list))
                yndialog.setStandardButtons(
                    QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Ok)
                yndialog.setDefaultButton(QtGui.QMessageBox.Cancel)
                ynresult = yndialog.exec_()
                if ynresult != QtGui.QMessageBox.Ok:
                    self.ignore_input = False
                    self.statusBar.showMessage("Ready")
                    return
            # assume at this point all-clear to delete old data list
            self.clear_all_scandata(suppress_update=True)
            for scandata in new_scandata_list:
                newscandata = dcparsing.move_scandata_attribute_to_front(
                    scandata, newfield)
                self.add_scandata_to_list(newscandata)
            self.list_scandata.setCurrentRow(0)
            self.temporary_fitdata_list = [None for x in new_scandata_list]
            self.current_scan_list_changed_since_last_update = True
            self.update_state()
            self.ignore_input = False
            self.statusBar.showMessage("Ready")

    def callback_sort_scan_list(self):
        # TODO: delegate to another file, e.g. dataclassparsing
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Processing...")
        primarykey = self.cmb_primarysort.currentText()
        secondarykey = self.cmb_secondarysort.currentText()
        field_index = 0  # since rearranges scandata anyway...
        oldscanlist = self.current_scan_list[:]
        oldtempfitlist = self.temporary_fitdata_list[:]
        self.clear_all_scandata(suppress_update=True)

        def scandatasortfcn(scandata_tempfitdata_tuple):
            scandata, _ = scandata_tempfitdata_tuple
            scaninfo = scandata.scaninfo_list[field_index]
            try:
                return (float(scaninfo[primarykey]),
                        float(scaninfo[secondarykey]))
            except KeyError:
                try:
                    return (float(scaninfo[primarykey]),
                            99999999)
                except KeyError:
                    try:
                        return (99999999,
                                float(scaninfo[secondarykey]))
                    except KeyError:
                        return (99999999, 99999999)
                    except ValueError:
                        print("Numerical sort keys only!")
                        return (99999999, 99999999)
                except ValueError:
                    print("Numerical sort keys only!")
                    return (99999999, 99999999)
            except ValueError:
                print("Numerical sort keys only!")
                return (99999999, 99999999)
            except AttributeError:
                print("Scandata expected as list element.")
                return (99999999, 99999999)

        newscanlist, newtempfitlist = zip(*sorted(
                                            zip(oldscanlist, oldtempfitlist),
                                            key=scandatasortfcn))
        for scandata in newscanlist:
            self.add_scandata_to_list(scandata)
        if len(self.current_scan_list) > 0:
            self.list_scandata.setCurrentRow(0)
        self.current_scan_list_changed_since_last_update = True
        self.update_state()
        self.ignore_input = False
        self.statusBar.showMessage("Ready")

    def callback_perform_fit(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Fitting...")
        fit_all = self.radio_fitall.isChecked()
        # TODO: REAL CODE
        if fit_all:
            self.temporary_fitdata_list = \
                [scandata.fitdata_list[0]
                 for scandata in self.current_scan_list]
        else:
            newtempfitdata = self.get_active_scandata().fitdata_list[0]
            self.replace_active_temp_fitdata(newtempfitdata)
        self.current_scan_list_changed_since_last_update = True
        self.update_state()
        self.ignore_input = False
        self.statusBar.showMessage("Ready")

    def callback_erase_temp_fit(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Erasing Temporary Fit(s)...")
        fit_all = self.radio_fitall.isChecked()
        if fit_all:
            self.temporary_fitdata_list = \
                [None for scandata in self.current_scan_list]
        else:
            self.replace_active_temp_fitdata(None)
        self.current_scan_list_changed_since_last_update = True
        self.update_state()
        self.ignore_input = False
        self.statusBar.showMessage("Ready")

    def callback_store_fit(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Storing Fit(s)...")
        fit_all = self.radio_fitall.isChecked()
        if fit_all:
            new_scan_list = []
            for ind, scandata in enumerate(self.current_scan_list):
                newfitdata_list = [fitdata
                                   for fitdata in scandata.fitdata_list]
                newfitdata_list[0] = self.temporary_fitdata_list[ind]
                self.temporary_fitdata_list[ind] = None
                scandata = ScanData(scandata.fields_list,
                                    scandata.scaninfo_list,
                                    scandata.dataseries_list,
                                    scandata.error_dataseries_list,
                                    newfitdata_list)
                new_scan_list.append(scandata)
            self.current_scan_list = new_scan_list
        else:
            scandata = self.get_active_scandata()
            newfitdata_list = [fitdata for fitdata in scandata.fitdata_list]
            newfitdata_list[0] = self.get_active_temp_fitdata()
            self.replace_active_temp_fitdata(None)
            scandata = ScanData(scandata.fields_list,
                                scandata.scaninfo_list,
                                scandata.dataseries_list,
                                scandata.error_dataseries_list,
                                newfitdata_list)
            self.replace_active_scandata(scandata)
        self.current_scan_list_changed_since_last_update = True
        self.update_state()
        self.ignore_input = False
        self.statusBar.showMessage("Ready")

    def callback_delete_fit(self):
        if self.ignore_input:
            return
        self.ignore_input = True
        self.statusBar.showMessage("Deleting Fit(s)...")
        fit_all = self.radio_fitall.isChecked()
        if fit_all:
            new_scan_list = []
            for scandata in self.current_scan_list:
                newfitdata_list = [None for x in range(len(scandata.fields))]
                scandata = ScanData(scandata.fields_list,
                                    scandata.scaninfo_list,
                                    scandata.dataseries_list,
                                    scandata.error_dataseries_list,
                                    newfitdata_list)
                new_scan_list.append(scandata)
            self.current_scan_list = new_scan_list
        else:
            scandata = self.get_active_scandata()
            newfitdata_list = [None for x in range(len(scandata.fields))]
            scandata = ScanData(scandata.fields_list,
                                scandata.scaninfo_list,
                                scandata.dataseries_list,
                                scandata.error_dataseries_list,
                                newfitdata_list)
            self.replace_active_scandata(scandata)
        self.current_scan_list_changed_since_last_update = True
        self.update_state()
        self.ignore_input = False
        self.statusBar.showMessage("Ready")

    def update_datatype_box_target(self, scandata):
        last_datatype = self.cmb_datatype.currentText()
        self.cmb_datatype.clear()
        for field in scandata.fields:
            self.cmb_datatype.addItem(field)
        if last_datatype in scandata.fields:
            for ind, field in enumerate(scandata.fields):
                if field == last_datatype:
                    self.cmb_datatype.setCurrentIndex(ind)

    def clear_all_scandata(self, suppress_update=False):
        while self.current_scan_list:  # while scandata remain
            self.list_scandata.takeItem(0)
            self.current_scan_list.pop(0)
        self.current_scan_list_changed_since_last_update = True
        if not suppress_update:
            self.update_state()

    def add_scandata_to_list(self, scandata_to_add):
        try:
            midtype_str = scandata_to_add.scaninfo_list[0]['MiddleScanType']
        except (KeyError, IndexError):
            midtype_str = "[Unknown]"
        try:
            midval_str = \
                str(scandata_to_add.scaninfo_list[0]['MiddleScanCoord'])
        except (KeyError, IndexError):
            midval_str = "[Unknown]"
        try:
            fasttype_str = scandata_to_add.scaninfo_list[0]['FastScanType']
        except (KeyError, IndexError):
            fasttype_str = "[Unknown]"
        try:
            start_str = str(scandata_to_add.dataseries_list[0].xvals()[0])
        except (KeyError, IndexError):
            start_str = "[Unknown]"
        try:
            stop_str = str(scandata_to_add.dataseries_list[0].xvals()[-1])
        except (KeyError, IndexError):
            stop_str = "[Unknown]"
        scandata_string = \
            "{midtype}: {midval}, {fasttype}: {start} to {stop}".format(
                midtype=midtype_str, midval=midval_str,
                fasttype=fasttype_str, start=start_str, stop=stop_str)
        self.list_scandata.addItem(scandata_string)
        self.current_scan_list.append(scandata_to_add)
        self.temporary_fitdata_list.append(None)
        self.current_scan_list_changed_since_last_update = True

    def get_active_scandata(self):
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                return self.current_scan_list[rownum]
        return None

    def get_active_temp_fitdata(self):
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                return self.temporary_fitdata_list[rownum]
        return None

    def replace_active_scandata(self, newscandata):
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                self.current_scan_list[rownum] = newscandata

    def replace_active_temp_fitdata(self, newtempfitdata):
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                self.temporary_fitdata_list[rownum] = newtempfitdata

    def update_state(self, new_current_selection=False):
        """
        Updates the current listboxes/plots to reflect the current state.
        """
        plot_active = False
        current_scandata = self.get_active_scandata()
        # Handle current selection - DataType, 1D ScanInfo, and plotting
        if current_scandata is None:
            self.cmb_datatype.clear()
            self.list_scaninfo.clear()
        else:
            if self.current_scan_list_changed_since_last_update or \
                                                        new_current_selection:
                self.update_datatype_box_target(current_scandata)
            self.display_scaninfo(current_scandata.scaninfo_list[0])
            if self.btn_plot1d.isChecked():
                if self.current_scan_list_changed_since_last_update:
                    self.plot_scandata(current_scandata)
                plot_active = True
        # 2D data plotting
        if self.btn_plot2d.isChecked():
            if self.current_scan_list_changed_since_last_update:
                self.plot_all_scandata(self.current_scan_list)
            plot_active = True
        # 2D Fit Parameter plotting
        if self.btn_plotfitparam.isChecked():
            plot_active = True
        # Always at end:
        if not plot_active:
            self.canvas.wipe()
        self.current_scan_list_changed_since_last_update = False

    def display_scaninfo(self, scaninfo):
        if scaninfo is not None:
            if scaninfo != self.last_scaninfo:
                self.list_scaninfo.clear()
                for key, val in scaninfo.items():
                    infostr = "{}: {}".format(key, val)
                    if len(infostr) > 50:
                        infostr = infostr[:50] + "..."
                    self.list_scaninfo.addItem(infostr)
                self.last_scaninfo = scaninfo.copy()

    def plot_scandata(self, scandata):
        """Plot a 1D scandata's primary dataset"""
        if scandata is not None:
            ind = self.cmb_datatype.currentIndex()
            self.statusBar.showMessage("Plotting...")
            self.canvas.wipe()
            dcgraphing.plot_scandata(scandata,
                                     field_index=ind,
                                     title=None,
                                     datatype=scandata.fields[0],
                                     axes=self.canvas.axes,
                                     plot_options={'suppress_legend': True})
            self.canvas.axes.set_aspect('auto')
            self.canvas.figure.tight_layout()
            self.canvas.draw()
            self.statusBar.showMessage("Ready")

    def plot_all_scandata(self, scandata_list):
        """Plot a 2D scandata list"""
        if scandata_list:
            self.statusBar.showMessage("Plotting...")
            self.canvas.wipe()
            # determine data type to plot, and get 1st scandata to be plotted:
            ref_scandata = self.get_active_scandata()
            if ref_scandata is None:
                ref_scandata = scandata_list[0]
            plotfield = self.cmb_datatype.currentText()
            if plotfield not in ref_scandata.fields:
                plotfield = ref_scandata.fields[0]
            # get first row of plot's length to ensure all rows match
            ind = ref_scandata.fields.index(plotfield)
            refdatayvals = \
                ref_scandata.dataseries_list[ind].yvals(unfiltered=True)
            plotdatalength = len(refdatayvals)
            # TODO: delegate to plot_scandata_2d from dataclassgraphing
            data2d = []
            for scandata in scandata_list:
                for ind, field in enumerate(scandata.fields):
                    if field == plotfield:
                        datayvals = scandata.dataseries_list[ind].yvals(
                                                            unfiltered=True)
                        if len(datayvals) == plotdatalength:
                            data2d.append(datayvals)
            if data2d:
                imageplot = self.canvas.axes.imshow(data2d,
                                                    interpolation="nearest")
                self.canvas.axes.set_aspect('auto')
                self.canvas.figure.colorbar(imageplot)
                self.canvas.figure.tight_layout()
                self.canvas.draw()
            self.statusBar.showMessage("Ready")

    @classmethod
    def launch_with_output(cls, app_saved_state=None):
        """
        Class method that creates a DataBrowserWindow instance and returns
        it and a linked AppOutput object
        """
        appoutput = AppOutput()
        apphandle = DataBrowserWindow(appoutput, app_saved_state)
        return apphandle, appoutput



# %%
class MPLCanvas(FigureCanvas):
    """Matplotlib figure that can be added as a widget."""
    def __init__(self, parent=None, width=4, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        super(MPLCanvas, self).__init__(fig)
        self.axes = fig.add_subplot(111)
        self.figure = self.axes.get_figure()
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def wipe(self):
        if len(self.figure.axes) == 1:
            self.axes.cla()
        else:
            self.figure.clf()
            self.axes = self.figure.add_subplot(111)
            self.axes.hold(False)
        self.draw()


# %%
class AppOutput(QtCore.QObject):
    """
    Simple class that exists only to catch the output of a Qt application.
    Retreive with appoutputobject.output after application has exited.
    """
    def __init__(self):
        super().__init__()
        self.output = None
        self.has_output = False

    def receive_output(self, outputdata):
        """
        Used by Qt application to send output data to this object
        """
        self.output = outputdata
        self.has_output = True

