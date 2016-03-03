# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:27:08 2016

@author: Michael
"""


from matplotlib.backends.backend_qt4agg \
    import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt4 import QtCore, QtGui, uic

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
import experimentdataanalysis.analysis.dataclassfitting as dcfitting
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
        uic.loadUi('experimentdataanalysis\\guis\\ui_databrowser.ui', self)
        self.fig = MPLCanvas(dpi=90)
        self.figure_container.addWidget(self.fig)
        # Set up application signals
        self.btn_loadfile.pressed.connect(self.load_data_file)
        self.btn_loaddir.pressed.connect(self.load_data_dir)
# FROM SERIESEDITOR - FOR REFERENCE ONLY
#        self.btn_seriesedit.pressed.connect(self.edit_scandata)
#        self.btn_seriesdelete.pressed.connect(self.delete_scandata)
#        self.btn_seriesadd.pressed.connect(self.add_scandata)
#        self.radio_linear.toggled.connect(self.update_settings)
#        self.radio_log.toggled.connect(self.update_settings)
#        self.edit_min.editingFinished.connect(self.update_settings)
#        self.edit_max.editingFinished.connect(self.update_settings)
#        self.edit_numvalues.editingFinished.connect(self.update_settings)
#        self.check_combinescandata.toggled.connect(self.update_settings)
#        self.check_breaksubsections.toggled.connect(self.update_settings)
#        self.slider_numsubsections.valueChanged.connect(self.update_settings)
#        self.check_shuffle.toggled.connect(self.update_settings)
#        self.btn_finish.pressed.connect(self.close)
        self.connect(self,  # emitter
                     QtCore.SIGNAL('output_signal(PyQt_PyObject)'),  # signal
                     appoutput.receive_output)  # receiver
        self.current_scan_list = []
        self.last_scaninfo = None
        if app_saved_state is not None:
            for scandata in app_saved_state['current_scan_list']:
                self.add_scandata_to_list(scandata)
        self.update_settings()
        self.show()

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

    def load_data_dir(self):
        try:
            scandata_list = \
                list(dcparsing.fetch_dir_as_unfit_scandata_iterator())
        # TODO: catch exceptions!
        except:
            raise
        for scandata in scandata_list:
            self.add_scandata_to_list(scandata)

    def load_data_file(self):
        try:
            scandata = dcparsing.fetch_csv_as_unfit_scandata()
        # TODO: catch exceptions!
        except:
            raise
        self.add_scandata_to_list(scandata)

    def add_scandata_to_list(self, scandata_to_add):
        self.current_scan_list.append(scandata_to_add)
        try:
            midtype_str = scandata_to_add.scaninfo['MiddleScanType']
        except KeyError:
            midtype_str = "[Unknown]"
        try:
            midval_str = str(scandata_to_add.scaninfo['MiddleScanCoord'])
        except KeyError:
            midval_str = "[Unknown]"
        try:
            fasttype_str = scandata_to_add.scaninfo['FastScanType']
        except KeyError:
            fasttype_str = "[Unknown]"
        try:
            start_str = str(scandata_to_add.dataseries.xvals()[0])
        except KeyError:
            start_str = "[Unknown]"
        try:
            stop_str = str(scandata_to_add.dataseries.xvals()[-1])
        except KeyError:
            stop_str = "[Unknown]"
        scandata_string = \
            "{midtype}: {midval}, {fasttype}: {start} to {stop}".format(
                midtype=midtype_str, midval=midval_str,
                fasttype=fasttype_str, start=start_str, stop=stop_str)
        self.list_currentscan.addItem(scandata_string)
        self.update_settings()

    def get_active_scandata(self):
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                return self.current_scan_list[rownum]
        return None

    def exclude_scandata(self):
        if self.current_scan_list:  # if any scandata have been saved
            rownum = self.list_scandata.currentRow()
            if rownum >= 0:  # if a row has been selected at all
                self.list_scandata.takeItem(rownum)
                self.current_scan_list.pop(rownum)
                self.update_settings()

    def display_scaninfo(self, scaninfo):
        if scaninfo != self.last_scaninfo:
            # INSERT DISPLAYING MOST DICT ENTRIES
            self.last_scaninfo = scaninfo.copy()

    def update_settings(self):
        """
        Updates the current listboxes/plots to reflect the data inputs,
        and calls are_settings_valid() to ensure data inputs are valid.
        Calls update_plot() if so.
        """
        if self.are_settings_valid():
            if self.btn_plot1d.isChecked():
                pass
            if self.btn_plot2d.isChecked():
                pass
            if self.btn_plotfitparam.isChecked():
                pass
            self.update_plot()
        else:
            self.active_series = None
            self.fig.clf()

    def are_settings_valid(self):
        """Returns False if any window settings are invalid."""
        valid = True
        # TODO: ACTUAL CHECK
        return valid

    def update_plot(self):
        """Update plot to match current settings"""
        scandata = self.get_active_scandata()
        if scandata is not None:
            self.fig.axes.plot(scandata.xvals(), scandata.yvals(), '.')
            self.fig.draw()

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
    def __init__(self, parent=None, width=3, height=3, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        super(MPLCanvas, self).__init__(fig)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


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


# %%
# For debugging only, likely won't play nice with other gui windows
if __name__ == '__main__':
    qapp = QApplicationStarter()
    # qapp.exec_()  # moved to QApplicationStarter
    window, windowoutput = DataBrowserWindow.launch_with_output()
