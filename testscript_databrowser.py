

from matplotlib.backends.backend_qt4agg \
    import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from PyQt4 import QtCore, QtGui, uic

from experimentdataanalysis.analysis.dataclasses \
    import FitData, FitFunc, ScanData, DataSeries
import experimentdataanalysis.analysis.dataclassfitting as dcfitting
import experimentdataanalysis.analysis.dataclassgraphing as dcgraphing
import experimentdataanalysis.parsing.dataclassparsing as dcparsing
import experimentdataanalysis.guis.databrowser as databrowser
from experimentdataanalysis.guis.guistarter import QApplicationStarter


# %%
if __name__ == "__main__":
    qapp = QApplicationStarter()
#    scandata_list = list(dcparsing.fetch_dir_as_unfit_scandata_iterator(
#       )
#        directorypath="C:\\Data\\febdata\\Experiment - Channel 2"))
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelCenter_200mT_Channel1_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime"))
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelCenter_200mT_Channel2_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime"))
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelCenter_200mT_Channel2_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime_run2"))
#        directorypath="C:\\Data\\160306\\DelayScan_OnChannelCenter_200mT_Channel3_033XT-B11_819.0nm_30K_2Dscan_Voltage_DelayTime"))
#        directorypath="C:\Data\decdata\Channel 3 Run 1"))
#        directorypath="C:\\Data\\decdata\\representative"))

#    scandata_list = list(dcfitting.fit_scandata_iterable(
#        scandata_list,
#        dataseriesfitfunction=None,
#        dataseriesfitfunction=dcfitting.fit_dataseries_with_one_decaying_cos,
#        dataseriesfitfunction=dcfitting.fit_dataseries_with_two_decaying_cos,
#        fit_drift=True, multiprocessing=True))

    app_saved_state = {'current_scan_list': scandata_list}
    window, windowoutput = \
        databrowser.DataBrowserWindow.launch_with_output(app_saved_state)
#    window, windowoutput = \
#        databrowser.DataBrowserWindow.launch_with_output()
    qapp.exec_()

#    scandata_list = windowoutput.output[0]["current_scan_list"]
