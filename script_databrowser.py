# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:27:08 2016

@author: Michael
"""

import experimentdataanalysis.guis.databrowser as databrowser
from experimentdataanalysis.guis.guistarter import QApplicationStarter


# %%
if __name__ == "__main__":

# %%    
    scandata_list = []
# %%
    qapp = QApplicationStarter()
    app_saved_state = {'current_scan_list': scandata_list}
    window, windowoutput = \
        databrowser.DataBrowserWindow.launch_with_output(app_saved_state)
#    window, windowoutput = \
#        databrowser.DataBrowserWindow.launch_with_output()
    qapp.exec_()

#    scandata_list = windowoutput.output[0]["current_scan_list"]
