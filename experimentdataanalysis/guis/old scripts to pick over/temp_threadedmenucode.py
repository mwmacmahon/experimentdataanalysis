# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:51:31 2016

@author: Michael
"""

import os
import sys

from PyQt4 import QtGui
from PyQt4.QtCore import QThread

import experimentdataanalysis.guis.dummyfilemenu as directoryselectgui


# %%
class DirectoryBrowseThread(QThread):

    def __init__(self, startdir):
        """
        TODO: this docstring
        """
        QThread.__init__(self)  # TODO: try switching to super()
        self.startdir = startdir

    def __del__(self):
        self.wait()

    def _browse_folder(self):
        """
        TODO: this docstring
        """
        self.listWidget.clear()
        getdirstr = "Pick a folder"
        directory = QtGui.QFileDialog.getExistingDirectory(self, getdirstr)
        if directory:
            for file_name in os.listdir(directory):
                self.listWidget.addItem(file_name)

    def run(self):
        """
        TODO: this docstring
        """
        self.browse_folder()
        self.done()
        # self.sleep(2)  # use this for self blocking


# %%
class DirectorySelectMenu(QtGui.QMainWindow,
                          directoryselectgui.Ui_MainWindow):
    def __init__(self):
        """
        TODO: this docstring
        """
        super().__init__()
        self.setupUi(self)
        self.btnBrowse.clicked.connect(self.browse_folder)

    def done(self):
        """
        TODO: this docstring
        """
        QtGui.QMessageBox.information(self, "Done",
                                      "Browse thread has finished processing.")

    def _start_browse_folder(self):
        """
        TODO: this docstring
        """
        self.dirbrowse_thread = DirectoryBrowseThread("C:\\pythonprojects\\")
        self.connect(self.dirbrowse_thread, SIGNAL("finished()"), self.done)
        self.dirbrowse_thread.start()


# %%
def main():
    """
    TODO: this docstring
    """
    # following https://nikolak.com/pyqt-qt-designer-getting-started/
    app = QtGui.QApplication(sys.argv)
    form = DirectorySelectMenu()
    form.show()
    app.exec_()


# %%
if __name__ == "__main__":
    main()
