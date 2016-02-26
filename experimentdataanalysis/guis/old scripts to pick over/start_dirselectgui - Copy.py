# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:51:31 2016

@author: Michael
"""

import os
import sys

from PyQt4 import QtCore, QtGui

import experimentdataanalysis.guis.dirselectgui as dirselectgui


# %%
class MenuOutput(QtCore.QObject):
    outputsignal = QtCore.pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.output = None
        self.outputsignal.connect(self.receive_output)

    def receive_output(self, outputdata):
        self.output = outputdata


# %%
class DirectorySelectMenu(QtGui.QMainWindow, dirselectgui.Ui_MainWindow):
    def __init__(self, defaultdir="C:\\"):
        """
        TODO: this docstring
        """
        super().__init__()
        self.setupUi(self)
        self.defaultdir = defaultdir  # TODO: check that it is a valid path
        self.currentdir = defaultdir
        #self.menuoutput = menuoutput
        self.btnBrowse.clicked.connect(self._browse_directories)
        self.btnGo.clicked.connect(self._return_directory)

    def _browse_directories(self):
        """
        TODO: this docstring
        """
        self.listWidget.clear()
        getdirstr = "Pick a folder"
        directory = QtGui.QFileDialog.getExistingDirectory(
                        self, getdirstr, directory=self.defaultdir)
        if directory:
            self.currentdir = directory
            for file_name in os.listdir(directory):
                self.listWidget.addItem(file_name)

    def _return_directory(self):
        """
        TODO: this docstring
        """
        output = self.currentdir
        #self.menuoutput.outputsignal.emit(output)
        self.close()


# %%
def testmain():
    """
    TODO: this docstring
    """
    qapp = QtGui.QApplication(sys.argv)
    menuoutput = MenuOutput()
    menu = DirectorySelectMenu()
    menu.show()
    qapp.exec_()
    dirpath = menuoutput.output
    print(dirpath)


# %%
if __name__ == "__main__":
    testmain()
