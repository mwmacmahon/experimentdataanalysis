# -*- coding: utf-8 -*-
"""
Contains the necessary QApplicationStarter() singleton class to start using
small Qt windows. Also contains several convenience functions to bring up
common basic gui windows like file/directory selection.

Created on Wed Jan 20 01:36:24 2016

@author: Michael
"""

from PyQt4 import QtGui


# %%
class QApplicationStarter:
    """
    Starts a QApplication. Since this code should only be run once, it is
    run as a singleton. It would be better to force the program to only
    execute it once, but the current paradigm involves multiple entrance
    functions that might be called in any order, multiple times, hence
    the need for a singleton.

    Ensures only one version of this object exists using the 'Borg'
    singleton technique. See http://www.aleax.it/Python/5ep.html
    """
    _shared_state = {"qapp": None}

    def __init__(self):
        self.__dict__ = self._shared_state
        if self.qapp is None:
            self.qapp = QtGui.QApplication([])

    def __getattr__(self, name):
        return getattr(self.qapp, name)  # never called on qapp, already exists

    def __setattr__(self, name, value):
        if name in ["__dict__", "qapp"]:
            object.__setattr__(self, name, value)
        else:
            setattr(self.qapp, name, value)


# %%
def get_dir_dialog(*, defaultpath=None):
    """
    Opens a directory browser dialog to return a directory
    """
    qapp = QApplicationStarter()
    kwargs = {'directory': defaultpath}
    directory = QtGui.QFileDialog.getExistingDirectory(
                        None, "Choose a directory", **kwargs)
    return directory


# %%
def get_file_dialog(*, defaultpath=None, extensionfilter=None):
    """
    Opens a directory browser dialog to return a single file.
    The optional 'extensionfilter' keyword argument should be formatted like:
    "Images (*.png *.xpm *.jpg)"
    or
    "Images (*.png *.xpm *.jpg);;Text files (*.txt)"
    """
    qapp = QApplicationStarter()
    kwargs = {'directory': defaultpath,
              'filter': extensionfilter}
    directory = QtGui.QFileDialog.getOpenFileName(
                        None, "Choose a file", **kwargs)
    return directory
