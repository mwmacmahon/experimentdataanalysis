# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:07:59 2016

@author: Michael
"""

import numpy as np

from collections import Iterable, namedtuple
from collections.abc import Sequence


# %%
# returned from curvefitting.py functions
FitData = namedtuple("FitData", ["fitparams", "fitparamstds",
                                 "fitparamstring", "fitdataseries",
                                 "meansquarederror"])
# used to fit scandata in dataclassfitting.py functions [DEPRECATED]
FitFunc = namedtuple("FitFunc", ["description", "fitfunction",
                                 "fitparamlist", "fitargumentlist"])


# %%
# default way of storing each csv file's data for active use
# "scaninfo" is a dict containing scan parameters, e.g. "Voltage: 5"
# note "dataseries" and "fitdata" are PLURAL - a tuple of entries is
# expected, and fields[i], dataseries[i], and fitdata[i] are correlated
ScanData = namedtuple("ScanData", ["fields",
                                   "scaninfo_list",
                                   "dataseries_list",
                                   "error_dataseries_list",
                                   "fitdata_list"])
# old:
# ScanData = namedtuple("ScanData", ["filepath", "scaninfo", "fields",
#                                    "dataseries", "fitdata"])

# %% Helper fcns used by DataSeries equality operator:
# NEEDS DESCRIPTION, TEST, SPHINX DOCUMENTATION
def approx_equal(x, y, tol):
    return abs(x - y) <= 0.5*(abs(x) + abs(y))*tol


# NEEDS DESCRIPTION, TEST, SPHINX DOCUMENTATION
def approx_equal_lists(x_list, y_list, tol):
    x_list = list(x_list)
    y_list = list(y_list)
    return all([approx_equal(x, y, tol) for (x, y) in zip(x_list, y_list)])


# %% NEEDS TEST, SPHINX DOCUMENTATION
class DataSeries(Sequence):
    """
    Defines an immutable data structure that contains two correlated
    iterables, xvals and yvals. When retreiving them with xvals() or
    yvals(), "unsorted=True" can be used to retreive them in the order
    provided upon construction. Adding or subtracting DataSeries
    instances ensures proper matching of yvals corresponding to the
    same xval. Scalar values and lists of matching length may also be
    added, subtracted, multiplied, and divided.

    Calling this object is equivalent to using its data() method
    """
    def __init__(self, xvals_or_datatuples, yvals=None):
        """
        Initializes a data structure that correlates a series of xvals
        with a series of yvals, and when adding or subtracting ensures
        they are properly paired up xvalwise regardless of order of
        data points. All values should be numerical.

        method 1:
        param 1: xvals: iterable of xvals
        param 2: yvals: iterable of yvals
            e.g. DataSeries([1,2,3], [4,5,6])

        method 2:
        param 1: datatuples: iterable of 2-tuples in form (x, y)
            e.g. DataSeries(data), where data = [(x1, y1),(x2, y2), ...]
                 DataSeries(dataseries(unsorted=True))
                 DataSeries((t,f(t)) for t in xvals), where f(t) is unary
        """
        # NOTE: unpack operation forces traversal of RHS. Important below...
        xvals_or_datatuples = list(xvals_or_datatuples)
        if len(xvals_or_datatuples) == 0:
            raise ValueError("DataSeries cannot be initialized with no data")
        if yvals is None:  # handle list of data points
            try:
                xvals, yvals = zip(*xvals_or_datatuples)  # outputs tuples
            except ValueError:
                raise ValueError("if DataSeries is initialized with a list " +
                                 "of data elements, all elements must be " +
                                 "in form (xval, yval)")
        else:
            xvals = tuple(xvals_or_datatuples)
            yvals = tuple(yvals)
            if len(xvals) != len(yvals):
                raise ValueError("uneven number of x and y values given " +
                                 "to DataSeries constructor")
        try:  # ensure numerical values
            xvals = np.array([float(x) for x in xvals])
            yvals = np.array([float(y) for y in yvals])
        except ValueError:
            raise ValueError("all values in DataSeries must be numeric!")
        # sort xy pairs by x and store sorted lists:
        self.map_to_sorted = xvals.argsort()
        self.map_to_unsorted = self.map_to_sorted.argsort()
        self._xvals = np.array(xvals[self.map_to_sorted])
        self._yvals = np.array(yvals[self.map_to_sorted])
        # set all immutable!
        self.map_to_sorted.setflags(write=False)
        self.map_to_unsorted.setflags(write=False)
        self._xvals.setflags(write=False)
        self._yvals.setflags(write=False)

    def xvals(self, *, unsorted=False, raw=False):
        """
        Return this object's xvals as a numpy array, optionally in the
        original order used when this DataSeries was created.
        
        Optional keyword params:
        unsorted (default: False) - if True, return in original ordering
        raw (default: False) - (sets all of the above keywords, legacy keyword)
        """
        if raw:
            unsorted = True
        xvals, _ = self.data(unsorted=unsorted)
        return xvals

    def yvals(self, *, unsorted=False, raw=False):
        """
        Return this object's yvals as a numpy array, optionally in the
        original order used when this DataSeries was created.
        
        Optional keyword params:
        unsorted (default: False) - if True, return in original ordering
        raw (default: False) - (sets all of the above keywords, legacy keyword)
        """
        if raw:
            unsorted = True
        _, yvals = self.data(unsorted=unsorted)
        return yvals

    def datatuples(self, *, unsorted=False, raw=False):
        """
        Return this object's data as a numpy array of (x, y) pairs, optionally
        in the original order used when this DataSeries was created.
        
        Optional keyword params:
        unsorted (default: False) - if True, return in original ordering
        raw (default: False) - (sets all of the above keywords, legacy keyword)
        """
        if raw:
            unsorted = True
        return np.vstack(self.data(unsorted=unsorted)).transpose()

    def data(self, *, unsorted=False, raw=False):
        """
        Return this object's data as two numpy arrays of xvals and yvals,
        optionally in the original order used when this DataSeries was created.
        Format is (xvals, yvals), identical to (xvals(), yvals())
        Consider using "pyplot.plot(*dataseries.datalists())" for plots
        
        Optional keyword params:
        unsorted (default: False) - if True, return in original ordering
        raw (default: False) - (sets all of the above keywords, legacy keyword)
        """
        if raw:
            unsorted = True
        if unsorted:
            xvals = self._xvals[self.map_to_unsorted]
            yvals = self._yvals[self.map_to_unsorted]
            return xvals, yvals
        else:
            return self._xvals, self._yvals

    def copy(self):
        """
        Simple copy constructor making a deep copy of another DataSeries
        object. Takes a single DataSeries instance as a parameter.
        """
        return self.__class__(*self.data(raw=True))

    def copy_subset(self, index_list):
        if any(ind > len(self) for ind in index_list):
            print("Warning: DataSeries.copy_subset(...) index_list " +
                  "parameter contains out of bounds indices, ignoring.")
        return self.__class__(
            [xypair for ind, xypair in enumerate(self.datatuples(raw=True))
             if ind in index_list])

    def __call__(self, *args, **kwargs):
        return self.data(*args, **kwargs)

    def __iter__(self):
        for index in range(len(self)):
            yield self[index]

    def __getitem__(self, index):
        if index >= len(self):
            raise IndexError('Index out of range')
        return self._xvals[index], self._yvals[index]

    def __len__(self):
        return len(self._xvals)

    def __repr__(self):
        xvals, yvals = self.data(raw=True)
        outputstr = "DataSeries({0}, {1})".format(repr(xvals),
                                                  repr(yvals))
        return outputstr

    def __str__(self):
        outputstr = ''
        for index, (xval, yval) in enumerate(zip(self._xvals, self._yvals)):
            if index > 0:
                outputstr += "\n"
            outputstr += "(x = {xval}, y = {yval})".format(xval=xval,
                                                           yval=yval)
        return outputstr

    def general_numeric_fcn(self, other, function, fcn_verb):
        fcn_verb = str(fcn_verb)
        map_to_unsorted = self.map_to_unsorted
        if isinstance(other, self.__class__):  # is a DataSeries
            try:
                otherxvals = other._xvals
                otheryvals = other._yvals
                othermap = other.map_to_unsorted
                if not np.array_equal(self._xvals, otherxvals):
                    raise TypeError('DataSeries instances have different ' +
                                    'xvalues, so failed to ' + fcn_verb)
                if not np.array_equal(self.map_to_unsorted, othermap):
                    print('Warning: combining two DataSeries instances ' +
                          'with different sort orders, resulting sort ' +
                          'order uncertain.')
                    if np.array_equal(map_to_unsorted, np.arange(len(self))):
                        map_to_unsorted = othermap  # if self-sorted, other map
                newyvals = np.array(function(self._yvals, otheryvals))
                newyvals = newyvals[map_to_unsorted]
            except AttributeError:
                print('Warning: failed to combine DataSeries with apparent ' +
                      'DataSeries, treating latter as generic list/array')
                newyvals = function(self.yvals(raw=True), other)
        else:  # not a DataSeries: let numpy handle it
            newyvals = function(self.yvals(raw=True), other)
        newxvals = self._xvals[map_to_unsorted]
        return self.__class__(newxvals, newyvals)

    def __add__(self, other):
        fcn_verb = "add"
        def fcn(x, y):
            return x + y
        return self.general_numeric_fcn(other, fcn, fcn_verb)

    def __sub__(self, other):
        return self.__add__(-1*other)

    def __mul__(self, other):
        fcn_verb = "multiply"
        def fcn(x, y):
            return x * y
        return self.general_numeric_fcn(other, fcn, fcn_verb)

    def __div__(self, other):
        fcn_verb = "divide"
        def fcn(x, y):
            return x / y
        return self.general_numeric_fcn(other, fcn, fcn_verb)

    def __pow__(self, other):
        fcn_verb = "raise"
        def fcn(x, y):
            return x**y
        return self.general_numeric_fcn(other, fcn, fcn_verb)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return -1*self.__sub__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        return -1*self

    def __pos__(self):
        return self

    def __abs__(self):
        xvals, yvals = self.data(raw=True)
        return self.__class__(xvals, abs(yvals))

    def __eq__(self, other):
        try:
            return all(
                [np.array_equal(self._xvals, other._xvals),
                 np.array_equal(self._yvals, other._yvals)])
        except AttributeError:  # not a DataSeries
            return False
