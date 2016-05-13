# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:07:59 2016

@author: Michael
"""

from collections import Iterable, namedtuple
from collections.abc import Sequence


# %%
# returned from curvefitting.py functions
FitData = namedtuple("FitData", ["fitparams", "fitparamstds",
                                 "fitparamstring", "fitdataseries"])
# used to fit scandata in dataclassfitting.py functions [DEPRECATED]
FitFunc = namedtuple("FitFunc", ["description", "fitfunction",
                                 "fitparamlist", "fitargumentlist"])


# %%
# default way of storing each csv file's data for active use
# "scaninfo" is a dict containing scan parameters, e.g. "Voltage: 5"
# note "dataseries" and "fitdata" are PLURAL - a tuple of entries is
# expected, and fields[i], dataseries[i], and fitdata[i] are correlated
ScanData = namedtuple("ScanData", ["filepath", "scaninfo", "fields",
                                   "dataseries", "fitdata"])
## ALTERNATE: (just added a few convenience functions)
#class ScanData(namedtuple("ScanData", ["filepath", "scaninfo", "fields",
#                                       "dataseries", "fitdata"])):
#    def primary_field(self):
#        return self.fields[0]
#
#    def primary_dataseries(self):
#        return self.dataseries[0]
#
#    def primary_fitdata(self):
#        return self.fitdata[0]


# %%
class DataSeries(Sequence):
    """
    Defines an immutable data structure that contains two correlated
    iterables, xvals and yvals. When retreiving them with xvals() or
    yvals(), "unfiltered=False" can be used to retreive the full data
    set, and "unsorted=True" can be used to retreive them in the order
    provided upon construction. Adding or subtracting DataSeries
    instances ensures proper matching of yvals corresponding to the
    same xval. Scalar yvals and lists of matching length may also be
    added or subtracted.

    Calling this object is equivalent to using its datatuples() method
    """
    def __init__(self, datatuples, *, excluded_intervals=None):
        """
        Initializes a data structure that correlates a series of xval points
        with a set of data, and when adding or subtracting ensures data is
        properly paired up xvalwise regardless of order of data points.
        Also allows setting start/end ranges of "excluded" data points
        that are not returned by most calls.

        param 1: datatuples: iterable of 2-tuples in form (data, yval)
            e.g. DataSeries(data), where data = [(t1,v1),(t2,v2),...]
                 DataSeries(dataseries(raw=True))
                 DataSeries(zip(xvals,yvals)), where len(xvals)=len(yvals)
                 DataSeries((t,f(t)) for t in xvals), where f(t) is unary

        optional param 1: excluded_intervals=None: excluded xval ranges given
                              by iterable of 2-tuples in form (start, end)
            e.g. DataSeries(data, excluded_intervals=[(-10000,0),(10,10000)])
        """
        # NOTE: unpack operation forces traversal of RHS. Important below...
        # TIMES AND VALUES
        try:
            xvals, yvals = zip(*datatuples)
        except ValueError:
            raise ValueError("all data elements must be in form (xval, yval)")
        sortedxvallist = sorted(enumerate(xvals),
                                key=lambda tuple: tuple[1])
        self.map_to_sorted, self._xvals = zip(*sortedxvallist)  # unzip
        unsortedxvallist = sorted(enumerate(self.map_to_sorted),
                                  key=lambda tuple: tuple[1])
        self.map_to_unsorted, _ = zip(*unsortedxvallist)  # unzip
        self._yvals = tuple(yvals[index] for index in self.map_to_sorted)

        # EXCLUDED INTERVALS
        if excluded_intervals is not None:
            # may need to iterate multiple xvals to handle
            # None, empty list, single two-element pair, etc.
            excluded_intervals = list(excluded_intervals)
        try:
            if excluded_intervals is None or len(excluded_intervals) == 0:
                self._start_xvals = None
                self._end_xvals = None
            else:
                self._start_xvals, self._end_xvals = zip(*excluded_intervals)
        except TypeError as typeerror:
            # see if they just gave a single interval
            try:
                self._start_xvals, self._end_xvals = zip(excluded_intervals)
            except ValueError:
                raise ValueError("all xval intervals must be in form " +
                                 "(start_xval, end_xval)")
            except TypeError:
                raise(typeerror)
        except ValueError:
            raise ValueError("all xval intervals must be in form " +
                             "(start_xval, end_xval)")

    def xvals(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's xvals, possibly filtered or in the
        original sorting used when this instance was created.
        """
        if raw:
            unfiltered = True
            unsorted = True
        xvals, _ = self.datalists(unfiltered=unfiltered, unsorted=unsorted)
        return xvals

    def yvals(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's yvals, possibly filtered or in the
        original sorting used when this instance was created.
        """
        if raw:
            unfiltered = True
            unsorted = True
        _, yvals = self.datalists(unfiltered=unfiltered, unsorted=unsorted)
        return yvals

    def datatuples(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's xvals and associated yvals, possibly
        filtered or in the original sorting used when this instance
        was created. Note this is a completely lazy operation,
        consisting of chained generator expressions.
        Format is an iterator (t1,v1), (t2,v2), ..., identical to the
        input format.
        """
        if raw:
            unfiltered = True
            unsorted = True
        if unsorted:
            tuples = ((self._xvals[i], self._yvals[i])
                      for i in self.map_to_unsorted)
        else:
            tuples = zip(self._xvals, self._yvals)
        if not unfiltered and self._start_xvals is not None:
            for start, end in zip(self._start_xvals,
                                  self._end_xvals):
                tuples = ((xval, yval) for xval, yval in tuples
                          if xval < start or xval > end)
        for tup in tuples:
            yield tup

    def datalists(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's xvals and associated yvals, possibly
        filtered or in the original sorting used when this instance
        was created.
        Format is (xvals, yvals), identical to (xvals(), yvals())
        Consider using "pyplot.plot(*dataseries.datalists())" for plots
        """
        if raw:
            unfiltered = True
            unsorted = True
        try:
            xvals, yvals = zip(*self.datatuples(unfiltered=unfiltered,
                                                unsorted=unsorted))
            xvals = list(xvals)
            yvals = list(yvals)
        except ValueError:  # happens if all yvals are filtered out
            xvals = []
            yvals = []
        return xvals, yvals

    def excluded_intervals(self):
        """
        Return this instance's list of excluded xval intervals.
        Note this is a completely lazy operation, consisting of
        chained generator expressions.
        Format is an iterator (start, stop), (start,stop), ...,
        identical to the excluded_intervals input format.
        """
        if self._start_xvals is None:
            return
        else:
            for start, end in zip(self._start_xvals, self._end_xvals):
                yield start, end

    def is_excluded(self, xval):
        if self._start_xvals is not None:
            for start, end in zip(self._start_xvals,
                                  self._end_xvals):
                if xval >= start and xval <= end:
                    return True
        return False

    def copy(self):
        """
        Simple copy constructor making a deep copy of another DataSeries
        object. Takes a single DataSeries instance as a parameter.
        """
        return self.__class__(zip(self.xvals(raw=True),
                              self.yvals(raw=True)),
                              excluded_intervals=self.excluded_intervals())

    def __call__(self, *args, **kwargs):
        return self.datatuples(*args, **kwargs)

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
        outputstr = 'DataSeries(zip('
        outputstr += repr(self.xvals(raw=True)) + ', '
        outputstr += repr(self.yvals(raw=True)) + ')'
        if self._start_xvals is not None:
            outputstr += ", excluded_intervals=["
            for interval in list(self.excluded_intervals()):
                outputstr += str(interval) + ", "
            outputstr += "])"
        else:
            outputstr += ")"
        return outputstr

    def __str__(self):
        outputstr = ''
        for index, (xval, yval) in enumerate(zip(self._xvals, self._yvals)):
            if index > 0:
                outputstr += "\n"
            outputstr += "(x = {xval}, y = {yval})".format(xval=xval,
                                                           yval=yval)
            if self.is_excluded(xval):
                outputstr += " [EXCLUDED]"
        return outputstr

    def __add__(self, other):
        excluded_intervals = []
        map_to_unsorted = self.map_to_unsorted
        if isinstance(other, Iterable):
            try:
                otherlist = other._yvals
                otherxvals = other._xvals
                other_intervals = other.excluded_intervals()
            except AttributeError:  # not a DataSeries
                otherlist = list(other)
                if len(otherlist) is not len(self):
                    raise TypeError('attemped to add list of non-matching' +
                                    ' length to DataSeries instance')
            else:  # is a DataSeries
                othermap = other.map_to_unsorted
                if self._xvals != otherxvals:
                    raise TypeError('DataSeries instances have different' +
                                    ' xvalues, so they cannot be added')
                if self.map_to_unsorted != othermap:
                    print('Warning: combining two DataSeries instances' +
                          ' with different sort orders, resulting sort' +
                          ' order uncertain.')
                    if self.map_to_unsorted == tuple(range(len(self))):
                        map_to_unsorted = othermap
                # Combine filter lists, excluding duplicates
                excluded_intervals = list(self.excluded_intervals())
                excluded_intervals.extend(
                    interval for interval in other_intervals
                    if interval not in self.excluded_intervals())
            newyvals = [yval + otherlist[index]
                        for index, yval in enumerate(self._yvals)]
            # Only unsort after adding sorted so 1:1 correspondence between
            # identical xvals. This might be odd when adding a naked list,
            # but better practice is to covert that list to DataSeries first.
            # Don't filter yvals at all, as we are setting filters above.
            newyvals = [newyvals[i] for i in map_to_unsorted]
        else:
            newyvals = [yval + other for yval in self.yvals(unsorted=True,
                                                            unfiltered=True)]
        newxvals = [self._xvals[i] for i in map_to_unsorted]
        return self.__class__(zip(newxvals, newyvals),
                              excluded_intervals=excluded_intervals)

    def __sub__(self, other):
        return self.__add__(-1*other)

    def __mul__(self, other):
        if isinstance(other, Iterable):
            raise TypeError('Can only multiply DataSeries by a scalar value')
        newyvals = [yval * other for yval in self.yvals(raw=True)]
        return self.__class__(zip(self.xvals(raw=True), newyvals),
                              excluded_intervals=self.excluded_intervals())

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
        newyvals = [abs(yval) for yval in self.yvals(raw=True)]
        return self.__class__(zip(self.xvals(raw=True), newyvals),
                              excluded_intervals=self.excluded_intervals())
