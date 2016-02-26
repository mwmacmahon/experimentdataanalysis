# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:07:59 2016

@author: Michael
"""

from collections import Iterable, namedtuple
from collections.abc import Sequence


# %%
# returned from curvefitting.py functions
FitData = namedtuple("FitData", ["fitparams", "fiterror",
                                 "fitparamstring", "fittimeseries"])
# default way of storing each csv file's data for active use
# "scaninfo" is a dict containing scan parameters, e.g. "Voltage: 5"
ScanData = namedtuple("ScanData", ["filepath", "scaninfo",
                                   "timeseries", "fitdata"])


# %%
class TimeSeries(Sequence):
    """
    Defines an immutable data structure that contains two correlated
    iterables, times and values. When retreiving them with times() or
    values(), "unfiltered=False" can be used to retreive the full data
    set, and "unsorted=True" can be used to retreive them in the order
    provided upon construction. Adding or subtracting TimeSeries
    instances ensures proper matching of values corresponding to the
    same time. Scalar values and lists of matching length may also be
    added or subtracted.

    Calling this object is equivalent to using its datatuples() method
    """
    def __init__(self, datatuples, *, excluded_intervals=None):
        """
        Initializes a data structure that correlates a series of time points
        with a set of data, and when adding or subtracting ensures data is
        properly paired up timewise regardless of order of data points.
        Also allows setting start/end ranges of "excluded" data points
        that are not returned by most calls.

        param 1: datatuples: iterable of 2-tuples in form (data, value)
            e.g. TimeSeries(data), where data = [(t1,v1),(t2,v2),...]
                 TimeSeries(timeseries(raw=True))
                 TimeSeries(zip(times,values)), where len(times)=len(values)
                 TimeSeries((t,f(t)) for t in times), where f(t) is unary

        optional param 1: excluded_intervals=None: excluded time ranges given
                              by iterable of 2-tuples in form (start, end)
            e.g. TimeSeries(data, excluded_intervals=[(-10000,0),(10,10000)])
        """
        # NOTE: unpack operation forces traversal of RHS. Important below...
        # TIMES AND VALUES
        try:
            times, values = zip(*datatuples)
        except ValueError:
            raise ValueError("all data elements must be in form (time, value)")
        sortedtimelist = sorted(enumerate(times),
                                key=lambda tuple: tuple[1])
        self.map_to_sorted, self._times = zip(*sortedtimelist)  # unzip
        unsortedtimelist = sorted(enumerate(self.map_to_sorted),
                                  key=lambda tuple: tuple[1])
        self.map_to_unsorted, _ = zip(*unsortedtimelist)  # unzip
        self._values = tuple(values[index] for index in self.map_to_sorted)

        # EXCLUDED INTERVALS
        if excluded_intervals is not None:
            # may need to iterate multiple times to handle
            # None, empty list, single two-element pair, etc.
            excluded_intervals = list(excluded_intervals)
        try:
            if excluded_intervals is None or len(excluded_intervals) == 0:
                self._start_times = None
                self._end_times = None
            else:
                self._start_times, self._end_times = zip(*excluded_intervals)
        except TypeError as typeerror:
            # see if they just gave a single interval
            try:
                self._start_times, self._end_times = zip(excluded_intervals)
            except ValueError:
                raise ValueError("all time intervals must be in form " +
                                 "(start_time, end_time)")
            except TypeError:
                raise(typeerror)
        except ValueError:
            raise ValueError("all time intervals must be in form " +
                             "(start_time, end_time)")

    def times(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's times, possibly filtered or in the
        original sorting used when this instance was created.
        """
        if raw:
            unfiltered = True
            unsorted = True
        times, _ = self.datalists(unfiltered=unfiltered, unsorted=unsorted)
        return times

    def values(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's values, possibly filtered or in the
        original sorting used when this instance was created.
        """
        if raw:
            unfiltered = True
            unsorted = True
        _, values = self.datalists(unfiltered=unfiltered, unsorted=unsorted)
        return values

    def datatuples(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's times and associated values, possibly
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
            tuples = ((self._times[i], self._values[i])
                      for i in self.map_to_unsorted)
        else:
            tuples = zip(self._times, self._values)
        if not unfiltered and self._start_times is not None:
            for start, end in zip(self._start_times,
                                  self._end_times):
                tuples = ((time, value) for time, value in tuples
                          if time < start or time > end)
        for tup in tuples:
            yield tup

    def datalists(self, *, unfiltered=False, unsorted=False, raw=False):
        """
        Return this instance's times and associated values, possibly
        filtered or in the original sorting used when this instance
        was created.
        Format is (times, values), identical to (times(), values())
        Consider using "pyplot.plot(*timeseries.datalists())" for plots
        """
        if raw:
            unfiltered = True
            unsorted = True
        try:
            times, values = zip(*self.datatuples(unfiltered=unfiltered,
                                                 unsorted=unsorted))
            times = list(times)
            values = list(values)
        except ValueError:  # happens if all values are filtered out
            times = []
            values = []
        return times, values

    def excluded_intervals(self):
        """
        Return this instance's list of excluded time intervals.
        Note this is a completely lazy operation, consisting of
        chained generator expressions.
        Format is an iterator (start, stop), (start,stop), ...,
        identical to the excluded_intervals input format.
        """
        if self._start_times is None:
            return
        else:
            for start, end in zip(self._start_times, self._end_times):
                yield start, end

    def is_excluded(self, time):
        if self._start_times is not None:
            for start, end in zip(self._start_times,
                                  self._end_times):
                if time >= start and time <= end:
                    return True
        return False

    def copy(self):
        """
        Simple copy constructor making a deep copy of another TimeSeries
        object. Takes a single TimeSeries instance as a parameter.
        """
        return self.__class__(zip(self.times(raw=True),
                              self.values(raw=True)),
                              excluded_intervals=self.excluded_intervals())

    def __call__(self, *args, **kwargs):
        return self.datatuples(*args, **kwargs)

    def __iter__(self):
        for index in range(len(self)):
            yield self[index]

    def __getitem__(self, index):
        if index >= len(self):
            raise IndexError('Index out of range')
        return self._times[index], self._values[index]

    def __len__(self):
        return len(self._times)

    def __repr__(self):
        outputstr = 'TimeSeries(zip('
        outputstr += repr(self.times(raw=True)) + ', '
        outputstr += repr(self.values(raw=True)) + ')'
        if self._start_times is not None:
            outputstr += ", excluded_intervals=["
            for interval in list(self.excluded_intervals()):
                outputstr += str(interval) + ", "
            outputstr += "])"
        else:
            outputstr += ")"
        return outputstr

    def __str__(self):
        outputstr = ''
        for index, (time, value) in enumerate(zip(self._times, self._values)):
            if index > 0:
                outputstr += "\n"
            outputstr += "t = {time}: {value}".format(
                         time=time, value=value)
            if self.is_excluded(time):
                outputstr += " [EXCLUDED]"
        return outputstr

    def __add__(self, other):
        excluded_intervals = []
        map_to_unsorted = self.map_to_unsorted
        if isinstance(other, Iterable):
            if isinstance(other, self.__class__):
                if self._times != other._times:
                    raise TypeError('TimeSeries instances have different' +
                                    ' time values, so they cannot be added')
                if self.map_to_unsorted != other.map_to_unsorted:
                    print('Warning: combining two TimeSeries instances' +
                          ' with different sort orders, resulting sort' +
                          ' order uncertain.')
                    if self.map_to_unsorted == tuple(range(len(self))):
                        map_to_unsorted = other.map_to_unsorted
                otherlist = other._values
                # Combine filter lists, excluding duplicates
                excluded_intervals = list(self.excluded_intervals())
                excluded_intervals.extend(
                    interval for interval in list(other.excluded_intervals())
                    if interval not in self.excluded_intervals())
#                if other.start_time is not None:
#                    if self.start_time is not None:
#                        start_time = max(self.start_time, other.start_time)
#                    else:
#                        start_time = other.start_time
#                if other.end_time is not None:
#                    if self.end_time is not None:
#                        end_time = min(self.end_time, other.end_time)
#                    else:
#                        end_time = other.end_time
            else:
                otherlist = list(other)
                if len(otherlist) is not len(self):
                    raise TypeError('attemped to add list of non-matching' +
                                    ' length to TimeSeries instance')
            newvalues = [value + otherlist[index]
                         for index, value in enumerate(self._values)]
            # Only unsort after adding sorted so 1:1 correspondence between
            # identical times. This might be odd when adding a naked list,
            # but better practice is to covert that list to TimeSeries first.
            # Don't filter values at all, as we are setting filters above.
            newvalues = [newvalues[i] for i in map_to_unsorted]
        else:
            newvalues = [value + other
                         for value in self.values(unsorted=True,
                                                  unfiltered=True)]
        newtimes = [self._times[i] for i in map_to_unsorted]
        return self.__class__(zip(newtimes, newvalues),
                              excluded_intervals=excluded_intervals)

    def __sub__(self, other):
        return self.__add__(-1*other)

    def __mul__(self, other):
        if isinstance(other, Iterable):
            raise TypeError('Can only multiply TimeSeries by a scalar value')
        newvalues = [value * other
                     for value in self.values(raw=True)]
        return self.__class__(zip(self.times(raw=True),
                              newvalues),
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
        newvalues = [abs(value)
                     for value in self.values(raw=True)]
        return self.__class__(zip(self.times(raw=True),
                              newvalues),
                              excluded_intervals=self.excluded_intervals())
