# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:07:59 2016

@author: Michael
"""

import numpy as np
import re  # regex commands for fixing invalid attribute names

from collections import Iterable, namedtuple
from collections.abc import Sequence
from copy import deepcopy


# %%
# returned from curvefitting.py functions
FitData = namedtuple("FitData", ["fitparams", "fitparamstds",
                                 "fitparamstring", "fityvals",
                                 "freeparamindices", "covariancematrix",
                                 "meansquarederror"])
## used to fit scandata in dataclassfitting.py functions [DEPRECATED]
#FitFunc = namedtuple("FitFunc", ["description", "fitfunction",
#                                 "fitparamlist", "fitargumentlist"])


# %% SCANDATA
class ScanDataReservedAlias():
    def __init__(self):
        self.alias = None  # set to be x, y, yerr, etc. by metaclass
        self.storage_alias = None  # same, but _x, _y, etc. Normally unused.

    def __get__(self, scandata_instance, scandata_class):
        if scandata_instance is None: return self
        fetch_fcn = self.get_other  # default behavior for unknown alias
        if self.alias == 'x': fetch_fcn = self.get_x
        if self.alias == 'y': fetch_fcn = self.get_y
        if self.alias == 'yerr': fetch_fcn = self.get_yerr
        if self.alias == 'xy': fetch_fcn = self.get_xy
        if self.alias == 'xyerr': fetch_fcn = self.get_xyerr
        if self.alias == 'yyerr': fetch_fcn = self.get_yyerr
        if self.alias == 'xyyerr': fetch_fcn = self.get_xyyerr
        if self.alias == 'fitdata': fetch_fcn = self.get_fitdata
        return fetch_fcn(scandata_instance)

    def __set__(self, scandata_instance, value):
        error_msg = ("Cannot set attribute for multiple-target " +
                     "reserved alias '{}'".format(self.alias))
        # if one of these aliases w/ multiple field connections, throw error
        if self.alias == 'xy': raise AttributeError(error_msg)
        if self.alias == 'xyerr': raise AttributeError(error_msg)
        if self.alias == 'yyerr': raise AttributeError(error_msg)
        if self.alias == 'xyyerr': raise AttributeError(error_msg)
        # otherwise actually set value
        store_fcn = self.set_other
        if self.alias == 'x': store_fcn = self.set_x
        if self.alias == 'y': store_fcn = self.set_y
        if self.alias == 'yerr': store_fcn = self.set_yerr
        if self.alias == 'fitdata': store_fcn = self.set_fitdata
        store_fcn(scandata_instance, value)

    def get_x(self, scandata_instance):  # throw AttributeError if not found
        return getattr(scandata_instance, scandata_instance.xfield)

    def get_y(self, scandata_instance):  # throw AttributeError if not found
        return getattr(scandata_instance, scandata_instance.yfield)

    def get_yerr(self, scandata_instance):  # default to None if not found
        return getattr(scandata_instance,
                       scandata_instance.yfield + '_error', None)

    def get_xy(self, scandata_instance):
        x = self.get_x(scandata_instance)
        y = self.get_y(scandata_instance)
        return x, y

    def get_xyerr(self, scandata_instance):
        x = self.get_x(scandata_instance)
        yerr = self.get_yerr(scandata_instance)
        return x, yerr

    def get_yyerr(self, scandata_instance):
        y = self.get_y(scandata_instance)
        yerr = self.get_yerr(scandata_instance)
        return y, yerr

    def get_xyyerr(self, scandata_instance):
        x, y = self.get_xy(scandata_instance)
        yerr = self.get_yerr(scandata_instance)
        return x, y, yerr

    def get_fitdata(self, scandata_instance):
        info_dict = scandata_instance.info
        return info_dict.get('fitdata_' + scandata_instance.yfield, None)

    def get_other(self, scandata_instance):  # return stored or throw error
        return getattr(scandata_instance, self.storage_alias)

    def set_array_check(self, scandata_instance, value):
        new_array = np.array(value)
        if new_array.shape == scandata_instance.shape:
            return new_array
        else:
            raise ValueError("attemping to replace a ScanData field " +
                             "with one of different shape " +
                             "({} -> {})".format(scandata_instance.shape,
                                                 new_array.shape))

    def set_x(self, scandata_instance, value):
        value = self.set_array_check(scandata_instance, value)
        setattr(scandata_instance, scandata_instance.xfield, value)

    def set_y(self, scandata_instance, value):
        value = self.set_array_check(scandata_instance, value)
        setattr(scandata_instance, scandata_instance.yfield, value)

    def set_yerr(self, scandata_instance, value):
        value = self.set_array_check(scandata_instance, value)
        setattr(scandata_instance,
                scandata_instance.yfield + '_error', value)

    def set_fitdata(self, scandata_instance, value):
        info_dict = scandata_instance.info
        info_dict['fitdata_' + scandata_instance.yfield] = value

    def set_other(self, scandata_instance, value):
        setattr(scandata_instance, self.storage_alias, value)


class SupportsReservedAliases(type):  # note: run only when _class_ is defined
    def __new__(meta, name, bases, class_dict):
        for key, value in class_dict.items():  # for each invoked descriptor
            if isinstance(value, ScanDataReservedAlias):
                value.alias = key
                value.storage_alias = '_' + key  # storage key inside ScanData
                class_dict['reserved_names'].append(key)
        modified_class = type.__new__(meta, name, bases, class_dict)
        return modified_class


def make_into_attribute_name(original_name):
    # grabbed some regular expressions to fix names off stackoverflow:
    new_name = original_name
    new_name = re.sub('[^0-9a-zA-Z_]', '', new_name)  # no invalid characters
    new_name = re.sub('^[^a-zA-Z_]+', '', new_name)  # first char _ or letter
    if new_name == "":
        raise ValueError("No valid attribute name can " +
                         "be derived from '{}'".format(original_name))
    is_changed = not (new_name==original_name)
    return is_changed, new_name


class ScanData(object, metaclass=SupportsReservedAliases):
    """
    Data object meant to store the equivalent of a table of data and retreive
    columns (as numpy arrays) via attributes: e.g. scandata.delaytimes.
    Reads in a list of field names and field arrays (or any iterables, but all
    iterables must be same shape and convertible to numpy arrays), and
    arrays will be saved as an attribute under the associated field name.
    A dictionary containing info from file header/filename or other sources
    is also provided as scandata.info.

    Field arrays may be multi-dimensional, but only if passed as numpy arrays.

    Note the "x" column will be the first column by default, and "y" the
    second. These can be overwritten at instantiation or by changing the
    attributes 'xfield' and 'yfield'. The following shortcut attributes and
    methods can be used to get one or more fields with a single, short call:

    Shortcut attributes:
    scandata.x -> array saved under attribute name given by scandata.xfield
    scandata.y -> array saved under attribute name given by scandata.yfield
    scandata.yerr -> array saved under attribute name "[scandata.yfield]_error"
    scandata.xy -> returns 2-tuple (scandata.x, scandata.y)
    scandata.xyerr -> returns 2-tuple (scandata.x, scandata.yerr)
    scandata.yyerr -> returns 2-tuple (scandata.y, scandata.yerr)
    scandata.xyyerr -> returns 3-tuple (scandata.x, scandata.y, scandata.yerr)
    scandata.fitdata -> returns FitData associated with scandata.y, or None

    Shortcut methods:
    scandata.get_field_y(field_name) -> same as scandata.y, but using
                                        field_name instead of scandata.yfield
    scandata.get_field_yerr(field_name) -> same, but scandata.yerr
    scandata.get_field_xy(field_name) -> same, but scandata.xy
    scandata.get_field_xyerr(field_name) -> same, but scandata.xyerr
    scandata.get_field_yyerr(field_name) -> same, but scandata.yyerr
    scandata.get_field_xyyerr(field_name) -> same, but scandata.xyyerr
    scandata.get_field_fitdata(field_name) -> same, but scandata.fitdata

    Other methods:
    
    """
    reserved_names = ['reserved_names', 'info', 'shape',
                      'fields', 'xfield', 'yfield']
    # un-settable descriptor attributes that return notable fields:
    x = ScanDataReservedAlias()  # note: these are all added to above list
    y = ScanDataReservedAlias()
    yerr = ScanDataReservedAlias()
    xy = ScanDataReservedAlias()  # not to be confused with get_field_xy(field)
    xyerr = ScanDataReservedAlias()  # or with get_field_xyerr(field)
    yyerr = ScanDataReservedAlias()  # or with get_field_yyerr(field)
    xyyerr = ScanDataReservedAlias()  # or with get_field_xyyerr(field)
    fitdata = ScanDataReservedAlias()  # or with get_field_fitdata(field)

    def __init__(self, field_names, field_arrays, info_dict={},
                 x_field_name=None, y_field_name=None):
        self.info = info_dict.copy()  # in case initialized from another!
        field_name_changed = False
        if x_field_name:
            is_changed, self.xfield = make_into_attribute_name(x_field_name)
            field_name_changed = field_name_changed or is_changed
        if y_field_name:
            is_changed, self.yfield = make_into_attribute_name(y_field_name)
            field_name_changed = field_name_changed or is_changed
        self.fields = []
        for raw_fieldname, field_array in zip(field_names, field_arrays):
            is_changed, fieldname = make_into_attribute_name(raw_fieldname)
            field_name_changed = field_name_changed or is_changed
            if fieldname in self.__class__.reserved_names:
                print('ScanData name conflict, renaming {} to {}'.format(
                        fieldname, fieldname + '_data'))
                fieldname = fieldname + '_data'
                field_name_changed = True
            field_array = np.array(field_array)
            field_array.flags.writeable = False
            if "shape" not in self.__dict__:
                self.shape = field_array.shape
                if "xfield" not in self.__dict__ : # first col is x by default
                    self.xfield = fieldname
            else:  # fields 2+
                if "yfield" not in self.__dict__ : # second col is y by default
                    self.yfield = fieldname
                if self.shape != field_array.shape:  # all fields same shape!
                    errstr = "ScanData field arrays not of " + \
                             "matching shape! shapes: "
                    for field in self.fields:
                        errstr += str(self.get(field).shape) + ", "
                    errstr += str(self.shape)
                    raise ValueError(errstr)
            setattr(self, fieldname, field_array)
            self.fields.append(fieldname)

        # check for issues during initialization, warn or raise error as needed
        if len(self.fields) < 2:  # not enough columns for x, y!
            raise ValueError("ScanData requires at least two fields!")
        if self.xfield not in self.fields:  # never found x field!
            print("ScanData: field {} not found, x field set to {}".format(
                    self.xfield, self.fields[0]))
            self.xfield = self.fields[0]
        if self.yfield not in self.fields:  # never found y field!
            print("ScanData: field {} not found, y field set to {}".format(
                    self.yfield, self.fields[1]))
            self.yfield = self.fields[1]
        if field_name_changed:
            print("ScanData: one or more field names changed, " +
                  "new field names: {}".format(self.fields))

    def get_field_y(self, field_name):
        y = getattr(self, field_name, None)
        return y

    def get_field_yerr(self, field_name):
        yerr = getattr(self, field_name + '_error', None)
        return yerr

    def get_field_xy(self, field_name):
        x = self.x
        y = self.get_field_y(field_name)
        return x, y

    def get_field_xyerr(self, field_name):
        x = self.x
        yerr = self.get_field_yerr(field_name)
        return x, yerr

    def get_field_yyerr(self, field_name):
        y = self.get_field_y(field_name)
        yerr = self.get_field_yerr(field_name)
        return y, yerr

    def get_field_xyyerr(self, field_name):
        x, y = self.get_field_xy(field_name)
        yerr = self.get_field_yerr(field_name)
        return x, y, yerr

    def get_field_fitdata(self, field_name):
        return self.info.get('fitdata_' + field_name, None)

    def copy(self):
        new_fields = self.fields[:]
        new_fields += [field_name + '_error'
                       for field_name in self.fields]
        new_field_arrays = [np.array(getattr(self, field_name))
                            for field_name in self.fields]
        new_field_arrays += [np.array(self.get_field_yerr(field_name))
                             for field_name in self.fields
                             if field_name + '_error' in self.__dict__]
        new_info = deepcopy(self.info)
        new_scandata = self.__class__(new_fields, new_field_arrays,
                                      new_info, self.xfield, self.yfield)
        for key, value in self.__dict__.items():  # other stuff in dicts
            if key not in self.fields:
                if key not in self.__class__.reserved_names:
                    setattr(new_scandata, key, value)
        return new_scandata


# %% OLD SCANDATA DEFINITION
# default way of storing each csv file's data for active use
# "scaninfo" is a dict containing scan parameters, e.g. "Voltage: 5"
# note "dataseries" and "fitdata" are PLURAL - a tuple of entries is
# expected, and fields[i], dataseries[i], and fitdata[i] are correlated
#ScanData = namedtuple("ScanData", ["fields",
#                                   "scaninfo_list",
#                                   "dataseries_list",
#                                   "error_dataseries_list",
#                                   "fitdata_list"])
# older:
# ScanData = namedtuple("ScanData", ["filepath", "scaninfo", "fields",
#                                    "dataseries", "fitdata"])


# %% Helper fcns used by DataSeries equality operator:
## NEEDS DESCRIPTION, TEST, SPHINX DOCUMENTATION
#def approx_equal(x, y, tol):
#    return abs(x - y) <= 0.5*(abs(x) + abs(y))*tol
#
#
## NEEDS DESCRIPTION, TEST, SPHINX DOCUMENTATION
#def approx_equal_lists(x_list, y_list, tol):
#    x_list = list(x_list)
#    y_list = list(y_list)
#    return all([approx_equal(x, y, tol) for (x, y) in zip(x_list, y_list)])


# %% NEEDS TEST, SPHINX DOCUMENTATION
#class DataSeries(Sequence):
#    """
#    Defines an immutable data structure that contains two correlated
#    iterables, xvals and yvals. When retreiving them with xvals() or
#    yvals(), "unsorted=True" can be used to retreive them in the order
#    provided upon construction. Adding or subtracting DataSeries
#    instances ensures proper matching of yvals corresponding to the
#    same xval. Scalar values and lists of matching length may also be
#    added, subtracted, multiplied, and divided. In the latter case,
#    operations are done element-by-element in sorted order.
#
#    Calling this object is equivalent to using its data() method
#    """
#    def __init__(self, xvals_or_datatuples, yvals=None):
#        """
#        Initializes a data structure that correlates a series of xvals
#        with a series of yvals, and when adding or subtracting ensures
#        they are properly paired up xvalwise regardless of order of
#        data points. All values should be numerical.
#
#        method 1:
#        param 1: xvals: iterable of xvals
#        param 2: yvals: iterable of yvals
#            e.g. DataSeries([1,2,3], [4,5,6])
#
#        method 2:
#        param 1: datatuples: iterable of 2-tuples in form (x, y)
#            e.g. DataSeries(data), where data = [(x1, y1),(x2, y2), ...]
#                 DataSeries(dataseries(unsorted=True))
#                 DataSeries((t,f(t)) for t in xvals), where f(t) is unary
#        """
#        # NOTE: unpack operation forces traversal of RHS. Important below...
#        xvals_or_datatuples = list(xvals_or_datatuples)
#        if len(xvals_or_datatuples) == 0:
#            raise ValueError("DataSeries cannot be initialized with no data")
#        if yvals is None:  # handle list of data points
#            try:
#                xvals, yvals = zip(*xvals_or_datatuples)  # outputs tuples
#            except ValueError:
#                raise ValueError("if DataSeries is initialized with a list " +
#                                 "of data elements, all elements must be " +
#                                 "in form (xval, yval)")
#        else:
#            xvals = tuple(xvals_or_datatuples)
#            yvals = tuple(yvals)
#            if len(xvals) != len(yvals):
#                raise ValueError("uneven number of x and y values given " +
#                                 "to DataSeries constructor")
#        try:  # ensure numerical values
#            xvals = np.array([float(x) for x in xvals])
#            yvals = np.array([float(y) for y in yvals])
#        except ValueError:
#            raise ValueError("all values in DataSeries must be numeric!")
#        # sort xy pairs by x and store sorted lists:
#        self.map_to_sorted = xvals.argsort()
#        self.map_to_unsorted = self.map_to_sorted.argsort()
#        self._xvals = np.array(xvals[self.map_to_sorted])
#        self._yvals = np.array(yvals[self.map_to_sorted])
#        # set all immutable!
#        self.map_to_sorted.setflags(write=False)
#        self.map_to_unsorted.setflags(write=False)
#        self._xvals.setflags(write=False)
#        self._yvals.setflags(write=False)
#
#    def xvals(self, *, unsorted=False, raw=False):
#        """
#        Return this object's xvals as a numpy array, optionally in the
#        original order used when this DataSeries was created.
#        
#        Optional keyword params:
#        unsorted (default: False) - if True, return in original ordering
#        raw (default: False) - (sets all of the above keywords, legacy keyword)
#        """
#        if raw:
#            unsorted = True
#        xvals, _ = self.data(unsorted=unsorted)
#        return xvals
#
#    def yvals(self, *, unsorted=False, raw=False):
#        """
#        Return this object's yvals as a numpy array, optionally in the
#        original order used when this DataSeries was created.
#        
#        Optional keyword params:
#        unsorted (default: False) - if True, return in original ordering
#        raw (default: False) - (sets all of the above keywords, legacy keyword)
#        """
#        if raw:
#            unsorted = True
#        _, yvals = self.data(unsorted=unsorted)
#        return yvals
#
#    def datatuples(self, *, unsorted=False, raw=False):
#        """
#        Return this object's data as a numpy array of (x, y) pairs, optionally
#        in the original order used when this DataSeries was created.
#        
#        Optional keyword params:
#        unsorted (default: False) - if True, return in original ordering
#        raw (default: False) - (sets all of the above keywords, legacy keyword)
#        """
#        if raw:
#            unsorted = True
#        return np.vstack(self.data(unsorted=unsorted)).transpose()
#
#    def data(self, *, unsorted=False, raw=False):
#        """
#        Return this object's data as two numpy arrays of xvals and yvals,
#        optionally in the original order used when this DataSeries was created.
#        Format is (xvals, yvals), identical to (xvals(), yvals())
#        Consider using "pyplot.plot(*dataseries.datalists())" for plots
#        
#        Optional keyword params:
#        unsorted (default: False) - if True, return in original ordering
#        raw (default: False) - (sets all of the above keywords, legacy keyword)
#        """
#        if raw:
#            unsorted = True
#        if unsorted:
#            xvals = self._xvals[self.map_to_unsorted]
#            yvals = self._yvals[self.map_to_unsorted]
#            return xvals, yvals
#        else:
#            return self._xvals, self._yvals
#
#    def copy(self):
#        """
#        Simple copy constructor making a deep copy of another DataSeries
#        object. Takes a single DataSeries instance as a parameter.
#        """
#        return self.__class__(*self.data(raw=True))
#
#    def copy_subset(self, index_list):
#        index_list = list(index_list)
#        if any(ind > len(self) for ind in index_list):
#            print("Warning: DataSeries.copy_subset(...) index_list " +
#                  "parameter contains out of bounds indices, ignoring.")
#        return self.__class__(
#            [xypair for ind, xypair in enumerate(self.datatuples(raw=True))
#             if ind in index_list])
#
#    def __call__(self, *args, **kwargs):
#        return self.data(*args, **kwargs)
#
#    def __iter__(self):
#        for index in range(len(self)):
#            yield self[index]
#
#    def __getitem__(self, indices):
#        if isinstance(indices, slice):
#            if indices.start is not None:
#                start = indices.start
#            else:
#                start = 0
#            if indices.stop is not None:
#                stop = indices.stop
#            else:
#                stop = len(self)
#            if indices.step is not None:
#                step = indices.step
#            else:
#                step = 1
#            return self.copy_subset(range(start, stop, step))
#        else:
#            if isinstance(indices, Iterable):
#                raise IndexError('DataSeries.__getitem__: ' +
#                                 'Index cannot be iterable.')
#            if indices >= len(self):
#                raise IndexError('DataSeries.__getitem__: Index out of range')
#            return self._xvals[indices], self._yvals[indices]
#
#    def __len__(self):
#        return len(self._xvals)
#
#    def __repr__(self):
#        xvals, yvals = self.data(raw=True)
#        outputstr = "DataSeries({0}, {1})".format(repr(xvals),
#                                                  repr(yvals))
#        return outputstr
#
#    def __str__(self):
#        outputstr = ''
#        for index, (xval, yval) in enumerate(zip(self._xvals, self._yvals)):
#            if index > 0:
#                outputstr += "\n"
#            outputstr += "(x = {xval}, y = {yval})".format(xval=xval,
#                                                           yval=yval)
#        return outputstr
#
#    def general_numeric_fcn(self, other, function, fcn_verb):
#        fcn_verb = str(fcn_verb)
#        map_to_unsorted = self.map_to_unsorted
#        if isinstance(other, self.__class__):  # is a DataSeries
#            try:
#                otherxvals = other._xvals
#                otheryvals = other._yvals
#                othermap = other.map_to_unsorted
#                if not np.array_equal(self._xvals, otherxvals):
#                    raise TypeError('DataSeries instances have different ' +
#                                    'xvalues, so failed to ' + fcn_verb)
#                if not np.array_equal(self.map_to_unsorted, othermap):
#                    print('Warning: combining two DataSeries instances ' +
#                          'with different sort orders, resulting sort ' +
#                          'order uncertain.')
#                    if np.array_equal(map_to_unsorted, np.arange(len(self))):
#                        map_to_unsorted = othermap  # if self-sorted, other map
#                newyvals = np.array(function(self._yvals, otheryvals))
#            except AttributeError:
#                print('Warning: failed to combine DataSeries with apparent ' +
#                      'DataSeries, treating latter as generic list/array')
#                newyvals = function(self.yvals(unsorted=False), other)
#        else:  # not a DataSeries: let numpy handle it
#            newyvals = function(self.yvals(unsorted=False), other)
#        newxvals = self._xvals[map_to_unsorted]
#        newyvals = newyvals[map_to_unsorted]
#        return self.__class__(newxvals, newyvals)
#
#    def __add__(self, other):
#        fcn_verb = "add"
#        def fcn(x, y):
#            return x + y
#        return self.general_numeric_fcn(other, fcn, fcn_verb)
#
#    def __sub__(self, other):
#        return self.__add__(-1*other)
#
#    def __mul__(self, other):
#        fcn_verb = "multiply"
#        def fcn(x, y):
#            return x * y
#        return self.general_numeric_fcn(other, fcn, fcn_verb)
#
#    def __div__(self, other):
#        fcn_verb = "divide"
#        def fcn(x, y):
#            return x / y
#        return self.general_numeric_fcn(other, fcn, fcn_verb)
#
#    def __truediv__(self, other):
#        fcn_verb = "divide"
#        def fcn(x, y):
#            return x / y
#        return self.general_numeric_fcn(other, fcn, fcn_verb)
#
#    def __pow__(self, other):
#        fcn_verb = "raise"
#        def fcn(x, y):
#            return x**y
#        return self.general_numeric_fcn(other, fcn, fcn_verb)
#
#    def __radd__(self, other):
#        return self.__add__(other)
#
#    def __rsub__(self, other):
#        return -1*self.__sub__(other)
#
#    def __rmul__(self, other):
#        return self.__mul__(other)
#
#    def __neg__(self):
#        return -1*self
#
#    def __pos__(self):
#        return self
#
#    def __abs__(self):
#        xvals, yvals = self.data(raw=True)
#        return self.__class__(xvals, abs(yvals))
#
#    def __eq__(self, other):
#        try:
#            return all(
#                [np.array_equal(self._xvals, other._xvals),
#                 np.array_equal(self._yvals, other._yvals)])
#        except AttributeError:  # not a DataSeries
#            return False
