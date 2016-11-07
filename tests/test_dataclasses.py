#! python
# -*- coding: utf-8 -*-
"""
Designed to test the DataSeries class.

Run from command line (recommended!) by inputting one of the following:
py.test                                                [this runs all tests]
python setup.py test                                   [this runs all tests]
py.test "tests/[this_test_name].py"                    [this test only]
python setup.py test -a "tests/[this_test_name].py"    [this test only]

Created on Mon Jan  4 03:45:43 2016

@author: Michael
"""

import numpy as np
import pytest
from experimentdataanalysis.analysis.dataclasses import DataSeries


# %% NON FIXTURE HELPER FUNCTIONS


# %% TEST FIXTURES


# JUST HERE FOR SYNTAX REFERENCE
#@pytest.fixture()
#def requestfactory_testcounterfcn():
#    def factory():
#        def testcounterfcn(actor):
#            actor.properties.setdefault('testcounter', 0)
#            actor.properties['testcounter'] += 1
#        return request.Request(testcounterfcn)
#    return factory


# %% TESTS
# don't forget to make column default names proper attribute names!!!!!
# tests needed for:
# -less than two fields
# -repeating field names
# -invalid field names, both reserved names and ones with spaces, etc.
# -can set x_field_name, y_field_name, even if original names invalid
# -rejects field arrays of different lengths
# -accepts various kinds of field array inputs, e.g. lists, tuples
# -properly forwards info keys (if feature not dropped)
# -properly handles info key name conflicts
# -proper functioning of descriptors and get_ equivalents w/[out] _error fields

#def test_xvalsareimmutable():
#    xvals = [1, 2, 3]
#    yvals = [10, 20, 30]
#    data = DataSeries(zip(xvals, yvals))
#    with pytest.raises(ValueError):
#        data._xvals[0] = 0
#    with pytest.raises(ValueError):
#        data._yvals[0] = 0
#    xvals[0] = 0
#    yvals[0] = 0
#    assert not np.array_equal(data[0], (0, 0))
#
#
#def test_sortingflagworkscorrectly():
#    data = DataSeries(zip([3, 1, 2], [30, 10, 20]))
#    assert np.array_equal(  data.yvals(),              [10, 20, 30])
#    assert np.array_equal(  data.yvals(raw=True),      [30, 10, 20])
#    assert np.array_equal(  data.yvals(unsorted=True), [30, 10, 20])
#
#
#def test_usesmorecomplexorderingonaddition():
#    series1 = DataSeries(zip([1, 2, 3], [10, 20, 30]))
#    series2 = DataSeries(zip([3, 1, 2], [10, 20, 30]))
#    series3 = series1 + series2
#    series4 = series2 + series1
#    assert np.array_equal(  series3.xvals(unsorted=True), [3, 1, 2])
#    assert np.array_equal(  series4.xvals(unsorted=True), [3, 1, 2])
#
#
#def test_basicmathoperationswork():
#    series1 = DataSeries(zip([1, 2, 3], [10, 20, 30]))
#    series2 = DataSeries(zip([3, 2, 1], [10, 20, 30]))
#    assert np.array_equal(  (series1 + series2).yvals(),   [40, 40, 40])
#    assert np.array_equal(  (series1 - series2).yvals(),   [-20, 0, 20])
#    assert np.array_equal(  (series1 + [1, 2, 3]).yvals(), [11, 22, 33])
#    assert np.array_equal(  (series1 - 10).yvals(),        [0, 10, 20])
#    assert np.array_equal(  (series1 * 10).yvals(),        [100, 200, 300])
#    assert np.array_equal(  (series1 * 0 + 1).yvals(),     [1, 1, 1])


# JUST HERE FOR SYNTAX REFERENCE
#def test_reactproperlycallsexternalfcn():
#    actor = reactactor()
#    request_count = 10
#    execution_count = 0
#
#    def testcounterfcn(actor):
#        nonlocal execution_count
#        execution_count += 1
#        actor.properties.setdefault('testcounter', 0)
#        actor.properties['testcounter'] += 1
#
#    for count in range(request_count):
#        newrequest = request.Request(testcounterfcn)
#        actor.request_queue.put((1, newrequest))
#    assert actor.properties.get('testcounter', 0) == 0
#    assert execution_count == 0
#    actor.react()
#    assert actor.properties.get('testcounter', 0) == request_count
#    assert execution_count == request_count
