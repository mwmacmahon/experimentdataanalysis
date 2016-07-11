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
def test_xvalsareimmutable():
    xvals = [1, 2, 3]
    yvals = [10, 20, 30]
    data = DataSeries(zip(xvals, yvals))
    with pytest.raises(TypeError):
        data._xvals[0] = 0
    with pytest.raises(TypeError):
        data._yvals[0] = 0
    xvals[0] = 0
    yvals[0] = 0
    assert data[0] != (0, 0)


def test_filtersareimmutable():
    data = DataSeries(zip([1, 2, 3], [10, 20, 30]),
                      excluded_intervals=[(1, 1)])
    with pytest.raises(TypeError):
        data._start_xvals[0] = 2
    with pytest.raises(TypeError):
        data._end_xvals[0] = 1


def test_filtersworkcorrectly():
    data = DataSeries(zip([1, 2, 3], [10, 20, 30]))
    assert data.yvals() == [10, 20, 30]

    data = DataSeries(data(raw=True), excluded_intervals=[(1, 1)])
    assert data.yvals() == [20, 30]
    data = DataSeries(data(raw=True), excluded_intervals=[(2, 2)])
    assert data.yvals() == [10, 30]
    data = DataSeries(data(raw=True), excluded_intervals=[(3, 3)])
    assert data.yvals() == [10, 20]
    data = DataSeries(data(raw=True), excluded_intervals=[(1, 2)])
    assert data.yvals() == [30]
    data = DataSeries(data(raw=True), excluded_intervals=[(2, 3)])
    assert data.yvals() == [10]
    data = DataSeries(data(raw=True), excluded_intervals=[(1, 3)])
    assert data.yvals() == []
    assert data.yvals(unfiltered=True) == [10, 20, 30]


def test_sortingflagworkscorrectly():
    data = DataSeries(zip([3, 1, 2], [30, 10, 20]))
    assert data.yvals() == [10, 20, 30]
    assert data.yvals(unsorted=True) == [30, 10, 20]


def test_usesmorecomplexorderingonaddition():
    series1 = DataSeries(zip([1, 2, 3], [10, 20, 30]))
    series2 = DataSeries(zip([3, 1, 2], [10, 20, 30]))
    series3 = series1 + series2
    series4 = series2 + series1
    assert series3.xvals(unsorted=True) == [3, 1, 2]
    assert series4.xvals(unsorted=True) == [3, 1, 2]


def test_basicmathoperationswork():
    series1 = DataSeries(zip([1, 2, 3], [10, 20, 30]))
    series2 = DataSeries(zip([3, 2, 1], [10, 20, 30]))
    assert (series1 + series2).yvals() == [40, 40, 40]
    assert (series1 - series2).yvals() == [-20, 0, 20]
    assert (series1 + [1, 2, 3]).yvals() == [11, 22, 33]
    assert (series1 - 10).yvals() == [0, 10, 20]
    assert (series1 * 10).yvals() == [100, 200, 300]
    assert (series1 * 0 + 1).yvals() == [1, 1, 1]


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
