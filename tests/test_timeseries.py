# -*- coding: utf-8 -*-
"""
Designed to test the TimeSeries class.

Run by typing the following in python while in the project directory:
import pytest
pytest.main(args=['-s'])

Created on Mon Jan  4 03:45:43 2016

@author: Michael
"""

import pytest
from experimentdataanalysis.analysis.dataclasses import TimeSeries


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
def test_timesareimmutable():
    times = [1, 2, 3]
    values = [10, 20, 30]
    data = TimeSeries(zip(times, values))
    with pytest.raises(TypeError):
        data._times[0] = 0
    with pytest.raises(TypeError):
        data._values[0] = 0
    times[0] = 0
    values[0] = 0
    assert data[0] != (0, 0)


def test_filtersareimmutable():
    data = TimeSeries(zip([1, 2, 3], [10, 20, 30]),
                      excluded_intervals=[(1, 1)])
    with pytest.raises(TypeError):
        data._start_times[0] = 2
    with pytest.raises(TypeError):
        data._end_times[0] = 1


def test_filtersworkcorrectly():
    data = TimeSeries(zip([1, 2, 3], [10, 20, 30]))
    assert data.values() == [10, 20, 30]

    data = TimeSeries(data(raw=True), excluded_intervals=[(1, 1)])
    assert data.values() == [20, 30]
    data = TimeSeries(data(raw=True), excluded_intervals=[(2, 2)])
    assert data.values() == [10, 30]
    data = TimeSeries(data(raw=True), excluded_intervals=[(3, 3)])
    assert data.values() == [10, 20]
    data = TimeSeries(data(raw=True), excluded_intervals=[(1, 2)])
    assert data.values() == [30]
    data = TimeSeries(data(raw=True), excluded_intervals=[(2, 3)])
    assert data.values() == [10]
    data = TimeSeries(data(raw=True), excluded_intervals=[(1, 3)])
    assert data.values() == []
    assert data.values(unfiltered=True) == [10, 20, 30]


def test_sortingflagworkscorrectly():
    data = TimeSeries(zip([3, 1, 2], [30, 10, 20]))
    assert data.values() == [10, 20, 30]
    assert data.values(unsorted=True) == [30, 10, 20]


def test_usesmorecomplexorderingonaddition():
    series1 = TimeSeries(zip([1, 2, 3], [10, 20, 30]))
    series2 = TimeSeries(zip([3, 1, 2], [10, 20, 30]))
    series3 = series1 + series2
    series4 = series2 + series1
    assert series3.times(unsorted=True) == [3, 1, 2]
    assert series4.times(unsorted=True) == [3, 1, 2]


def test_basicmathoperationswork():
    series1 = TimeSeries(zip([1, 2, 3], [10, 20, 30]))
    series2 = TimeSeries(zip([3, 2, 1], [10, 20, 30]))
    assert (series1 + series2).values() == [40, 40, 40]
    assert (series1 - series2).values() == [-20, 0, 20]
    assert (series1 + [1, 2, 3]).values() == [11, 22, 33]
    assert (series1 - 10).values() == [0, 10, 20]
    assert (series1 * 10).values() == [100, 200, 300]
    assert (series1 * 0 + 1).values() == [1, 1, 1]


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
