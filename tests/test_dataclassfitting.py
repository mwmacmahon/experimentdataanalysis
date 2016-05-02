#! python
# -*- coding: utf-8 -*-
"""
Designed to test the dataclassfitting module.

Run by typing the following in python while in the project directory:
import pytest
pytest.main(args=['-s'])

Created on Fri Feb 26 06:44:22 2016

@author: Michael
"""

import pytest

from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
import experimentdataanalysis.parsing.dataclassparsing as dcparsing
import experimentdataanalysis.analysis.dataclassfitting as dcfitting


# %% NON FIXTURE HELPER FUNCTIONS


# %% TEST FIXTURES
@pytest.fixture(scope="module")
def loadscandatalist():
    test_dir_path = __file__[:__file__.rfind("\\")]
    filepath = (test_dir_path + "\\representativetwocosdata")
    scandatalist = list(
                    dcparsing.fetch_dir_as_unfit_scandata_iterator(filepath))
    return scandatalist


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
def test_no_multiprocessing_fit_scandata_list_twocos(loadscandatalist):
    fitfunc = dcfitting.fit_dataseries_with_two_decaying_cos
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=True,
                                             multiprocessing=False)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=False,
                                             multiprocessing=False)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000


def test_no_multiprocessing_fit_scandata_list_onecos(loadscandatalist):
    fitfunc = dcfitting.fit_dataseries_with_one_decaying_cos
    output = dcfitting.fit_scandata_iterable([loadscandatalist[1]],
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=True,
                                             multiprocessing=False)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000
    output = dcfitting.fit_scandata_iterable([loadscandatalist[1]],
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=False,
                                             multiprocessing=False)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000


def test_multiprocessing_fit_scandata_list_twocos(loadscandatalist):
    fitfunc = dcfitting.fit_dataseries_with_two_decaying_cos
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=True,
                                             multiprocessing=True)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=False,
                                             multiprocessing=True)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000


def test_multiprocessing_fit_scandata_list_onecos(loadscandatalist):
    fitfunc = dcfitting.fit_dataseries_with_one_decaying_cos
    output = dcfitting.fit_scandata_iterable([loadscandatalist[1]],
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=True,
                                             multiprocessing=True)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000
    output = dcfitting.fit_scandata_iterable([loadscandatalist[1]],
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=False,
                                             multiprocessing=True)
    assert abs(next(output).fitdata[0].fitparamstds[0]) < 100000


def test_multiprocessing_fit_no_fitfunc(loadscandatalist):
    fitfunc = None
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=True,
                                             multiprocessing=True)
    assert next(output).fitdata[0] is None
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=False,
                                             multiprocessing=True)
    assert next(output).fitdata[0] is None


def test_no_multiprocessing_fit_no_fitfunc(loadscandatalist):
    fitfunc = None
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=True,
                                             multiprocessing=False)
    assert next(output).fitdata[0] is None
    output = dcfitting.fit_scandata_iterable(loadscandatalist,
                                             dataseriesfitfunction=fitfunc,
                                             fit_drift=False,
                                             multiprocessing=False)
    assert next(output).fitdata[0] is None


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
