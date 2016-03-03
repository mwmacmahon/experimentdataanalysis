# -*- coding: utf-8 -*-
"""
Designed to test the curvefitting module.

Run by typing the following in python while in the project directory:
import pytest
pytest.main(args=['-s'])

Created on Fri Feb 26 16:23:00 2016

@author: Michael
"""

import pytest

import experimentdataanalysis.analysis.curvefitting as curvefitting
from experimentdataanalysis.analysis.dataclasses \
    import FitData, ScanData, DataSeries
import experimentdataanalysis.parsing.dataclassparsing as dcparsing


# %% NON FIXTURE HELPER FUNCTIONS


# %% TEST FIXTURES
@pytest.fixture(scope="module")
def loadscandatalist():
    filepath = ("C:\\pythonprojects\\experimentdataanalysis_project\\" +
                "tests\\representativetwocosdata")
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
