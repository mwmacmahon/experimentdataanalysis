#! python
# -*- coding: utf-8 -*-
"""
Designed to test the dataclassfitting module.

Run from command line (recommended!) by inputting one of the following:
py.test                                                [this runs all tests]
python setup.py test                                   [this runs all tests]
py.test "tests/[this_test_name].py"                    [this test only]
python setup.py test -a "tests/[this_test_name].py"    [this test only]

Created on Fri Feb 26 06:44:22 2016

@author: Michael
"""

import numpy as np
import pytest

from experimentdataanalysis.analysis.multidataseriesprocessing \
    import scandata_iterable_fit
import experimentdataanalysis.analysis.fitfunctions as fitfcns
import experimentdataanalysis.parsing.dataseriesparsing as dsparsing


# %% NON FIXTURE HELPER FUNCTIONS


# %% TEST FIXTURES
@pytest.fixture(scope="module")
def loadscandatalist():
    test_dir_path = __file__[:__file__.rfind("\\")]
    filepath = (test_dir_path + "\\representative3ddata")
    scandatalist = list(
                    dsparsing.fetch_dir_as_unfit_scandata_iterator(filepath))
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
def test_fit_scandata_list(loadscandatalist):
    output_list = scandata_iterable_fit(
                    loadscandatalist,
                    field_index=0,
                    fitfunction=fitfcns.fitfcn_two_exp_sin_decay,
                    free_params = [True, True, True, True,
                                   True, True, True, True],
                    initial_params = [0.002, 50, 0.002, 1000, 810, 0, 0, 0],
                    param_bounds = [(-1, 1), (1, 200), (-1, 1), (10, 1e6),
                                    (810, 830), (0, 2*np.pi),
                                    (-1e-6, 1e-6), (-0.01, 0.01)],
                    multiprocessing=False)
    assert abs(output_list[0].fitdata_list[0].fitparamstds[0]) < 100000
    output_list = scandata_iterable_fit(
                    loadscandatalist,
                    field_index=0,
                    fitfunction=fitfcns.fitfcn_two_exp_sin_decay,
                    free_params = [True, True, True, True,
                                   True, True, True, True],
                    initial_params = [0.002, 50, 0.002, 1000, 810, 0, 0, 0],
                    param_bounds = [(-1, 1), (1, 200), (-1, 1), (10, 1e6),
                                    (810, 830), (0, 2*np.pi),
                                    (-1e-6, 1e-6), (-0.01, 0.01)],
                    multiprocessing=True)
    assert abs(output_list[0].fitdata_list[0].fitparamstds[0]) < 100000


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
