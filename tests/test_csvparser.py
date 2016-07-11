#! python
# -*- coding: utf-8 -*-
"""
Designed to test the csvparser module.

Run from command line (recommended!) by inputting one of the following:
py.test                                                [this runs all tests]
python setup.py test                                   [this runs all tests]
py.test "tests/[this_test_name].py"                    [this test only]
python setup.py test -a "tests/[this_test_name].py"    [this test only]

Created on Fri Mar  4 18:46:09 2016

@author: Michael
"""

import pytest

import experimentdataanalysis.parsing.csvparser as csvparser

# %% NON FIXTURE HELPER FUNCTIONS


# %% TEST FIXTURES
@pytest.fixture(scope="module")
def test_csvs_dir_path():
    test_dir_path = __file__[:__file__.rfind("\\")]
    return test_dir_path + "\\csvs_to_test\\"

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
def test_extrarowsatend(test_csvs_dir_path):
    filepath, rawdata = csvparser.parse_csv(
        test_csvs_dir_path + "extra_tabs.dat")
    fields, columndata = rawdata
    assert len(fields) == 3
    assert len(columndata) == 3
    assert len(columndata[0]) == len(columndata[1])
    assert len(columndata[1]) == len(columndata[2])


def test_headered_csv(test_csvs_dir_path):
    filepath, rawdata = csvparser.parse_csv(
        test_csvs_dir_path + "headers.dat")
    fields, columndata = rawdata
    assert fields == tuple(["col1", "col2", "col3"])
    assert columndata[0] == tuple([1.0, 2.0, 3.0, 4.0, 5.0])
    assert columndata[1] == tuple([10, 20, 30, 40, 50])
    assert columndata[2] == tuple([1000, 2000, 3000, 4000, 5000])


def test_nonheadered_csv(test_csvs_dir_path):
    filepath, rawdata = csvparser.parse_csv(
        test_csvs_dir_path + "noheaders.dat")
    fields, columndata = rawdata
    assert fields == tuple(["Column 1", "Column 2", "Column 3"])
    assert columndata[0] == tuple([1.0, 2.0, 3.0, 4.0, 5.0])
    assert columndata[1] == tuple([10, 20, 30, 40, 50])
    assert columndata[2] == tuple([1000, 2000, 3000, 4000, 5000])


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
