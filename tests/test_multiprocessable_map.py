#! python
# -*- coding: utf-8 -*-
"""
Designed to test the multiprocessable_map function in generalutilities.py

Run from command line (recommended!) by inputting one of the following:
py.test                                                [this runs all tests]
python setup.py test                                   [this runs all tests]
py.test "tests/[this_test_name].py"                    [this test only]
python setup.py test -a "tests/[this_test_name].py"    [this test only]

WARNING: takes over twice as long as executing the same code manually...

Created on Fri Feb 26 05:46:52 2016

@author: Michael
"""

from functools import partial

import pytest

from experimentdataanalysis.analysis.generalutilities \
    import multiprocessable_map


# %% NON FIXTURE HELPER FUNCTIONS
def simple_fcn(x):
    return x*2


def simple_fcn_two_args(x=1, y=1):
    return x*y


def simple_fcn_with_nested_partial_fcn(x):
    partialfcn = partial(simple_fcn_two_args, y=10)
    return partialfcn(x)  # should be x*10


def fcn_that_curries_a_given_fcn(fcn):
    partialfcn = partial(fcn, 10)
    return partialfcn(2)


def fcn_with_nested_fucntional_fcn(fcn, x):
    def curry(fcn, y):
        return fcn(10, y)
    return curry(fcn, x)


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
def test_can_mmap_simple_helper_fcn():
    input_args_list = [[1], [2], [3]]
    output = multiprocessable_map(simple_fcn, input_args_list,
                                  multiprocessing=False)
    assert list(output) == [2, 4, 6]
    output = multiprocessable_map(simple_fcn, input_args_list,
                                  multiprocessing=True)
    assert list(output) == [2, 4, 6]


def test_can_mmap_variable_num_args():
    input_args_list = [[1], [5], [3, 5]]
    output = multiprocessable_map(simple_fcn_two_args, input_args_list,
                                  multiprocessing=False)
    assert list(output) == [1, 5, 15]
    output = multiprocessable_map(simple_fcn_two_args, input_args_list,
                                  multiprocessing=True)
    assert list(output) == [1, 5, 15]


def test_can_mmap_fcn_with_partial_subfcn():
    input_args_list = [[1], [2], [3]]
    output = multiprocessable_map(simple_fcn_with_nested_partial_fcn,
                                  input_args_list, multiprocessing=False)
    assert list(output) == [10, 20, 30]
    output = multiprocessable_map(simple_fcn_with_nested_partial_fcn,
                                  input_args_list, multiprocessing=True)
    assert list(output) == [10, 20, 30]


def test_can_mmap_fcn_that_curries_given_fcn():
    input_args_list = [[simple_fcn_two_args], [simple_fcn_two_args]]
    output = multiprocessable_map(fcn_that_curries_a_given_fcn,
                                  input_args_list, multiprocessing=False)
    assert list(output) == [20, 20]
    output = multiprocessable_map(fcn_that_curries_a_given_fcn,
                                  input_args_list, multiprocessing=True)
    assert list(output) == [20, 20]


def test_can_mmap_fcn_containing_functional_subfcn():
    input_args_list = [[simple_fcn_two_args, 1],
                       [simple_fcn_two_args, 2]]
    output = multiprocessable_map(fcn_with_nested_fucntional_fcn,
                                  input_args_list, multiprocessing=False)
    assert list(output) == [10, 20]
    output = multiprocessable_map(fcn_with_nested_fucntional_fcn,
                                  input_args_list, multiprocessing=True)
    assert list(output) == [10, 20]


# NOPE, hangs forever
#def test_can_mmap_simple_nested_fcn():
#    def simple_nested_fcn(x):
#        return x+1
#    output = multiprocessable_map(simple_nested_fcn, [1, 2, 3],
#                                  multiprocessing=True)
#    assert list(output) == [2, 3, 4]


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
