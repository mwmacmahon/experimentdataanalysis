# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:20:02 2016

@author: vsih-lab
"""

from functools import partial


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
