# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:02:49 2016

@author: vsih-lab
"""
from experimentdataanalysis.analysis.generalutilities \
    import multiprocessable_map
from testscript_temp import \
    simple_fcn, simple_fcn_two_args, simple_fcn_with_nested_partial_fcn, \
    fcn_that_curries_a_given_fcn, fcn_with_nested_fucntional_fcn


# %%
if __name__ == "__main__":
    input_args_list = [[simple_fcn_two_args, 1],
                       [simple_fcn_two_args, 2]]
    output = multiprocessable_map(fcn_with_nested_fucntional_fcn,
                                  input_args_list, multiprocessing=False)
    print(list(output))
    output = multiprocessable_map(fcn_with_nested_fucntional_fcn,
                                  input_args_list, multiprocessing=True)
    print(list(output))
