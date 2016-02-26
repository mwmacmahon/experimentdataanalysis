# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:22:19 2016

@author: Michael
"""

from concurrent.futures import ProcessPoolExecutor
import time


# %%
def multiprocessable_map(processfunction, input_iterable,
                         multiprocessing=False):
    """
    multiprocessing = False:
    identical to map(processfunction, input_iterable)
    ---
    multiprocessing = True:
    Generic ProcessPoolExecutor loop. Supports any arbitrary
    function, whether processes or threads or faster may depend
    on the function and computer used.
    WARNING: pickle fails unless the input/output object classes
        are defined via IMPORT statement. Cannot define in module!
    WARNING: do not call even indirectly from interpreter, run
        script with main() method!
    WARNING: numpy's curve_map is not threadsafe in a way that may
        cause issues. Use multiprocessing instead (faster anyway)
    WARNING: pickle has to send the code to be run, and it looks
        like it cannot handle sending functions that dynamically
        create local subfunctions or are local subfunctions
        themselves. Avoid that, I guess.
            e.g. sending this instead of processfunction failed:
                def safer_map_function(timeseries):
                    try:
                        return processfunction(timeseries)
                        except <stuff>...

    Returns iterable pointing to the outputs of the processfunction
    acting on each input value.

    Positional arguments:
        processfunction -- should take a timeseries and return in form
        input_iterable -- iterable pointing to TimeSeries
                               objects containing data to map.
                               should not be reused elsewhere.
    """
    if multiprocessing:
        start_processing_time = time.time()
        input_iterable = list(input_iterable)
#        processfunction_iterable = [processfunction]*len(input_iterable)
#        with ThreadPoolExecutor(max_workers=None) as executor:
        with ProcessPoolExecutor(max_workers=None) as executor:
            output_iter = executor.map(processfunction, input_iterable,
                                       timeout=30, chunksize=1)
        elapsed_processing_time = time.time() - start_processing_time
        print('{} seconds elapsed during multiprocessing'.format(
                                                    elapsed_processing_time))
        return output_iter
    else:
        return map(processfunction, input_iterable)


# %%
def mltprctstfcn1(x=0):
    return x+2
def mltprctstfcn2(x=0, stuff=False):
    if stuff:
        return x+2
    else:
        return (x, x+2)
def testprintoutput(a):
    print(next(a))
    print(next(a))
    print(next(a))
    try:
        print(next(a))
    except StopIteration:
        print('iterable exhausted')
if __name__ == "__main__":
    a = multiprocessable_map(mltprctstfcn1, [1,2,3], multiprocessing=False)
    testprintoutput(a)
    a = multiprocessable_map(mltprctstfcn2, [1,2,3], multiprocessing=False)
    testprintoutput(a)
    a = multiprocessable_map(mltprctstfcn1, [1,2,3], multiprocessing=True)
    testprintoutput(a)
    a = multiprocessable_map(mltprctstfcn2, [1,2,3], multiprocessing=True)
    testprintoutput(a)
    # result: OK to have multiple args in function given
    # to multiprocessable_map as long as rest have defaults
