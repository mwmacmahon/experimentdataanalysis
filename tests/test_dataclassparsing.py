# -*- coding: utf-8 -*-
"""
Designed to test the dataclassparsing module.

Run by typing the following in python while in the project directory:
import pytest
pytest.main(args=['-s'])

Created on Thu Feb 25 14:21:43 2016

@author: Michael
"""

import pytest

import experimentdataanalysis.parsing.csvparser as csvparser
import experimentdataanalysis.parsing.dataclassparsing as dcparsing


# %% NON FIXTURE HELPER FUNCTIONS
def check_scandata_veracity(scandata, scancoord_index):
    filepathinfo = dcparsing.analyze_scan_filepath(scandata.filepath)
    assert scandata.scaninfo['Voltage'] == filepathinfo['Voltage']
    scancoord_dataseries = scandata.dataseries[scancoord_index]
    for xval, yval in scancoord_dataseries.datatuples():
        assert abs(xval - yval) == 0  # should default to scancoord


# %% TEST FIXTURES
@pytest.fixture(scope="module")
def loadcsvdir():
    filepath = ("C:\\pythonprojects\\experimentdataanalysis_project\\" +
                "tests\\representativetwocosdata")
    csvdirdata = csvparser.parse_csv_directory(filepath, delimiter='\t')
    filepath_list, rawcsvdata_list = zip(*csvdirdata)
    return filepath_list, rawcsvdata_list


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
def test_extract_scandata_iter_from_filepath(loadcsvdir):
    filepath = ("C:\\pythonprojects\\experimentdataanalysis_project\\" +
                "tests\\representativetwocosdata")
    scandatalist = dcparsing.fetch_dir_as_unfit_scandata_iterator(filepath)
    for scandata in scandatalist:
        check_scandata_veracity(scandata, scancoord_index=1)


def test_extract_scandata_from_filepath(loadcsvdir):
    filepath_list, rawcsvdata_list = loadcsvdir
    for filepath in filepath_list:
        scandata = dcparsing.fetch_csv_as_unfit_scandata(
                                    filepath, 'scancoord')
        check_scandata_veracity(scandata, scancoord_index=0)


def test_tryinvalidattribute_singlefile(loadcsvdir):
    filepath_list, rawcsvdata_list = loadcsvdir
    for filepath in filepath_list:
        scandata = dcparsing.fetch_csv_as_unfit_scandata(
                                    filepath, 'invalid')
        check_scandata_veracity(scandata, scancoord_index=0)


def test_filenameparsing():
    filepath = ("C:\\pythonprojects\\experimentdataanalysis_project\\" +
                "tests\\csvs_to_test\\" +
                "Experiment_Channel3_033XT-B11" +
                "_819.0nm_30K_2Dscan_Voltage_DelayTime_run2\\" +
                "Ind_1_DelayTime -400_to_6100 Voltage 2.5x.dat")
    scaninfo = dcparsing.analyze_scan_filepath(filepath)
    print(scaninfo)
    assert scaninfo["Voltage"] == 2.5
    assert scaninfo["Channel 3"] == True
    assert scaninfo["Wavelength"] == 819.0
    assert scaninfo["SetTemperature"] == 30.0
    assert scaninfo["FastScanIndex"] == 1
    assert scaninfo["FastScanType"] == "DelayTime"
    assert scaninfo["MiddleScanCoord"] == 2.5
    assert scaninfo["MiddleScanType"] == "Voltage"


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
