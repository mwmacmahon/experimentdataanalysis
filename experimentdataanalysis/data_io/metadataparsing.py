# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:07:53 2016

@author: Michael
"""

import os.path
import time

import numpy as np


# %%
def analyze_scan_filepath(filepath, scaninfo=None, keywordlists=None):
    """
    Scans the filepath (including filename) of a scan, and returns
    a dict containing all the info and tags it can match.

    Custom keyword lists can be passed by the keywordlists keyword
    argument, where None for a keywordlist type means use defaults
    """
    if scaninfo is None:
        scaninfo = {}  # NECESSARY, if defined in header, shared between calls!
    scaninfo['Filepath'] = filepath
    try:
        scaninfo['File Last Modified'] = time.ctime(os.path.getmtime(filepath))
    except OSError:
        pass
    except FileNotFoundError:
        pass
    if keywordlists is not None:
        this_element_keyword_list, \
            next_element_keyword_list, \
            inside_this_element_keyword_list = keywordlists
    else:
        # DEFAULT SEARCH TERMS AND SEARCH RULES:
        # 1. If first string found, register second string as
        #    tag containing third string/value
        #        e.g. ("warmup", "Warmup?", "Yes")
        #             result: "..._warmup_..." -> {"Warmup?": "Yes"}
        this_element_keyword_list = [
            ("ExampleTargetKey", "ExampleStorageKey", "ExampleStoredValue"),
        ]
        # 2. Grab next element(s) if this one CONTAINS first string,
        #    tag next element(s) as second string(s)
        #        e.g. ("Ind", "FastScanIndex")
        #             result: "..._Ind_3_..." -> {"FastScanIndex": 3}
        #        e.g. ("2DScan", ["SecondScanType", "FirstScanType"])
        #             result: "..._2Dscan_MirrorY_MirrorZ_..."
        #                         -> {"SecondScanType": "MirrorY",
        #                             "FirstScanType": "MirrorZ"}
        next_element_keyword_list = [
            ("ExampleTargetKey1", "ExampleStorageKey1"),
            ("ExampleTargetKey2", ["ExampleStorageKey2",
                                   "ExampleStorageKey3"]),
        ]
        # 3. Grab this element if it CONTAINS first string,
        #    tag remainder as second string
        #        e.g. ("K", "SetTemperature")
        #             result: "..._30K_..." -> {"SetTemperature": 30}
        inside_this_element_keyword_list = [
            ("ExampleTargetKey", "ExampleStorageKey"),
        ]

    # get rid of idiosyncratic delimiters by swapping with _
    filepath = filepath.replace("\\", "__")
    filepath = filepath.replace(" ", "_")
    filepath = filepath.replace(".dat", "_")
    next_element_tags = []
    for element in filepath.split("_"):
        if len(next_element_tags) > 0:
            try:
                value = float(element)
            except ValueError:  # if not a numeric value:
                # ignore trailing keywords, e.g. units
                value = element
                if len(inside_this_element_keyword_list) > 0:
                    for matchstr, _ in inside_this_element_keyword_list:
                        if element.endswith(matchstr):
                            value = element.replace(matchstr, "")
                    try:
                        value = float(value)
                    except ValueError:  # still non-numeric
                        value = element
            scaninfo[next_element_tags.pop(0)] = value
        elif len(next_element_keyword_list) > 0:
            for matchstr, tags in next_element_keyword_list:
                if element == matchstr:
                    if isinstance(tags, str):  # only one string
                        next_element_tags = [tags]
                    else:
                        next_element_tags = list(tags)  # copy old list!
        for matchstr, tag, value in this_element_keyword_list:
            if element == matchstr:
                scaninfo[tag] = value
        # a little trickier, must use _best_ keyword match or none
        # e.g. "0Vcm": '0' units 'Vcm', NOT ALSO '0cm' units 'V'
        if len(inside_this_element_keyword_list) > 0:
            inside_this_element_keys, inside_this_element_tags = \
                zip(*inside_this_element_keyword_list)
            keymatches = [key in element for key in inside_this_element_keys]
            keylengths = [len(key) for key in inside_this_element_keys]
            best_match_index = \
                np.argmax([matched * length  # use bool -> 0, 1 conversion
                           for matched, length in zip(keymatches, keylengths)])
            matchstr, tag = inside_this_element_keyword_list[best_match_index]
            if matchstr in element:  # avoid likely case of no match
                value = element.replace(matchstr, "")
                try:
                    value = float(value)
                    scaninfo[tag] = value
                except ValueError:
                    pass  # by this rule, only take numerical values
    return scaninfo


# %%
def analyze_string_for_dict_pairs(infostr, scaninfo=None):
    """
    Currently looks for key, value pairs in form "key: value" in the
    provided strings and adds them to the dict given (or otherwise
    creates a new dict).
    """
    if scaninfo is None:
        scaninfo = {}  # NECESSARY, if defined in header, shared between calls!
    strrows = infostr.splitlines()
    for row in strrows:
        key, value = "", ""
        try:
            key, value = row.split(":")
        except ValueError:
            pass
        if key:
            key = key.strip()
            value = value.strip()
            scaninfo[key] = value
    return scaninfo
