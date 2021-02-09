#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
import dmsky.targets

def test_create_target():
    # Build the library of pre-defined rosters
    target_library = dmsky.targets.TargetLibrary()

    # Ursa Minor using implicit defaults
    umi = target_library.create_target('ursa_minor')
    assert umi.ra == 227.2854
    assert umi.dec == 67.2225
    assert umi.distance == 76.0

    # Check that implicit and explicit defaults the same
    umi2015 = target_library.create_target('ursa_minor',default='default')
    assert umi2015.ra == umi.ra
    assert umi2015.dec == umi.dec
    assert umi2015.distance == umi.distance

    umi2020 = target_library.create_target('ursa_minor',default='default2020')
    assert umi2020.ra == 227.242
    assert umi2020.dec == 67.2221
    assert umi2020.distance == 76.2

    #print(umi2020)

test_create_target()
