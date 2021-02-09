#!/usr/bin/env python
"""
Test the creation of rosters
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np

import dmsky.roster

def test_create_roster():
    # Build the library of pre-defined rosters
    roster_library = dmsky.roster.RosterLibrary()

    ackermann2015 = roster_library.create_roster('ackermann2015_dsphs')
    assert len(ackermann2015) == 15
    assert 'grus_II' not in ackermann2015

    import pdb; pdb.set_trace()

    albert2017 = roster_library.create_roster('albert2017_dsphs')
    assert len(albert2017) == 41
    assert 'grus_II' in albert2017

    # Create test2015
    test2015 = roster_library.create_roster('test2015')
    assert len(test2015) == 16
    assert 'segue_1' in test2015
    assert 'bootes_I' in test2015
    assert 'bootes_II' in test2015
    assert 'grus_II' not in test2015

    # Create the test with updated cooredinates from 2020
    test2020 = roster_library.create_roster('test2015',default='default2020')
    assert len(test2020) == 16

test_create_roster()
