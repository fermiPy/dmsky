#!/usr/bin/env python
"""
Test for target creation with defaults.
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

    umi2022 = target_library.create_target('ursa_minor',default='default2022')
    assert umi2022.ra == 227.242
    assert umi2022.dec == 67.2221
    assert umi2022.distance == 76.2

    #print(umi2020)

def test_target_jfactor():
    # Test jfactor consistency

    # Build the library of pre-defined rosters
    target_library = dmsky.targets.TargetLibrary()

    target = target_library.create_target('ursa_minor:geringer-sameth2015_nfw',default='default')
    np.testing.assert_allclose(target.j_photo(), 5.038e18, rtol=1e-3)
    np.testing.assert_allclose(target.j_integ, 8.51e18)
    np.testing.assert_allclose(target.j_sigma, 0.23)

    target = target_library.create_target('ursa_minor:geringer-sameth2015_nfw',default='default2022')
    np.testing.assert_allclose(target.j_photo(), 3.498e18, rtol=1e-3)
    np.testing.assert_allclose(target.j_integ, 8.51e18)
    np.testing.assert_allclose(target.j_sigma, 0.23)

    target = target_library.create_target('ursa_minor:pace2019_nfw',default='default2022')
    np.testing.assert_allclose(target.j_photo(), 3.498e18, rtol=1e-3)
    np.testing.assert_allclose(target.j_integ, 5.62e18)
    np.testing.assert_allclose(target.j_sigma, 0.12)


if __name__ == '__main__':
    test_create_target()
    test_target_jfactor()
