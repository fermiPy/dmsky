#!/usr/bin/env python
"""
Simple module to deal with unit conversion
"""
__all__ = ['Units']

import numpy as np

class Units(object):
    # Could be replaced by astropy...

    # length
    cm = 1
    m  = 1e2           # m to cm
    km = m*1e3         # km to cm
    pc = 3.08568e18   # pc to cm
    kpc = pc*1e3      # kpc to cm
    m2 = 1e4

    # mass
    g = 1.0
    msun = 1.98892e33 # solar mass to g
    gev = 1.78266e-24 # gev to g

    # density
    msun_pc3 = msun*np.power(pc,-3) 
    msun_kpc3 = msun*np.power(kpc,-3)
    msun2_pc5 = np.power(msun,2)*np.power(pc,-5)
    msun2_kpc5 = np.power(msun,2)*np.power(kpc,-5)
    gev2_cm5 = np.power(gev,2)
    gev_cm3 = np.power(gev,1)
    gev_cm2 = np.power(gev,1)
    g_cm3 = 1.0
    cm3_s = 1.0

    # random
    hr = 3600.
    deg2 = np.power(np.pi/180.,2)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
