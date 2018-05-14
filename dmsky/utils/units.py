#!/usr/bin/env python
"""
Simple module to deal with unit conversion
"""
__all__ = ['Units']

import numpy as np


class Units(object):
    """Simple class to keep track of unit conversions"""
    # Could be replaced by astropy...

    # length
    cm = 1
    m = 1e2           # m to cm
    km = m * 1e3         # km to cm
    pc = 3.08568e18   # pc to cm
    kpc = pc * 1e3      # kpc to cm
    m2 = 1e4

    # mass
    g = 1.0
    msun = 1.98892e33  # solar mass to g
    gev = 1.78266e-24  # gev to g

    # density
    msun_pc3 = msun * np.power(pc, -3)
    msun_kpc3 = msun * np.power(kpc, -3)
    msun2_pc5 = np.power(msun, 2) * np.power(pc, -5)
    msun2_kpc5 = np.power(msun, 2) * np.power(kpc, -5)
    gev2_cm5 = np.power(gev, 2)
    gev_cm3 = np.power(gev, 1)
    gev_cm2 = np.power(gev, 1)
    g_cm3 = 1.0
    cm3_s = 1.0

    # random
    hr = 3600.
    deg2 = np.power(np.pi / 180., 2)

    # This is to convert stuff to strings that astropy units understands
    map_to_astropy = {'gev2_cm5': 'GeV2 / cm5',
                      'gev_cm3': 'GeV / cm3',
                      'gev_cm2': 'GeV / cm2',
                      'g_cm3': 'g / cm3',
                      'cm3_s': 'cm3 / s'}
    map_from_astropy = {}
    for k, v in map_to_astropy.items():
        map_from_astropy[v] = k

    @staticmethod
    def get_value(key):
        """ Get a conversion value based on a key

        This is here to make it easy to automate unit conversion

        Parameters
        ----------
        key   : str, a key corresponding to one of the globals defined above

        Returns
        -------
        the conversion constant
        """
        if key is None:
            return None
        try:
            return getattr(Units, key)
        except AttributeError:
            try:
                newkey = Units.map_from_astropy[key]
                return getattr(Units, newkey)
            except:
                raise KeyError("Did not recoginze units %s" % key)

    @staticmethod
    def convert_to(value, key):
        """ Convert from cgs units to a different type of units

        Parameters
        ----------
        value : scalar or array-like, the input value(s)
        key   : str, a key corresponding to one of the globals defined above

        Returns
        -------
        the input values, converted to requested units
        """
        conv = Units.get_value(key)
        if conv is None:
            return value
        return value / conv

    @staticmethod
    def convert_from(value, key):
        """ Convert to cgs units from a different type of units

        Parameters
        ----------
        value : scalar or array-like, the input value(s)
        key   : str, a key corresponding to one of the globals defined above

        Returns
        -------
        the input values, converted to cgs units
        """
        conv = Units.get_value(key)
        if conv is None:
            return value
        return value * conv


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
