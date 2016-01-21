#!/usr/bin/env python
"""
Utilities for working with coordinates.
"""
import numpy as np
from  astropy.coordinates import SkyCoord

def gal2cel(glon,glat):
    coords = SkyCoord(glon,glat,frame='galactic',unit='deg')
    return coords.icrs.ra.deg,coords.icrs.dec.deg

def cel2gal(ra,dec):
    coords = SkyCoord(ra,dec,frame='icrs',unit='deg')
    return coords.galactic.l.deg,coords.galactic.b.deg

def angsep(lon1,lat1,lon2,lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Faster than creating astropy coordinate objects.

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    lon1,lat1 = np.radians([lon1,lat1])
    lon2,lat2 = np.radians([lon2,lat2])
    
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1,num2), denominator))

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
