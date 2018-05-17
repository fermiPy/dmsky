#!/usr/bin/env python
"""
Utilities for working with coordinates.
"""
import numpy as np


def gal2cel(glon, glat):
    """
    Convert Galactic coordinates to celestial coordinates
    (J2000). (Much faster than astropy for small arrays.)

    Parameters
    ----------
    glon : `numpy.array`
    glat : `numpy.array`
        Galactic Coordinates (in degrees)

    Returns
    -------
    ra : `numpy.array`
    dec : `numpy.array`
        Celestical Coordinates (in degrees)

    """
    glat = np.radians(glat)
    sin_glat = np.sin(glat)
    cos_glat = np.cos(glat)

    glon = np.radians(glon)
    ra_gp = np.radians(192.85948)
    de_gp = np.radians(27.12825)
    lcp = np.radians(122.932)

    sin_lcp_glon = np.sin(lcp - glon)
    cos_lcp_glon = np.cos(lcp - glon)

    sin_d = (np.sin(de_gp) * sin_glat) \
        + (np.cos(de_gp) * cos_glat * cos_lcp_glon)
    ramragp = np.arctan2(cos_glat * sin_lcp_glon,
                         (np.cos(de_gp) * sin_glat)
                         - (np.sin(de_gp) * cos_glat * cos_lcp_glon))
    dec = np.arcsin(sin_d)
    ra = (ramragp + ra_gp + (2. * np.pi)) % (2. * np.pi)
    return np.degrees(ra), np.degrees(dec)


def cel2gal(ra, dec):
    """
    Convert celestial coordinates (J2000) to Galactic
    coordinates. (Much faster than astropy for small arrays.)

    Parameters
    ----------
    ra : `numpy.array`
    dec : `numpy.array`
        Celestical Coordinates (in degrees)


    Returns
    -------
    glon : `numpy.array`
    glat : `numpy.array`
        Galactic Coordinates (in degrees)

    """
    dec = np.radians(dec)
    sin_dec = np.sin(dec)
    cos_dec = np.cos(dec)

    ra = np.radians(ra)
    ra_gp = np.radians(192.85948)
    de_gp = np.radians(27.12825)

    sin_ra_gp = np.sin(ra - ra_gp)
    cos_ra_gp = np.cos(ra - ra_gp)

    lcp = np.radians(122.932)
    sin_b = (np.sin(de_gp) * sin_dec) + (np.cos(de_gp) * cos_dec * cos_ra_gp)

    lcpml = np.arctan2(cos_dec * sin_ra_gp,
                       (np.cos(de_gp) * sin_dec)
                       - (np.sin(de_gp) * cos_dec * cos_ra_gp))

    glat = np.arcsin(sin_b)
    glon = (lcp - lcpml + (2. * np.pi)) % (2. * np.pi)
    return np.degrees(glon), np.degrees(glat)


def angsep(lon1, lat1, lon2, lat2):
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
    lon1, lat1 = np.radians([lon1, lat1])
    lon2, lat2 = np.radians([lon2, lat2])

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1, num2), denominator))
