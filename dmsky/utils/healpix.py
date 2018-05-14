#!/usr/bin/env python
"""
Utilities for working with healpix. Originally branched from ugali.
"""
import numpy as np
import healpy

############################################################


def superpixel(subpix, nside_subpix, nside_superpix):
    """
    Return the indices of the super-pixels which contain each of the sub-pixels.
    """
    if nside_subpix == nside_superpix:
        return subpix
    theta, phi = healpy.pix2ang(nside_subpix, subpix)
    return healpy.ang2pix(nside_superpix, theta, phi)


def subpixel(superpix, nside_superpix, nside_subpix):
    """
    Return the indices of sub-pixels (resolution nside_subpix) within
    the super-pixel with (resolution nside_superpix).
    """
    if nside_superpix == nside_subpix:
        return superpix
    vec = healpy.pix2vec(nside_superpix, superpix)
    radius = np.degrees(2. * healpy.max_pixrad(nside_superpix))
    subpix = query_disc(nside_subpix, vec, radius=radius)
    pix_for_subpix = superpixel(subpix, nside_subpix, nside_superpix)
    # Might be able to speed up array indexing...
    return subpix[pix_for_subpix == superpix]

############################################################


def pix2ang(nside, pix):
    """
    Return (lon, lat) in degrees instead of (theta, phi) in radians
    """
    theta, phi = healpy.pix2ang(nside, pix)
    lon = np.degrees(phi)
    lat = 90. - np.degrees(theta)
    return lon, lat


def ang2pix(nside, lon, lat):
    """
    Input (lon, lat) in degrees instead of (theta, phi) in radians
    """
    theta = np.radians(90. - lat)
    phi = np.radians(lon)
    return healpy.ang2pix(nside, theta, phi)


def ang2vec(lon, lat):
    """
    Convert (lon, lat) to a 3-vector
    """
    theta = np.radians(90. - lat)
    phi = np.radians(lon)
    vec = healpy.ang2vec(theta, phi)
    return vec

pixToAng = pix2ang
angToPix = ang2pix
angToVec = ang2vec

############################################################


def healpixMap(nside, lon, lat, fill_value=0.):
    """
    Input (lon, lat) in degrees instead of (theta, phi) in radians.
    Returns HEALPix map at the desired resolution
    """
    pix = angToPix(nside, lon, lat)
    m = np.histogram(pix, np.arange(healpy.nside2npix(nside) + 1))[0].astype(float)
    if fill_value != 0.:
        m[m == 0.] = fill_value
    return m

############################################################


def in_pixels(lon, lat, pixels, nside):
    """
    Check if (lon,lat) in pixel list.
    """
    pix = ang2pix(nside, lon, lat)
    return np.in1d(pix, pixels)


def index_pixels(lon, lat, pixels, nside):
    """
    Find the index for object amoung a subset of healpix pixels.
    Set index of objects outside the pixel subset to -1

    # ADW: Not really safe to set index = -1 (accesses last entry);
    # -np.inf would be better, but breaks other code...
    """
    pix = ang2pix(nside, lon, lat)
    # pixels should be pre-sorted, otherwise...???
    index = np.searchsorted(pixels, pix)
    if np.isscalar(index):
        if not np.in1d(pix, pixels).any():
            index = -1
    else:
        # Find objects that are outside the roi
        #index[np.take(pixels,index,mode='clip')!=pix] = -1
        index[~np.in1d(pix, pixels)] = -1
    return index

############################################################


def query_disc(nside, lon, lat, radius, **kwargs):
    """
    Wrapper around healpy.query_disc to deal with old healpy implementation.

    nside : int
      The nside of the Healpix map.
    vec : float, sequence of 3 elements
      The coordinates of unit vector defining the disk center.
    radius : float
      The radius (in degrees) of the disk
    inclusive : bool, optional
      If False, return the exact set of pixels whose pixel centers lie
      within the disk; if True, return all pixels that overlap with the disk,
      and maybe a few more. Default: False
    fact : int, optional
      Only used when inclusive=True. The overlapping test will be done at
      the resolution fact*nside. For NESTED ordering, fact must be a power of 2,
      else it can be any positive integer. Default: 4.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

    """
    vec = ang2vec(lon, lat)
    return healpy.query_disc(nside, vec, np.radians(radius), **kwargs)
