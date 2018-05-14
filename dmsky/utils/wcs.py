#!/usr/bin/env python
"""
Interface with wcs.

Adapted from fermipy.skymap
"""
__author__ = "Alex Drlica-Wagner"
import numpy as np

from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord


def create_wcs(skydir, coordsys='CEL', projection='AIT',
               cdelt=1.0, crpix=1., naxis=2, energies=None):
    """Create a WCS object.
    Parameters
    ----------
    skydir : `~astropy.coordinates.SkyCoord`
        Sky coordinate of the WCS reference point.
    coordsys : str
    projection : str
    cdelt : float
    crpix : float or (float,float)
        In the first case the same value is used for x and y axes
    naxis : {2, 3}
       Number of dimensions of the projection.
    energies : array-like
       Array of energies that defines the third dimension if naxis=3.
    """

    w = WCS(naxis=naxis)

    if coordsys == 'CEL':
        w.wcs.ctype[0] = 'RA---%s' % (projection)
        w.wcs.ctype[1] = 'DEC--%s' % (projection)
        w.wcs.crval[0] = skydir.icrs.ra.deg
        w.wcs.crval[1] = skydir.icrs.dec.deg
    elif coordsys == 'GAL':
        w.wcs.ctype[0] = 'GLON-%s' % (projection)
        w.wcs.ctype[1] = 'GLAT-%s' % (projection)
        w.wcs.crval[0] = skydir.galactic.l.deg
        w.wcs.crval[1] = skydir.galactic.b.deg
    else:
        raise Exception('Unrecognized coordinate system.')

    try:
        w.wcs.crpix[0] = crpix[0]
        w.wcs.crpix[1] = crpix[1]
    except IndexError:
        w.wcs.crpix[0] = crpix
        w.wcs.crpix[1] = crpix
    w.wcs.cdelt[0] = -cdelt
    w.wcs.cdelt[1] = cdelt

    w = WCS(w.to_header())
    if naxis == 3 and energies is not None:
        w.wcs.crpix[2] = 1
        w.wcs.crval[2] = energies[0]
        w.wcs.cdelt[2] = energies[1] - energies[0]
        w.wcs.ctype[2] = 'Energy'
        w.wcs.cunit[2] = 'MeV'

    return w


def create_image_wcs(skydir, cdelt, npix, coordsys='CEL', projection='AIT'):
    """Create a blank image and associated WCS object
    """
    if np.isscalar(npix):
        npix = [npix, npix]
    crpix = np.array([n / 2. + 0.5 for n in npix])
    wcs = create_wcs(skydir, coordsys, projection, cdelt, crpix)
    return np.zeros(npix).T, wcs


def get_pixel_skydirs(npix, wcs):
    """Get a list of sky coordinates for the centers of every pixel.
    """
    if np.isscalar(npix):
        npix = [npix, npix]
    xpix = np.linspace(0, npix[0] - 1., npix[0])
    ypix = np.linspace(0, npix[1] - 1., npix[1])
    xypix = np.meshgrid(xpix, ypix, indexing='ij')
    return SkyCoord.from_pixel(np.ravel(xypix[0]),
                               np.ravel(xypix[1]), wcs)


def create_image_hdu(data, wcs, name=None):
    """Create a `astropy.io.fits.ImageHDU` object
    """
    if name is None:
        return fits.PrimaryHDU(data, header=wcs.to_header())
    return fits.ImageHDU(data, header=wcs.to_header(),
                         name=name)


def write_image_hdu(filename, data, wcs, name=None, clobber=False):
    """Write an image to a file
    """
    hdu = create_image_hdu(data, wcs, name)
    hdu.writeto(filename, clobber=clobber)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
