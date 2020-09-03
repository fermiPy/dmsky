#!/usr/bin/env python
"""
Module for creating HEALPix skymaps.
"""
from __future__ import absolute_import, division, print_function

import healpy as hp
import numpy as np

from dmsky.utils.healpix import pix2ang, query_disc, nside_to_order
from dmsky.utils.coords import gal2cel
from dmsky.roster import Roster

from astropy.io import fits

class Skymap(object):
    """Object to fill a HEALPix sky-map with J values from a `Roster`
    """

    def __init__(self, *targets, **kwargs):
        """C'tor

        Parameters
        ----------

        targets : list
            List of objects to include in this map


        Keyword arguments
        -----------------

        nside : int
            HEALPix n-side parameter [512]

        coord : str
            Coordinate system 'cel' or 'gal' ['gal']

        """
        self.filename = kwargs.get('filename', None)
        if self.filename is None:
            self._construct(*targets, **kwargs)
        else:
            self._read_fits(filename)


    def _construct(self, *targets, **kwargs):        
        self.nside = kwargs.get('nside', 512)
        self.coord = kwargs.get('coord', 'gal')
        self.ann = kwargs.get('ann', True)
        if self.ann:
            self.units = 'GeV^2 cm^5 s^-2'
        else:
            self.units = 'GeV cm^3'
        self.values = np.zeros(hp.nside2npix(self.nside))
        self.pix = np.arange(hp.nside2npix(self.nside))
        self.lon, self.lat = pix2ang(self.nside, self.pix)
        self.roster = Roster(*targets)
        self._fill()


    def _fill(self):
        """Fills map values """
        for target in self.roster.values():
            print(target)
            if self.coord == 'gal':
                idx = query_disc(self.nside, target.glon, target.glat, target.psi_max)
                lon, lat = gal2cel(self.lon[idx], self.lat[idx])
            else:
                idx = query_disc(self.nside, target.ra, target.dec, target.psi_max)
                lon, lat = self.lon[idx], self.lat[idx]

            if self.ann:
                self.values[idx] += target.jvalue(lon, lat)
            else:
                self.values[idx] += target.dvalue(lon, lat)


    def write_fits(self, filename, overwrite=True):
        """Write this Skymap to a FITS file"""
        col = fits.Column(name='Bin0', format='1D', array=self.values)
        header = fits.Header()
        header['COORDSYS'] = self.coord
        header['COORDTYPE'] = self.coord
        header['ORDERING'] = 'RING'
        header['INDXSCHM'] = 'IMPLICIT'
        header['NSIDE'] = self.nside
        header['ORDER'] = nside_to_order(self.nside)
        header['FIRSTPIX'] = 0
        header['LASTPIX'] = len(self.values) - 1
        header['HPX_CONV'] = 'GALPROP'
        header['BANDSHDU'] = ''

        primhdu = fits.PrimaryHDU()
        skyhdu = fits.BinTableHDU.from_columns([col], header=header, name='SKYMAP')

        hdulist = fits.HDUList([primhdu, skyhdu])
        hdulist.writeto(filename, overwrite=overwrite)

    def _assert_header_key(self, header, key, value):
        """Internal method to verify that a FITS header key matches the expected value"""
        check = header.get(key)
        if check != value:
            raise ValueError("Header card %s=%s it should be %s" % (key, check, value))

    def _read_fits(self, filename):
        """Internal method to read this from a FITS file"""
        f = fits.open(filename)
        header = f[1].header
        self._assert_header_key(header, 'HPX_CONV', 'GALPROP')
        self._assert_header_key(header, 'ORDERING', 'RING')
        self._assert_header_key(header, 'INDXSCHM', 'IMPLICIT')

        self.nside = header['NSIDE']
        self.coord = header['COORDSYS']
        self.values = np.zeros(hp.nside2npix(self.nside))
        self.pix = np.arange(hp.nside2npix(self.nside))
        self.lon, self.lat = pix2ang(self.nside, self.pix)
        self.values = f[2].data.field('Bin0')
