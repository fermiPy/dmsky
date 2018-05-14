#!/usr/bin/env python
"""
Module for creating HEALPix skymaps.
"""
from __future__ import absolute_import, division, print_function

import healpy as hp
import numpy as np

from dmsky.utils.healpix import pix2ang, query_disc
from dmsky.utils.coords import gal2cel
from dmsky.roster import Roster


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
        self.nside = kwargs.get('nside', 512)
        self.coord = kwargs.get('coord', 'gal')
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

            self.values[idx] += target.jvalue(lon, lat)
