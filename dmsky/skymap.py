#!/usr/bin/env python
"""
Module for creating HEALPix skymaps.
"""
import healpy as hp
import numpy as np

from dmsky.utils.healpix import pix2ang, query_disc
from dmsky.utils.coords import gal2cel
from dmsky.roster import Roster

class Skymap(object):
    _defaults = (
        # 512 has a pixel size of ~0.11 deg
        ('nside',    128, 'Healpix nside'          ),
        ('coord',  'gal', 'Coordinate frame'       ),
    )

    def __init__(self, *args, **kwargs):
        self._load(**kwargs)
        self.roster = Roster(*args)
    
    def _load(self, **kwargs):
        kw = dict([(d[0],d[1]) for d in self._defaults])
        kw.update(kwargs)
        self.__dict__.update(kw)

    def fill(self):
        self.values = np.zeros(hp.nside2npix(self.nside))
        self.pix = np.arange(hp.nside2npix(self.nside))
        self.lon,self.lat = pix2ang(self.nside,self.pix)
        
        for target in self.roster.values():
            print target
            if self.coord == 'gal':
                idx = query_disc(self.nside,target.glon,target.glat,target.psi_max)
                lon,lat = gal2cel(self.lon[idx],self.lat[idx])
            else:
                idx = query_disc(self.nside,target.ra,target.dec,target.psi_max)
                lon,lat = self.lon[idx],self.lat[idx]

            self.values[idx] += target.jvalue(lon,lat)
            
        
if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
