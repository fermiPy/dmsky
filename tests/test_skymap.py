#!/usr/bin/env python
import os

import healpy as hp
import numpy as np
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')
import pylab as plt

import dmsky.targets
import dmsky.roster
import dmsky.skymap
from dmsky.utils.healpix import pix2ang,query_disc
from dmsky.utils.units import Units


def test_skymap_from_roster():

    targets = dmsky.targets.TargetLibrary()
    rosters = dmsky.roster.RosterLibrary()
    gc = rosters.create_roster('galactic:ackermann2013_nfw')
    dsphs = rosters.create_roster('martinez2015_nfw')
    
    r = gc + dsphs

    nside = 128
    skymap = dmsky.skymap.Skymap(r, nside=nside)

    gc_skymap = dmsky.skymap.Skymap(gc, nside=nside)

    dsphs_skymap = dmsky.skymap.Skymap(dsphs, nside=nside)

    hp.mollview(np.log10(skymap.values))
    plt.savefig('dsphs_gc_skymap.png',bbox_inches='tight')

    mask = np.zeros_like(skymap.values,dtype=bool)
    mask[query_disc(nside,gc[0].glon,gc[0].glat,radius=41)] = True
    mask[(np.abs(skymap.lon-360.*(skymap.lon>180))>6) & (np.abs(skymap.lat)<5)] = False

    hp.mollview(np.log10(skymap.values*mask))
    plt.savefig('gc_region.png',bbox_inches='tight')

    jfactor = gc_skymap.values[~mask].sum() * hp.nside2pixarea(nside) * Units.msun2_kpc5/Units.gev2_cm5


if __name__ == '__main__':
    test_skymap_from_roster()
    plt.ion()
    plt.show()
