#!/usr/bin/env python
import os

import healpy as hp
import numpy as np
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')
import pylab as plt

import dmsky.targets
import dmsky.roster
import dmsky.skymap
from dmsky.utils.healpix import pix2ang,query_disc
from dmsky.utils.units import Units

targets = dmsky.targets.TargetLibrary()
#draco = targets.create_target('draco')

rosters = dmsky.roster.RosterLibrary()
#test = rosters.create_roster('test2015')
#dsphs_gc = rosters.create_roster('dsphs_gc')
gc = rosters.create_roster('galactic:ackermann2013_nfw')
dsphs = rosters.create_roster('martinez2015_nfw')

#psi = np.deg2rad(np.linspace(0.1,179,500))
#gc[0].jlosfn(psi)

#target = dmsky.targets.Target(profile=dict(rhos=1.0,rs=1.0,type='NFW'),distance=10,ra=0,dec=0)

r = gc + dsphs
#r = dsphs

nside = 128
skymap = dmsky.skymap.Skymap(r, nside=nside)
skymap.fill()

gc_skymap = dmsky.skymap.Skymap(gc, nside=nside)
gc_skymap.fill()

dsphs_skymap = dmsky.skymap.Skymap(dsphs, nside=nside)
dsphs_skymap.fill()

hp.mollview(np.log10(skymap.values))
plt.savefig('dsphs_gc_skymap.png',bbox_inches='tight')
#hp.mollview(np.log10(dsphs_skymap.values))
#hp.mollview(np.log10(gc_skymap.values))
#hp.mollview(np.log10(dsphs_skymap.values+gc_skymap.values))
#hp.mollview(np.log10(50*dsphs_skymap.values+gc_skymap.values))

mask = np.zeros_like(skymap.values,dtype=bool)
mask[query_disc(nside,gc[0].glon,gc[0].glat,radius=41)] = True
mask[(np.abs(skymap.lon-360.*(skymap.lon>180))>6) & (np.abs(skymap.lat)<5)] = False

hp.mollview(np.log10(skymap.values*mask))
plt.savefig('gc_region.png',bbox_inches='tight')

jfactor = gc_skymap.values[~mask].sum() * hp.nside2pixarea(nside) * Units.msun2_kpc5/Units.gev2_cm5

plt.ion()
plt.show()
