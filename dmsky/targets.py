#!/usr/bin/env python
"""
Module for target objects. A target carries the physical information
about a target (ra, dec, distance, density profile, etc.) and
interfaces with the `jcalc` module to calculate the l.o.s. integral at
a given sky position. Classes inherit from the `Target` baseclass and
can be created with the `factory` function.
"""
import sys
import os, os.path
from os.path import abspath, dirname, join

import numpy as np

from dmsky.jcalc import LoSIntegral, LoSIntegralFast, LoSIntegralInterp
from dmsky.utils import coords
from dmsky.utils.units import Units
from dmsky.utils.tools import update_dict, merge_dict, yaml_load, get_items, item_version
from dmsky.library import ObjectLibrary
from dmsky.density import factory as density_factory

class Target(object):
    defaults = (
        ('title',     'Target', 'Human-readable name'          ),
        ('name',      'target', 'Machine-readable name'        ),
        ('abbr',      'Tar',    'Title abbreviation'           ),
        ('altnames',   [],      'Alternative names'            ),
        ('ra',         0.0,     'Right Ascension (deg)'        ),
        ('dec',        0.0,     'Declination (deg)'            ),
        ('distance',   0.0,     'Distance (kpc)'               ),
        ('profile',    {},      'Density Profile (see `jcalc`)'),
        ('references', [],      'Literature references'        ),
        ('color',      'k',     'Plotting color'               ),
    )

    def __init__(self, **kwargs):
        self._load(**kwargs)

    def _load(self, **kwargs):
        kw = dict([(d[0],d[1]) for d in self.defaults])
        kw.update(kwargs)

        self.__dict__.update(kw)
        self._load_profile()

    def _load_profile(self, mode='interp'):
        # Convert profile units to msun/kpc3
        units = self.profile.pop('units',None)
        if units:
            density,distance = units.rsplit('_',1)
            self.profile['rhos'] *= getattr(Units,density)/Units.msun_kpc3
            self.profile['rs']   *= getattr(Units,distance)/Units.kpc
        self.density = density_factory(**self.profile)

        distance = self.distance
        if mode == 'interp':
            self.jlosfn = LoSIntegralInterp(self.density, distance, ann=True)
            self.dlosfn = LoSIntegralInterp(self.density, distance, ann=False)
        elif mode == 'fast':
            self.jlosfn = LoSIntegralFast(self.density, distance, ann=True)
            self.dlosfn = LoSIntegralFast(self.density, distance, ann=False)
        else:
            self.jlosfn = LoSIntegral(self.density, distance, ann=True)
            self.dlosfn = LoSIntegral(self.density, distance, ann=False)

    def __str__(self):
        ret = self.__class__.__name__
        for k in ['name','ra','dec','distance','density']:
            v = self.__dict__[k]
            ret += '\n  %-15s: %s'%(k,v)
        return ret

    @property
    def glon(self):
        return coords.cel2gal(self.ra,self.dec)[0]

    @property
    def glat(self):
        return coords.cel2gal(self.ra,self.dec)[1]

    def jvalue(self,ra,dec):
        sep = coords.angsep(self.ra,self.dec,ra,dec)
        return self.jlosfn(np.radians(sep))

    def jsigma(self,ra,dec):
        raise Exception('Not implemented')

    def dvalue(self,ra,dec):
        sep = coords.angsep(self.ra,self.dec,ra,dec)
        return self.dlosfn(np.radians(sep))
    
    def dsigma(self,ra,dec):
        raise Exception('Not implemented')

class Galactic(Target): pass
class Dwarf(Target): pass
class Galaxy(Target): pass
class Cluster(Target): pass
class Isotropic(Target): pass

def factory(type, **kwargs):
    import dmsky.factory
    return dmsky.factory.factory(type, module=__name__,**kwargs)

class TargetLibrary(ObjectLibrary):
    _defaults = (
        ('path', join(dirname(abspath(__file__)),'data/targets')),
    )

    def get_target_dict(self, name, version=None, **kwargs):
        """ Step through the various levels of dependencies to get the
        full dictionary for a target.

        target: version -> ... -> target: default -> default: type
        """
        n,v = item_version(name)
        if version is not None and v is not None:
            msg = "Version specified twice: %s, %s"%(name,version)
            raise ValueError(msg)

        if v is not None:   version = v
        if version is None: version = 'default'
        name = n

        # Start with the target:version requested
        ret = self.library[name][version]

        # Walk down the chain until we either return None or the
        # 'default' version
        while (version is not None) and (version != 'default'):
            version = ret.get('base','default')
            ret = merge_dict(self.library[name][version], ret)
        ret['version'] = version

        # And finally, overwrite with kwargs
        update_dict(ret,kwargs)
        return ret

    def create_target(self, name, version=None, **kwargs):
        kw = self.get_target_dict(name,version,**kwargs)
        return factory(**kw)

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

