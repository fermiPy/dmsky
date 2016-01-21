#!/usr/bin/env python
"""
Module for target objects. A target carries the physical information
about a target (ra, dec, distance, density profile, etc.) and
interfaces with the `jcalc` module to calculate the l.o.s. integral at
a given sky position. Classes inherit from the `Target` baseclass and
can be created with the `factory` function.
"""
import sys
import copy

from astropy.coordinates import SkyCoord

from jcalc import DensityProfile, LoSIntegralFn
from utils import coords

class Target(object):
    defaults = (
        ('title','Target', 'Human-readable name'),
        ('name', 'target', 'Machine-readable name'),
        ('ra',        0.0, 'Right Ascension (deg)'),
        ('dec',       0.0, 'Declination (deg)'),
        ('distance',  0.0, 'Distance (kpc)'),
        ('profile',    {}, 'Density Profile (see `jcalc`)')
        )

    def __init__(self, **kwargs):
        self._load(**kwargs)

    def _load(self, **kwargs):
        kw = dict([ (d[0],d[1]) for d in self.defaults])
        kw.update(kwargs)

        ra,dec = kw.pop('ra'),kw.pop('dec')
        self.coords = SkyCoord(ra,dec,frame='icrs',unit='deg')

        self.__dict__.update(kwargs)
        self.density = DensityProfile.create(**vars(self).get('profile',{}))
        self.jlosfn = LoSIntegralFn(self.density, self.distance, ann=True)
        self.dlosfn = LoSIntegralFn(self.density, self.distance, ann=False)

    @property
    def ra(self):
        return self.coords.ra.deg

    @property
    def dec(self):
        return self.coords.dec.deg

    @property
    def glon(self):
        return self.coords.galactic.l.deg

    @property
    def glat(self):
        return self.coords.galactic.b.deg

    def jfactor(ra,dec):
        sep = coords.angsep(self.ra,self.dec,ra,dec)
        return self.jlosfn(np.radians(sep))

    def jsigma(ra,dec):
        raise Exception('Not implemented')

    def dfactor(ra,dec):
        sep = coords.angsep(self.ra,self.dec,ra,dec)
        return self.dlosfn(np.radians(sep))
    
    def dsigma(ra,dec):
        raise Exception('Not implemented')

class Dwarf(Target):
    pass

class Galaxy(Target):
    pass

class Galactic(Target):
    pass

class Cluster(Target):
    pass

class Isotropic(Target):
    pass

def factory(type, **kwargs):
    """
    Factory for creating objects. Arguments are passed directly to the
    constructor of the chosen class.
    """
    from collections import OrderedDict as odict
    import inspect

    cls = type
    fn = lambda member: inspect.isclass(member) and member.__module__==__name__
    classes = odict(inspect.getmembers(sys.modules[__name__], fn))
    members = odict([(k.lower(),v) for k,v in classes.items()])
    
    lower = cls.lower()
    if lower not in members.keys():
        msg = "%s not found in kernels:\n %s"%(cls,classes.keys())
        #logger.error(msg)
        print msg
        msg = "Unrecognized kernel: %s"%cls
        raise Exception(msg)
 
    return members[lower](**kwargs)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
























