#!/usr/bin/env python
"""
Module for target objects. A target carries the physical information
about a target (ra, dec, distance, density profile, etc.) and
interfaces with the `jcalc` module to calculate the l.o.s. integral at
a given sky position. Classes inherit from the `Target` baseclass and
can be created with the `factory` function.
"""

import copy

from jcalc import DensityProfile
from utils import coords

class Target(object):
    def __init__(self, *args, **kwargs):
        self._load(**kwargs)

    def _load(self, **kwargs):
        self.params = copy.deepcopy(kwargs)
        self.density = DensityProfile.create(**self.params.get('profile',None))
        self.jlosfn = LoSIntegralFn(self.density, self.distance, ann=True)
        self.dlosfn = LoSIntegralFn(self.density, self.distance, ann=False)

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

def factory(name, **kwargs):
    """
    Factory for creating objects. Arguments are passed directly to the
    constructor of the chosen class.
    """
    fn = lambda member: inspect.isclass(member) and member.__module__==__name__
    classes = odict(inspect.getmembers(sys.modules[__name__], fn))
    members = odict([(k.lower(),v) for k,v in classes.items()])
    
    namelower = name.lower()
    if namelower not in members.keys():
        msg = "%s not found in kernels:\n %s"%(name,classes.keys())
        #logger.error(msg)
        print msg
        msg = "Unrecognized kernel: %s"%name
        raise Exception(msg)
 
    return members[namelower](**kwargs)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
