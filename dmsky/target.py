#!/usr/bin/env python
"""
Module for building targets.
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

class GalacticCenter(Target):
    pass

class Cluster(Target):
    pass

class Isotropic(Target):
    pass

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
