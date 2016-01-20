#!/usr/bin/env python
"""
Module for building targets.
"""
import copy

class Target(object):
    def __init__(self, *args, **kwargs):
        self._load(**kwargs)

    def _load(**kwargs):
        self.params = copy.deepcopy(kwargs)

    def jfactor(ra,dec):
        pass

    def jsigma(ra,dec):
        pass

    def dfactor(ra,dec):
        pass
    
    def dsigma(ra,dec):
        pass

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
