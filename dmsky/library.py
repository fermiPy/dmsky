#!/usr/bin/env python
"""
Module for performing target lookup.
"""
import os, os.path
from os.path import abspath, dirname, join
import glob
import copy

from collections import defaultdict as ddict
from dmsky.utils.tools import update_dict, merge_dict, yaml_load
import dmsky.targets

class TargetLibrary(object):

    def __init__(self, path=None):
        self.paths = [] if path is None else [path]
        if os.getenv('DMSKY_PATH') is not None:
            self.paths += os.getenv('DMSKY_PATH').split(':')
        self.paths += [join(dirname(abspath(__file__)),'data/targets')]
        self.defaults = self.load_defaults(self.paths)
        self.targets = self.load_targets(self.paths)

    @classmethod
    def load_defaults(cls, paths):
        defaults = dict()
        for path in paths[::-1]:
            infiles = glob.glob(path+'/defaults.yaml')
            for f in infiles:
                d = yaml_load(f)
                update_dict(defaults,d)
        return defaults

    @classmethod
    def load_targets(cls, paths):
        targets = dict()
        for path in paths[::-1]:
            subdirs = os.walk(path).next()[1]
            for subdir in subdirs:
                infiles = glob.glob(join(path,subdir)+'/*.yaml')
                for f in infiles:
                    d = yaml_load(f)
                    update_dict(targets,d)
        return targets

    def get_target_dict(self, name, version='default', **kwargs):
        """ Step through the various levels of dependencies to get the
        full dictionary for a target.

        target: version -> ... -> target: default -> default: type
        """
        
        # Start with the target:version requested
        ret = self.targets[name][version]

        # Walk down the chain until we either return None or the
        # 'default' version
        while (version is not None) and (version != 'default'):
            version = ret.get('base','default')
            ret = merge_dict(self.targets[name][version], ret)

        # Now update with the default type dictionary
        ret = merge_dict(self.defaults.get(ret.get('type','target')),ret)

        # And finally, overwrite with kwargs
        update_dict(ret,kwargs)
        return ret

    def create_target(self, name, version='default', **kwargs):
        kw = self.get_target_dict(name,version,**kwargs)
        return dmsky.targets.factory(**kw)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
