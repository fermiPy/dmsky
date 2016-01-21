#!/usr/bin/env python
"""
Module for performing target lookup.
"""
import os, os.path
from os.path import abspath, dirname, join
import glob
import copy
import logging

from collections import defaultdict as ddict, OrderedDict as odict

from dmsky.utils.tools import update_dict, merge_dict, yaml_load, get_items
import dmsky.targets

class ObjectLibrary(object):
    _defaults = (
        ('path', join(dirname(abspath(__file__)),'data')),
    )

    def __init__(self, path=None):
        self.defaults = odict([(d[0],d[1]) for d in self._defaults])

        self.paths = [self.defaults['path']]
        if os.getenv('DMSKY_PATH') is not None:
            self.paths += os.getenv('DMSKY_PATH').split(':')[::-1]
        if path is not None: self.paths += [path]
        logging.debug('PATHS: %s'%self.paths)

        self.library = self.load_library(self.paths)

    @classmethod
    def load_library(cls, paths):
        library = dict()
        for path in paths:
            subdirs = [path] + os.walk(path).next()[1]
            for subdir in subdirs:
                infiles = glob.glob(join(path,subdir)+'/*.yaml')
                for f in infiles:
                    d = yaml_load(f)
                    update_dict(library,d)
        return library

        

class TargetLibrary(ObjectLibrary):
    _defaults = (
        ('path', join(dirname(abspath(__file__)),'data/targets')),
    )


    def get_target_dict(self, name, version='default', **kwargs):
        """ Step through the various levels of dependencies to get the
        full dictionary for a target.

        target: version -> ... -> target: default -> default: type
        """
        
        # Start with the target:version requested
        ret = self.library[name][version]

        # Walk down the chain until we either return None or the
        # 'default' version
        while (version is not None) and (version != 'default'):
            version = ret.get('base','default')
            ret = merge_dict(self.library[name][version], ret)

        # And finally, overwrite with kwargs
        update_dict(ret,kwargs)
        return ret

    def create_target(self, name, version='default', **kwargs):
        kw = self.get_target_dict(name,version,**kwargs)
        return dmsky.targets.factory(**kw)


class RosterLibrary(ObjectLibrary):
    _defaults = (
        ('path', join(dirname(abspath(__file__)),'data/rosters')),
    )

    def get_roster_list(self, name, *args):
        items = [name] + list(args)
        ret = get_items(items,self.library)

        if len(ret) != len(set(ret)):
            msg = "Duplicate entries"
            logging.warning(msg)

        return ret

    def create_roster(self, name, *args):
        ar = self.get_roster_list(name,*args)
        return dmsky.rosters.factory(*ar)


def TargetFactory(name, version='default', **kwargs):
    targets = TargetLibrary()


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
