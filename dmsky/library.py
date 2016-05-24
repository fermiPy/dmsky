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
            print path
            subdirs = [path] + os.walk(path).next()[1]
            for subdir in subdirs:
                infiles = glob.glob(join(path,subdir)+'/*.yaml')
                for f in infiles:
                    d = yaml_load(f)
                    update_dict(library,d)
        return library

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    lib = ObjectLibrary()
