#!/usr/bin/env python
"""
Module for performing target lookup.
"""
from __future__ import absolute_import, division, print_function

import os
import os.path
from os.path import abspath, dirname, join
import glob
import logging

from collections import OrderedDict as odict

from dmsky.utils.tools import update_dict, yaml_load


class ObjectLibrary(object):
    """Library that keeps track of object we are building
    """

    _suffix = ''

    _defaults = (
        ('path', join(dirname(abspath(__file__)), 'data')),
    )

    def __init__(self, path=None):
        self.defaults = odict([(d[0], d[1]) for d in self._defaults])
        self.paths = [self.defaults['path']]
        if os.getenv('DMSKY_PATH') is not None:
            self.paths += [join(p, self._suffix) for p in os.getenv('DMSKY_PATH').split(':')]
        if path is not None:
            self.paths += [path]
        logging.debug('PATHS: %s' % self.paths)

        self.library = self.load_library(self.paths)

    @classmethod
    def load_library(cls, paths):
        """Build the library by walking through paths
        """
        library = dict()
        for path in paths:
            # Should use logging
            print("Using %s for %s" % (path, cls.__name__))
            subdirs = [path] + next(os.walk(path))[1]
            for subdir in subdirs:
                infiles = glob.glob(join(path, subdir) + '/*.yaml')
                for f in infiles:
                    d = yaml_load(f)
                    update_dict(library, d)
        return library

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    lib = ObjectLibrary()
