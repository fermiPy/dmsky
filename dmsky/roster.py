#!/usr/bin/env python
"""
Module for `Rosters` = collections of `Targets`.
"""
from collections import OrderedDict as odict

from dmsky.targets import Target

class Roster(odict):
    def __init__(self, *args, **kwargs):
        super(self,Roster).__init__(*args, **kwargs)
    
    def __add__(self, other):
        """ Add two Rosters.
        """
        if isinstance(other,Target):
            items = self.list + [other]
        else:
            items = self.items() + other.items()
        return self.__class__(*items)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
