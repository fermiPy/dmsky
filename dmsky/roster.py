#!/usr/bin/env python
"""
Module for `Rosters` (collections of `Targets`).
"""
import os, os.path
from os.path import abspath, dirname, join
from collections import OrderedDict as odict

from dmsky.utils.tools import update_dict, merge_dict, yaml_load, get_items
from dmsky.library import ObjectLibrary
from dmsky.targets import Target, TargetLibrary

global targetlib
targetlib = TargetLibrary()

class Roster(odict):

    def __init__(self, *args, **kwargs):
        """
        Targets can be passed either:
        1) As string elements of `args` (looked up in target library)
        2) As object elements of `args` (name taken from the target)
        3) As a list in args
        4) As an existing Roster
        """
        super(Roster,self).__init__()
        self._parse_args(args)
        
    def _parse_args(self, args):
        for arg in args:
            if isinstance(arg,basestring):
                target = targetlib.create_target(arg)
                self[target.name] = target
            elif isinstance(arg,Target):
                self[arg.name] = arg
            elif isinstance(arg,Roster):
                self._parse_args(arg.values())
            elif hasattr(arg,'__iter__'):
                self._parse_args(arg)
            else:
                msg = "Unrecognized argument type: %s"%(type(arg))
                raise ValueError(msg)
            
    def __add__(self, other):
        """ Add two Rosters.
        """
        if isinstance(other,Target):
            items = self.values() + [other]
        elif isinstance(other,Roster):
            items = self.values() + other.values()
        return self.__class__(*items)

    def __getitem__(self, key):
        try:
            return super(Roster,self).__getitem__(key)
        except KeyError, e:
            if isinstance(key,int):
                return self.values()[key]
            else:
                raise e

def factory(*args, **kwargs):
    """Factory for creating a roster. Arguments are passed directly to
    the Roster constructor.
    """
    return Roster(*args, **kwargs)

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
        return factory(*ar)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    rl = RosterLibrary()
    ro = rl.create_roster('martinez2015_nfw')

    draco = ro['draco']
    
    colNames = ['name','ra','dec','distance','proftype','glat','glon','j_integ','d_integ']

    from dmsky.file_io import table as table_io    
    tab = table_io.make_table_for_roster(colNames,ro)
