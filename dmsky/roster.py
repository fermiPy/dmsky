#!/usr/bin/env python
"""
Module for `Rosters` (collections of `Targets`).
"""
from __future__ import absolute_import, division, print_function

from os.path import abspath, dirname, join
from collections import OrderedDict as odict

from dmsky.utils.tools import get_items
from dmsky.library import ObjectLibrary
from dmsky.targets import Target, TargetLibrary

global targetlib
targetlib = TargetLibrary()


class Roster(odict):
    """A self-consistent set of search targets, typically with exactly one version
    for any given target.

    """
    def __init__(self, *targets):
        """C'tor

        Targets can be passed either:
        1) As string elements of `targets` (looked up in target library)
        2) As object elements of `targets` (name taken from the target)
        3) As a list in args
        4) As an existing Roster
        """
        super(Roster, self).__init__()
        self._parse_args(targets)

    def _parse_args(self, targets):
        """Internal function to parse all the various possible
        type of inputs to the c'tor of this class
        """
        for arg in targets:
            if isinstance(arg, str):
                target = targetlib.create_target(arg)
                self[target.name] = target
            elif isinstance(arg, Target):
                self[arg.name] = arg
            elif isinstance(arg, Roster):
                self._parse_args(arg.values())
            elif hasattr(arg, '__iter__'):
                self._parse_args(arg)
            else:
                msg = "Unrecognized argument type: %s" % (type(arg))
                raise ValueError(msg)

    def __add__(self, other):
        """Add two Rosters.
        """
        if isinstance(other, Target):
            items = self.values() + [other]
        elif isinstance(other, Roster):
            items = []
            items += self.values()
            items += other.values()
        return self.__class__(*items)

    def __getitem__(self, key):
        """Get a member of this Roster by name
        """
        try:
            return super(Roster, self).__getitem__(key)
        except KeyError as e:
            if isinstance(key, int):
                return list(self.values())[key]
            else:
                raise e


def factory(*targets):
    """Factory for creating a roster. Arguments are passed directly to
    the Roster constructor.
    """
    return Roster(*targets)


class RosterLibrary(ObjectLibrary):
    """A top-level object, keeping track of all the `Roster` objects
    that we have created
    """

    _suffix = 'rosters'

    _defaults = (
        ('path', join(dirname(abspath(__file__)), 'data', _suffix)),
    )

    def get_roster_list(self, name, *targets):
        """Get the list of objects on a roster,
        creating them in needed.

        Parameters
        ----------

        name : str
            A name for the list

        targets : list
            Objects that will go on the list


        Returns
        -------

        items : list
            List of items needed to create a `Roster`

        """

        items = [name] + list(targets)
        ret = get_items(items, self.library)

        if len(ret) != len(set(ret)):
            msg = "Duplicate entries"
            print ("Warning: %s" % msg)
            #logging.warning(msg)

        return ret

    def create_roster(self, name, *targets):
        """Create a roster

        Parameters
        ----------

        name : str
            A name for the Roster

        targets : list
            Objects that will go on the Roster


        Returns
        -------

        roster : `Roster`
            The newly created Roster

        """
        ar = self.get_roster_list(name, *targets)
        return factory(*ar)
