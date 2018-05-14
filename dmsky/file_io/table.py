#!/usr/bin/env python
"""
Module file IO using astropy.table class.
"""

from astropy import table

from pymodeler.parameter import Parameter, Derived, Property
from dmsky.targets import Target

NullTarget = Target(title='Null', name='Null', abbr='Nul', profile={})


def get_column_kwargs(name, prop):
    """Get keyword arguments needed to build a column for a `Property` object
    """
    opt_keys = ['data', 'unit', 'shape', 'length']
    d = dict(data=None,
             name=name,
             dtype=prop.dtype,
             shape=(),
             length=0,
             description=prop.help,
             unit=None,
             format=prop.format)
    for k in opt_keys:
        try:
            d[k] = getattr(prop, k)
        except AttributeError:
            pass
    # FIXME, force strings to 20 characters
    if isinstance("", d['dtype']):
        d['dtype'] = 'S20'
    return d


def get_row_values(model, keys):
    """Get the values needed to fill a row in the table
    """
    d = {}
    for k in keys:
        d[k] = getattr(model, k)
    return d


def columns_from_property(name, prop):
    """Build a set of `Columns` objects for a `Property` object
    """
    d = get_column_kwargs(name, prop)
    return [table.Column(**d)]


def columns_from_parameter(name, prop):
    """Build a set of `Column` objects for a `Parameter` object
    """
    d = get_column_kwargs(name, prop)
    c_val = table.Column(**d)
    d['name'] += '_err'
    c_err = table.Column(**d)
    d['name'] += '_errp'
    c_errp = table.Column(**d)
    d['name'] += '_errn'
    c_errn = table.Column(**d)
    return [c_val, c_err, c_errp, c_errn]


def columns_from_derived(name, prop):
    """Build a set of `Column` objects for a `Derived` object
    """
    d = get_column_kwargs(name, prop)
    return [table.Column(**d)]


def make_columns_for_prop(name, prop):
    """Generic function to make a set of `Column` objects
    for any `Property` sub-class
    """
    if isinstance(prop, Derived):
        return columns_from_derived(name, prop)
    elif isinstance(prop, Parameter):
        return columns_from_derived(name, prop)
    elif isinstance(prop, Property):
        return columns_from_property(name, prop)
    else:
        raise TypeError("Can not make columns for %s" % type(prop))
    return None


def make_columns_for_model(names, model):
    """Make a set of `Column` objects needed to
    describe a `Model` object.
    """
    clist = []
    for n in names:
        p = model.getp(n)
        clist += make_columns_for_prop(n, p)
    return clist


def fill_table_from_targetlist(tab, targetList):
    """Fill a table for a set of `Target` objects
    """
    cnames = tab.colnames
    for t in targetList:
        tab.add_row(get_row_values(t, cnames))
    return


def make_table_for_targetlist(colNames, targetList):
    """Build a table for a set of `Target` objects
    """
    clist = make_columns_for_model(colNames, NullTarget)
    tab = table.Table(data=clist)
    fill_table_from_targetlist(tab, targetList)
    return tab


def make_table_for_roster(colNames, roster):
    """Make a table for a `Roster` object
    """
    clist = make_columns_for_model(colNames, NullTarget)
    tab = table.Table(data=clist)
    fill_table_from_targetlist(tab, roster.values())
    return tab


def row_to_dict(row):
    """Convert a `Table` row into a dict
    """
    d = {}
    for c in row.colnames:
        d[c] = row[c]
    return d


def make_target_from_row(row):
    """Make a `Target` object from a `Table` row
    """
    d = row_to_dict(row)
    d['abbr'] = row['name'][0:3]
    d['title'] = row['name']
    d['profile'] = {}
    t = Target(**d)
    return t
