#!/usr/bin/env python
"""
Module that handles file IO using `astropy.table.Table` class.
"""

from astropy import table

from pymodeler.parameter import Parameter, Derived, Property
from dmsky.targets import Target

NullTarget = Target(title='Null', name='Null', abbr='Nul', profile={})


def get_column_kwargs(name, prop):
    """Get keyword arguments needed to build a column for a `dmsky.Property` object

    Parameters
    ----------

    name : str
        Name of the Column we will build

    prop : `dmsky.Property`
        Property object we are building the column for


    Returns
    -------
    
    opt_keys : dict
        A dictionary we can use to constuct a `astropy.table.Column`

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

    Parameters
    ----------

    model : `dmsky.Model`
        Model that we are getting the data from

    keys : list
        Names of the properties we are reading


    Returns
    -------
    
    vals : dict
        A dictionary we can use to fill an `astropy.table.Column` row

    """
    d = {}
    for k in keys:
        d[k] = getattr(model, k)
    return d


def columns_from_property(name, prop):
    """Build some of `astropy.table.Column` objects for a `dmsky.Property` object

    In this case we only build a single column

    Parameters
    ----------

    name : str
        Name of the Column we will build

    prop : `dmsky.Property`
        Property object we are building the column for
   
    Returns
    -------
    
    cols : list
        A list of `astropy.table.Column` 

    """
    d = get_column_kwargs(name, prop)
    return [table.Column(**d)]


def columns_from_parameter(name, prop):
    """Build a set of `astropy.table.Column` objects for a `dmsky.Parameter` object

    In this case we build a several columns

    Parameters
    ----------

    name : str
        Base of the names of the Column we will build

    prop : `dmsky.Property`
        Property object we are building the column for
   
    Returns
    -------
    
    cols : list
        A list of `astropy.table.Column` 

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

    In this case we only build a single column

    Parameters
    ----------

    name : str
        Name of the Column we will build

    prop : `dmsky.Property`
        Property object we are building the column for
   
    Returns
    -------
    
    cols : list
        A list of `astropy.table.Column` 

    """
    d = get_column_kwargs(name, prop)
    return [table.Column(**d)]


def make_columns_for_prop(name, prop):
    """Generic function to make a set of `astropy.table.Column` objects
    for any `dmsky.Property` sub-class

    Parameters
    ----------

    name : str
        Name of the Column we will build

    prop : `dmsky.Property`
        Property object we are building the column for
   
    Returns
    -------
    
    cols : list
        A list of `astropy.table.Column` 

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

    Parameters
    ----------

    names : list
        List of the names of the properties to convert

    model : `dmsky.Model`
        Model that we are getting the data from


    Returns
    -------
    
    clist : list
        A list of `astropy.table.Column` 

    """
    clist = []
    for n in names:
        p = model.getp(n)
        clist += make_columns_for_prop(n, p)
    return clist


def fill_table_from_targetlist(tab, targetList):
    """Fill a table for a set of `dmsky.Target` objects

    Parameters
    ----------

    tab : `astropy.table.Table`
        The table we are filling

    targetList : list
        List of 'dmsky.Target' object used to fill the table

    """
    cnames = tab.colnames
    for t in targetList:
        tab.add_row(get_row_values(t, cnames))



def make_table_for_targetlist(colNames, targetList):
    """Build a table for a set of `Target` objects

    Parameters
    ----------

    colNames : list
        The names of the properties to include in the table

    targetList : list
        List of 'dmsky.Target' object used to fill the table


    Returns
    -------

    table : `astropy.table.Table`
        A table with the data from those targets

    """
    clist = make_columns_for_model(colNames, NullTarget)
    tab = table.Table(data=clist)
    fill_table_from_targetlist(tab, targetList)
    return tab


def make_table_for_roster(colNames, roster):
    """Make a table for a `Roster` object

    Parameters
    ----------

    colNames : list
        The names of the properties to include in the table

    roster : `dmsky.Roster`
        Roster used to fill the table


    Returns
    -------

    table : `astropy.table.Table`
        A table with the data from that Roster


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
