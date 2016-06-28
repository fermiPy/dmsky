#!/usr/bin/env python
"""
Module file IO using astropy.table class.
"""

from pymodeler.parameter import *
from astropy import table
from dmsky.targets import Target

NullTarget = Target(title='Null',name='Null',abbr='Nul',profile={})

def get_column_kwargs(name,prop):
    opt_keys = ['data','unit','shape','length']
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
            d[k] = getattr(prop,k)
        except AttributeError:
            pass
    # FIXME, force strings to 16 characters
    if d['dtype'] == type(""):
        d['dtype'] =  'S16' 
    return d

def get_row_values(model,keys):
    d = {}
    for k in keys:
        d[k] = getattr(model,k)
    return d

def columns_from_property(name,prop):
    d = get_column_kwargs(name,prop)
    return [table.Column(**d)]

def columns_from_parameter(name,prop):
    d = get_column_kwargs(name,prop)
    c_val = table.Column(**d)
    d['name'] += '_err'
    c_err = table.Column(**d)
    d['name'] += '_errp'
    c_errp = table.Column(**d)
    d['name'] += '_errn'
    c_errn = table.Column(**d)
    return [c_val,c_err,c_errp,c_errn]

def columns_from_derived(name,prop):
    d = get_column_kwargs(name,prop)
    return [table.Column(**d)]

def make_columns_for_prop(name,prop):
    if isinstance(prop,Derived):
        return columns_from_derived(name,prop)
    elif isinstance(prop,Parameter):
        return columns_from_derived(name,prop)
    elif isinstance(prop,Property):
        return columns_from_property(name,prop)

def make_columns_for_model(names,model):
    clist = []
    for n in names:
        p = model.getp(n)
        clist += make_columns_for_prop(n,p)
    return clist

def fill_table_from_targetlist(tab,targetList):
    cnames = tab.colnames
    for t in targetList:
        tab.add_row( get_row_values(t,cnames) )
    return 


def make_table_for_targetlist(colNames,targetList):
    clist = make_columns_for_model(colNames,NullTarget)
    tab = table.Table(data=clist)
    fill_table_from_targetlist(tab,targetList)
    return tab

def make_table_for_roster(colNames,roster):
    clist = make_columns_for_model(colNames,NullTarget)
    tab = table.Table(data=clist)
    fill_table_from_targetlist(tab,roster.values())
    return tab


def row_to_dict(row):
    d = {}
    for c in row.colnames:
        d[c] = row[c]
    return d


def make_target_from_row(row):
    d = row_to_dict(row)
    d['abbr'] = row['name'][0:3]
    d['title'] = row['name']
    d['profile'] = {}
    t = Target(**d)
    return t

