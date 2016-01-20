#!/usr/bin/env python
"""
"""
import os
import yaml
import numpy as np
import copy

def yaml_load(filename):
    """ Load a yaml file (use libyaml when available) """
    if not os.path.exists(filename):
        raise Exception('File does not exist: %s'%(filename))
    try: 
        ret = yaml.load(open(filename),Loader=yaml.CLoader)
    except:
        ret = yaml.load(open(filename),Loader=yaml.Loader)
    return ret
    
def yaml_dump(x, filename):    
    """ Dump object to a yaml file (use libyaml when available) 
        x        : output to dump to the file
        filename : output file (can be file-type or path string)
    """
    if isinstance(filename, basestring):
        out = open(filename,'w')
    elif isinstance(filename, file):
        out = filename
    else:
        raise Exception("Unrecognized file: ",filename)
 
    try:
        out.write(yaml.dump(x,Dumper=yaml.CDumper))
    except:
        out.write(yaml.dump(x,Dumper=yaml.Dumper))
    out.close()

def update_dict(d0,d1,add_keys=True,append=False):
    """Recursively update the contents of dictionary d0 with
    the contents of python dictionary d1.

    Input:
    ------
    d0: Base dictionary
    d1: Update dictionary
    
    Returns:
    -------
    None
    """
    # Also see:
    # http://stackoverflow.com/questions/3232943/

    if d1 is None: return
    if d0 is None:
        d0 = copy.copy(d1)
        return
    
    for k, v in d0.iteritems():
        if not k in d1: continue

        if isinstance(v,dict) and isinstance(d1[k],dict):
            update_dict(d0[k],d1[k],add_keys,append)
        elif isinstance(v,list) and isinstance(d1[k],str):
            d0[k] = d1[k].split(',')            
        elif isinstance(v,np.ndarray) and append:
            d0[k] = np.concatenate((v,d1[k]))
        else: d0[k] = d1[k]

    if add_keys:
        for k, v in d1.iteritems(): 
            if not k in d0: d0[k] = d1[k]

def merge_dict(d0,d1,add_keys=True,append=False):
    """ 
    Merge two target dicts into a new dict.

    Input:
    ------
    d0: Base dictionary
    d1: Update dictionary

    return:
    d: A new dictionary merging d and u
    """
    d = copy.copy(d0)
    update_dict(d,d1,add_keys,append)
    return d

def getnest(d, *keys):
    """ Function to return nested keys.
    """
    ret = copy.copy(d)
    for arg in args:
        ret = ret[key]
    return ret

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
