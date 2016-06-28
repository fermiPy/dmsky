#!/usr/bin/env python
"""
Factory for generating instances of classes
"""
import sys
from collections import OrderedDict as odict
import inspect

def factory(type, module=None, **kwargs):
    """
    Factory for creating objects. Arguments are passed directly to the
    constructor of the chosen class.
    """
    cls = type
    if module is None: module = __name__
    fn = lambda member: inspect.isclass(member) and member.__module__==module
    classes = odict(inspect.getmembers(sys.modules[module], fn))
    members = odict([(k.lower(),v) for k,v in classes.items()])
    
    lower = cls.lower()
    if lower not in members.keys():
        #msg = "%s not found in:\n %s"%(cls,classes.keys())
        #logging.error(msg)
        msg = "Unrecognized class: %s"%cls
        raise Exception(msg)

    return members[lower](**kwargs)


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
