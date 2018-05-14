#!/usr/bin/env python
"""
Utilities for speed testing
"""
import time


def speedtest(func, *args, **kwargs):
    """ Test the speed of a function. """
    n = 100
    start = time.time()
    i = 0
    while i < n:
        i += 1
        func(*args, **kwargs)
    end = time.time()
    return (end - start) / n

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage, description=description)
    (opts, args) = parser.parse_args()
