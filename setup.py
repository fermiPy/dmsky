import sys
import os
try: 
    from setuptools import setup, find_packages
except ImportError: 
    from distutils.core import setup
    def find_packages():
        return ['dmsky']

import versioneer

NAME = 'dmsky'
HERE = os.path.abspath(os.path.dirname(__file__))
CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Python
Natural Language :: English
Topic :: Scientific/Engineering
"""
URL = 'https://github.com/kadrlica/dmsky'
DESC = "Dark matter skymaps."
LONG_DESC = "See %s"%URL

def read(filename):
    return open(os.path.join(HERE,filename)).read()

def find_data():
    # Copy the data files
    datadir = os.path.join('dmsky','data')
    datafiles = [(d, [os.path.join(d,f) for f in files])
                 for d, folders, files in os.walk(datadir)]
    return datafiles

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url=URL,
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = [],
    install_requires=[
        'python >= 2.7.0',
        'numpy >= 1.9.0',
        'scipy >= 0.14.0',
        'healpy >= 1.6.0',
        'pyyaml >= 3.10',
        'pymodeler >= 0.1.0',
    ],
    packages=find_packages(),
    data_files=find_data(),
    description=DESC,
    long_description=LONG_DESC,
    platforms='any',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f]
)
