import sys
import os
try: 
    from setuptools import setup, find_packages
except ImportError: 
    from distutils.core import setup
    def find_packages():
        return ['dmsky']

from dmsky.get_version import get_version, write_version_py

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
VERSION = get_version()
write_version_py(version=VERSION)

def read(filename):
    return open(os.path.join(HERE,filename)).read()

setup(
    name=NAME,
    version=VERSION,
    url='',
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = [],
    install_requires=[
        'python >= 2.7.0',
        'numpy >= 1.9.0',
        'scipy >= 0.14.0',
        'healpy >= 1.6.0',
        'pyfits >= 3.1',
        'pyyaml >= 3.10',
    ],
    packages=find_packages(),
    package_data={'dmsky': ['data']},
    description="Dark matter skymaps.",
    long_description=read('README.md'),
    platforms='any',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f]
)