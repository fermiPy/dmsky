import sys
import os

from setuptools import setup, find_packages
import versioneer

NAME = 'dmsky'
CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 2.7
Natural Language :: English
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Astronomy
Operating System :: MacOS
Operating System :: POSIX
License :: OSI Approved :: MIT License
"""
URL = 'https://github.com/kadrlica/dmsky'
DESC = "Dark matter skymaps."
LONG_DESC = "See %s"%URL

if sys.version_info[:2] < (2, 7):
    raise RuntimeError("Python version >= 2.7 required.")

def find_data_files(datadir=None):
    """
    http://stackoverflow.com/a/13629066/4075339
    """
    # Copy the data files
    if not datadir: datadir = os.path.join(NAME,'data')
    data_files = [(d, [os.path.join(d,f) for f in files])
                 for d, folders, files in os.walk(datadir)]
    return data_files

def find_package_data(datadir=None):
    """
    http://stackoverflow.com/a/36693250/4075339
    """
    if not datadir: datadir = os.path.join(NAME,'data')
    paths = []
    for (d, folders, files) in os.walk(datadir):
        for f in files:
            paths.append(os.path.join('..', d, f))
    package_data = {'':paths}
    return package_data

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url=URL,
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = [],
    install_requires=[
        'setuptools',
        'numpy >= 1.9.0',
        'scipy >= 0.14.0',
        'healpy >= 1.6.0',
        'pyyaml >= 3.10',
        'pymodeler >= 0.1.0',
    ],
    packages=find_packages(),
    #data_files=find_data_files(),
    package_data=find_package_data(),
    include_package_data=True,
    description=DESC,
    long_description=LONG_DESC,
    platforms='any',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f]
)
