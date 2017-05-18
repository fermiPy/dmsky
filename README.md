Dark Matter Sky Maps
====================

[![Build](https://img.shields.io/travis/fermiPy/dmsky.svg)](https://travis-ci.org/fermiPy/dmsky)
[![Release](https://img.shields.io/github/release/fermiPy/dmsky.svg)](../../releases)
[![PyPI](https://img.shields.io/pypi/v/dmsky.svg)](https://pypi.python.org/pypi/dmsky)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](../../)

Introduction
------------
Tool for generating dark matter skymaps.

Installation
------------

The easiest way to install ``dmsky`` is with ``pip``:

```
# For an initial installation
pip install dmsky

# update just dmsky
pip install dmsky --no-deps --upgrade --ignore-installed

# update dmsky and all dependencies
pip install dmsky --upgrade
```

You can also download the latest release or the bleeding edge of the source code from github:

```
git clone https://github.com/kadrlica/dmsky.git
cd dmsky
python setup.py install
```

Usage
------------
The hope is that ``dmsky`` will provide a simple and transparent interface towards making dark matter sky maps. Functionality is provided through the Python API, though there should eventually be a simple command-line interface. Several examples can be found in the [examples](examples/) and [tests](tests/) directories.

Want to Help?
-------------

Help is always appreciated and contributions will be accepted at all levels (just like NPR). Some potential  ways to get involved:

* Use the code and provide feedback
* Contribute a configuration file of your favorite dark matter source
* Help code up some missing functionality
