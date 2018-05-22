.. _install:

Installation
============

.. _pipinstall:

Installing with pip
-------------------

These instructions cover installation with the ``pip`` package
management tool.  This will install fermipy and its dependencies into
your python distribution.

Before starting the installation process, you will need to determine
whether you have setuptools and pip installed in your local python
environment.  You may need to install these packages if you are
running with the binary version of the Fermi Science Tools distributed
by the FSSC.  The following command will install both packages in your
local environment:

.. code-block:: bash

   $ curl https://bootstrap.pypa.io/get-pip.py | python -

Check if pip is correctly installed:

.. code-block:: bash

   $ which pip

Once again, if this isn't the pip in your python environment something went wrong.
Now install fermipy by running:

.. code-block:: bash

   $ pip install dmsky

To run the ipython notebook examples you will also need to install
jupyter notebook:
   
.. code-block:: bash

   $ pip install jupyter

.. Running pip and setup.py with the ``user`` flag is recommended if you do not
.. have write access to your python installation (for instance if you are
.. running in a UNIX/Linux environment with a shared python
.. installation).  To install fermipy into the common package directory
.. of your python installation the ``user`` flag should be ommitted.

Finally, check that dmsky imports:

.. code-block:: bash

   $ python
   Python 2.7.8 (default, Aug 20 2015, 11:36:15)
   [GCC 4.2.1 Compatible Apple LLVM 6.0 (clang-600.0.56)] on darwin
   Type "help", "copyright", "credits" or "license" for more information. 
   >>> import dmsky
   >>> dmsky.__file__


The instructions describe how to install development versions of
Fermipy.  Before installing a development version we recommend first
installing a tagged release following the :ref:`pipinstall` 
instructions above.

The development version of Fermipy can be installed by running ``pip
install`` with the URL of the git repository:

.. code-block:: bash
                
   $ pip install git+https://github.com/fermiPy/dmsky.git

This will install the most recent commit on the master branch.  Note
that care should be taken when using development versions as
features/APIs under active development may change in subsequent
versions without notice.
   
   
Upgrading
---------

By default installing fermipy with ``pip`` or ``conda`` will get the latest tagged
released available on the `PyPi <https://pypi.python.org/pypi>`_
package respository.  You can check your currently installed version
of fermipy with ``pip show``:

.. code-block:: bash

   $ pip show fermipy

or ``conda info``:

.. code-block:: bash

   $ conda info fermipy
   
To upgrade your fermipy installation to the latest version run the pip
installation command with ``--upgrade --no-deps`` (remember to also
include the ``--user`` option if you're running at SLAC):
   
.. code-block:: bash
   
   $ pip install fermipy --upgrade --no-deps
   Collecting fermipy
   Installing collected packages: fermipy
     Found existing installation: fermipy 0.6.6
       Uninstalling fermipy-0.6.6:
         Successfully uninstalled fermipy-0.6.6
   Successfully installed fermipy-0.6.7

If you installed fermipy with ``conda`` the equivalent command is:

.. code-block:: bash

   $ conda update fermipy
   
   
.. _gitinstall:
   
Developer Installation
----------------------

These instructions describe how to install fermipy from its git source
code repository using the ``setup.py`` script.  Installing from source
can be useful if you want to make your own modifications to the
fermipy source code.  Note that non-developers are recommended to
install a tagged release of fermipy following the :ref:`pipinstall` or
instructions above.

First clone the fermipy git repository and cd to the root directory of
the repository:

.. code-block:: bash

   $ git clone https://github.com/fermiPy/fermipy.git
   $ cd fermipy
   
To install the latest commit in the master branch run ``setup.py
install`` from the root directory:

.. code-block:: bash

   # Install the latest commit
   $ git checkout master
   $ python setup.py install --user 

A useful option if you are doing active code development is to install
your working copy of the package.  This will create an installation in
your python distribution that is linked to the copy of the code in
your local repository.  This allows you to run with any local
modifications without having to reinstall the package each time you
make a change.  To install your working copy of fermipy run with the
``develop`` argument:

.. code-block:: bash

   # Install a link to your source code installation
   $ python setup.py develop --user 

You can later remove the link to your working copy by running the same
command with the ``--uninstall`` flag:

.. code-block:: bash

   # Install a link to your source code installation
   $ python setup.py develop --user --uninstall
   

Specific release tags can be installed by running ``git checkout``
before running the installation command:
   
.. code-block:: bash
   
   # Checkout a specific release tag
   $ git checkout X.X.X 
   $ python setup.py install --user 

To see the list of available release tags run ``git tag``.
   
Issues
------

If you get an error about importing matplotlib (specifically something
about the macosx backend) you might change your default backend to get
it working.  The `customizing matplotlib page
<http://matplotlib.org/users/customizing.html>`_ details the
instructions to modify your default matplotlibrc file (you can pick
GTK or WX as an alternative).  Specifically the ``TkAgg`` and
``macosx`` backends currently do not work on OSX if you upgrade
matplotlib to the version required by fermipy.  To get around this
issue you can switch to the ``Agg`` backend at runtime before
importing fermipy:

.. code-block:: bash

   >>> import matplotlib
   >>> matplotlib.use('Agg')

However note that this backend does not support interactive plotting.

If you are running OSX El Capitan or newer you may see errors like the following:

.. code-block:: bash
                
   dyld: Library not loaded

In this case you will need to disable the System Integrity Protections
(SIP).  See `here
<http://www.macworld.com/article/2986118/security/how-to-modify-system-integrity-protection-in-el-capitan.html>`_
for instructions on disabling SIP on your machine.

In some cases the setup.py script will fail to properly install the
fermipy package dependecies.  If installation fails you can try
running a forced upgrade of these packages with ``pip install --upgrade``:

.. code-block:: bash

   $ pip install --upgrade --user numpy matplotlib scipy astropy pyyaml healpy wcsaxes ipython jupyter
