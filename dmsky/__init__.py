"""
Dark matter skymaps.
"""
__author__ = "Alex Drlica-Wagner"
__email__ = "kadrlica@fnal.gov"

# Automatically grab the version from the git tag to ensure that the
# tag and the code agree. For more details, see `python-versioneer`

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
