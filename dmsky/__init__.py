"""
Dark matter skymaps.
"""
__author__ = "Alex Drlica-Wagner"
__email__ = "kadrlica@fnal.gov"


# Automatically grab the version from the git tag to ensure that the
# tag and the code agree. For more details, see `get_version.py`
try:
    from .version import __version__
except ImportError:
    from .get_version import get_version
    __version__ = get_version()
