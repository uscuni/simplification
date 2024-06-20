import contextlib
from importlib.metadata import PackageNotFoundError, version

from . import geometry, stats, utils, viz

with contextlib.suppress(PackageNotFoundError):
    __version__ = version("core")
