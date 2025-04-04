import contextlib
from importlib.metadata import PackageNotFoundError, version

from . import protocol, stats, utils, viz

with contextlib.suppress(PackageNotFoundError):
    __version__ = version("core")
