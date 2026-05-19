"""Python helpers for loading cfDNAlab output files."""

from .ends import (
    EndMotifCounts,
    GlobalEndMotifCounts,
    GroupedEndMotifCounts,
    WindowedEndMotifCounts,
    read_end_motifs,
)
from .lengths import (
    GlobalLengthCounts,
    GroupedLengthCounts,
    LengthCounts,
    WindowedLengthCounts,
    read_lengths,
)
from .midpoints import MidpointProfiles, read_midpoints


def get_version():
    import importlib.metadata

    return importlib.metadata.version("cfdnalab")


__version__ = get_version()


__all__ = [
    "EndMotifCounts",
    "GlobalEndMotifCounts",
    "GlobalLengthCounts",
    "GroupedEndMotifCounts",
    "GroupedLengthCounts",
    "LengthCounts",
    "MidpointProfiles",
    "WindowedEndMotifCounts",
    "WindowedLengthCounts",
    "__version__",
    "read_end_motifs",
    "read_lengths",
    "read_midpoints",
]
