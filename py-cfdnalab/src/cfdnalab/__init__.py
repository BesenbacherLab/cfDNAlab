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
from .ref_kmers import (
    GlobalRefKmerFrequencies,
    GroupedRefKmerFrequencies,
    RefKmerFrequencies,
    WindowedRefKmerFrequencies,
    read_ref_kmers,
)


def get_version():
    import importlib.metadata

    return importlib.metadata.version("cfdnalab")


__version__ = get_version()


__all__ = [
    "EndMotifCounts",
    "GlobalEndMotifCounts",
    "GlobalLengthCounts",
    "GlobalRefKmerFrequencies",
    "GroupedEndMotifCounts",
    "GroupedLengthCounts",
    "GroupedRefKmerFrequencies",
    "LengthCounts",
    "MidpointProfiles",
    "RefKmerFrequencies",
    "WindowedEndMotifCounts",
    "WindowedLengthCounts",
    "WindowedRefKmerFrequencies",
    "__version__",
    "read_end_motifs",
    "read_lengths",
    "read_midpoints",
    "read_ref_kmers",
]
