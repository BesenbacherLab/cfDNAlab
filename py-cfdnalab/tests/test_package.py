from importlib.metadata import version

import cfdnalab


def test_package_version_is_exported() -> None:
    assert cfdnalab.__version__ == version("cfdnalab")


def test_public_midpoint_loader_is_exported() -> None:
    assert callable(cfdnalab.read_midpoints)


def test_public_end_motif_loader_is_exported() -> None:
    assert callable(cfdnalab.read_end_motifs)


def test_public_end_motif_classes_are_exported() -> None:
    assert cfdnalab.GlobalEndMotifCounts
    assert cfdnalab.WindowedEndMotifCounts
    assert cfdnalab.GroupedEndMotifCounts
