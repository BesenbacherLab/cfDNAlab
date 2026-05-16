import cfdnalab


def test_public_midpoint_loader_is_exported() -> None:
    assert callable(cfdnalab.load_midpoints)
