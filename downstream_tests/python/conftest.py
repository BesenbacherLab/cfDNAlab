from __future__ import annotations

import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def midpoint_zarr_path() -> Path:
    path = Path(
        os.environ.get(
            "CFDNALAB_MIDPOINT_ZARR",
            Path(__file__).parents[1] / "tmp" / "tiny.midpoint_profiles.zarr",
        )
    )
    if not path.exists():
        pytest.fail(
            "Missing cfDNAlab-generated midpoint Zarr fixture at "
            f"{path}. Generate it with the ignored "
            "`generate_midpoint_zarr_fixture_with_cfdnalab` integration test."
        )
    return path
