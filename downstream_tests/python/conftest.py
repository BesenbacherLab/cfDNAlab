from __future__ import annotations

import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def midpoint_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_MIDPOINT_ZARR",
        "tiny.midpoint_profiles.zarr",
        "midpoint",
    )


@pytest.fixture(scope="session")
def dense_global_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_DENSE_GLOBAL_ZARR",
        "tiny_dense_global.end_motifs.zarr",
        "dense global end-motif",
    )


@pytest.fixture(scope="session")
def sparse_windowed_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_SPARSE_WINDOWED_ZARR",
        "tiny_sparse_windowed.end_motifs.zarr",
        "sparse windowed end-motif",
    )


@pytest.fixture(scope="session")
def sparse_grouped_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_SPARSE_GROUPED_ZARR",
        "tiny_sparse_grouped.end_motifs.zarr",
        "sparse grouped end-motif",
    )


@pytest.fixture(scope="session")
def global_length_counts_path() -> Path:
    return _fixture_path(
        "CFDNALAB_LENGTHS_GLOBAL_TSV",
        "tiny_lengths_global.length_counts.tsv.zst",
        "global length-count",
    )


@pytest.fixture(scope="session")
def windowed_length_counts_path() -> Path:
    return _fixture_path(
        "CFDNALAB_LENGTHS_WINDOWED_TSV",
        "tiny_lengths_windowed.length_counts.tsv.zst",
        "windowed length-count",
    )


@pytest.fixture(scope="session")
def grouped_length_counts_path() -> Path:
    return _fixture_path(
        "CFDNALAB_LENGTHS_GROUPED_TSV",
        "tiny_lengths_grouped.length_counts.tsv.zst",
        "grouped length-count",
    )


@pytest.fixture(scope="session")
def windowed_length_counts_no_blacklist_path() -> Path:
    return _fixture_path(
        "CFDNALAB_LENGTHS_WINDOWED_NO_BLACKLIST_TSV",
        "tiny_lengths_windowed_no_blacklist.length_counts.tsv.zst",
        "windowed length-count without blacklist",
    )


@pytest.fixture(scope="session")
def grouped_length_counts_no_blacklist_path() -> Path:
    return _fixture_path(
        "CFDNALAB_LENGTHS_GROUPED_NO_BLACKLIST_TSV",
        "tiny_lengths_grouped_no_blacklist.length_counts.tsv.zst",
        "grouped length-count without blacklist",
    )


def _fixture_path(env_var: str, default_name: str, label: str) -> Path:
    path = Path(
        os.environ.get(env_var, Path(__file__).parents[1] / "tmp" / default_name)
    )
    if not path.exists():
        pytest.fail(
            f"Missing cfDNAlab-generated {label} fixture at {path}. "
            "Generate it with the ignored downstream fixture integration tests."
        )
    return path
