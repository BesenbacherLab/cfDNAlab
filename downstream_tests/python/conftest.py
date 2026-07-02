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
def sparse_windowed_selected_motifs_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_SPARSE_WINDOWED_SELECTED_MOTIFS_ZARR",
        "tiny_sparse_windowed_selected_motifs.end_motifs.zarr",
        "sparse windowed selected-motifs end-motif",
    )


@pytest.fixture(scope="session")
def sparse_grouped_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_SPARSE_GROUPED_ZARR",
        "tiny_sparse_grouped.end_motifs.zarr",
        "sparse grouped end-motif",
    )


@pytest.fixture(scope="session")
def sparse_grouped_motif_group_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_SPARSE_GROUPED_MOTIF_GROUP_ZARR",
        "tiny_sparse_grouped_motif_groups.end_motifs.zarr",
        "sparse grouped motif-group end-motif",
    )


@pytest.fixture(scope="session")
def sparse_grouped_wide_motif_group_end_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_ENDS_SPARSE_GROUPED_WIDE_MOTIF_GROUP_ZARR",
        "tiny_sparse_grouped_wide_motif_groups.end_motifs.zarr",
        "sparse grouped wide motif-group end-motif",
    )


@pytest.fixture(scope="session")
def dense_global_ref_kmer_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_REF_KMERS_DENSE_GLOBAL_ZARR",
        "tiny_ref_kmers_dense_global.ref_kmer_counts.zarr",
        "dense global reference k-mer",
    )


@pytest.fixture(scope="session")
def sparse_windowed_ref_kmer_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_REF_KMERS_SPARSE_WINDOWED_ZARR",
        "tiny_ref_kmers_sparse_windowed.ref_kmer_counts.zarr",
        "sparse windowed reference k-mer",
    )


@pytest.fixture(scope="session")
def sparse_grouped_ref_kmer_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_REF_KMERS_SPARSE_GROUPED_ZARR",
        "tiny_ref_kmers_sparse_grouped.ref_kmer_counts.zarr",
        "sparse grouped reference k-mer",
    )


@pytest.fixture(scope="session")
def dense_grouped_motif_group_ref_kmer_zarr_path() -> Path:
    return _fixture_path(
        "CFDNALAB_REF_KMERS_DENSE_GROUPED_MOTIF_GROUP_ZARR",
        "tiny_ref_kmers_dense_grouped_motif_groups.ref_kmer_counts.zarr",
        "dense grouped motif-group reference k-mer",
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
