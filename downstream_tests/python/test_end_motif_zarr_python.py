from __future__ import annotations

from pathlib import Path

import numpy as np
import zarr


def decode_motifs(store: zarr.Group) -> list[str]:
    motif_ascii = store["motif_ascii"][:]
    return [bytes(row).decode("ascii") for row in motif_ascii]


def test_python_zarr_reads_dense_global_end_motif_schema(
    dense_global_end_zarr_path: Path,
) -> None:
    store = zarr.open_group(str(dense_global_end_zarr_path), mode="r", zarr_format=3)

    assert store.attrs["cfdnalab_schema"] == "end_motif_counts"
    assert store.attrs["storage_mode"] == "dense"
    assert store.attrs["row_mode"] == "global"
    assert store.attrs["primary_array"] == "counts"
    assert store.attrs["primary_group"] is None
    assert tuple(store["counts"].metadata.dimension_names) == ("row", "motif")
    assert decode_motifs(store) == ["_A", "_C", "_G", "_T"]
    np.testing.assert_allclose(
        store["counts"][:],
        np.array([[1.0, 0.0, 1.0, 0.0]], dtype=np.float64),
    )
    assert store["row"].attrs["labels"] == ["global"]


def test_python_zarr_reads_sparse_grouped_end_motif_schema(
    sparse_grouped_end_zarr_path: Path,
) -> None:
    store = zarr.open_group(str(sparse_grouped_end_zarr_path), mode="r", zarr_format=3)

    assert store.attrs["cfdnalab_schema"] == "end_motif_counts"
    assert store.attrs["storage_mode"] == "sparse_coo"
    assert store.attrs["row_mode"] == "grouped_bed"
    assert store.attrs["primary_array"] is None
    assert store.attrs["primary_group"] == "sparse"
    assert store["group"].attrs["labels"] == ["beta", "alpha", "gamma"]
    assert decode_motifs(store) == ["_A", "_G"]

    sparse_group = store["sparse"]
    assert tuple(sparse_group["row"].metadata.dimension_names) == ("nnz",)
    assert tuple(sparse_group["motif"].metadata.dimension_names) == ("nnz",)
    assert tuple(sparse_group["count"].metadata.dimension_names) == ("nnz",)
    np.testing.assert_array_equal(sparse_group["shape"][:], np.array([3, 2]))
    np.testing.assert_array_equal(sparse_group["row"][:], np.array([0, 0, 1]))
    np.testing.assert_array_equal(sparse_group["motif"][:], np.array([0, 1, 0]))
    np.testing.assert_allclose(sparse_group["count"][:], np.array([1.0, 2.0, 1.0]))

    dense = np.zeros(tuple(sparse_group["shape"][:].astype(int)), dtype=np.float64)
    dense[sparse_group["row"][:].astype(int), sparse_group["motif"][:].astype(int)] = (
        sparse_group["count"][:]
    )
    np.testing.assert_allclose(
        dense,
        np.array([[1.0, 2.0], [1.0, 0.0], [0.0, 0.0]], dtype=np.float64),
    )
