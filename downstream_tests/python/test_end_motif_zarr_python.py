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

    sparse_row = store["sparse/row"]
    sparse_motif = store["sparse/motif"]
    sparse_count = store["sparse/count"]
    sparse_shape = store["sparse/shape"]
    assert tuple(sparse_row.metadata.dimension_names) == ("nnz",)
    assert tuple(sparse_motif.metadata.dimension_names) == ("nnz",)
    assert tuple(sparse_count.metadata.dimension_names) == ("nnz",)
    assert tuple(sparse_shape.metadata.dimension_names) == ("sparse_dimension",)
    np.testing.assert_array_equal(sparse_shape[:], np.array([3, 2], dtype=np.int32))
    np.testing.assert_array_equal(sparse_row[:], np.array([0, 0, 1], dtype=np.int32))
    np.testing.assert_array_equal(sparse_motif[:], np.array([0, 1, 0], dtype=np.int32))
    np.testing.assert_allclose(sparse_count[:], np.array([1.0, 2.0, 1.0]))

    dense = np.zeros(tuple(sparse_shape[:].astype(int)), dtype=np.float64)
    dense[sparse_row[:].astype(int), sparse_motif[:].astype(int)] = sparse_count[:]
    np.testing.assert_allclose(
        dense,
        np.array([[1.0, 2.0], [1.0, 0.0], [0.0, 0.0]], dtype=np.float64),
    )
