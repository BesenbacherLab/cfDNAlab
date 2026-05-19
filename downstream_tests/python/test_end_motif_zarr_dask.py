from __future__ import annotations

from pathlib import Path

import dask.array as da
import numpy as np


def test_dask_reads_dense_end_motif_count_matrix(
    dense_global_end_zarr_path: Path,
) -> None:
    counts = da.from_zarr(str(dense_global_end_zarr_path), component="counts")

    assert counts.shape == (1, 4)
    assert counts.dtype == np.dtype("float64")
    np.testing.assert_allclose(
        counts.compute(),
        np.array([[1.0, 0.0, 1.0, 0.0]], dtype=np.float64),
    )


def test_dask_reads_sparse_end_motif_coordinate_arrays(
    sparse_windowed_end_zarr_path: Path,
) -> None:
    sparse_shape = da.from_zarr(
        str(sparse_windowed_end_zarr_path),
        component="sparse/shape",
    )
    sparse_row = da.from_zarr(
        str(sparse_windowed_end_zarr_path),
        component="sparse/row",
    )
    sparse_motif = da.from_zarr(
        str(sparse_windowed_end_zarr_path),
        component="sparse/motif",
    )
    sparse_count = da.from_zarr(
        str(sparse_windowed_end_zarr_path),
        component="sparse/count",
    )

    np.testing.assert_array_equal(
        sparse_shape.compute(),
        np.array([3, 2], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        sparse_row.compute(),
        np.array([0, 1, 2], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        sparse_motif.compute(),
        np.array([1, 0, 1], dtype=np.int32),
    )
    np.testing.assert_allclose(
        sparse_count.compute(),
        np.array([1.0, 1.0, 1.0], dtype=np.float64),
    )
