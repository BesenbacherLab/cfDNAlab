from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def test_xarray_opens_dense_end_motif_counts_with_named_dimensions(
    dense_global_end_zarr_path: Path,
) -> None:
    dataset = xr.open_zarr(str(dense_global_end_zarr_path), consolidated=False)

    assert dataset.attrs["cfdnalab_schema"] == "end_motif_counts"
    assert dataset.attrs["primary_array"] == "counts"
    assert dataset["counts"].dims == ("row", "motif")
    assert dataset["counts"].shape == (1, 4)
    assert "row" in dataset
    assert "motif_index" in dataset
    np.testing.assert_allclose(
        dataset["counts"].values,
        np.array([[1.0, 0.0, 1.0, 0.0]], dtype=np.float64),
    )


def test_xarray_builds_dataframe_from_dense_end_motif_counts(
    dense_global_end_zarr_path: Path,
) -> None:
    dataset = xr.open_zarr(str(dense_global_end_zarr_path), consolidated=False)
    motif_labels = [
        bytes(row).decode("ascii")
        for row in dataset["motif_ascii"].values.astype(np.uint8)
    ]

    frame = pd.DataFrame(
        {
            "row_label": dataset["row"].attrs["labels"][0],
            "motif_index": dataset["motif_index"].values,
            "motif": motif_labels,
            "count": dataset["counts"].isel(row=0).values,
        }
    )

    pd.testing.assert_frame_equal(
        frame,
        pd.DataFrame(
            {
                "row_label": np.array(["global", "global", "global", "global"], dtype=object),
                "motif_index": np.array([0, 1, 2, 3], dtype=np.int32),
                "motif": np.array(["_A", "_C", "_G", "_T"], dtype=object),
                "count": np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64),
            }
        ),
    )


def test_xarray_reads_selected_end_motif_labels_and_sparse_coordinates(
    sparse_windowed_selected_motifs_end_zarr_path: Path,
) -> None:
    dataset = xr.open_zarr(
        str(sparse_windowed_selected_motifs_end_zarr_path),
        consolidated=False,
    )
    sparse = xr.open_zarr(
        str(sparse_windowed_selected_motifs_end_zarr_path),
        group="sparse",
        consolidated=False,
    )
    motif_labels = [
        bytes(row).decode("ascii")
        for row in dataset["motif_ascii"].values.astype(np.uint8)
    ]

    assert dataset.attrs["motif_axis_kind"] == "motif"
    assert motif_labels == ["GT_AC", "AC_GT"]
    assert dataset["motif_ascii"].dims == ("motif", "motif_byte")
    np.testing.assert_array_equal(
        sparse["shape"].values,
        np.array([3, 2], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        sparse["motif"].values,
        np.array([1, 0, 1], dtype=np.int32),
    )
