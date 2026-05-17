from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def test_xarray_opens_counts_with_named_dimensions(midpoint_zarr_path: Path) -> None:
    dataset = xr.open_zarr(str(midpoint_zarr_path), consolidated=False)

    assert dataset["counts"].dims == ("group", "length_bin", "position")
    assert dataset["counts"].shape == (3, 3, 5)
    assert "group" in dataset
    assert "length_bin" in dataset
    assert "position" in dataset


def test_xarray_builds_dataframe_for_one_group(midpoint_zarr_path: Path) -> None:
    dataset = xr.open_zarr(str(midpoint_zarr_path), consolidated=False)
    group_index = 0
    group_names = dataset["group"].attrs["labels"]
    profile = dataset["counts"].isel(group=group_index).values
    length_index, position_index = np.indices(profile.shape)

    frame = pd.DataFrame(
        {
            "group_name": group_names[group_index],
            "eligible_intervals": int(
                dataset["eligible_intervals"].isel(group=group_index).values
            ),
            "length_start_bp": dataset["length_start_bp"].values[length_index.ravel()],
            "length_end_bp": dataset["length_end_bp"].values[length_index.ravel()],
            "position_bin_start_bp": dataset["position_bin_start_bp"].values[
                position_index.ravel()
            ],
            "position_bin_end_bp": dataset["position_bin_end_bp"].values[
                position_index.ravel()
            ],
            "count": profile.ravel(),
        }
    )

    assert frame.shape == (15, 7)
    assert frame["group_name"].unique().tolist() == ["LYL1"]
    assert frame["count"].tolist() == [
        1.0,
        0.5,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.5,
        0.5,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        2.0,
    ]
