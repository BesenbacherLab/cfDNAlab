from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import zarr


EXPECTED_COUNTS = np.array(
    [
        [
            [1.0, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.5, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 2.0],
        ],
        [
            [0.5, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.5, 0.0, 0.5],
            [0.0, 0.5, 0.0, 1.0, 0.0],
        ],
        [
            [0.0, 0.0, 2.5, 0.0, 0.0],
            [0.5, 0.0, 0.0, 1.5, 0.0],
            [0.0, 0.0, 0.5, 0.0, 1.0],
        ],
    ],
    dtype=np.float32,
)
EXPECTED_ARRAYS = {
    "counts",
    "group",
    "eligible_intervals",
    "length_bin",
    "length_start_bp",
    "length_end_bp",
    "position",
    "position_bin_start_bp",
    "position_bin_end_bp",
}


def test_python_zarr_reads_midpoint_profile_schema(midpoint_zarr_path: Path) -> None:
    store = zarr.open_group(str(midpoint_zarr_path), mode="r", zarr_format=3)

    assert store.attrs["cfdnalab_schema"] == "midpoint_profiles"
    assert store.attrs["cfdnalab_schema_version"] == 1
    assert store.attrs["primary_array"] == "counts"
    assert store.attrs["count_units"] == "weighted_midpoint_count"
    assert set(store.array_keys()) == EXPECTED_ARRAYS
    assert store["counts"].shape == (3, 3, 5)
    assert store["counts"].dtype == np.dtype("float32")
    assert tuple(store["counts"].metadata.dimension_names) == (
        "group",
        "length_bin",
        "position",
    )
    assert store["group"].attrs["label_field"] == "group_name"
    assert store["group"].attrs["labels"] == ["LYL1", "beta-site", "gamma_long"]
    np.testing.assert_allclose(store["counts"][:], EXPECTED_COUNTS)
    np.testing.assert_array_equal(store["group"][:], np.array([0, 1, 2], dtype=np.int32))
    np.testing.assert_array_equal(
        store["length_bin"][:],
        np.array([0, 1, 2], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        store["position"][:],
        np.array([0, 1, 2, 3, 4], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        store["eligible_intervals"][:],
        np.array([2, 2, 2], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        store["length_start_bp"][:],
        np.array([30, 50, 70], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        store["length_end_bp"][:],
        np.array([50, 70, 100], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        store["position_bin_start_bp"][:],
        np.array([0, 2, 4, 6, 8], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        store["position_bin_end_bp"][:],
        np.array([2, 4, 6, 8, 10], dtype=np.int32),
    )


def test_python_zarr_builds_plotting_dataframe_for_one_group(
    midpoint_zarr_path: Path,
) -> None:
    store = zarr.open_group(str(midpoint_zarr_path), mode="r", zarr_format=3)
    group_index = 1
    group_names = store["group"].attrs["labels"]
    profile = store["counts"][group_index, :, :]
    length_index, position_index = np.indices(profile.shape)

    frame = pd.DataFrame(
        {
            "group_name": group_names[group_index],
            "eligible_intervals": int(store["eligible_intervals"][group_index]),
            "length_start_bp": store["length_start_bp"][:][length_index.ravel()],
            "length_end_bp": store["length_end_bp"][:][length_index.ravel()],
            "position_bin_start_bp": store["position_bin_start_bp"][:][
                position_index.ravel()
            ],
            "position_bin_end_bp": store["position_bin_end_bp"][:][
                position_index.ravel()
            ],
            "count": profile.ravel(),
        }
    )

    assert list(frame.columns) == [
        "group_name",
        "eligible_intervals",
        "length_start_bp",
        "length_end_bp",
        "position_bin_start_bp",
        "position_bin_end_bp",
        "count",
    ]
    assert frame.shape == (15, 7)
    assert frame["count"].tolist() == [
        0.5,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.5,
        0.0,
        0.5,
        0.0,
        0.5,
        0.0,
        1.0,
        0.0,
    ]
