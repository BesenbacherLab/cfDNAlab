from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import cfdnalab

from test_midpoint_zarr_python import EXPECTED_COUNTS


def test_cfdnalab_package_reads_midpoint_fixture_metadata(
    midpoint_zarr_path: Path,
) -> None:
    profiles = cfdnalab.read_midpoints(midpoint_zarr_path)

    pd.testing.assert_frame_equal(
        profiles.group_metadata(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2], dtype=np.int32),
                "group_name": np.array(["LYL1", "beta-site", "gamma_long"], dtype=str),
                "eligible_intervals": np.array([2, 2, 2], dtype=np.int32),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        profiles.length_bins(),
        pd.DataFrame(
            {
                "length_bin": np.array([0, 1, 2], dtype=np.int32),
                "length_start_bp": np.array([30, 50, 70], dtype=np.int32),
                "length_end_bp": np.array([50, 70, 100], dtype=np.int32),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        profiles.positions(),
        pd.DataFrame(
            {
                "position": np.array([0, 1, 2, 3, 4], dtype=np.int32),
                "position_bin_start_bp": np.array([0, 2, 4, 6, 8], dtype=np.int32),
                "position_bin_end_bp": np.array([2, 4, 6, 8, 10], dtype=np.int32),
            }
        ),
    )
    np.testing.assert_allclose(profiles.counts_array(), EXPECTED_COUNTS)


def test_cfdnalab_package_resolves_fixture_indices_and_slices(
    midpoint_zarr_path: Path,
) -> None:
    profiles = cfdnalab.read_midpoints(midpoint_zarr_path)

    assert profiles.group_idx("beta-site") == 1
    assert profiles.length_bin_idx(49) == 0
    assert profiles.length_bin_idx(50) == 1
    assert profiles.length_bin_idx(70) == 2

    np.testing.assert_allclose(
        profiles.counts_array(group_idxs=1, length_bin_idxs=1),
        np.array([[[0.0, 0.0, 1.5, 0.0, 0.5]]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        profiles.counts_array(groups="LYL1"),
        EXPECTED_COUNTS[0:1],
    )
    np.testing.assert_allclose(
        profiles.counts_array(with_lengths=70),
        EXPECTED_COUNTS[:, 2:3, :],
    )
    np.testing.assert_allclose(
        profiles.counts_array(
            groups=["gamma_long", "LYL1"],
            with_length_range=(50, 100),
        ),
        np.array(
            [
                [
                    [0.5, 0.0, 0.0, 1.5, 0.0],
                    [0.0, 0.0, 0.5, 0.0, 1.0],
                ],
                [
                    [0.0, 0.0, 1.5, 0.5, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 2.0],
                ],
            ],
            dtype=np.float32,
        ),
    )


def test_cfdnalab_package_builds_profile_dataframe_from_fixture(
    midpoint_zarr_path: Path,
) -> None:
    profiles = cfdnalab.read_midpoints(midpoint_zarr_path)

    frame = profiles.data_frame(group_idxs=0, length_bin_idxs=1)
    range_frame = profiles.data_frame(groups="LYL1", with_length_range=(50, 100))

    pd.testing.assert_frame_equal(
        frame,
        pd.DataFrame(
            {
                "group_idx": [0, 0, 0, 0, 0],
                "group_name": ["LYL1"] * 5,
                "eligible_intervals": [2, 2, 2, 2, 2],
                "length_bin": [1, 1, 1, 1, 1],
                "length_start_bp": [50, 50, 50, 50, 50],
                "length_end_bp": [70, 70, 70, 70, 70],
                "position": np.array([0, 1, 2, 3, 4], dtype=np.int32),
                "position_bin_start_bp": np.array([0, 2, 4, 6, 8], dtype=np.int32),
                "position_bin_end_bp": np.array([2, 4, 6, 8, 10], dtype=np.int32),
                "count": np.array([0.0, 0.0, 1.5, 0.5, 0.0], dtype=np.float32),
            }
        ),
    )
    assert range_frame["group_name"].unique().tolist() == ["LYL1"]
    assert range_frame["length_bin"].tolist() == [1] * 5 + [2] * 5
    assert range_frame["count"].tolist() == [
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


def test_cfdnalab_package_builds_group_and_length_dataframes_from_fixture(
    midpoint_zarr_path: Path,
) -> None:
    profiles = cfdnalab.read_midpoints(midpoint_zarr_path)

    group_frame = profiles.data_frame(groups="beta-site")
    length_frame = profiles.data_frame(length_bin_idxs=2)

    assert group_frame.shape == (15, 10)
    assert group_frame["group_name"].unique().tolist() == ["beta-site"]
    assert group_frame["length_bin"].tolist() == [0] * 5 + [1] * 5 + [2] * 5
    assert group_frame["position"].tolist() == [0, 1, 2, 3, 4] * 3
    assert group_frame["count"].tolist() == [
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

    assert length_frame.shape == (15, 10)
    assert length_frame["length_bin"].unique().tolist() == [2]
    assert length_frame["group_name"].tolist() == (
        ["LYL1"] * 5 + ["beta-site"] * 5 + ["gamma_long"] * 5
    )
    assert length_frame["position"].tolist() == [0, 1, 2, 3, 4] * 3
    assert length_frame["count"].tolist() == [
        0.0,
        0.0,
        0.0,
        0.0,
        2.0,
        0.0,
        0.5,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.5,
        0.0,
        1.0,
    ]
