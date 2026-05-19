from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cfdnalab


def test_cfdnalab_package_reads_global_length_counts(
    global_length_counts_path: Path,
) -> None:
    lengths = cfdnalab.read_lengths(global_length_counts_path)

    assert isinstance(lengths, cfdnalab.GlobalLengthCounts)
    pd.testing.assert_frame_equal(
        lengths.length_bins(),
        pd.DataFrame(
            {
                "length_bin": np.array([0, 1, 2], dtype=np.int32),
                "length_start_bp": np.array([30, 50, 70], dtype=np.int64),
                "length_end_bp": np.array([50, 70, 100], dtype=np.int64),
                "length_midpoint_bp": np.array([40.0, 60.0, 85.0], dtype=np.float64),
                "length_width_bp": np.array([20, 20, 30], dtype=np.int64),
            }
        ),
    )
    np.testing.assert_allclose(lengths.counts_vec(), np.array([3.0, 2.0, 1.0]))
    np.testing.assert_allclose(
        lengths.data_frame(value="density")["density"].to_numpy(),
        np.array([0.025, 1 / 60, 1 / 180], dtype=np.float64),
    )


def test_cfdnalab_package_reads_windowed_length_counts(
    windowed_length_counts_path: Path,
) -> None:
    lengths = cfdnalab.read_lengths(windowed_length_counts_path)

    assert isinstance(lengths, cfdnalab.WindowedLengthCounts)
    pd.testing.assert_frame_equal(
        lengths.windows(),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1, 2, 3], dtype=np.int32),
                "chrom": np.array(["chr1", "chr1", "chr1", "chr1"], dtype=object),
                "start": np.array([0, 100, 200, 300], dtype=np.int64),
                "end": np.array([100, 200, 300, 360], dtype=np.int64),
                "blacklisted_fraction": np.array([0.04, 0.05, 0.1, 0.25]),
            }
        ),
    )
    np.testing.assert_allclose(
        lengths.counts_matrix(),
        np.array(
            [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]
        ),
    )
    np.testing.assert_allclose(lengths.counts_for_window(3), np.array([1.0, 0.0, 0.0]))

    selected = lengths.data_frame(window_idxs=[1, 3], value="fraction", keep_wide=True)
    assert selected["window_idx"].tolist() == [1, 3]
    np.testing.assert_allclose(selected["fraction_30_50"], np.array([0.0, 1.0]))
    np.testing.assert_allclose(selected["fraction_50_70"], np.array([1.0, 0.0]))
    np.testing.assert_allclose(selected["fraction_70_100"], np.array([0.0, 0.0]))

    filtered = lengths.data_frame(max_blacklisted_fraction=0.05)
    assert filtered["window_idx"].unique().tolist() == [0, 1]


def test_cfdnalab_package_reads_grouped_length_counts(
    grouped_length_counts_path: Path,
) -> None:
    lengths = cfdnalab.read_lengths(grouped_length_counts_path)

    assert isinstance(lengths, cfdnalab.GroupedLengthCounts)
    assert lengths.group_idx("gamma") == 2
    pd.testing.assert_frame_equal(
        lengths.groups(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2, 3], dtype=np.int32),
                "group_name": np.array(["beta", "alpha", "gamma", "zero"], dtype=object),
                "eligible_windows": np.array([2, 1, 1, 1], dtype=np.int64),
                "blacklisted_fraction": np.array([0.07, 0.05, 0.25, 0.333]),
            }
        ),
    )
    beta = lengths.data_frame(groups="beta")
    assert beta["count"].tolist() == [2.0, 0.0, 1.0]
    np.testing.assert_allclose(
        lengths.counts_for_group("beta"), np.array([2.0, 0.0, 1.0])
    )

    wide_density = lengths.data_frame(
        groups=["alpha", "zero"],
        value="density",
        keep_wide=True,
    )
    np.testing.assert_allclose(wide_density["density_30_50"], np.array([0.0, np.nan]))
    np.testing.assert_allclose(wide_density["density_50_70"], np.array([1 / 20, np.nan]))
    np.testing.assert_allclose(wide_density["density_70_100"], np.array([0.0, np.nan]))


def test_cfdnalab_package_reads_no_blacklist_length_counts(
    windowed_length_counts_no_blacklist_path: Path,
    grouped_length_counts_no_blacklist_path: Path,
) -> None:
    windowed = cfdnalab.read_lengths(windowed_length_counts_no_blacklist_path)
    grouped = cfdnalab.read_lengths(grouped_length_counts_no_blacklist_path)

    assert "blacklisted_fraction" not in windowed.windows().columns
    assert "blacklisted_fraction" not in grouped.groups().columns
    assert windowed.data_frame(max_blacklisted_fraction=1.0)["count"].tolist() == [
        2.0,
        0.0,
        0.0,
        0.0,
        2.0,
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        0.0,
        0.0,
    ]
    assert grouped.data_frame(groups="beta")["count"].tolist() == [2.0, 0.0, 1.0]
    with pytest.raises(ValueError, match="has no blacklisted_fraction column"):
        windowed.data_frame(max_blacklisted_fraction=0.5)
    with pytest.raises(ValueError, match="has no blacklisted_fraction column"):
        grouped.data_frame(max_blacklisted_fraction=0.5)
