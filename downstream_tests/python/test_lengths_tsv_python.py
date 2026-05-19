from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ZSTD_MAGIC = bytes([0x28, 0xB5, 0x2F, 0xFD])


def test_length_count_tsv_fixtures_are_zstd_frames(
    global_length_counts_path: Path,
    windowed_length_counts_path: Path,
    grouped_length_counts_path: Path,
    windowed_length_counts_no_blacklist_path: Path,
    grouped_length_counts_no_blacklist_path: Path,
) -> None:
    for path in (
        global_length_counts_path,
        windowed_length_counts_path,
        grouped_length_counts_path,
        windowed_length_counts_no_blacklist_path,
        grouped_length_counts_no_blacklist_path,
    ):
        with path.open("rb") as file_handle:
            assert file_handle.read(4) == ZSTD_MAGIC


def test_pandas_reads_global_length_count_tsv(
    global_length_counts_path: Path,
) -> None:
    table = pd.read_csv(global_length_counts_path, sep="\t", compression="infer")

    pd.testing.assert_frame_equal(
        table,
        pd.DataFrame(
            {
                "count_30_50": np.array([3], dtype=np.int64),
                "count_50_70": np.array([2], dtype=np.int64),
                "count_70_100": np.array([1], dtype=np.int64),
            }
        ),
    )


def test_pandas_reads_windowed_length_count_tsv(
    windowed_length_counts_path: Path,
    windowed_length_counts_no_blacklist_path: Path,
) -> None:
    table = pd.read_csv(windowed_length_counts_path, sep="\t", compression="infer")
    no_blacklist = pd.read_csv(
        windowed_length_counts_no_blacklist_path,
        sep="\t",
        compression="infer",
    )

    assert table.columns.tolist() == [
        "chrom",
        "start",
        "end",
        "blacklisted_fraction",
        "count_30_50",
        "count_50_70",
        "count_70_100",
    ]
    assert no_blacklist.columns.tolist() == [
        "chrom",
        "start",
        "end",
        "count_30_50",
        "count_50_70",
        "count_70_100",
    ]
    assert table["chrom"].tolist() == ["chr1", "chr1", "chr1", "chr1"]
    assert table["start"].tolist() == [0, 100, 200, 300]
    assert table["end"].tolist() == [100, 200, 300, 360]
    np.testing.assert_allclose(
        table["blacklisted_fraction"].to_numpy(),
        np.array([0.04, 0.05, 0.1, 0.25], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        table[["count_30_50", "count_50_70", "count_70_100"]].to_numpy(),
        np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64),
    )


def test_pandas_reads_grouped_length_count_tsv(
    grouped_length_counts_path: Path,
    grouped_length_counts_no_blacklist_path: Path,
) -> None:
    table = pd.read_csv(grouped_length_counts_path, sep="\t", compression="infer")
    no_blacklist = pd.read_csv(
        grouped_length_counts_no_blacklist_path,
        sep="\t",
        compression="infer",
    )

    assert table.columns.tolist() == [
        "group_name",
        "eligible_windows",
        "blacklisted_fraction",
        "count_30_50",
        "count_50_70",
        "count_70_100",
    ]
    assert no_blacklist.columns.tolist() == [
        "group_name",
        "eligible_windows",
        "count_30_50",
        "count_50_70",
        "count_70_100",
    ]
    assert table["group_name"].tolist() == ["beta", "alpha", "gamma", "zero"]
    assert table["eligible_windows"].tolist() == [2, 1, 1, 1]
    np.testing.assert_allclose(
        table["blacklisted_fraction"].to_numpy(),
        np.array([0.07, 0.05, 0.25, 0.333], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        table[["count_30_50", "count_50_70", "count_70_100"]].to_numpy(),
        np.array([[2, 0, 1], [0, 2, 0], [1, 0, 0], [0, 0, 0]], dtype=np.int64),
    )
