from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import zstandard as zstd

import cfdnalab


def _write_length_tsv(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_length_tsv_zst(path: Path, lines: list[str]) -> Path:
    tsv_bytes = ("\n".join(lines) + "\n").encode("utf-8")
    path.write_bytes(zstd.ZstdCompressor().compress(tsv_bytes))
    return path


def test_global_length_counts_expose_bins_counts_and_data_frames(tmp_path: Path) -> None:
    path = _write_length_tsv(
        tmp_path / "sample.length_counts.tsv",
        [
            "count_30\tcount_31_40",
            "12\t3",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    assert isinstance(lengths, cfdnalab.GlobalLengthCounts)
    assert repr(lengths) == (
        "GlobalLengthCounts("
        f"path={str(path)!r}, "
        "mode='global', "
        "shape=(1, 2)"
        ")"
    )
    pd.testing.assert_frame_equal(
        lengths.length_bins(),
        pd.DataFrame(
            {
                "length_bin": np.array([0, 1], dtype=np.int32),
                "length_start_bp": np.array([30, 31], dtype=np.int64),
                "length_end_bp": np.array([31, 40], dtype=np.int64),
                "length_midpoint_bp": np.array([30.5, 35.5], dtype=np.float64),
                "length_width_bp": np.array([1, 9], dtype=np.int64),
            }
        ),
    )
    assert lengths.length_bin_idx(30) == 0
    assert lengths.length_bin_idx(39) == 1
    with pytest.raises(TypeError, match="Fragment length must be an integer"):
        lengths.length_bin_idx(True)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Fragment length must be an integer"):
        lengths.length_bin_idx(30.5)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Fragment length must be an integer"):
        lengths.length_bin_idx("30")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="Fragment length must be non-negative"):
        lengths.length_bin_idx(-1)
    with pytest.raises(KeyError, match="No length-count bin contains length 40"):
        lengths.length_bin_idx(40)
    np.testing.assert_allclose(lengths.counts_matrix(), np.array([[12.0, 3.0]]))
    np.testing.assert_allclose(lengths.counts_vec(), np.array([12.0, 3.0]))

    pd.testing.assert_frame_equal(
        lengths.data_frame(value="fraction"),
        pd.DataFrame(
            {
                "length_bin": np.array([0, 1], dtype=np.int32),
                "length_start_bp": np.array([30, 31], dtype=np.int64),
                "length_end_bp": np.array([31, 40], dtype=np.int64),
                "length_midpoint_bp": np.array([30.5, 35.5], dtype=np.float64),
                "length_width_bp": np.array([1, 9], dtype=np.int64),
                "fraction": np.array([0.8, 0.2], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        lengths.data_frame(value="density", keep_wide=True),
        pd.DataFrame(
            {
                "density_30": np.array([0.8], dtype=np.float64),
                "density_31_40": np.array([0.2 / 9], dtype=np.float64),
            }
        ),
    )
    with pytest.raises(TypeError, match="keep_wide must be a bool"):
        lengths.data_frame(keep_wide="yes")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        lengths.data_frame(max_blacklisted_fraction=1.0)  # type: ignore[call-arg]


def test_length_wide_data_frames_preserve_source_count_column_labels(
    tmp_path: Path,
) -> None:
    path = _write_length_tsv(
        tmp_path / "padded.length_counts.tsv",
        [
            "count_030\tcount_031_040",
            "12\t3",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    assert "count_column" not in lengths.length_bins().columns
    assert lengths.data_frame(keep_wide=True).columns.tolist() == [
        "count_030",
        "count_031_040",
    ]
    assert lengths.data_frame(value="fraction", keep_wide=True).columns.tolist() == [
        "fraction_030",
        "fraction_031_040",
    ]


def test_windowed_length_counts_filter_blacklisted_rows(tmp_path: Path) -> None:
    path = _write_length_tsv(
        tmp_path / "windowed.length_counts.tsv",
        [
            "chrom\tstart\tend\tblacklisted_fraction\tcount_30\tcount_31_40",
            "chr1\t10\t20\t0.25\t12\t3",
            "chr2\t30\t45\t0\t0\t5",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    assert isinstance(lengths, cfdnalab.WindowedLengthCounts)
    assert repr(lengths) == (
        "WindowedLengthCounts("
        f"path={str(path)!r}, "
        "mode='windowed', "
        "shape=(2, 2)"
        ")"
    )
    pd.testing.assert_frame_equal(
        lengths.windows(),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chrom": np.array(["chr1", "chr2"], dtype=object),
                "start": np.array([10, 30], dtype=np.int64),
                "end": np.array([20, 45], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.0], dtype=np.float64),
            }
        ),
    )
    selected = lengths.data_frame(window_idxs=1)
    assert selected["window_idx"].tolist() == [1, 1]
    assert selected["count"].tolist() == [0.0, 5.0]
    np.testing.assert_allclose(lengths.counts_for_window(1), np.array([0.0, 5.0]))
    with pytest.raises(TypeError, match="window_idx must be an integer"):
        lengths.counts_for_window(True)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="window_idx must be an integer"):
        lengths.data_frame(window_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(IndexError, match="window_idx 2 is outside 0..1"):
        lengths.data_frame(window_idxs=2)
    with pytest.raises(ValueError, match="window_idxs contains duplicate values"):
        lengths.data_frame(window_idxs=[1, 1])
    assert lengths.data_frame(window_idxs=[]).empty

    filtered = lengths.data_frame(max_blacklisted_fraction=0.1)
    assert filtered["window_idx"].unique().tolist() == [1]
    assert lengths.data_frame(window_idxs=0, max_blacklisted_fraction=0.1).empty
    with pytest.raises(ValueError, match="max_blacklisted_fraction must be"):
        lengths.data_frame(max_blacklisted_fraction=1.1)


def test_grouped_length_counts_support_group_selectors_and_zero_rows(
    tmp_path: Path,
) -> None:
    path = _write_length_tsv(
        tmp_path / "grouped.length_counts.tsv",
        [
            "group_name\teligible_windows\tblacklisted_fraction\tcount_30\tcount_31_40",
            "alpha\t2\t0.25\t12\t3",
            "beta\t1\t0\t0\t5",
            "zero\t0\t0\t0\t0",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    assert isinstance(lengths, cfdnalab.GroupedLengthCounts)
    assert repr(lengths) == (
        "GroupedLengthCounts("
        f"path={str(path)!r}, "
        "mode='grouped', "
        "shape=(3, 2)"
        ")"
    )
    assert lengths.group_idx("beta") == 1
    pd.testing.assert_frame_equal(
        lengths.groups(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2], dtype=np.int32),
                "group_name": np.array(["alpha", "beta", "zero"], dtype=object),
                "eligible_windows": np.array([2, 1, 0], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.0, 0.0], dtype=np.float64),
            }
        ),
    )

    beta = lengths.data_frame(groups="beta")
    assert beta["group_idx"].tolist() == [1, 1]
    assert beta["count"].tolist() == [0.0, 5.0]
    np.testing.assert_allclose(lengths.counts_for_group("beta"), np.array([0.0, 5.0]))
    np.testing.assert_allclose(lengths.counts_for_group(1), np.array([0.0, 5.0]))
    with pytest.raises(TypeError, match="group_idx must be an integer"):
        lengths.counts_for_group(True)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="group_idx must be an integer"):
        lengths.data_frame(group_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(IndexError, match="group_idx 3 is outside 0..2"):
        lengths.data_frame(group_idxs=3)
    with pytest.raises(TypeError, match="groups must contain strings"):
        lengths.data_frame(groups=[0, "beta"])  # type: ignore[list-item]
    with pytest.raises(ValueError, match="groups contains duplicate values"):
        lengths.data_frame(groups=["beta", "beta"])
    with pytest.raises(ValueError, match="group_idxs contains duplicate values"):
        lengths.data_frame(group_idxs=[1, 1])
    with pytest.raises(ValueError, match="value must be one of"):
        lengths.data_frame(value="invalid")
    assert lengths.data_frame(groups=[]).empty
    assert lengths.data_frame(group_idxs=[]).empty
    assert lengths.data_frame(groups=["beta", "alpha"])["group_idx"].tolist() == [
        1,
        1,
        0,
        0,
    ]

    filtered_wide = lengths.data_frame(
        groups=["alpha", "zero"],
        value="density",
        keep_wide=True,
        max_blacklisted_fraction=0.1,
    )
    pd.testing.assert_frame_equal(
        filtered_wide,
        pd.DataFrame(
            {
                "group_idx": np.array([2], dtype=np.int32),
                "group_name": np.array(["zero"], dtype=object),
                "eligible_windows": np.array([0], dtype=np.int64),
                "blacklisted_fraction": np.array([0.0], dtype=np.float64),
                "density_30": np.array([np.nan], dtype=np.float64),
                "density_31_40": np.array([np.nan], dtype=np.float64),
            }
        ),
    )


def test_no_blacklist_length_counts_keep_all_by_default_and_error_for_cutoff(
    tmp_path: Path,
) -> None:
    path = _write_length_tsv(
        tmp_path / "no_blacklist.length_counts.tsv",
        [
            "chrom\tstart\tend\tcount_30",
            "chr1\t10\t20\t1",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    pd.testing.assert_frame_equal(
        lengths.windows(),
        pd.DataFrame(
            {
                "window_idx": np.array([0], dtype=np.int32),
                "chrom": np.array(["chr1"], dtype=object),
                "start": np.array([10], dtype=np.int64),
                "end": np.array([20], dtype=np.int64),
            }
        ),
    )
    assert lengths.data_frame(max_blacklisted_fraction=1.0)["count"].tolist() == [1.0]
    with pytest.raises(ValueError, match="has no blacklisted_fraction column"):
        lengths.data_frame(max_blacklisted_fraction=0.5)


def test_length_tsv_validation_rejects_ambiguous_shapes(tmp_path: Path) -> None:
    missing = tmp_path / "missing.length_counts.tsv"
    with pytest.raises(FileNotFoundError, match="Length-count TSV does not exist"):
        cfdnalab.read_lengths(missing)

    directory = tmp_path / "length_counts.tsv"
    directory.mkdir()
    with pytest.raises(IsADirectoryError, match="exists but is a directory"):
        cfdnalab.read_lengths(directory)

    wrong_extension = _write_length_tsv(
        tmp_path / "wrong.length_counts.csv",
        [
            "count_30",
            "1",
        ],
    )
    with pytest.raises(ValueError, match="must end in '.tsv' or '.tsv.zst'"):
        cfdnalab.read_lengths(wrong_extension)

    wrong_gzip_extension = _write_length_tsv(
        tmp_path / "wrong.length_counts.tsv.gz",
        [
            "count_30",
            "1",
        ],
    )
    with pytest.raises(ValueError, match="must end in '.tsv' or '.tsv.zst'"):
        cfdnalab.read_lengths(wrong_gzip_extension)

    with pytest.raises(TypeError, match="path-like string"):
        cfdnalab.read_lengths(123)  # type: ignore[arg-type]

    bad_count = _write_length_tsv(
        tmp_path / "bad_count.length_counts.tsv",
        [
            "count_30\tother",
            "1\t2",
        ],
    )
    with pytest.raises(ValueError, match="count columns must be contiguous"):
        cfdnalab.read_lengths(bad_count)

    multiple_global_rows = _write_length_tsv(
        tmp_path / "multiple_global.length_counts.tsv",
        [
            "count_30",
            "1",
            "2",
        ],
    )
    with pytest.raises(ValueError, match="Global length-count output must contain exactly one row"):
        cfdnalab.read_lengths(multiple_global_rows)

    duplicate_bins = _write_length_tsv(
        tmp_path / "duplicate_bins.length_counts.tsv",
        [
            "count_30\tcount_30_31",
            "1\t2",
        ],
    )
    with pytest.raises(ValueError, match="duplicate length bins"):
        cfdnalab.read_lengths(duplicate_bins)

    duplicate_columns = _write_length_tsv(
        tmp_path / "duplicate_columns.length_counts.tsv",
        [
            "count_30\tcount_30",
            "1\t2",
        ],
    )
    with pytest.raises(ValueError, match="column names must be unique"):
        cfdnalab.read_lengths(duplicate_columns)

    unsupported_metadata = _write_length_tsv(
        tmp_path / "unsupported.length_counts.tsv",
        [
            "window_idx\tchrom\tstart\tend\tcount_30",
            "0\tchr1\t10\t20\t1",
        ],
    )
    with pytest.raises(ValueError, match="Could not infer length-count output mode"):
        cfdnalab.read_lengths(unsupported_metadata)

    negative_count = _write_length_tsv(
        tmp_path / "negative_count.length_counts.tsv",
        [
            "count_30",
            "-1",
        ],
    )
    with pytest.raises(ValueError, match="count_30 must contain non-negative values"):
        cfdnalab.read_lengths(negative_count)

    nonfinite_count = _write_length_tsv(
        tmp_path / "nonfinite_count.length_counts.tsv",
        [
            "count_30",
            "Inf",
        ],
    )
    with pytest.raises(ValueError, match="count_30 must contain finite values"):
        cfdnalab.read_lengths(nonfinite_count)

    negative_window_start = _write_length_tsv(
        tmp_path / "negative_window_start.length_counts.tsv",
        [
            "chrom\tstart\tend\tcount_30",
            "chr1\t-1\t20\t1",
        ],
    )
    with pytest.raises(ValueError, match="start must contain non-negative"):
        cfdnalab.read_lengths(negative_window_start)


def test_length_tsv_reads_zstandard_and_numeric_chromosomes(tmp_path: Path) -> None:
    path = _write_length_tsv_zst(
        tmp_path / "windowed.length_counts.tsv.zst",
        [
            "chrom\tstart\tend\tcount_30",
            "1\t10\t20\t5",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    pd.testing.assert_frame_equal(
        lengths.windows(),
        pd.DataFrame(
            {
                "window_idx": np.array([0], dtype=np.int32),
                "chrom": np.array(["1"], dtype=object),
                "start": np.array([10], dtype=np.int64),
                "end": np.array([20], dtype=np.int64),
            }
        ),
    )


def test_length_bin_lookup_rejects_gaps_and_overlaps(tmp_path: Path) -> None:
    gapped = cfdnalab.read_lengths(
        _write_length_tsv(
            tmp_path / "gapped.length_counts.tsv",
            [
                "count_30_40\tcount_50_60",
                "1\t2",
            ],
        )
    )
    with pytest.raises(KeyError, match="No length-count bin contains length 45"):
        gapped.length_bin_idx(45)

    overlapping = cfdnalab.read_lengths(
        _write_length_tsv(
            tmp_path / "overlapping.length_counts.tsv",
            [
                "count_30_50\tcount_40_60",
                "1\t2",
            ],
        )
    )
    with pytest.raises(ValueError, match="Multiple length-count bins contain length 45"):
        overlapping.length_bin_idx(45)
