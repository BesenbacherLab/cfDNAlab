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


def test_global_length_counts_expose_bins_counts_and_data_frames(
    tmp_path: Path,
) -> None:
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
        f"GlobalLengthCounts(path={str(path)!r}, mode='global', shape=(1, 2))"
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
    np.testing.assert_allclose(lengths.counts_array(), np.array([[12.0, 3.0]]))
    np.testing.assert_allclose(lengths.counts_array(with_lengths=39), np.array([[3.0]]))
    np.testing.assert_allclose(
        lengths.counts_array(with_length_range=(31, 35)), np.array([[3.0]])
    )

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
        lengths.data_frame(with_length_range=(31, 35), value="fraction"),
        pd.DataFrame(
            {
                "length_bin": np.array([1], dtype=np.int32),
                "length_start_bp": np.array([31], dtype=np.int64),
                "length_end_bp": np.array([40], dtype=np.int64),
                "length_midpoint_bp": np.array([35.5], dtype=np.float64),
                "length_width_bp": np.array([9], dtype=np.int64),
                "fraction": np.array([0.2], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        lengths.data_frame(
            with_length_range=(31, 35),
            value="fraction",
            denominator="selected_bins",
        ),
        pd.DataFrame(
            {
                "length_bin": np.array([1], dtype=np.int32),
                "length_start_bp": np.array([31], dtype=np.int64),
                "length_end_bp": np.array([40], dtype=np.int64),
                "length_midpoint_bp": np.array([35.5], dtype=np.float64),
                "length_width_bp": np.array([9], dtype=np.int64),
                "fraction": np.array([1.0], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        lengths.data_frame(
            with_length_range=(31, 35),
            value="density",
            keep_wide=True,
        ),
        pd.DataFrame(
            {
                "density_31_40": np.array([0.2 / 9], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        lengths.data_frame(
            with_length_range=(31, 35),
            value="density",
            denominator="selected_bins",
            keep_wide=True,
        ),
        pd.DataFrame(
            {
                "density_31_40": np.array([1.0 / 9], dtype=np.float64),
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


def test_length_bin_selectors_preserve_order_boundaries_and_totals(
    tmp_path: Path,
) -> None:
    path = _write_length_tsv(
        tmp_path / "ranges.length_counts.tsv",
        [
            "count_30_50\tcount_50_70\tcount_70_100",
            "3\t2\t1",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    np.testing.assert_allclose(
        lengths.counts_array(with_lengths=[75, 30]),
        np.array([[1.0, 3.0]]),
    )
    np.testing.assert_allclose(
        lengths.counts_array(length_bin_idxs=[2, 0]),
        np.array([[1.0, 3.0]]),
    )
    empty_counts = lengths.counts_array(with_lengths=[])
    assert empty_counts.shape == (1, 0)

    exact_range = lengths.data_frame(with_length_range=(50, 70), value="fraction")
    assert exact_range["length_bin"].tolist() == [1]
    assert exact_range["fraction"].tolist() == [2 / 6]
    exact_range_selected = lengths.data_frame(
        with_length_range=(50, 70),
        value="fraction",
        denominator="selected_bins",
    )
    assert exact_range_selected["length_bin"].tolist() == [1]
    assert exact_range_selected["fraction"].tolist() == [1.0]

    left_edge_range = lengths.data_frame(with_length_range=(49, 50))
    assert left_edge_range["length_bin"].tolist() == [0]
    assert left_edge_range["count"].tolist() == [3.0]

    right_edge_range = lengths.data_frame(with_length_range=(70, 100))
    assert right_edge_range["length_bin"].tolist() == [2]
    assert right_edge_range["count"].tolist() == [1.0]

    overlap_range = lengths.data_frame(with_length_range=(49, 71))
    assert overlap_range["length_bin"].tolist() == [0, 1, 2]
    assert overlap_range["count"].tolist() == [3.0, 2.0, 1.0]

    wide_selected = lengths.data_frame(length_bin_idxs=[2, 0], keep_wide=True)
    pd.testing.assert_frame_equal(
        wide_selected,
        pd.DataFrame(
            {
                "count_70_100": np.array([1.0], dtype=np.float64),
                "count_30_50": np.array([3.0], dtype=np.float64),
            }
        ),
    )
    wide_selected_fraction = lengths.data_frame(
        length_bin_idxs=[2, 0],
        value="fraction",
        denominator="selected_bins",
        keep_wide=True,
    )
    pd.testing.assert_frame_equal(
        wide_selected_fraction,
        pd.DataFrame(
            {
                "fraction_70_100": np.array([0.25], dtype=np.float64),
                "fraction_30_50": np.array([0.75], dtype=np.float64),
            }
        ),
    )
    empty_wide = lengths.data_frame(length_bin_idxs=[], keep_wide=True)
    assert empty_wide.shape == (1, 0)
    assert empty_wide.columns.tolist() == []
    assert lengths.data_frame(with_lengths=[]).empty
    assert lengths.data_frame(length_bin_idxs=[], value="density").empty


def test_length_bin_selectors_reject_invalid_values(tmp_path: Path) -> None:
    path = _write_length_tsv(
        tmp_path / "selector_errors.length_counts.tsv",
        [
            "count_30_50\tcount_50_70",
            "3\t2",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    with pytest.raises(TypeError, match="Fragment length must be an integer"):
        lengths.counts_array(with_lengths=True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="with_lengths contains duplicate values"):
        lengths.counts_array(with_lengths=[35, 35])
    with pytest.raises(TypeError, match="length_bin_idxs must be an integer"):
        lengths.counts_array(length_bin_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="length_bin_idxs contains duplicate values"):
        lengths.counts_array(length_bin_idxs=[0, 0])
    with pytest.raises(TypeError, match="pair of non-negative integer bp bounds"):
        lengths.data_frame(with_length_range="30-50")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="exactly two bounds"):
        lengths.data_frame(with_length_range=(30, 50, 70))
    with pytest.raises(TypeError, match="Fragment length must be an integer"):
        lengths.data_frame(with_length_range=(30.5, 50))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="Fragment length must be non-negative"):
        lengths.data_frame(with_length_range=(-1, 50))
    with pytest.raises(ValueError, match="start must be smaller than end"):
        lengths.data_frame(with_length_range=(50, 50))
    with pytest.raises(ValueError, match="denominator must be one of"):
        lengths.data_frame(value="fraction", denominator="selected")


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
        f"WindowedLengthCounts(path={str(path)!r}, mode='windowed', shape=(2, 2))"
    )
    pd.testing.assert_frame_equal(
        lengths.window_metadata(),
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
    selected_range = lengths.data_frame(window_idxs=1, with_length_range=(31, 35))
    assert selected_range["window_idx"].tolist() == [1]
    assert selected_range["length_bin"].tolist() == [1]
    assert selected_range["count"].tolist() == [5.0]
    np.testing.assert_allclose(
        lengths.counts_array(window_idxs=1), np.array([[0.0, 5.0]])
    )
    np.testing.assert_allclose(
        lengths.counts_array(window_idxs=1, length_bin_idxs=1), np.array([[5.0]])
    )
    with pytest.raises(TypeError, match="window_idxs must be an integer"):
        lengths.counts_array(window_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="window_idxs must be an integer"):
        lengths.data_frame(window_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(IndexError, match="window_idx 2 is outside 0..1"):
        lengths.data_frame(window_idxs=2)
    with pytest.raises(ValueError, match="window_idxs contains duplicate values"):
        lengths.data_frame(window_idxs=[1, 1])
    assert lengths.data_frame(window_idxs=[]).empty
    assert lengths.data_frame(length_bin_idxs=[]).empty
    metadata_only = lengths.data_frame(
        window_idxs=1,
        length_bin_idxs=[],
        keep_wide=True,
    )
    pd.testing.assert_frame_equal(
        metadata_only,
        pd.DataFrame(
            {
                "window_idx": np.array([1], dtype=np.int32),
                "chrom": np.array(["chr2"], dtype=object),
                "start": np.array([30], dtype=np.int64),
                "end": np.array([45], dtype=np.int64),
                "blacklisted_fraction": np.array([0.0], dtype=np.float64),
            }
        ),
    )

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
        f"GroupedLengthCounts(path={str(path)!r}, mode='grouped', shape=(3, 2))"
    )
    assert lengths.group_idx("beta") == 1
    pd.testing.assert_frame_equal(
        lengths.group_metadata(),
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
    beta_length = lengths.data_frame(groups="beta", with_lengths=39)
    assert beta_length["group_idx"].tolist() == [1]
    assert beta_length["length_bin"].tolist() == [1]
    assert beta_length["count"].tolist() == [5.0]
    selected_fraction = lengths.data_frame(
        groups=["alpha", "beta", "zero"],
        with_length_range=(31, 35),
        value="fraction",
        denominator="selected_bins",
    )
    assert selected_fraction["group_idx"].tolist() == [0, 1, 2]
    assert selected_fraction["fraction"].tolist()[:2] == [1.0, 1.0]
    assert np.isnan(selected_fraction["fraction"].iloc[2])
    np.testing.assert_allclose(
        lengths.counts_array(groups="beta"), np.array([[0.0, 5.0]])
    )
    np.testing.assert_allclose(
        lengths.counts_array(group_idxs=1), np.array([[0.0, 5.0]])
    )
    with pytest.raises(TypeError, match="group_idxs must be an integer"):
        lengths.counts_array(group_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="group_idxs must be an integer"):
        lengths.data_frame(group_idxs=True)  # type: ignore[arg-type]
    with pytest.raises(IndexError, match="group_idx 3 is outside 0..2"):
        lengths.data_frame(group_idxs=3)
    with pytest.raises(TypeError, match="groups must contain strings"):
        lengths.data_frame(groups=[0, "beta"])  # type: ignore[list-item]
    with pytest.raises(ValueError, match="groups contains duplicate values"):
        lengths.data_frame(groups=["beta", "beta"])
    with pytest.raises(ValueError, match="group_idxs contains duplicate values"):
        lengths.data_frame(group_idxs=[1, 1])
    with pytest.raises(ValueError, match="Use only one of with_lengths"):
        lengths.data_frame(with_lengths=39, with_length_range=(31, 35))
    with pytest.raises(
        ValueError,
        match=(
            "with_lengths values must resolve to distinct length bins; "
            "31 and 39 both resolve to length_bin_idx 1"
        ),
    ):
        lengths.data_frame(with_lengths=[31, 39])
    with pytest.raises(ValueError, match="with_lengths contains duplicate values"):
        lengths.data_frame(with_lengths=[39, 39])
    with pytest.raises(ValueError, match="does not overlap any length bins"):
        lengths.data_frame(with_length_range=(40, 45))
    with pytest.raises(ValueError, match="value must be one of"):
        lengths.data_frame(value="invalid")
    assert lengths.data_frame(groups=[]).empty
    assert lengths.data_frame(group_idxs=[]).empty
    metadata_only = lengths.data_frame(
        groups="beta",
        length_bin_idxs=[],
        keep_wide=True,
    )
    pd.testing.assert_frame_equal(
        metadata_only,
        pd.DataFrame(
            {
                "group_idx": np.array([1], dtype=np.int32),
                "group_name": np.array(["beta"], dtype=object),
                "eligible_windows": np.array([1], dtype=np.int64),
                "blacklisted_fraction": np.array([0.0], dtype=np.float64),
            }
        ),
    )
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


def test_grouped_length_selectors_preserve_requested_row_and_bin_order(
    tmp_path: Path,
) -> None:
    path = _write_length_tsv(
        tmp_path / "group_order.length_counts.tsv",
        [
            "group_name\teligible_windows\tcount_30\tcount_31_40",
            "alpha\t2\t1\t2",
            "beta\t3\t3\t4",
            "gamma\t5\t5\t6",
        ],
    )

    lengths = cfdnalab.read_lengths(path)

    np.testing.assert_allclose(
        lengths.counts_array(groups=["gamma", "alpha"], length_bin_idxs=[1, 0]),
        np.array([[6.0, 5.0], [2.0, 1.0]]),
    )
    selected_wide = lengths.data_frame(
        groups=["gamma", "alpha"],
        length_bin_idxs=[1, 0],
        keep_wide=True,
    )
    pd.testing.assert_frame_equal(
        selected_wide,
        pd.DataFrame(
            {
                "group_idx": np.array([2, 0], dtype=np.int32),
                "group_name": np.array(["gamma", "alpha"], dtype=object),
                "eligible_windows": np.array([5, 2], dtype=np.int64),
                "count_31_40": np.array([6.0, 2.0], dtype=np.float64),
                "count_30": np.array([5.0, 1.0], dtype=np.float64),
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
        lengths.window_metadata(),
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
    with pytest.raises(
        ValueError, match="Global length-count output must contain exactly one row"
    ):
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
        lengths.window_metadata(),
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
    with pytest.raises(
        ValueError, match="Multiple length-count bins contain length 45"
    ):
        overlapping.length_bin_idx(45)
