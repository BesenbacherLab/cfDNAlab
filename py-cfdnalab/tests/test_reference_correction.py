from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import zarr

import cfdnalab
from test_ends import (
    _create_array,
    _create_labeled_axis,
    _create_motif_axis,
    _write_dense_global_store,
    _write_dense_global_motif_group_store,
    _write_dense_grouped_store,
    _write_dense_window_store,
    _write_reference_correction_ref_kmer_store,
    _write_sparse_window_store,
)
from test_ref_kmers import (
    _write_dense_global_motif_group_store as _write_dense_global_motif_group_ref_kmer_store,
)


def test_reference_corrected_end_motifs_keep_count_scale(tmp_path: Path) -> None:
    # Arrange: There are three reference motifs, so each reference frequency is
    # multiplied by 3 before correcting counts. A uniform reference would have
    # scale 1 and leave counts unchanged.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(ref_kmers=ref_kmers)

    # Assert
    assert "corrected_count" in corrected
    assert "corrected_frequency" in corrected
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 5.0 / 3.0, 1.0 / 3.0, 16.0 / 3.0, 0.0]),
    )


def test_two_sided_correction_requires_reference_kmers(tmp_path: Path) -> None:
    end_path = _write_dense_global_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(end_path)

    with pytest.raises(ValueError, match="two_sided_correction requires ref_kmers"):
        ends.data_frame(two_sided_correction="joint")


def test_motif_group_reference_correction_rejects_two_sided_mode(
    tmp_path: Path,
) -> None:
    end_path = _write_dense_global_motif_group_store(
        tmp_path / "sample.end_motifs.zarr"
    )
    ref_path = _write_dense_global_motif_group_ref_kmer_store(
        tmp_path / "reference.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    for mode in ["joint", "split", "outside", "inside"]:
        with pytest.raises(ValueError, match="Motif-group"):
            ends.data_frame(ref_kmers=ref_kmers, two_sided_correction=mode)


def test_reference_corrected_end_motifs_rejects_canonical_reference(
    tmp_path: Path,
) -> None:
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "reference.ref_kmers.zarr",
        motifs=np.array(["AA", "CC"], dtype=object),
        frequencies=np.array([[0.5, 0.5], [0.5, 0.5]], dtype=np.float64),
    )
    ref_root = zarr.open_group(str(ref_path), mode="a", zarr_format=3)
    ref_root.attrs["canonical"] = True
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    with pytest.raises(ValueError, match="non-canonical"):
        ends.data_frame(ref_kmers=ref_kmers)


def test_one_sided_reference_correction_rejects_explicit_mode(
    tmp_path: Path,
) -> None:
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "reference.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    with pytest.raises(ValueError, match="One-sided"):
        ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="joint")


def test_reference_corrected_end_motifs_selectors_match_filtered_full_correction(
    tmp_path: Path,
) -> None:
    # Arrange: Motif selection should be an output filter, not a request to
    # renormalize the reference over only the selected motif. The selected row
    # should therefore keep the same scale as it has in the full correction.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)
    full_correction = ends.data_frame(ref_kmers=ref_kmers)

    # Act
    selected = ends.data_frame(
        ref_kmers=ref_kmers,
        window_idxs=1,
        motifs="_CC",
    )

    # Assert
    expected = full_correction.loc[
        (full_correction["window_idx"] == 1) & (full_correction["motif"] == "_CC")
    ].reset_index(drop=True)
    pd.testing.assert_frame_equal(selected, expected)
    np.testing.assert_allclose(
        selected["corrected_count"].to_numpy(),
        np.array([16.0 / 3.0]),
    )


def test_reference_corrected_end_motifs_blacklist_filter_uses_end_rows(
    tmp_path: Path,
) -> None:
    # Arrange: The end output keeps window 1 at max_blacklisted_fraction=0.1.
    # The reference output has the opposite blacklist fractions, so applying
    # the reference blacklist independently would select the wrong row.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        blacklisted_fraction=np.array([0.0, 0.25], dtype=np.float64),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(
        ref_kmers=ref_kmers,
        max_blacklisted_fraction=0.1,
    )

    # Assert
    np.testing.assert_array_equal(
        corrected["window_idx"].to_numpy(),
        np.array([1, 1, 1], dtype=np.int32),
    )
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0 / 3.0, 16.0 / 3.0, 0.0], dtype=np.float64),
    )


def test_reference_corrected_end_motif_count_extractors_keep_selector_shape(
    tmp_path: Path,
) -> None:
    # Arrange: Corrected arrays and sparse matrices should use the selected
    # row and motif order, while keeping the same values as the corrected data
    # frame.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)
    expected = np.array(
        [[16.0 / 3.0, 1.0 / 3.0], [0.0, 1.0]],
        dtype=np.float64,
    )

    # Act
    corrected_array = ends.corrected_counts_array(
        ref_kmers,
        window_idxs=[1, 0],
        motifs=["_CC", "_AA"],
    )
    corrected_matrix = ends.sparse_corrected_counts_matrix(
        ref_kmers,
        window_idxs=[1, 0],
        motifs=["_CC", "_AA"],
    )

    # Assert
    np.testing.assert_allclose(corrected_array, expected)
    np.testing.assert_allclose(corrected_matrix.toarray(), expected)


def test_reference_corrected_global_end_motifs_keep_count_scale(
    tmp_path: Path,
) -> None:
    # Arrange: A matched global reference uses the same per-row scaling rule as
    # matched windowed and grouped references.
    end_path = _write_dense_global_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array([[1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0]], dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(ref_kmers=ref_kmers)

    # Assert
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 5.0 / 3.0], dtype=np.float64),
    )


def test_sparse_reference_corrected_end_motif_matrix_uses_sparse_input(
    tmp_path: Path,
) -> None:
    # Arrange: The sparse corrected matrix should not require densifying the
    # end-motif store. It should still return the full selected matrix shape.
    end_path = _write_sparse_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected_matrix = ends.sparse_corrected_counts_matrix(ref_kmers)

    # Assert
    np.testing.assert_allclose(
        corrected_matrix.toarray(),
        np.array([[1.0, 0.0, 5.0 / 3.0], [0.0, 16.0 / 3.0, 0.0]]),
    )
    with pytest.raises(ValueError, match="sparse_corrected_counts_matrix"):
        ends.corrected_counts_array(ref_kmers)


def test_sparse_reference_corrected_end_motif_matrix_uses_sparse_reference_support(
    tmp_path: Path,
) -> None:
    # Arrange: Both inputs are sparse. Each reference row has two supported
    # motifs, so a uniform reference over that row's sparse support keeps count
    # scale unchanged for supported motifs.
    end_path = _write_sparse_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array(
            [[1.0 / 2.0, 0.0, 1.0 / 2.0], [1.0 / 2.0, 1.0 / 2.0, 0.0]],
            dtype=np.float64,
        ),
        storage_mode="sparse_coo",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected_matrix = ends.sparse_corrected_counts_matrix(ref_kmers)

    # Assert
    np.testing.assert_allclose(
        corrected_matrix.toarray(),
        np.array([[1.0, 0.0, 2.5], [0.0, 4.0, 0.0]], dtype=np.float64),
    )


def test_reference_corrected_end_motifs_uses_row_sparse_reference_support(
    tmp_path: Path,
) -> None:
    # Arrange: The sparse reference omits CC in the first window, so that row
    # has two positive reference motifs. The second row has all three. The
    # correction motif count must follow the row's reference support, because a
    # uniform reference over that support should leave counts unchanged.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array(
            [[1.0 / 2.0, 0.0, 1.0 / 2.0], [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]],
            dtype=np.float64,
        ),
        storage_mode="sparse_coo",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(ref_kmers=ref_kmers)

    # Assert
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 2.5, 0.5, 4.0, 0.0]),
    )


def test_reference_corrected_end_motifs_rejects_positive_count_at_zero_reference(
    tmp_path: Path,
) -> None:
    # Arrange: CC has count 4 in the second window, but the matching reference
    # frequency is zero. Dividing by zero would make the correction undefined.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array(
            [[1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0], [1.0 / 2.0, 0.0, 1.0 / 2.0]],
            dtype=np.float64,
        ),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(
        ValueError,
        match="Positive-count end motifs have no positive reference-based correction factor",
    ):
        ends.data_frame(ref_kmers=ref_kmers)


def test_reference_corrected_end_motifs_rejects_non_finite_corrected_counts(
    tmp_path: Path,
) -> None:
    # Arrange: AA has count 1 and the smallest positive float as its reference
    # frequency. Multiplying by the three-motif support still leaves a
    # denominator so small that dividing 1 by it overflows.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    smallest_positive_float = np.nextafter(0.0, 1.0)
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array(
            [
                [smallest_positive_float, 0.5, 0.5],
                [smallest_positive_float, 0.5, 0.5],
            ],
            dtype=np.float64,
        ),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(ValueError, match="non-finite corrected counts.*_AA"):
        ends.data_frame(ref_kmers=ref_kmers)


def test_reference_corrected_end_motifs_rejects_missing_reference_motif_by_default(
    tmp_path: Path,
) -> None:
    # Arrange: CC is observed in the end-motif output, but it is not present in
    # the reference motif axis. There is no reference frequency to divide by.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        motifs=np.array(["AA", "GG"], dtype=object),
        frequencies=np.array([[0.5, 0.5], [0.5, 0.5]], dtype=np.float64),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(
        ValueError,
        match="Pass unsupported_motifs='drop'",
    ):
        ends.data_frame(ref_kmers=ref_kmers)


def test_reference_corrected_end_motifs_can_drop_unsupported_motifs(
    tmp_path: Path,
) -> None:
    # Arrange: CC is outside the reference motif axis. The drop policy returns
    # only motifs with a positive reference frequency in the matched row.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        motifs=np.array(["AA", "GG"], dtype=object),
        frequencies=np.array([[0.5, 0.5], [0.5, 0.5]], dtype=np.float64),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(
        ref_kmers=ref_kmers,
        unsupported_motifs="drop",
    )

    # Assert
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0, 2.5, 0.5, 0.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["corrected_frequency"].to_numpy(),
        np.array([2.0 / 7.0, 5.0 / 7.0, 1.0, 0.0], dtype=np.float64),
    )
    assert "_CC" not in corrected["motif"].tolist()


def test_reference_corrected_end_motif_matrix_extractors_reject_drop_policy(
    tmp_path: Path,
) -> None:
    # Arrange: Dropping unsupported motifs changes the motif axis length, which
    # fixed-shape arrays and sparse matrices cannot represent.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        motifs=np.array(["AA", "GG"], dtype=object),
        frequencies=np.array([[0.5, 0.5], [0.5, 0.5]], dtype=np.float64),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(ValueError, match="fixed-shape corrected_counts_array"):
        ends.corrected_counts_array(ref_kmers, unsupported_motifs="drop")
    with pytest.raises(ValueError, match="fixed-shape sparse_corrected_counts_matrix"):
        ends.sparse_corrected_counts_matrix(ref_kmers, unsupported_motifs="drop")


def test_reference_corrected_end_motifs_can_keep_unsupported_motifs_as_nan(
    tmp_path: Path,
) -> None:
    # Arrange: CC has no reference frequency. Zero-count unsupported rows stay
    # zero, while positive unsupported rows are marked as undefined.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        motifs=np.array(["AA", "GG"], dtype=object),
        frequencies=np.array([[0.5, 0.5], [0.5, 0.5]], dtype=np.float64),
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(
        ref_kmers=ref_kmers,
        unsupported_motifs="keep_na",
    )

    # Assert
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 2.5, 0.5, np.nan, 0.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["corrected_frequency"].to_numpy(),
        np.array([2.0 / 7.0, 0.0, 5.0 / 7.0, np.nan, np.nan, np.nan]),
    )


def test_reference_corrected_end_motifs_requires_opt_in_for_global_reference(
    tmp_path: Path,
) -> None:
    # Arrange
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array([[1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0]], dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(ValueError, match="Pass use_global_bias=True"):
        ends.data_frame(ref_kmers=ref_kmers)


def test_sparse_reference_correction_validates_rows_for_empty_motif_selection(
    tmp_path: Path,
) -> None:
    # Arrange: The first reference interval starts at 11 instead of 10, so the
    # full end and reference row axes do not match. An empty motif selection
    # still has to validate the inputs before returning a 2-by-0 matrix.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ref_root = zarr.open_group(str(ref_path), mode="a", zarr_format=3)
    ref_root["row_start_bp"][:] = np.array([11, 40], dtype=np.int64)
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(ValueError, match="rows do not match"):
        ends.sparse_corrected_counts_matrix(ref_kmers, motifs=[])


def test_sparse_reference_correction_validates_global_bias_for_empty_selection(
    tmp_path: Path,
) -> None:
    # Arrange
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(TypeError, match="use_global_bias must be a bool"):
        ends.sparse_corrected_counts_matrix(
            ref_kmers,
            motifs=[],
            use_global_bias="yes",  # type: ignore[arg-type]
        )


def test_reference_corrected_end_motifs_can_use_global_reference_bias(
    tmp_path: Path,
) -> None:
    # Arrange: The same global reference frequencies are applied to both
    # windows. With three reference motifs, the scales are 1, 0.5, and 1.5.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array([[1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0]], dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(
        ref_kmers=ref_kmers,
        use_global_bias=True,
    )
    selected = ends.data_frame(
        ref_kmers=ref_kmers,
        window_idxs=1,
        motifs="_CC",
        use_global_bias=True,
    )
    filtered = ends.data_frame(
        ref_kmers=ref_kmers,
        max_blacklisted_fraction=0.1,
        use_global_bias=True,
    )

    # Assert
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 5.0 / 3.0, 0.5, 8.0, 0.0]),
    )
    expected = corrected.loc[
        (corrected["window_idx"] == 1) & (corrected["motif"] == "_CC")
    ].reset_index(drop=True)
    pd.testing.assert_frame_equal(selected, expected)
    expected_filtered = corrected.loc[corrected["window_idx"] == 1].reset_index(
        drop=True
    )
    pd.testing.assert_frame_equal(filtered, expected_filtered)
    with pytest.raises(TypeError, match="unexpected keyword argument 'groups'"):
        ends.data_frame(  # type: ignore[call-arg]
            ref_kmers=ref_kmers,
            groups="A",
            use_global_bias=True,
        )


def test_reference_corrected_end_motifs_rejects_global_bias_for_matched_reference(
    tmp_path: Path,
) -> None:
    # Arrange: The reference already has rows matching the end output. The
    # global-bias flag would be misleading because there is no global reference
    # row to broadcast.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(
        ValueError,
        match="use_global_bias=True requires a global reference k-mer output",
    ):
        ends.data_frame(ref_kmers=ref_kmers, use_global_bias=True)


def test_reference_corrected_group_indices_match_reference_rows_by_group_name(
    tmp_path: Path,
) -> None:
    # Arrange: `group_idxs=1` selects `long_group` in the end output. The
    # reference has the same group names in a different order, so correction
    # must map by group name before reading reference motif rows.
    end_path = _write_dense_grouped_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmers.zarr",
        frequencies=np.array(
            [
                [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                [0.5, 0.25, 0.25],
                [1.0, 0.0, 0.0],
            ],
            dtype=np.float64,
        ),
        group_names=np.array(["long_group", "A", "mid"], dtype=object),
        row_mode="grouped_bed",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(
        ref_kmers=ref_kmers,
        group_idxs=1,
        motifs="_CC",
    )

    # Assert
    assert corrected["group_name"].tolist() == ["long_group"]
    np.testing.assert_allclose(
        corrected["corrected_count"].to_numpy(),
        np.array([4.25], dtype=np.float64),
    )


def test_two_sided_reference_correction_modes_use_mode_axis(tmp_path: Path) -> None:
    # Arrange: This fixture is shared with the Rust and R correction tests.
    # Every joint and side reference frequency is positive and differs from its
    # uniform baseline, so every correction mode has a visible effect.
    end_path = _write_two_sided_global_end_motif_store(
        tmp_path / "sample.end_motifs.zarr"
    )
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "reference.ref_kmers.zarr",
        motifs=np.array(["AC", "AG", "TC", "TG"], dtype=object),
        frequencies=np.array(
            [[1.0 / 8.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 2.0]],
            dtype=np.float64,
        ),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # A two-sided motif axis has no single implied correction. The caller must
    # choose joint motifs, split side factors, or a derived side axis.
    with pytest.raises(ValueError, match="two-sided"):
        ends.data_frame(ref_kmers=ref_kmers)

    joint = ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="joint")
    # Four positive reference motifs make the uniform frequency 1/4. Relative
    # to uniform, frequencies [1/8, 1/8, 1/4, 1/2] give correction factors
    # [1/2, 1/2, 1, 2] for [AC, AG, TC, TG]. Dividing original counts
    # [2, 4, 6, 8] by those factors gives [4, 8, 6, 4]. Their total is 22, so
    # dividing each corrected count by 22 gives [2/11, 4/11, 3/11, 2/11].
    assert joint["motif"].tolist() == ["A_C", "A_G", "T_C", "T_G"]
    np.testing.assert_allclose(joint["count"], [2.0, 4.0, 6.0, 8.0])
    np.testing.assert_allclose(joint["corrected_count"], [4.0, 8.0, 6.0, 4.0])
    np.testing.assert_allclose(
        joint["corrected_frequency"],
        [2.0 / 11.0, 4.0 / 11.0, 3.0 / 11.0, 2.0 / 11.0],
    )

    split = ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="split")
    assert split["motif"].tolist() == ["A_C", "A_G", "T_C", "T_G"]
    assert "corrected_count" in split
    assert "corrected_frequency" in split
    assert "reference_frequency" not in split
    np.testing.assert_allclose(split["count"], [2.0, 4.0, 6.0, 8.0])
    # Two positive labels on each side make each side's uniform frequency 1/2.
    # Outside frequencies A=1/4 and T=3/4 give factors 1/2 and 3/2. Inside
    # frequencies C=3/8 and G=5/8 give factors 3/4 and 5/4. Multiplying matching
    # side factors gives [3/8, 5/8, 9/8, 15/8] for [A_C, A_G, T_C, T_G].
    # Dividing original counts [2, 4, 6, 8] by those factors gives
    # [16/3, 32/5, 16/3, 64/15]. The corrected counts total
    # 64/3, so normalization gives [1/4, 3/10, 1/4, 1/5].
    np.testing.assert_allclose(
        split["corrected_count"].to_numpy(),
        np.array(
            [16.0 / 3.0, 32.0 / 5.0, 16.0 / 3.0, 64.0 / 15.0],
            dtype=np.float64,
        ),
    )
    np.testing.assert_allclose(
        split["corrected_frequency"].to_numpy(),
        np.array([1.0 / 4.0, 3.0 / 10.0, 1.0 / 4.0, 1.0 / 5.0]),
    )

    outside = ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="outside")
    assert outside["motif"].tolist() == ["A_", "T_"]
    # Counts aggregate to [6, 14]. Relative to the uniform outside frequency
    # 1/2, reference frequencies [1/4, 3/4] give factors [1/2, 3/2]. Dividing
    # the aggregated counts by them gives [12, 28/3]. These total 64/3, so
    # normalization gives frequencies [9/16, 7/16].
    np.testing.assert_allclose(
        outside["count"].to_numpy(),
        np.array([6.0, 14.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        outside["corrected_count"].to_numpy(),
        np.array([12.0, 28.0 / 3.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        outside["corrected_frequency"].to_numpy(),
        np.array([9.0 / 16.0, 7.0 / 16.0], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.corrected_motifs_metadata(
            ref_kmers,
            two_sided_correction="outside",
        ),
        pd.DataFrame(
            {
                "matrix_column": [0, 1],
                "motif_index": [0, 1],
                "motif": ["A_", "T_"],
            }
        ),
    )
    np.testing.assert_allclose(
        ends.corrected_counts_array(ref_kmers, two_sided_correction="outside"),
        np.array([[12.0, 28.0 / 3.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        ends.sparse_corrected_counts_matrix(
            ref_kmers,
            two_sided_correction="outside",
        ).toarray(),
        np.array([[12.0, 28.0 / 3.0]], dtype=np.float64),
    )

    inside = ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="inside")
    # Counts aggregate to [8, 12]. Relative to the uniform inside frequency
    # 1/2, reference frequencies [3/8, 5/8] give factors [3/4, 5/4]. Dividing
    # the aggregated counts by them gives [32/3, 48/5]. These total 304/15, so
    # normalization gives frequencies [10/19, 9/19].
    assert inside["motif"].tolist() == ["_C", "_G"]
    np.testing.assert_allclose(inside["count"], [8.0, 12.0])
    np.testing.assert_allclose(
        inside["corrected_count"],
        np.array([32.0 / 3.0, 48.0 / 5.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        inside["corrected_frequency"],
        np.array([10.0 / 19.0, 9.0 / 19.0], dtype=np.float64),
    )
    selected_inside = ends.data_frame(
        ref_kmers=ref_kmers,
        two_sided_correction="inside",
        motifs="_G",
    )
    expected_inside = inside.loc[inside["motif"] == "_G"].reset_index(drop=True)
    pd.testing.assert_frame_equal(selected_inside, expected_inside)
    np.testing.assert_allclose(
        selected_inside["count"].to_numpy(),
        np.array([12.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        selected_inside["corrected_count"].to_numpy(),
        np.array([48.0 / 5.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        selected_inside["corrected_frequency"].to_numpy(),
        np.array([9.0 / 19.0], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.corrected_motifs_metadata(
            ref_kmers,
            two_sided_correction="inside",
            motifs="_G",
        ),
        pd.DataFrame({"matrix_column": [0], "motif_index": [1], "motif": ["_G"]}),
    )
    with pytest.raises(ValueError, match="motif index selectors"):
        ends.data_frame(
            ref_kmers=ref_kmers,
            two_sided_correction="outside",
            motif_idxs=0,
        )
    with pytest.raises(ValueError, match="Side-mode motif axis"):
        ends.data_frame(
            ref_kmers=ref_kmers,
            two_sided_correction="outside",
            motifs="A_C",
        )


def test_side_correction_applies_unsupported_policy_after_aggregation(
    tmp_path: Path,
) -> None:
    # Arrange: The G side has a positive aggregated sample count but no
    # reference support. The C side remains supported.
    end_path = _write_two_sided_global_end_motif_store(
        tmp_path / "sample.end_motifs.zarr"
    )
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "reference.ref_kmers.zarr",
        motifs=np.array(["AC"], dtype=object),
        frequencies=np.array([[1.0]], dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # After inside aggregation, _G has a positive sample count but no positive
    # reference correction factor. The default policy must name it as an error.
    with pytest.raises(ValueError, match="_G"):
        ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="inside")

    dropped = ends.data_frame(
        ref_kmers=ref_kmers,
        two_sided_correction="inside",
        unsupported_motifs="drop",
    )
    assert dropped["motif"].tolist() == ["_C"]
    np.testing.assert_allclose(dropped["corrected_count"], np.array([8.0]))
    np.testing.assert_allclose(dropped["corrected_frequency"], np.array([1.0]))

    kept = ends.data_frame(
        ref_kmers=ref_kmers,
        two_sided_correction="inside",
        unsupported_motifs="keep_na",
    )
    assert kept["motif"].tolist() == ["_C", "_G"]
    np.testing.assert_allclose(kept["corrected_count"], np.array([8.0, np.nan]))
    assert kept["corrected_frequency"].isna().all()


def test_corrected_frequencies_are_zero_when_corrected_total_is_zero(
    tmp_path: Path,
) -> None:
    end_path = _write_two_sided_global_end_motif_store(
        tmp_path / "sample.end_motifs.zarr",
        counts=np.zeros((1, 4), dtype=np.float64),
    )
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "reference.ref_kmers.zarr",
        motifs=np.array(["AC", "AG", "TC", "TG"], dtype=object),
        frequencies=np.full((1, 4), 0.25, dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    corrected = ends.data_frame(
        ref_kmers=ref_kmers,
        two_sided_correction="split",
    )

    np.testing.assert_allclose(corrected["corrected_count"], np.zeros(4))
    np.testing.assert_allclose(corrected["corrected_frequency"], np.zeros(4))


def _write_two_sided_global_end_motif_store(
    path: Path,
    counts: np.ndarray | None = None,
) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 2
    root.attrs["storage_mode"] = "dense"
    root.attrs["row_mode"] = "global"
    root.attrs["motif_axis_kind"] = "motif"

    _create_motif_axis(root, np.array(["A_C", "A_G", "T_C", "T_G"], dtype=object))
    _create_labeled_axis(
        root,
        "row",
        np.array([0], dtype=np.int32),
        "row_label",
        np.array(["global"], dtype=object),
    )
    _create_array(
        root,
        "counts",
        (
            np.array([[2.0, 4.0, 6.0, 8.0]], dtype=np.float64)
            if counts is None
            else counts
        ),
        chunks=(1, 4),
        dimension_names=("row", "motif"),
    )
    return path
