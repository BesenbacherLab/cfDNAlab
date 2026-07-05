from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cfdnalab
from test_ends import (
    _write_dense_global_store,
    _write_dense_grouped_store,
    _write_dense_window_store,
    _write_reference_correction_ref_kmer_store,
    _write_sparse_window_store,
)


def test_reference_corrected_end_motifs_keep_count_scale(tmp_path: Path) -> None:
    # Arrange: There are three reference motifs, so each reference frequency is
    # multiplied by 3 before correcting counts. A uniform reference would have
    # scale 1 and leave counts unchanged.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr"
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(ref_kmers=ref_kmers)

    # Assert
    assert corrected["reference_motif"].tolist() == ["AA", "CC", "GG"] * 2
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([3, 3, 3, 3, 3, 3], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_scale"].to_numpy(),
        np.array([1.0, 0.5, 1.5, 1.5, 0.75, 0.75], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 5.0 / 3.0, 1.0 / 3.0, 16.0 / 3.0, 0.0]),
    )


def test_reference_corrected_end_motifs_selectors_match_filtered_full_correction(
    tmp_path: Path,
) -> None:
    # Arrange: Motif selection should be an output filter, not a request to
    # renormalize the reference over only the selected motif. The selected row
    # should therefore keep the same scale as it has in the full correction.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr"
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
    np.testing.assert_array_equal(
        selected["correction_motif_count"].to_numpy(),
        np.array([3], dtype=np.int64),
    )
    np.testing.assert_allclose(
        selected["reference_corrected_count"].to_numpy(),
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
        corrected["reference_corrected_count"].to_numpy(),
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
        tmp_path / "hg38.ref_kmer_counts.zarr"
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
        frequencies=np.array([[1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0]], dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act
    corrected = ends.data_frame(ref_kmers=ref_kmers)

    # Assert
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([3, 3, 3], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_scale"].to_numpy(),
        np.array([1.0, 0.5, 1.5], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 5.0 / 3.0], dtype=np.float64),
    )


def test_sparse_reference_corrected_end_motif_matrix_uses_sparse_input(
    tmp_path: Path,
) -> None:
    # Arrange: The sparse corrected matrix should not require densifying the
    # end-motif store. It should still return the full selected matrix shape.
    end_path = _write_sparse_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr"
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([2, 2, 2, 3, 3, 3], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_scale"].to_numpy(),
        np.array([1.0, 0.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 2.5, 0.5, 4.0, 0.0]),
    )


def test_reference_corrected_end_motifs_rejects_positive_count_at_zero_reference(
    tmp_path: Path,
) -> None:
    # Arrange: CC has count 4 in the second window, but the matching reference
    # frequency is zero. Dividing by zero would make the correction undefined.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
        match="Positive-count end motifs have no positive reference frequency",
    ):
        ends.data_frame(ref_kmers=ref_kmers)


def test_reference_corrected_end_motifs_rejects_missing_reference_motif_by_default(
    tmp_path: Path,
) -> None:
    # Arrange: CC is observed in the end-motif output, but it is not present in
    # the reference motif axis. There is no reference frequency to divide by.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
    assert corrected["reference_motif"].tolist() == ["AA", "GG", "AA", "GG"]
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([2, 2, 2, 2], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_scale"].to_numpy(),
        np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
        np.array([1.0, 2.5, 0.5, 0.0], dtype=np.float64),
    )


def test_reference_corrected_end_motif_matrix_extractors_reject_drop_policy(
    tmp_path: Path,
) -> None:
    # Arrange: Dropping unsupported motifs changes the motif axis length, which
    # fixed-shape arrays and sparse matrices cannot represent.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([2, 2, 2, 2, 2, 2], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_scale"].to_numpy(),
        np.array([1.0, 0.0, 1.0, 1.0, 0.0, 1.0], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
        np.array([1.0, 0.0, 2.5, 0.5, np.nan, 0.0], dtype=np.float64),
    )


def test_reference_corrected_end_motifs_requires_opt_in_for_global_reference(
    tmp_path: Path,
) -> None:
    # Arrange
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr",
        frequencies=np.array([[1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0]], dtype=np.float64),
        row_mode="global",
    )
    ends = cfdnalab.read_end_motifs(end_path)
    ref_kmers = cfdnalab.read_ref_kmers(ref_path)

    # Act / Assert
    with pytest.raises(ValueError, match="Pass use_global_bias=True"):
        ends.data_frame(ref_kmers=ref_kmers)


def test_reference_corrected_end_motifs_can_use_global_reference_bias(
    tmp_path: Path,
) -> None:
    # Arrange: The same global reference frequencies are applied to both
    # windows. With three reference motifs, the scales are 1, 0.5, and 1.5.
    end_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ref_path = _write_reference_correction_ref_kmer_store(
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([3, 3, 3, 3, 3, 3], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_scale"].to_numpy(),
        np.array([1.0, 0.5, 1.5, 1.0, 0.5, 1.5], dtype=np.float64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
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
        tmp_path / "hg38.ref_kmer_counts.zarr"
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
        tmp_path / "hg38.ref_kmer_counts.zarr",
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
        corrected["reference_frequency"].to_numpy(),
        np.array([1.0 / 3.0], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        corrected["correction_motif_count"].to_numpy(),
        np.array([3], dtype=np.int64),
    )
    np.testing.assert_allclose(
        corrected["reference_corrected_count"].to_numpy(),
        np.array([4.25], dtype=np.float64),
    )
