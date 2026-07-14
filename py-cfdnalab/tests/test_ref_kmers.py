from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sparse
import zarr

import cfdnalab


MOTIF_NAMES = np.array(["AA", "AC", "GT"], dtype=object)
MOTIF_GROUP_NAMES = np.array(["left", "right"], dtype=object)
REFERENCE_FOOTPRINT = [
    {"name": "chr2", "length": 100},
    {"name": "chr10", "length": 120},
]


def test_dense_windowed_ref_kmers_load_metadata_and_arrays(tmp_path: Path) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.ref_kmers.zarr")

    ref_kmers = cfdnalab.read_ref_kmers(store_path)

    assert ref_kmers.path == store_path
    assert isinstance(ref_kmers, cfdnalab.WindowedRefKmerFrequencies)
    assert repr(ref_kmers) == (
        "WindowedRefKmerFrequencies("
        f"path={str(store_path)!r}, "
        "schema_version=1, "
        "storage_mode='dense', "
        "row_mode='bed', "
        "motif_axis_kind='motif', "
        "shape=(2, 3)"
        ")"
    )
    assert not hasattr(ref_kmers, "group_metadata")
    assert ref_kmers.storage_mode() == "dense"
    assert ref_kmers.row_mode() == "bed"
    assert ref_kmers.motif_axis_kind() == "motif"
    assert ref_kmers.kmer_size() == 2
    assert not ref_kmers.canonical()
    assert not ref_kmers.all_motifs()
    assert ref_kmers.assign_by() == "count-overlap"
    assert ref_kmers.reference_contig_footprint() == REFERENCE_FOOTPRINT
    assert ref_kmers.has_motif("AC")
    assert not ref_kmers.has_motif("TT")
    assert ref_kmers.motif_idx("GT") == 2
    pd.testing.assert_frame_equal(
        ref_kmers.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": np.array([0, 1, 2], dtype=np.int32),
                "motif": MOTIF_NAMES,
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.window_metadata(),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chrom": np.array(["chr2", "chr10"], dtype=object),
                "start": np.array([10, 40], dtype=np.int64),
                "end": np.array([20, 60], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.0], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.row_scaling_factors(),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chrom": np.array(["chr2", "chr10"], dtype=object),
                "start": np.array([10, 40], dtype=np.int64),
                "end": np.array([20, 60], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.0], dtype=np.float64),
                "row_scaling_factor": np.array([4.0, 2.0], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_frequencies_array(),
        np.array([[0.25, 0.0, 0.75], [0.5, 0.5, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(),
        np.array([[1.0, 0.0, 3.0], [1.0, 1.0, 0.0]], dtype=np.float64),
    )
    assert ref_kmers.dense_frequencies_zarr_array().shape == (2, 3)
    assert sparse.isspmatrix_coo(ref_kmers.sparse_frequencies_matrix())


def test_dense_windowed_ref_kmer_selectors_return_expected_frames(
    tmp_path: Path,
) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.ref_kmers.zarr")
    ref_kmers = cfdnalab.read_ref_kmers(store_path)

    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(window_idxs=0),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 0, 0], dtype=np.int32),
                "chrom": np.array(["chr2", "chr2", "chr2"], dtype=object),
                "start": np.array([10, 10, 10], dtype=np.int64),
                "end": np.array([20, 20, 20], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.25, 0.25]),
                "motif_index": np.array([0, 1, 2], dtype=np.int32),
                "motif": MOTIF_NAMES,
                "frequency": np.array([0.25, 0.0, 0.75], dtype=np.float64),
                "count": np.array([1.0, 0.0, 3.0], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(motifs="AC"),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chrom": np.array(["chr2", "chr10"], dtype=object),
                "start": np.array([10, 40], dtype=np.int64),
                "end": np.array([20, 60], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.0]),
                "motif_index": np.array([1, 1], dtype=np.int32),
                "motif": ["AC", "AC"],
                "frequency": np.array([0.0, 0.5], dtype=np.float64),
                "count": np.array([0.0, 1.0], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_frequencies_array(window_idxs=1),
        np.array([[0.5, 0.5, 0.0]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(motifs="GT"),
        np.array([[3.0], [0.0]], dtype=np.float64),
    )
    filtered = ref_kmers.data_frame(motifs="AC", max_blacklisted_fraction=0.1)
    assert filtered["window_idx"].tolist() == [1]
    assert filtered["count"].tolist() == [1.0]
    with pytest.raises(
        ValueError, match="max_blacklisted_fraction must be a single finite fraction"
    ):
        ref_kmers.data_frame(window_idxs=0, max_blacklisted_fraction=1.1)
    with pytest.raises(ValueError, match="Use either motifs or motif_idxs"):
        ref_kmers.data_frame(motifs="AA", motif_idxs=0)
    with pytest.raises(ValueError, match="window_idxs contains duplicate values"):
        ref_kmers.data_frame(window_idxs=[0, 0])


def test_sparse_grouped_ref_kmers_reconstruct_dense_counts_and_metadata(
    tmp_path: Path,
) -> None:
    store_path = _write_sparse_grouped_store(tmp_path / "sample.ref_kmers.zarr")

    ref_kmers = cfdnalab.read_ref_kmers(store_path)

    assert isinstance(ref_kmers, cfdnalab.GroupedRefKmerFrequencies)
    assert not hasattr(ref_kmers, "window_metadata")
    assert ref_kmers.storage_mode() == "sparse_coo"
    assert ref_kmers.row_mode() == "grouped_bed"
    assert ref_kmers.group_idx("long_group") == 1
    pd.testing.assert_frame_equal(
        ref_kmers.group_metadata(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2], dtype=np.int32),
                "group_name": np.array(["A", "long_group", "empty"], dtype=object),
                "eligible_windows": np.array([1, 2, 0], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.125, 0.0], dtype=np.float64),
            }
        ),
    )
    with pytest.raises(ValueError, match="would turn sparse reference k-mer output"):
        ref_kmers.dense_frequencies_array()
    np.testing.assert_allclose(
        ref_kmers.dense_frequencies_array(allow_densify=True),
        np.array(
            [[0.25, 0.0, 0.75], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            dtype=np.float64,
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(allow_densify=True),
        np.array(
            [[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]],
            dtype=np.float64,
        ),
    )
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(groups="A", allow_densify=True),
        np.array([[1.0, 0.0, 3.0]], dtype=np.float64),
    )
    motif_coo = ref_kmers.sparse_counts_matrix(motifs="GT")
    assert motif_coo.shape == (3, 1)
    np.testing.assert_array_equal(motif_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(motif_coo.col, np.array([0], dtype=np.int32))
    np.testing.assert_allclose(motif_coo.data, np.array([3.0]))


def test_sparse_grouped_ref_kmer_frames_preserve_selected_order(
    tmp_path: Path,
) -> None:
    store_path = _write_sparse_grouped_store(tmp_path / "sample.ref_kmers.zarr")
    ref_kmers = cfdnalab.read_ref_kmers(store_path)

    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 0, 1], dtype=np.int32),
                "group_name": np.array(["A", "A", "long_group"], dtype=object),
                "eligible_windows": np.array([1, 1, 2], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.125]),
                "motif_index": np.array([0, 2, 1], dtype=np.int32),
                "motif": np.array(["AA", "GT", "AC"], dtype=object),
                "frequency": np.array([0.25, 0.75, 1.0], dtype=np.float64),
                "count": np.array([1.0, 3.0, 2.0], dtype=np.float64),
            }
        ),
    )
    selected = ref_kmers.data_frame(
        groups=["long_group", "A"],
        motifs=["GT", "AA"],
        densify=True,
    )
    assert selected["group_name"].tolist() == ["long_group", "long_group", "A", "A"]
    assert selected["motif"].tolist() == ["GT", "AA", "GT", "AA"]
    assert selected["count"].tolist() == [0.0, 0.0, 3.0, 1.0]
    empty_dense = ref_kmers.data_frame(groups="empty", densify=True)
    assert empty_dense["motif"].tolist() == ["AA", "AC", "GT"]
    assert empty_dense["frequency"].tolist() == [0.0, 0.0, 0.0]
    assert empty_dense["count"].tolist() == [0.0, 0.0, 0.0]


def test_dense_global_ref_kmers_support_motif_group_axis(tmp_path: Path) -> None:
    store_path = _write_dense_global_motif_group_store(
        tmp_path / "sample.ref_kmers.zarr"
    )

    ref_kmers = cfdnalab.read_ref_kmers(store_path)

    assert isinstance(ref_kmers, cfdnalab.GlobalRefKmerFrequencies)
    assert ref_kmers.storage_mode() == "dense"
    assert ref_kmers.row_mode() == "global"
    assert ref_kmers.motif_axis_kind() == "motif_group"
    pd.testing.assert_frame_equal(
        ref_kmers.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": np.array([0, 1], dtype=np.int32),
                "motif": MOTIF_GROUP_NAMES,
            }
        ),
    )
    assert ref_kmers.motif_idx("right") == 1
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(motifs="right"),
        np.array([[3.0]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ref_kmers.data_frame(),
        pd.DataFrame(
            {
                "row_label": np.array(["global", "global"], dtype=object),
                "motif_index": np.array([0, 1], dtype=np.int32),
                "motif": MOTIF_GROUP_NAMES,
                "frequency": np.array([0.25, 0.75], dtype=np.float64),
                "count": np.array([1.0, 3.0], dtype=np.float64),
            }
        ),
    )


def test_ref_kmer_loader_rejects_invalid_paths(tmp_path: Path) -> None:
    missing_path = tmp_path / "missing.ref_kmers.zarr"
    file_path = tmp_path / "not_a_directory.zarr"
    file_path.write_text("not a zarr store")
    wrong_suffix_path = tmp_path / "sample.not_zarr"
    wrong_suffix_path.mkdir()

    with pytest.raises(FileNotFoundError, match="does not exist"):
        cfdnalab.read_ref_kmers(missing_path)
    with pytest.raises(NotADirectoryError, match="not a directory"):
        cfdnalab.read_ref_kmers(file_path)
    with pytest.raises(ValueError, match="should end with '.zarr'"):
        cfdnalab.read_ref_kmers(wrong_suffix_path)


def test_ref_kmer_loader_rejects_schema_and_shape_problems(tmp_path: Path) -> None:
    wrong_schema = _write_dense_window_store(
        tmp_path / "wrong_schema.ref_kmers.zarr",
        schema="other",
    )
    wrong_version = _write_dense_window_store(
        tmp_path / "wrong_version.ref_kmers.zarr",
        schema_version=99,
    )
    missing_array = _write_dense_window_store(
        tmp_path / "missing_array.ref_kmers.zarr",
        omit={"row_scaling_factor"},
    )
    wrong_dimensions = _write_dense_window_store(
        tmp_path / "wrong_dimensions.ref_kmers.zarr",
        frequencies_dimension_names=("motif", "row"),
    )
    wrong_motif_axis_dimensions = _write_dense_window_store(
        tmp_path / "wrong_motif_axis_dimensions.ref_kmers.zarr",
    )
    _patch_array_metadata(
        wrong_motif_axis_dimensions,
        "motif_index",
        dimension_names=("row",),
    )
    wrong_window_metadata_dimensions = _write_dense_window_store(
        tmp_path / "wrong_window_metadata_dimensions.ref_kmers.zarr",
    )
    _patch_array_metadata(
        wrong_window_metadata_dimensions,
        "row_start_bp",
        dimension_names=("motif",),
    )
    shape_mismatch = _write_sparse_grouped_store(
        tmp_path / "shape_mismatch.ref_kmers.zarr",
        sparse_shape=np.array([3, 2], dtype=np.int32),
    )
    duplicate_sparse_coordinate = _write_sparse_grouped_store(
        tmp_path / "duplicate_sparse_coordinate.ref_kmers.zarr",
        sparse_row=np.array([0, 0, 0], dtype=np.int32),
        sparse_motif=np.array([0, 0, 2], dtype=np.int32),
    )
    unsorted_sparse_coordinate = _write_sparse_grouped_store(
        tmp_path / "unsorted_sparse_coordinate.ref_kmers.zarr",
        sparse_row=np.array([0, 0, 0], dtype=np.int32),
        sparse_motif=np.array([2, 0, 1], dtype=np.int32),
    )
    float_sparse_coordinate = _write_sparse_grouped_store(
        tmp_path / "float_sparse_coordinate.ref_kmers.zarr",
        sparse_row=np.array([0.0, 0.0, 1.0], dtype=np.float64),
    )
    wrong_sparse_dimension_labels = _write_sparse_grouped_store(
        tmp_path / "wrong_sparse_dimension_labels.ref_kmers.zarr",
        sparse_dimension_labels=np.array(["motif", "row"], dtype=object),
    )

    with pytest.raises(ValueError, match="Expected cfdnalab_schema"):
        cfdnalab.read_ref_kmers(wrong_schema)
    with pytest.raises(ValueError, match="Unsupported reference k-mer schema version"):
        cfdnalab.read_ref_kmers(wrong_version)
    with pytest.raises(ValueError, match="missing arrays: \\['row_scaling_factor'\\]"):
        cfdnalab.read_ref_kmers(missing_array)
    with pytest.raises(ValueError, match="dense frequencies dimensions must be"):
        cfdnalab.read_ref_kmers(wrong_dimensions)
    with pytest.raises(ValueError, match="motif_index dimensions must be"):
        cfdnalab.read_ref_kmers(wrong_motif_axis_dimensions)
    with pytest.raises(ValueError, match="row_start_bp dimensions must be"):
        cfdnalab.read_ref_kmers(wrong_window_metadata_dimensions)
    with pytest.raises(ValueError, match="sparse/shape does not match"):
        cfdnalab.read_ref_kmers(shape_mismatch)
    with pytest.raises(ValueError, match="sorted and unique"):
        cfdnalab.read_ref_kmers(duplicate_sparse_coordinate)
    with pytest.raises(ValueError, match="sorted and unique"):
        cfdnalab.read_ref_kmers(unsorted_sparse_coordinate)
    with pytest.raises(ValueError, match="sparse/row must contain integer values"):
        cfdnalab.read_ref_kmers(float_sparse_coordinate)
    with pytest.raises(ValueError, match="sparse_dimension labels must be"):
        cfdnalab.read_ref_kmers(wrong_sparse_dimension_labels)


def test_ref_kmer_loader_rejects_metadata_that_changes_count_meaning(
    tmp_path: Path,
) -> None:
    wrong_units = _write_dense_window_store(
        tmp_path / "wrong_units.ref_kmers.zarr",
        value_units="other",
    )
    wrong_count_reconstruction = _write_dense_window_store(
        tmp_path / "wrong_reconstruction.ref_kmers.zarr",
        count_reconstruction="count = frequency",
    )
    wrong_scaling_array = _write_dense_window_store(
        tmp_path / "wrong_scaling.ref_kmers.zarr",
        row_scaling_factor_array="other",
    )

    with pytest.raises(ValueError, match="value_units"):
        cfdnalab.read_ref_kmers(wrong_units)
    with pytest.raises(ValueError, match="count_reconstruction"):
        cfdnalab.read_ref_kmers(wrong_count_reconstruction)
    with pytest.raises(ValueError, match="row_scaling_factor_array"):
        cfdnalab.read_ref_kmers(wrong_scaling_array)


def test_ref_kmer_loader_rejects_invalid_values_and_motifs(tmp_path: Path) -> None:
    bad_sparse_frequency = _write_sparse_grouped_store(
        tmp_path / "bad_sparse_frequency.ref_kmers.zarr",
        sparse_frequency=np.array([0.25, 1.25, 1.0], dtype=np.float64),
    )
    bad_scaling = _write_dense_window_store(
        tmp_path / "bad_scaling.ref_kmers.zarr",
        row_scaling_factor=np.array([4.0, np.nan], dtype=np.float64),
    )
    invalid_base = _write_dense_window_store(
        tmp_path / "invalid_base.ref_kmers.zarr",
        motif_names=np.array(["AA", "AN"], dtype=object),
    )
    noncanonical = _write_dense_window_store(
        tmp_path / "noncanonical.ref_kmers.zarr",
        motif_names=np.array(["AA", "GT"], dtype=object),
        canonical=True,
    )
    bad_dense_frequency = _write_dense_window_store(
        tmp_path / "bad_dense_frequency.ref_kmers.zarr",
        frequencies=np.array([[0.25, 1.25, 0.0], [0.5, 0.5, 0.0]]),
    )

    with pytest.raises(ValueError, match="sparse/frequency"):
        cfdnalab.read_ref_kmers(bad_sparse_frequency)
    with pytest.raises(ValueError, match="row_scaling_factor"):
        cfdnalab.read_ref_kmers(bad_scaling)
    with pytest.raises(ValueError, match="invalid bases"):
        cfdnalab.read_ref_kmers(invalid_base)
    with pytest.raises(ValueError, match="canonical reference k-mer motif label"):
        cfdnalab.read_ref_kmers(noncanonical)

    ref_kmers = cfdnalab.read_ref_kmers(bad_dense_frequency)
    with pytest.raises(ValueError, match="selected frequencies"):
        ref_kmers.dense_frequencies_array()


def test_ref_kmer_loader_rejects_invalid_row_metadata(tmp_path: Path) -> None:
    bad_interval = _write_dense_window_store(
        tmp_path / "bad_interval.ref_kmers.zarr",
        row_start_bp=np.array([10, 60], dtype=np.int64),
        row_end_bp=np.array([20, 40], dtype=np.int64),
    )
    bad_window_fraction = _write_dense_window_store(
        tmp_path / "bad_window_fraction.ref_kmers.zarr",
        blacklisted_fraction=np.array([0.25, 1.25], dtype=np.float64),
    )
    bad_chromosome_index = _write_dense_window_store(
        tmp_path / "bad_chromosome_index.ref_kmers.zarr",
        row_chromosome=np.array([0, 2], dtype=np.int32),
    )
    bad_eligible_windows = _write_sparse_grouped_store(
        tmp_path / "bad_eligible_windows.ref_kmers.zarr",
        eligible_windows=np.array([1, -1, 0], dtype=np.int32),
    )
    bad_group_fraction = _write_sparse_grouped_store(
        tmp_path / "bad_group_fraction.ref_kmers.zarr",
        blacklisted_fraction=np.array([0.0, np.nan, 0.0], dtype=np.float64),
    )

    with pytest.raises(
        ValueError, match="row_start_bp must be smaller than row_end_bp"
    ):
        cfdnalab.read_ref_kmers(bad_interval)
    with pytest.raises(ValueError, match="blacklisted_fraction"):
        cfdnalab.read_ref_kmers(bad_window_fraction)
    with pytest.raises(ValueError, match="row_chromosome contains an index outside"):
        cfdnalab.read_ref_kmers(bad_chromosome_index)
    with pytest.raises(ValueError, match="eligible_windows"):
        cfdnalab.read_ref_kmers(bad_eligible_windows)
    with pytest.raises(ValueError, match="blacklisted_fraction"):
        cfdnalab.read_ref_kmers(bad_group_fraction)


def test_ref_kmer_loader_rejects_invalid_json_labels(tmp_path: Path) -> None:
    numeric_group_labels = _write_dense_global_motif_group_store(
        tmp_path / "numeric_group_labels.ref_kmers.zarr"
    )
    _patch_array_metadata(
        numeric_group_labels,
        "motif_index",
        attributes={"label_field": "motif_group", "labels": [1, 2]},
    )
    control_character_label = _write_sparse_grouped_store(
        tmp_path / "control_character_label.ref_kmers.zarr",
        group_names=["A", "bad\nlabel", "empty"],
    )

    with pytest.raises(ValueError, match="labels must be character strings"):
        cfdnalab.read_ref_kmers(numeric_group_labels)
    with pytest.raises(ValueError, match="labels must not contain control characters"):
        cfdnalab.read_ref_kmers(control_character_label)


def test_sparse_ref_kmers_allow_empty_stored_coordinates(tmp_path: Path) -> None:
    store_path = _write_sparse_grouped_store(
        tmp_path / "empty_sparse.ref_kmers.zarr",
        sparse_row=np.asarray([], dtype=np.int32),
        sparse_motif=np.asarray([], dtype=np.int32),
        sparse_frequency=np.asarray([], dtype=np.float64),
    )

    ref_kmers = cfdnalab.read_ref_kmers(store_path)

    stored = ref_kmers.data_frame()
    assert list(stored.columns) == [
        "group_idx",
        "group_name",
        "eligible_windows",
        "blacklisted_fraction",
        "motif_index",
        "motif",
        "frequency",
        "count",
    ]
    assert len(stored) == 0
    coo = ref_kmers.sparse_frequencies_matrix()
    assert coo.shape == (3, 3)
    assert coo.nnz == 0
    np.testing.assert_allclose(
        ref_kmers.dense_counts_array(allow_densify=True),
        np.zeros((3, 3), dtype=np.float64),
    )


def _write_dense_window_store(
    path: Path,
    *,
    schema: str = "ref_kmer_frequencies",
    schema_version: int = 1,
    value_units: str = "reference_kmer_frequency",
    count_reconstruction: str = (
        "reference_kmer_count = frequency * row_scaling_factor[row]"
    ),
    row_scaling_factor_array: str = "row_scaling_factor",
    motif_names: np.ndarray = MOTIF_NAMES,
    canonical: bool = False,
    frequencies: np.ndarray | None = None,
    row_scaling_factor: np.ndarray | None = None,
    row_chromosome: np.ndarray | None = None,
    row_start_bp: np.ndarray | None = None,
    row_end_bp: np.ndarray | None = None,
    blacklisted_fraction: np.ndarray | None = None,
    omit: set[str] | None = None,
    frequencies_dimension_names: tuple[str, str] = ("row", "motif"),
) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    frequency_values = (
        frequencies
        if frequencies is not None
        else _default_dense_frequencies(len(motif_names))
    )
    _set_ref_kmer_root_attrs(
        root,
        schema=schema,
        schema_version=schema_version,
        storage_mode="dense",
        row_mode="bed",
        motif_axis_kind="motif",
        value_units=value_units,
        count_reconstruction=count_reconstruction,
        row_scaling_factor_array=row_scaling_factor_array,
        kmer_size=len(motif_names[0]),
        canonical=canonical,
    )

    _create_motif_axis(root, motif_names)
    _create_array(
        root,
        "row",
        np.array([0, 1], dtype=np.int32),
        chunks=(2,),
        dimension_names=("row",),
    )
    _create_labeled_axis(
        root,
        "chromosome",
        np.array([0, 1], dtype=np.int32),
        "chromosome_name",
        np.array(["chr2", "chr10"], dtype=object),
        dimension_names=("chromosome",),
    )
    _create_array(
        root,
        "row_chromosome",
        (
            row_chromosome
            if row_chromosome is not None
            else np.array([0, 1], dtype=np.int32)
        ),
        chunks=(2,),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "row_start_bp",
        (
            row_start_bp
            if row_start_bp is not None
            else np.array([10, 40], dtype=np.int64)
        ),
        chunks=(2,),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "row_end_bp",
        (
            row_end_bp
            if row_end_bp is not None
            else np.array([20, 60], dtype=np.int64)
        ),
        chunks=(2,),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "blacklisted_fraction",
        (
            blacklisted_fraction
            if blacklisted_fraction is not None
            else np.array([0.25, 0.0], dtype=np.float64)
        ),
        chunks=(2,),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "row_scaling_factor",
        (
            row_scaling_factor
            if row_scaling_factor is not None
            else np.array([4.0, 2.0], dtype=np.float64)
        ),
        chunks=(2,),
        dimension_names=("row",),
    )
    _create_reference_footprint(root)
    _create_array(
        root,
        "frequencies",
        frequency_values,
        chunks=frequency_values.shape,
        dimension_names=frequencies_dimension_names,
    )

    for array_name in omit or set():
        _delete_array(root, array_name)

    return path


def _default_dense_frequencies(n_motifs: int) -> np.ndarray:
    frequencies = np.zeros((2, n_motifs), dtype=np.float64)
    if n_motifs >= 1:
        frequencies[0, 0] = 0.25
        frequencies[1, 0] = 0.5
    if n_motifs >= 2:
        frequencies[1, 1] = 0.5
    if n_motifs >= 3:
        frequencies[0, 2] = 0.75
    return frequencies


def _write_sparse_grouped_store(
    path: Path,
    *,
    sparse_row: np.ndarray = np.array([0, 0, 1], dtype=np.int32),
    sparse_motif: np.ndarray = np.array([0, 2, 1], dtype=np.int32),
    sparse_shape: np.ndarray = np.array([3, 3], dtype=np.int32),
    sparse_frequency: np.ndarray = np.array([0.25, 0.75, 1.0], dtype=np.float64),
    sparse_dimension_labels: np.ndarray = np.array(["row", "motif"], dtype=object),
    group_names: list[str] | None = None,
    eligible_windows: np.ndarray = np.array([1, 2, 0], dtype=np.int32),
    blacklisted_fraction: np.ndarray = np.array([0.0, 0.125, 0.0], dtype=np.float64),
) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    _set_ref_kmer_root_attrs(
        root,
        storage_mode="sparse_coo",
        row_mode="grouped_bed",
        motif_axis_kind="motif",
        kmer_size=2,
    )

    _create_motif_axis(root, MOTIF_NAMES)
    _create_array(
        root,
        "row",
        np.array([0, 1, 2], dtype=np.int32),
        chunks=(3,),
        dimension_names=("row",),
    )
    _create_labeled_axis(
        root,
        "group",
        np.array([0, 1, 2], dtype=np.int32),
        "group_name",
        group_names if group_names is not None else ["A", "long_group", "empty"],
        dimension_names=("row",),
    )
    _create_array(
        root,
        "eligible_windows",
        eligible_windows,
        chunks=(3,),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "blacklisted_fraction",
        blacklisted_fraction,
        chunks=(3,),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "row_scaling_factor",
        np.array([4.0, 2.0, 0.0], dtype=np.float64),
        chunks=(3,),
        dimension_names=("row",),
    )
    _create_reference_footprint(root)

    sparse_group = root.create_group("sparse")
    _create_array(
        sparse_group,
        "row",
        sparse_row,
        chunks=(max(len(sparse_row), 1),),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse_group,
        "motif",
        sparse_motif,
        chunks=(max(len(sparse_motif), 1),),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse_group,
        "frequency",
        sparse_frequency,
        chunks=(max(len(sparse_frequency), 1),),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse_group,
        "shape",
        sparse_shape,
        chunks=(2,),
        dimension_names=("sparse_dimension",),
    )
    _create_labeled_axis(
        sparse_group,
        "sparse_dimension",
        np.array([0, 1], dtype=np.int32),
        "sparse_dimension_name",
        sparse_dimension_labels,
        dimension_names=("sparse_dimension",),
    )

    return path


def _write_dense_global_motif_group_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    _set_ref_kmer_root_attrs(
        root,
        storage_mode="dense",
        row_mode="global",
        motif_axis_kind="motif_group",
        kmer_size=2,
    )

    _create_motif_group_axis(root, MOTIF_GROUP_NAMES)
    _create_labeled_axis(
        root,
        "row",
        np.array([0], dtype=np.int32),
        "row_label",
        np.array(["global"], dtype=object),
        dimension_names=("row",),
    )
    _create_array(
        root,
        "row_scaling_factor",
        np.array([4.0], dtype=np.float64),
        chunks=(1,),
        dimension_names=("row",),
    )
    _create_reference_footprint(root)
    _create_array(
        root,
        "frequencies",
        np.array([[0.25, 0.75]], dtype=np.float64),
        chunks=(1, 2),
        dimension_names=("row", "motif"),
    )

    return path


def _set_ref_kmer_root_attrs(
    root: zarr.Group,
    *,
    schema: str = "ref_kmer_frequencies",
    schema_version: int = 1,
    storage_mode: str,
    row_mode: str,
    motif_axis_kind: str,
    value_units: str = "reference_kmer_frequency",
    count_units: str = "reference_kmer_count",
    row_scaling_factor_array: str = "row_scaling_factor",
    count_reconstruction: str = (
        "reference_kmer_count = frequency * row_scaling_factor[row]"
    ),
    kmer_size: int,
    canonical: bool = False,
    all_motifs: bool = False,
    assign_by: str = "count-overlap",
) -> None:
    root.attrs["cfdnalab_schema"] = schema
    root.attrs["cfdnalab_schema_version"] = schema_version
    root.attrs["storage_mode"] = storage_mode
    root.attrs["row_mode"] = row_mode
    root.attrs["motif_axis_kind"] = motif_axis_kind
    root.attrs["value_units"] = value_units
    root.attrs["count_units"] = count_units
    root.attrs["row_scaling_factor_array"] = row_scaling_factor_array
    root.attrs["count_reconstruction"] = count_reconstruction
    root.attrs["kmer_size"] = kmer_size
    root.attrs["canonical"] = canonical
    root.attrs["all_motifs"] = all_motifs
    root.attrs["assign_by"] = assign_by
    if storage_mode == "dense":
        root.attrs["primary_array"] = "frequencies"
        root.attrs["primary_group"] = None
    else:
        root.attrs["primary_array"] = None
        root.attrs["primary_group"] = "sparse"
        root.attrs["sparse_format"] = "coo"
        root.attrs["sparse_indices_base"] = 0


def _create_labeled_axis(
    root: zarr.Group,
    name: str,
    values: np.ndarray,
    label_field: str,
    labels: list[str] | np.ndarray | None,
    dimension_names: tuple[str, ...] | None = None,
) -> zarr.Array:
    axis = _create_array(
        root,
        name,
        values,
        chunks=(max(len(values), 1),),
        dimension_names=dimension_names,
    )
    axis.attrs["label_field"] = label_field
    if labels is not None:
        axis.attrs["labels"] = (
            labels.tolist() if isinstance(labels, np.ndarray) else labels
        )
    return axis


def _create_motif_axis(root: zarr.Group, labels: np.ndarray) -> None:
    _create_array(
        root,
        "motif_index",
        np.arange(len(labels), dtype=np.int32),
        chunks=(max(len(labels), 1),),
        dimension_names=("motif",),
    )
    motif_width = len(labels[0]) if len(labels) else 0
    _create_array(
        root,
        "motif_byte",
        np.arange(motif_width, dtype=np.int32),
        chunks=(max(motif_width, 1),),
        dimension_names=("motif_byte",),
    )
    motif_ascii = np.frombuffer("".join(labels.tolist()).encode("ascii"), dtype=np.uint8)
    motif_ascii = motif_ascii.reshape((len(labels), motif_width))
    _create_array(
        root,
        "motif_ascii",
        motif_ascii,
        chunks=(max(len(labels), 1), max(motif_width, 1)),
        dimension_names=("motif", "motif_byte"),
    )


def _create_motif_group_axis(root: zarr.Group, labels: np.ndarray) -> None:
    _create_labeled_axis(
        root,
        "motif_index",
        np.arange(len(labels), dtype=np.int32),
        "motif_group",
        labels,
        dimension_names=("motif",),
    )


def _create_reference_footprint(root: zarr.Group) -> None:
    encoded = json.dumps(REFERENCE_FOOTPRINT)
    _create_array(
        root,
        "reference_contig_footprint_json",
        np.frombuffer(encoded.encode("utf-8"), dtype=np.uint8),
        chunks=(len(encoded),),
        dimension_names=("json_byte",),
    )


def _create_array(
    root: zarr.Group,
    name: str,
    values: np.ndarray,
    *,
    chunks: tuple[int, ...],
    dimension_names: tuple[str, ...] | None = None,
) -> zarr.Array:
    return root.create_array(
        name,
        data=values,
        chunks=chunks,
        dimension_names=dimension_names,
    )


def _patch_array_metadata(
    store_path: Path,
    array_name: str,
    *,
    dimension_names: tuple[str, ...] | None = None,
    attributes: dict[str, object] | None = None,
) -> None:
    metadata_path = store_path.joinpath(*array_name.split("/"), "zarr.json")
    metadata = json.loads(metadata_path.read_text())
    if dimension_names is not None:
        metadata["dimension_names"] = list(dimension_names)
    if attributes is not None:
        metadata["attributes"] = attributes
    metadata_path.write_text(json.dumps(metadata))


def _delete_array(root: zarr.Group, name: str) -> None:
    if "/" not in name:
        del root[name]
        return
    group_name, array_name = name.split("/", 1)
    del root[group_name][array_name]
