from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sparse
import zarr

import cfdnalab


MOTIF_INDEX = np.array([0, 1, 2], dtype=np.int32)
MOTIF_NAMES = np.array(["_AA", "_CC", "_GG"], dtype=object)
MOTIF_GROUP_INDEX = np.array([0, 1], dtype=np.int32)
MOTIF_GROUP_NAMES = np.array(["short", "group-two"], dtype=object)


def test_dense_windowed_end_motifs_load_metadata_and_arrays(tmp_path: Path) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert ends.path == store_path
    assert isinstance(ends, cfdnalab.WindowedEndMotifCounts)
    assert repr(ends) == (
        "WindowedEndMotifCounts("
        f"path={str(store_path)!r}, "
        "schema_version=1, "
        "storage_mode='dense', "
        "row_mode='bed', "
        "shape=(2, 3)"
        ")"
    )
    assert not hasattr(ends, "group_metadata")
    assert ends.storage_mode() == "dense"
    assert ends.row_mode() == "bed"
    assert ends.has_motif("_AA")
    assert not ends.has_motif("_TT")
    assert ends.motif_idx("_GG") == 2
    pd.testing.assert_frame_equal(
        ends.motifs_metadata(),
        pd.DataFrame({"motif_index": MOTIF_INDEX, "motif": MOTIF_NAMES}),
    )
    pd.testing.assert_frame_equal(
        ends.window_metadata(),
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
    np.testing.assert_array_equal(
        ends.dense_counts_array(),
        np.array([[1.0, 0.0, 2.5], [0.5, 4.0, 0.0]], dtype=np.float64),
    )
    assert ends.dense_counts_zarr_array().shape == (2, 3)
    assert sparse.isspmatrix_coo(ends.sparse_counts_matrix())


def test_dense_windowed_end_motif_slice_helpers_return_expected_frames(
    tmp_path: Path,
) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    window_frame = ends.data_frame(window_idxs=0)
    motif_frame = ends.data_frame(motifs="_CC")

    pd.testing.assert_frame_equal(
        window_frame,
        pd.DataFrame(
            {
                "window_idx": np.array([0, 0, 0], dtype=np.int32),
                "chrom": np.array(["chr2", "chr2", "chr2"], dtype=object),
                "start": np.array([10, 10, 10], dtype=np.int64),
                "end": np.array([20, 20, 20], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.25, 0.25]),
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([1.0, 0.0, 2.5], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        motif_frame,
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chrom": np.array(["chr2", "chr10"], dtype=object),
                "start": np.array([10, 40], dtype=np.int64),
                "end": np.array([20, 60], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.0], dtype=np.float64),
                "motif_index": np.array([1, 1], dtype=np.int32),
                "motif": ["_CC", "_CC"],
                "count": np.array([0.0, 4.0], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(window_idxs=1),
        np.array([[0.5, 4.0, 0.0]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(motif_idxs=2),
        np.array([[2.5], [0.0]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(motifs="_CC"),
        np.array([[0.0], [4.0]], dtype=np.float64),
    )
    filtered_motif_frame = ends.data_frame(motifs="_CC", max_blacklisted_fraction=0.1)
    pd.testing.assert_frame_equal(
        filtered_motif_frame,
        motif_frame.iloc[[1]].reset_index(drop=True),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(motif_idxs=1),
        motif_frame,
    )
    assert ends.data_frame(window_idxs=0, max_blacklisted_fraction=0.1).empty
    with pytest.raises(
        ValueError, match="max_blacklisted_fraction must be a single finite fraction"
    ):
        ends.data_frame(window_idxs=0, max_blacklisted_fraction=1.1)
    assert ends.sparse_counts_matrix(window_idxs=1).shape == (1, 3)


def test_sparse_grouped_end_motifs_reconstruct_dense_matrix_and_metadata(
    tmp_path: Path,
) -> None:
    store_path = _write_sparse_grouped_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GroupedEndMotifCounts)
    assert not hasattr(ends, "window_metadata")
    assert ends.storage_mode() == "sparse_coo"
    assert ends.row_mode() == "grouped_bed"
    pd.testing.assert_frame_equal(
        ends.group_metadata(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 1, 2], dtype=np.int32),
                "group_name": np.array(["A", "long_group", "mid"], dtype=object),
                "eligible_windows": np.array([1, 2, 0], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.125, 0.0], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(allow_densify=True),
        np.array(
            [[1.0, 0.0, 2.5], [0.0, 4.25, 0.0], [0.75, 0.0, 0.0]],
            dtype=np.float64,
        ),
    )
    with pytest.raises(ValueError, match="would densify a sparse end-motif store"):
        ends.dense_counts_array()
    np.testing.assert_array_equal(
        ends.dense_counts_array(groups="A", allow_densify=True),
        np.array([[1.0, 0.0, 2.5]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(motifs="_AA", allow_densify=True),
        np.array([[1.0], [0.0], [0.75]], dtype=np.float64),
    )
    motif_coo = ends.sparse_counts_matrix(motifs="_AA")
    assert motif_coo.shape == (3, 1)
    np.testing.assert_array_equal(motif_coo.row, np.array([0, 2], dtype=np.int32))
    np.testing.assert_array_equal(motif_coo.col, np.array([0, 0], dtype=np.int32))
    np.testing.assert_array_equal(motif_coo.data, np.array([1.0, 0.75]))
    motif_idx_coo = ends.sparse_counts_matrix(motif_idxs=1)
    assert motif_idx_coo.shape == (3, 1)
    np.testing.assert_array_equal(motif_idx_coo.row, np.array([1], dtype=np.int32))
    np.testing.assert_array_equal(motif_idx_coo.data, np.array([4.25]))
    assert ends.group_idx("long_group") == 1
    group_coo = ends.sparse_counts_matrix(groups="long_group")
    assert group_coo.shape == (1, 3)
    np.testing.assert_array_equal(group_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(group_coo.col, np.array([1], dtype=np.int32))
    np.testing.assert_array_equal(group_coo.data, np.array([4.25]))


def test_sparse_rows_preserve_selected_order_and_metadata(tmp_path: Path) -> None:
    store_path = _write_sparse_grouped_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    pd.testing.assert_frame_equal(
        ends.data_frame(),
        pd.DataFrame(
            {
                "group_idx": np.array([0, 0, 1, 2], dtype=np.int32),
                "group_name": np.array(["A", "A", "long_group", "mid"], dtype=object),
                "eligible_windows": np.array([1, 1, 2, 0], dtype=np.int32),
                "blacklisted_fraction": np.array([0.0, 0.0, 0.125, 0.0]),
                "motif_index": np.array([0, 2, 1, 0], dtype=np.int32),
                "motif": np.array(["_AA", "_GG", "_CC", "_AA"], dtype=object),
                "count": np.array([1.0, 2.5, 4.25, 0.75], dtype=np.float64),
            }
        ),
    )


def test_sparse_windowed_end_motifs_slice_window_without_dense_roundtrip(
    tmp_path: Path,
) -> None:
    store_path = _write_sparse_window_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.WindowedEndMotifCounts)
    assert ends.storage_mode() == "sparse_coo"
    assert ends.has_motif("_CC")
    assert not ends.has_motif("_TT")
    with pytest.raises(KeyError, match="Unknown end-motif label"):
        ends.dense_counts_array(motifs="_TT")
    window_coo = ends.sparse_counts_matrix(window_idxs=1)

    assert window_coo.shape == (1, 3)
    np.testing.assert_array_equal(window_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(window_coo.col, np.array([1], dtype=np.int32))
    np.testing.assert_array_equal(window_coo.data, np.array([4.0]))
    with pytest.raises(ValueError, match="dense_counts_zarr_array"):
        ends.dense_counts_zarr_array()
    pd.testing.assert_frame_equal(
        ends.data_frame(window_idxs=0),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 0], dtype=np.int32),
                "chrom": np.array(["chr2", "chr2"], dtype=object),
                "start": np.array([10, 10], dtype=np.int64),
                "end": np.array([20, 20], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.25]),
                "motif_index": np.array([0, 2], dtype=np.int32),
                "motif": np.array(["_AA", "_GG"], dtype=object),
                "count": np.array([1.0, 2.5], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(window_idxs=0, densify=True),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 0, 0], dtype=np.int32),
                "chrom": np.array(["chr2", "chr2", "chr2"], dtype=object),
                "start": np.array([10, 10, 10], dtype=np.int64),
                "end": np.array([20, 20, 20], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.25, 0.25]),
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([1.0, 0.0, 2.5], dtype=np.float64),
            }
        ),
    )


def test_size_mode_end_motifs_use_windowed_helpers(tmp_path: Path) -> None:
    store_path = _write_sparse_window_store(
        tmp_path / "sample.end_motifs.zarr",
        row_mode="size",
    )
    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.WindowedEndMotifCounts)
    assert ends.row_mode() == "size"
    assert ends.window_metadata().shape == (2, 5)
    np.testing.assert_array_equal(
        ends.dense_counts_array(window_idxs=0, allow_densify=True),
        np.array([[1.0, 0.0, 2.5]], dtype=np.float64),
    )


def test_end_motif_data_frame_rejects_invalid_selector_combinations(
    tmp_path: Path,
) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    with pytest.raises(TypeError, match="unexpected keyword argument 'groups'"):
        ends.data_frame(groups="A")
    with pytest.raises(ValueError, match="Use either motifs or motif_idxs"):
        ends.data_frame(motifs="_AA", motif_idxs=0)
    with pytest.raises(ValueError, match="window_idxs contains duplicate values"):
        ends.data_frame(window_idxs=[0, 0])
    with pytest.raises(ValueError, match="motifs contains duplicate values"):
        ends.data_frame(motifs=["_AA", "_AA"])


def test_global_end_motifs_read_row_label_from_json_attrs(tmp_path: Path) -> None:
    store_path = _write_sparse_global_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GlobalEndMotifCounts)
    assert not hasattr(ends, "window_metadata")
    assert not hasattr(ends, "group_metadata")
    assert ends.row_mode() == "global"
    pd.testing.assert_frame_equal(
        ends.data_frame(motifs="_AA", densify=True),
        pd.DataFrame(
            {
                "row_label": np.array(["global"], dtype=object),
                "motif_index": np.array([0], dtype=np.int32),
                "motif": ["_AA"],
                "count": np.array([1.25], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(allow_densify=True),
        np.array([[1.25, 0.0, 3.5]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(),
        pd.DataFrame(
            {
                "row_label": np.array(["global", "global"], dtype=object),
                "motif_index": np.array([0, 2], dtype=np.int32),
                "motif": np.array(["_AA", "_GG"], dtype=object),
                "count": np.array([1.25, 3.5], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(motif_idxs=1, allow_densify=True),
        np.array([[0.0]], dtype=np.float64),
    )
    assert ends.sparse_counts_matrix(motif_idxs=1).nnz == 0


def test_dense_global_end_motifs_load_without_sparse_arrays(tmp_path: Path) -> None:
    store_path = _write_dense_global_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GlobalEndMotifCounts)
    assert ends.storage_mode() == "dense"
    assert ends.row_mode() == "global"
    assert ends.dense_counts_zarr_array().shape == (1, 3)
    np.testing.assert_array_equal(
        ends.dense_counts_array(),
        np.array([[1.0, 0.0, 2.5]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(),
        pd.DataFrame(
            {
                "row_label": np.array(["global", "global", "global"], dtype=object),
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([1.0, 0.0, 2.5], dtype=np.float64),
            }
        ),
    )


def test_dense_global_motif_group_axis_uses_json_labels(tmp_path: Path) -> None:
    store_path = _write_dense_global_motif_group_store(
        tmp_path / "sample.end_motifs.zarr"
    )

    ends = cfdnalab.read_end_motifs(store_path)

    pd.testing.assert_frame_equal(
        ends.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": MOTIF_GROUP_INDEX,
                "motif": MOTIF_GROUP_NAMES,
            }
        ),
    )
    assert ends.motif_idx("group-two") == 1
    assert ends.has_motif("short")
    np.testing.assert_array_equal(
        ends.dense_counts_array(motifs="group-two"),
        np.array([[3.0]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(motifs="group-two"),
        pd.DataFrame(
            {
                "row_label": np.array(["global"], dtype=object),
                "motif_index": np.array([1], dtype=np.int32),
                "motif": np.array(["group-two"], dtype=object),
                "count": np.array([3.0], dtype=np.float64),
            }
        ),
    )


def test_sparse_global_motif_group_axis_converts_with_motif_columns(
    tmp_path: Path,
) -> None:
    store_path = _write_sparse_global_motif_group_store(
        tmp_path / "sample.end_motifs.zarr"
    )

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GlobalEndMotifCounts)
    assert ends.storage_mode() == "sparse_coo"
    pd.testing.assert_frame_equal(
        ends.motifs_metadata(),
        pd.DataFrame(
            {
                "motif_index": MOTIF_GROUP_INDEX,
                "motif": MOTIF_GROUP_NAMES,
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.sparse_counts_matrix().toarray(),
        np.array([[1.5, 0.0]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.sparse_counts_matrix(motifs=["group-two", "short"]).toarray(),
        np.array([[0.0, 1.5]], dtype=np.float64),
    )
    with pytest.raises(ValueError, match="would densify a sparse end-motif store"):
        ends.dense_counts_array()
    np.testing.assert_array_equal(
        ends.dense_counts_array(allow_densify=True),
        np.array([[1.5, 0.0]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(motif_idxs=[1, 0], allow_densify=True),
        np.array([[0.0, 1.5]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(motifs="group-two", allow_densify=True),
        np.array([[0.0]], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(),
        pd.DataFrame(
            {
                "row_label": np.array(["global"], dtype=object),
                "motif_index": np.array([0], dtype=np.int32),
                "motif": np.array(["short"], dtype=object),
                "count": np.array([1.5], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(densify=True),
        pd.DataFrame(
            {
                "row_label": np.array(["global", "global"], dtype=object),
                "motif_index": MOTIF_GROUP_INDEX,
                "motif": MOTIF_GROUP_NAMES,
                "count": np.array([1.5, 0.0], dtype=np.float64),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        ends.data_frame(motifs="group-two", densify=True),
        pd.DataFrame(
            {
                "row_label": np.array(["global"], dtype=object),
                "motif_index": np.array([1], dtype=np.int32),
                "motif": np.array(["group-two"], dtype=object),
                "count": np.array([0.0], dtype=np.float64),
            }
        ),
    )


def test_dense_grouped_end_motifs_use_group_helpers(tmp_path: Path) -> None:
    store_path = _write_dense_grouped_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GroupedEndMotifCounts)
    assert ends.storage_mode() == "dense"
    assert ends.row_mode() == "grouped_bed"
    assert ends.group_idx("long_group") == 1
    np.testing.assert_array_equal(
        ends.dense_counts_array(group_idxs=0),
        np.array([[1.0, 0.0, 2.5]], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_array(groups="long_group"),
        np.array([[0.0, 4.25, 0.0]], dtype=np.float64),
    )
    group_coo = ends.sparse_counts_matrix(groups="long_group")
    assert group_coo.shape == (1, 3)
    np.testing.assert_array_equal(group_coo.col, np.array([1], dtype=np.int32))
    pd.testing.assert_frame_equal(
        ends.data_frame(groups="long_group"),
        pd.DataFrame(
            {
                "group_idx": np.array([1, 1, 1], dtype=np.int32),
                "group_name": np.array(
                    ["long_group", "long_group", "long_group"],
                    dtype=object,
                ),
                "eligible_windows": np.array([2, 2, 2], dtype=np.int32),
                "blacklisted_fraction": np.array([0.125, 0.125, 0.125]),
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([0.0, 4.25, 0.0], dtype=np.float64),
            }
        ),
    )
    assert ends.data_frame(groups="long_group", max_blacklisted_fraction=0.1).empty
    filtered_motif_frame = ends.data_frame(motifs="_CC", max_blacklisted_fraction=0.1)
    assert filtered_motif_frame["group_name"].tolist() == ["A", "mid"]
    assert filtered_motif_frame["count"].tolist() == [0.0, 0.0]


def test_end_motif_loader_rejects_invalid_paths(tmp_path: Path) -> None:
    missing_path = tmp_path / "missing.end_motifs.zarr"
    file_path = tmp_path / "not_a_directory.zarr"
    file_path.write_text("not a zarr store")
    wrong_suffix_path = tmp_path / "sample.not_zarr"
    wrong_suffix_path.mkdir()

    with pytest.raises(FileNotFoundError, match="does not exist"):
        cfdnalab.read_end_motifs(missing_path)
    with pytest.raises(NotADirectoryError, match="not a directory"):
        cfdnalab.read_end_motifs(file_path)
    with pytest.raises(ValueError, match="should end with '.zarr'"):
        cfdnalab.read_end_motifs(wrong_suffix_path)


def test_end_motif_loader_rejects_schema_and_shape_problems(tmp_path: Path) -> None:
    wrong_schema = _write_dense_window_store(
        tmp_path / "wrong_schema.end_motifs.zarr",
        schema="other",
    )
    wrong_version = _write_dense_window_store(
        tmp_path / "wrong_version.end_motifs.zarr",
        schema_version=99,
    )
    missing_array = _write_dense_window_store(
        tmp_path / "missing_array.end_motifs.zarr",
        omit={"row_end_bp"},
    )
    wrong_dimensions = _write_dense_window_store(
        tmp_path / "wrong_dimensions.end_motifs.zarr",
        counts_dimension_names=("motif", "row"),
    )
    shape_mismatch = _write_sparse_grouped_store(
        tmp_path / "shape_mismatch.end_motifs.zarr",
        sparse_shape=np.array([3, 2], dtype=np.int32),
    )
    wrong_sparse_dimensions = _write_sparse_grouped_store(
        tmp_path / "wrong_sparse_dimensions.end_motifs.zarr",
        sparse_row_dimension_names=("row",),
    )
    wrong_sparse_dimension_labels = _write_sparse_grouped_store(
        tmp_path / "wrong_sparse_dimension_labels.end_motifs.zarr",
        sparse_dimension_labels=np.array(["motif", "row"], dtype=object),
    )
    empty_sparse_counts = _write_empty_sparse_global_store(
        tmp_path / "empty_sparse_counts.end_motifs.zarr",
    )
    missing_motif_ascii = _write_dense_window_store(
        tmp_path / "missing_motif_ascii.end_motifs.zarr",
        omit={"motif_ascii"},
    )

    with pytest.raises(ValueError, match="Expected cfdnalab_schema"):
        cfdnalab.read_end_motifs(wrong_schema)
    with pytest.raises(ValueError, match="Unsupported end-motif schema version"):
        cfdnalab.read_end_motifs(wrong_version)
    with pytest.raises(ValueError, match="missing arrays: \\['row_end_bp'\\]"):
        cfdnalab.read_end_motifs(missing_array)
    with pytest.raises(ValueError, match="dense counts dimensions must be"):
        cfdnalab.read_end_motifs(wrong_dimensions)
    with pytest.raises(ValueError, match="sparse/shape does not match"):
        cfdnalab.read_end_motifs(shape_mismatch)
    with pytest.raises(ValueError, match="sparse/row dimensions must be"):
        cfdnalab.read_end_motifs(wrong_sparse_dimensions)
    with pytest.raises(ValueError, match="sparse_dimension labels must be"):
        cfdnalab.read_end_motifs(wrong_sparse_dimension_labels)
    with pytest.raises(ValueError, match="No end-motif counts are available"):
        cfdnalab.read_end_motifs(empty_sparse_counts)
    with pytest.raises(ValueError, match="missing arrays: \\['motif_ascii'\\]"):
        cfdnalab.read_end_motifs(missing_motif_ascii)


def _write_dense_window_store(
    path: Path,
    *,
    schema: str = "end_motif_counts",
    schema_version: int = 1,
    omit: set[str] | None = None,
    counts_dimension_names: tuple[str, str] = ("row", "motif"),
) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = schema
    root.attrs["cfdnalab_schema_version"] = schema_version
    root.attrs["storage_mode"] = "dense"
    root.attrs["row_mode"] = "bed"

    _create_motif_axis(root, MOTIF_NAMES)
    _create_array(root, "row", np.array([0, 1], dtype=np.int32), chunks=(2,))
    _create_labeled_axis(
        root,
        "chromosome",
        np.array([0, 1], dtype=np.int32),
        "chromosome_name",
        np.array(["chr2", "chr10"], dtype=object),
    )
    _create_array(
        root,
        "row_chromosome",
        np.array([0, 1], dtype=np.int32),
        chunks=(2,),
    )
    _create_array(root, "row_start_bp", np.array([10, 40], dtype=np.int64), chunks=(2,))
    _create_array(root, "row_end_bp", np.array([20, 60], dtype=np.int64), chunks=(2,))
    _create_array(
        root,
        "blacklisted_fraction",
        np.array([0.25, 0.0], dtype=np.float64),
        chunks=(2,),
    )
    _create_array(
        root,
        "counts",
        np.array([[1.0, 0.0, 2.5], [0.5, 4.0, 0.0]], dtype=np.float64),
        chunks=(2, 3),
        dimension_names=counts_dimension_names,
    )

    for array_name in omit or set():
        del root[array_name]

    return path


def _write_reference_correction_ref_kmer_store(
    path: Path,
    *,
    motifs: np.ndarray | None = None,
    frequencies: np.ndarray | None = None,
    blacklisted_fraction: np.ndarray | None = None,
    storage_mode: str = "dense",
    row_mode: str = "bed",
) -> Path:
    if motifs is None:
        motifs = np.array(["AA", "CC", "GG"], dtype=object)
    motifs = np.asarray(motifs, dtype=object)
    if frequencies is None:
        if len(motifs) == 3:
            frequencies = np.array(
                [
                    [1.0 / 3.0, 1.0 / 6.0, 1.0 / 2.0],
                    [1.0 / 2.0, 1.0 / 4.0, 1.0 / 4.0],
                ],
                dtype=np.float64,
            )
        else:
            frequencies = np.full((2, len(motifs)), 1.0 / len(motifs))
    frequencies = np.asarray(frequencies, dtype=np.float64)
    number_of_rows = frequencies.shape[0]
    if blacklisted_fraction is None:
        blacklisted_fraction = np.zeros(number_of_rows, dtype=np.float64)
        if number_of_rows == 2:
            blacklisted_fraction = np.array([0.25, 0.0], dtype=np.float64)
    blacklisted_fraction = np.asarray(blacklisted_fraction, dtype=np.float64)
    if len(blacklisted_fraction) != number_of_rows:
        raise ValueError("blacklisted_fraction length must match frequencies rows")

    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "ref_kmer_frequencies"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = storage_mode
    root.attrs["row_mode"] = row_mode
    root.attrs["motif_axis_kind"] = "motif"
    root.attrs["value_units"] = "reference_kmer_frequency"
    root.attrs["count_units"] = "reference_kmer_count"
    root.attrs["row_scaling_factor_array"] = "row_scaling_factor"
    root.attrs["count_reconstruction"] = (
        "reference_kmer_count = frequency * row_scaling_factor[row]"
    )
    root.attrs["kmer_size"] = len(str(motifs[0]))
    root.attrs["canonical"] = False
    root.attrs["all_motifs"] = True
    root.attrs["assign_by"] = "count-overlap"
    if storage_mode == "dense":
        root.attrs["primary_array"] = "frequencies"
        root.attrs["primary_group"] = None
    else:
        root.attrs["primary_array"] = None
        root.attrs["primary_group"] = "sparse"
        root.attrs["sparse_format"] = "coo"
        root.attrs["sparse_indices_base"] = 0

    _create_motif_axis(root, motifs)
    if row_mode == "global":
        _create_labeled_axis(
            root,
            "row",
            np.array([0], dtype=np.int32),
            "row_label",
            np.array(["global"], dtype=object),
            dimension_names=("row",),
        )
    else:
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
            np.array([0, 1], dtype=np.int32),
            chunks=(2,),
            dimension_names=("row",),
        )
        _create_array(
            root,
            "row_start_bp",
            np.array([10, 40], dtype=np.int64),
            chunks=(2,),
            dimension_names=("row",),
        )
        _create_array(
            root,
            "row_end_bp",
            np.array([20, 60], dtype=np.int64),
            chunks=(2,),
            dimension_names=("row",),
        )
        _create_array(
            root,
            "blacklisted_fraction",
            blacklisted_fraction,
            chunks=(number_of_rows,),
            dimension_names=("row",),
        )
    _create_array(
        root,
        "row_scaling_factor",
        np.repeat(6.0, number_of_rows).astype(np.float64),
        chunks=(number_of_rows,),
        dimension_names=("row",),
    )
    encoded_footprint = json.dumps([{"name": "chr2", "length": 100}])
    _create_array(
        root,
        "reference_contig_footprint_json",
        np.frombuffer(encoded_footprint.encode("utf-8"), dtype=np.uint8),
        chunks=(len(encoded_footprint),),
        dimension_names=("json_byte",),
    )
    if storage_mode == "dense":
        _create_array(
            root,
            "frequencies",
            frequencies,
            chunks=frequencies.shape,
            dimension_names=("row", "motif"),
        )
    else:
        sparse_rows, sparse_motifs = np.nonzero(frequencies > 0.0)
        sparse_group = root.create_group("sparse")
        _create_array(
            sparse_group,
            "row",
            sparse_rows.astype(np.int32),
            chunks=(max(len(sparse_rows), 1),),
            dimension_names=("nnz",),
        )
        _create_array(
            sparse_group,
            "motif",
            sparse_motifs.astype(np.int32),
            chunks=(max(len(sparse_motifs), 1),),
            dimension_names=("nnz",),
        )
        _create_array(
            sparse_group,
            "frequency",
            frequencies[sparse_rows, sparse_motifs].astype(np.float64),
            chunks=(max(len(sparse_rows), 1),),
            dimension_names=("nnz",),
        )
        _create_array(
            sparse_group,
            "shape",
            np.array([number_of_rows, len(motifs)], dtype=np.int32),
            chunks=(2,),
            dimension_names=("sparse_dimension",),
        )
        _create_labeled_axis(
            sparse_group,
            "sparse_dimension",
            np.array([0, 1], dtype=np.int32),
            "sparse_dimension_name",
            np.array(["row", "motif"], dtype=object),
            dimension_names=("sparse_dimension",),
        )
    return path


def _write_dense_global_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = "dense"
    root.attrs["row_mode"] = "global"

    _create_motif_axis(root, MOTIF_NAMES)
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
        np.array([[1.0, 0.0, 2.5]], dtype=np.float64),
        chunks=(1, 3),
        dimension_names=("row", "motif"),
    )

    return path


def _write_dense_global_motif_group_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 2
    root.attrs["storage_mode"] = "dense"
    root.attrs["row_mode"] = "global"
    root.attrs["motif_axis_kind"] = "motif_group"

    _create_motif_group_axis(root, MOTIF_GROUP_NAMES)
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
        np.array([[1.5, 3.0]], dtype=np.float64),
        chunks=(1, 2),
        dimension_names=("row", "motif"),
    )

    return path


def _write_sparse_global_motif_group_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 2
    root.attrs["storage_mode"] = "sparse_coo"
    root.attrs["row_mode"] = "global"
    root.attrs["motif_axis_kind"] = "motif_group"

    _create_motif_group_axis(root, MOTIF_GROUP_NAMES)
    _create_labeled_axis(
        root,
        "row",
        np.array([0], dtype=np.int32),
        "row_label",
        np.array(["global"], dtype=object),
    )

    sparse_group = root.create_group("sparse")
    _create_array(
        sparse_group,
        "row",
        np.array([0], dtype=np.int32),
        chunks=(1,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse_group,
        "motif",
        np.array([0], dtype=np.int32),
        chunks=(1,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse_group,
        "count",
        np.array([1.5], dtype=np.float64),
        chunks=(1,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse_group,
        "shape",
        np.array([1, 2], dtype=np.int32),
        chunks=(2,),
        dimension_names=("sparse_dimension",),
    )
    _create_labeled_axis(
        sparse_group,
        "sparse_dimension",
        np.array([0, 1], dtype=np.int32),
        "sparse_dimension_name",
        np.array(["row", "motif"], dtype=object),
        dimension_names=("sparse_dimension",),
    )

    return path


def _write_sparse_window_store(path: Path, *, row_mode: str = "bed") -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = "sparse_coo"
    root.attrs["row_mode"] = row_mode

    _create_motif_axis(root, MOTIF_NAMES)
    _create_array(root, "row", np.array([0, 1], dtype=np.int32), chunks=(2,))
    _create_labeled_axis(
        root,
        "chromosome",
        np.array([0, 1], dtype=np.int32),
        "chromosome_name",
        np.array(["chr2", "chr10"], dtype=object),
    )
    _create_array(
        root,
        "row_chromosome",
        np.array([0, 1], dtype=np.int32),
        chunks=(2,),
    )
    _create_array(root, "row_start_bp", np.array([10, 40], dtype=np.int64), chunks=(2,))
    _create_array(root, "row_end_bp", np.array([20, 60], dtype=np.int64), chunks=(2,))
    _create_array(
        root,
        "blacklisted_fraction",
        np.array([0.25, 0.0], dtype=np.float64),
        chunks=(2,),
    )

    sparse = root.create_group("sparse")
    _create_array(
        sparse,
        "row",
        np.array([0, 0, 1], dtype=np.int32),
        chunks=(3,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "motif",
        np.array([0, 2, 1], dtype=np.int32),
        chunks=(3,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "count",
        np.array([1.0, 2.5, 4.0], dtype=np.float64),
        chunks=(3,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "shape",
        np.array([2, 3], dtype=np.int32),
        chunks=(2,),
        dimension_names=("sparse_dimension",),
    )
    _create_labeled_axis(
        sparse,
        "sparse_dimension",
        np.array([0, 1], dtype=np.int32),
        "sparse_dimension_name",
        np.array(["row", "motif"], dtype=object),
        dimension_names=("sparse_dimension",),
    )

    return path


def _write_dense_grouped_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = "dense"
    root.attrs["row_mode"] = "grouped_bed"

    _create_motif_axis(root, MOTIF_NAMES)
    _create_array(root, "row", np.array([0, 1, 2], dtype=np.int32), chunks=(3,))
    _create_labeled_axis(
        root,
        "group",
        np.array([0, 1, 2], dtype=np.int32),
        "group_name",
        np.array(["A", "long_group", "mid"], dtype=object),
    )
    _create_array(
        root,
        "eligible_windows",
        np.array([1, 2, 0], dtype=np.int32),
        chunks=(3,),
    )
    _create_array(
        root,
        "blacklisted_fraction",
        np.array([0.0, 0.125, 0.0], dtype=np.float64),
        chunks=(3,),
    )
    _create_array(
        root,
        "counts",
        np.array(
            [[1.0, 0.0, 2.5], [0.0, 4.25, 0.0], [0.75, 0.0, 0.0]],
            dtype=np.float64,
        ),
        chunks=(3, 3),
        dimension_names=("row", "motif"),
    )

    return path


def _write_sparse_grouped_store(
    path: Path,
    *,
    sparse_shape: np.ndarray = np.array([3, 3], dtype=np.int32),
    sparse_row_dimension_names: tuple[str, ...] = ("nnz",),
    sparse_dimension_labels: np.ndarray = np.array(["row", "motif"], dtype=object),
) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = "sparse_coo"
    root.attrs["row_mode"] = "grouped_bed"

    _create_motif_axis(root, MOTIF_NAMES)
    _create_array(root, "row", np.array([0, 1, 2], dtype=np.int32), chunks=(3,))
    _create_labeled_axis(
        root,
        "group",
        np.array([0, 1, 2], dtype=np.int32),
        "group_name",
        np.array(["A", "long_group", "mid"], dtype=object),
    )
    _create_array(
        root,
        "eligible_windows",
        np.array([1, 2, 0], dtype=np.int32),
        chunks=(3,),
    )
    _create_array(
        root,
        "blacklisted_fraction",
        np.array([0.0, 0.125, 0.0], dtype=np.float64),
        chunks=(3,),
    )

    sparse = root.create_group("sparse")
    _create_array(
        sparse,
        "row",
        np.array([0, 0, 1, 2], dtype=np.int32),
        chunks=(4,),
        dimension_names=sparse_row_dimension_names,
    )
    _create_array(
        sparse,
        "motif",
        np.array([0, 2, 1, 0], dtype=np.int32),
        chunks=(4,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "count",
        np.array([1.0, 2.5, 4.25, 0.75], dtype=np.float64),
        chunks=(4,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "shape",
        sparse_shape,
        chunks=(2,),
        dimension_names=("sparse_dimension",),
    )
    _create_labeled_axis(
        sparse,
        "sparse_dimension",
        np.array([0, 1], dtype=np.int32),
        "sparse_dimension_name",
        sparse_dimension_labels,
        dimension_names=("sparse_dimension",),
    )

    return path


def _write_sparse_global_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = "sparse_coo"
    root.attrs["row_mode"] = "global"

    _create_motif_axis(root, MOTIF_NAMES)
    _create_labeled_axis(
        root,
        "row",
        np.array([0], dtype=np.int32),
        "row_label",
        np.array(["global"], dtype=object),
    )

    sparse = root.create_group("sparse")
    _create_array(
        sparse,
        "row",
        np.array([0, 0], dtype=np.int32),
        chunks=(2,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "motif",
        np.array([0, 2], dtype=np.int32),
        chunks=(2,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "count",
        np.array([1.25, 3.5], dtype=np.float64),
        chunks=(2,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "shape",
        np.array([1, 3], dtype=np.int32),
        chunks=(2,),
        dimension_names=("sparse_dimension",),
    )
    _create_labeled_axis(
        sparse,
        "sparse_dimension",
        np.array([0, 1], dtype=np.int32),
        "sparse_dimension_name",
        np.array(["row", "motif"], dtype=object),
        dimension_names=("sparse_dimension",),
    )

    return path


def _write_empty_sparse_global_store(path: Path) -> Path:
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = "end_motif_counts"
    root.attrs["cfdnalab_schema_version"] = 1
    root.attrs["storage_mode"] = "sparse_coo"
    root.attrs["row_mode"] = "global"

    _create_motif_axis(root, MOTIF_NAMES)
    _create_labeled_axis(
        root,
        "row",
        np.array([0], dtype=np.int32),
        "row_label",
        np.array(["global"], dtype=object),
    )

    sparse = root.create_group("sparse")
    _create_array(
        sparse,
        "row",
        np.array([], dtype=np.int32),
        chunks=(1,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "motif",
        np.array([], dtype=np.int32),
        chunks=(1,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "count",
        np.array([], dtype=np.float64),
        chunks=(1,),
        dimension_names=("nnz",),
    )
    _create_array(
        sparse,
        "shape",
        np.array([1, 3], dtype=np.int32),
        chunks=(2,),
        dimension_names=("sparse_dimension",),
    )
    _create_labeled_axis(
        sparse,
        "sparse_dimension",
        np.array([0, 1], dtype=np.int32),
        "sparse_dimension_name",
        np.array(["row", "motif"], dtype=object),
        dimension_names=("sparse_dimension",),
    )

    return path


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
        chunks=(len(values),),
        dimension_names=dimension_names,
    )
    axis.attrs["label_field"] = label_field
    if labels is not None:
        if isinstance(labels, np.ndarray):
            axis.attrs["labels"] = labels.tolist()
        else:
            axis.attrs["labels"] = labels
    return axis


def _create_motif_axis(root: zarr.Group, labels: np.ndarray) -> None:
    _create_array(
        root,
        "motif_index",
        np.arange(len(labels), dtype=np.int32),
        chunks=(len(labels),),
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
        MOTIF_GROUP_INDEX,
        "motif_group",
        labels,
        dimension_names=("motif",),
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
