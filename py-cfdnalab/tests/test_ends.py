from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sparse
import zarr

import cfdnalab


MOTIF_INDEX = np.array([0, 1, 2], dtype=np.int32)
MOTIF_NAMES = np.array(["_AA", "_CC", "_GG"], dtype=object)


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
    assert not hasattr(ends, "groups")
    assert ends.storage_mode() == "dense"
    assert ends.row_mode() == "bed"
    assert ends.motifs() == ["_AA", "_CC", "_GG"]
    assert ends.has_motif("_AA")
    assert not ends.has_motif("_TT")
    assert ends.motif_idx("_GG") == 2
    pd.testing.assert_frame_equal(
        ends.motif_metadata(),
        pd.DataFrame({"motif_index": MOTIF_INDEX, "motif": MOTIF_NAMES}),
    )
    pd.testing.assert_frame_equal(
        ends.windows(),
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chromosome": np.array([0, 1], dtype=np.int32),
                "chromosome_name": np.array(["chr2", "chr10"], dtype=object),
                "window_start_bp": np.array([10, 40], dtype=np.uint64),
                "window_end_bp": np.array([20, 60], dtype=np.uint64),
                "blacklisted_fraction": np.array([0.25, 0.0], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_matrix(),
        np.array([[1.0, 0.0, 2.5], [0.5, 4.0, 0.0]], dtype=np.float64),
    )
    assert ends.dense_counts_zarr_array().shape == (2, 3)
    assert sparse.isspmatrix_coo(ends.sparse_coo())


def test_dense_windowed_end_motif_slice_helpers_return_expected_frames(
    tmp_path: Path,
) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    window_frame = ends.dense_data_frame_for_window(0)
    motif_frame = ends.dense_data_frame_for_motif("_CC")

    pd.testing.assert_frame_equal(
        window_frame,
        pd.DataFrame(
            {
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([1.0, 0.0, 2.5], dtype=np.float64),
                "window_idx": np.array([0, 0, 0], dtype=np.int64),
                "chromosome": np.array([0, 0, 0], dtype=np.int64),
                "chromosome_name": np.array(["chr2", "chr2", "chr2"], dtype=object),
                "window_start_bp": np.array([10, 10, 10], dtype=np.int64),
                "window_end_bp": np.array([20, 20, 20], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.25, 0.25]),
            }
        ),
    )
    pd.testing.assert_frame_equal(
        motif_frame,
        pd.DataFrame(
            {
                "window_idx": np.array([0, 1], dtype=np.int32),
                "chromosome": np.array([0, 1], dtype=np.int32),
                "chromosome_name": np.array(["chr2", "chr10"], dtype=object),
                "window_start_bp": np.array([10, 40], dtype=np.uint64),
                "window_end_bp": np.array([20, 60], dtype=np.uint64),
                "blacklisted_fraction": np.array([0.25, 0.0], dtype=np.float64),
                "motif_index": np.array([1, 1], dtype=np.int64),
                "motif": ["_CC", "_CC"],
                "count": np.array([0.0, 4.0], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_for_window(1),
        np.array([0.5, 4.0, 0.0], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_for_motif_idx(2), np.array([2.5, 0.0], dtype=np.float64)
    )
    np.testing.assert_array_equal(
        ends.dense_counts_for_motif("_CC"), np.array([0.0, 4.0], dtype=np.float64)
    )
    pd.testing.assert_frame_equal(
        ends.dense_data_frame_for_motif_idx(1),
        motif_frame,
    )
    assert ends.sparse_coo_for_window(1).shape == (1, 3)


def test_sparse_grouped_end_motifs_reconstruct_dense_matrix_and_metadata(
    tmp_path: Path,
) -> None:
    store_path = _write_sparse_grouped_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GroupedEndMotifCounts)
    assert not hasattr(ends, "windows")
    assert ends.storage_mode() == "sparse_coo"
    assert ends.row_mode() == "grouped_bed"
    pd.testing.assert_frame_equal(
        ends.groups(),
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
        ends.dense_counts_matrix(allow_densify=True),
        np.array(
            [[1.0, 0.0, 2.5], [0.0, 4.25, 0.0], [0.75, 0.0, 0.0]],
            dtype=np.float64,
        ),
    )
    with pytest.raises(ValueError, match="would densify a sparse end-motif store"):
        ends.dense_counts_matrix()
    np.testing.assert_array_equal(
        ends.dense_counts_for_group("A", allow_densify=True),
        np.array([1.0, 0.0, 2.5], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_for_motif("_AA", allow_densify=True),
        np.array([1.0, 0.0, 0.75], dtype=np.float64),
    )
    motif_coo = ends.sparse_coo_for_motif("_AA")
    assert motif_coo.shape == (3, 1)
    np.testing.assert_array_equal(motif_coo.row, np.array([0, 2], dtype=np.int32))
    np.testing.assert_array_equal(motif_coo.col, np.array([0, 0], dtype=np.int32))
    np.testing.assert_array_equal(motif_coo.data, np.array([1.0, 0.75]))
    motif_idx_coo = ends.sparse_coo_for_motif_idx(1)
    assert motif_idx_coo.shape == (3, 1)
    np.testing.assert_array_equal(motif_idx_coo.row, np.array([1], dtype=np.int32))
    np.testing.assert_array_equal(motif_idx_coo.data, np.array([4.25]))
    assert ends.group_idx("long_group") == 1
    group_coo = ends.sparse_coo_for_group("long_group")
    assert group_coo.shape == (1, 3)
    np.testing.assert_array_equal(group_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(group_coo.col, np.array([1], dtype=np.int32))
    np.testing.assert_array_equal(group_coo.data, np.array([4.25]))


def test_sparse_coo_data_frame_preserves_sorted_payload(tmp_path: Path) -> None:
    store_path = _write_sparse_grouped_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    pd.testing.assert_frame_equal(
        ends.sparse_coo_data_frame(),
        pd.DataFrame(
            {
                "row": np.array([0, 0, 1, 2], dtype=np.int32),
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
        ends.dense_counts_for_motif("_TT")
    window_coo = ends.sparse_coo_for_window(1)

    assert window_coo.shape == (1, 3)
    np.testing.assert_array_equal(window_coo.row, np.array([0], dtype=np.int32))
    np.testing.assert_array_equal(window_coo.col, np.array([1], dtype=np.int32))
    np.testing.assert_array_equal(window_coo.data, np.array([4.0]))
    with pytest.raises(ValueError, match="dense_counts_zarr_array"):
        ends.dense_counts_zarr_array()
    with pytest.raises(ValueError, match="would densify a sparse end-motif store"):
        ends.dense_data_frame_for_window(0)
    pd.testing.assert_frame_equal(
        ends.dense_data_frame_for_window(0, allow_densify=True),
        pd.DataFrame(
            {
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([1.0, 0.0, 2.5], dtype=np.float64),
                "window_idx": np.array([0, 0, 0], dtype=np.int64),
                "chromosome": np.array([0, 0, 0], dtype=np.int64),
                "chromosome_name": np.array(["chr2", "chr2", "chr2"], dtype=object),
                "window_start_bp": np.array([10, 10, 10], dtype=np.int64),
                "window_end_bp": np.array([20, 20, 20], dtype=np.int64),
                "blacklisted_fraction": np.array([0.25, 0.25, 0.25]),
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
    assert ends.windows().shape == (2, 6)
    np.testing.assert_array_equal(
        ends.dense_counts_for_window(0, allow_densify=True),
        np.array([1.0, 0.0, 2.5], dtype=np.float64),
    )


def test_sparse_coo_data_frame_rejects_dense_stores(tmp_path: Path) -> None:
    store_path = _write_dense_window_store(tmp_path / "sample.end_motifs.zarr")
    ends = cfdnalab.read_end_motifs(store_path)

    with pytest.raises(ValueError, match="only available for sparse_coo output"):
        ends.sparse_coo_data_frame()


def test_global_end_motifs_read_row_label_from_json_attrs(tmp_path: Path) -> None:
    store_path = _write_sparse_global_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GlobalEndMotifCounts)
    assert not hasattr(ends, "windows")
    assert not hasattr(ends, "groups")
    assert ends.row_mode() == "global"
    pd.testing.assert_frame_equal(
        ends.dense_data_frame_for_motif("_AA", allow_densify=True),
        pd.DataFrame(
            {
                "row_label": np.array(["global"], dtype=object),
                "motif_index": np.array([0], dtype=np.int64),
                "motif": ["_AA"],
                "count": np.array([1.25], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_vec(allow_densify=True),
        np.array([1.25, 0.0, 3.5], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.sparse_coo_data_frame(),
        pd.DataFrame(
            {
                "row": np.array([0, 0], dtype=np.int32),
                "motif_index": np.array([0, 2], dtype=np.int32),
                "motif": np.array(["_AA", "_GG"], dtype=object),
                "count": np.array([1.25, 3.5], dtype=np.float64),
            }
        ),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_for_motif_idx(1, allow_densify=True),
        np.array([0.0], dtype=np.float64),
    )
    assert ends.sparse_coo_for_motif_idx(1).nnz == 0


def test_dense_global_end_motifs_load_without_sparse_arrays(tmp_path: Path) -> None:
    store_path = _write_dense_global_store(tmp_path / "sample.end_motifs.zarr")

    ends = cfdnalab.read_end_motifs(store_path)

    assert isinstance(ends, cfdnalab.GlobalEndMotifCounts)
    assert ends.storage_mode() == "dense"
    assert ends.row_mode() == "global"
    assert ends.dense_counts_zarr_array().shape == (1, 3)
    np.testing.assert_array_equal(
        ends.dense_counts_vec(),
        np.array([1.0, 0.0, 2.5], dtype=np.float64),
    )
    pd.testing.assert_frame_equal(
        ends.dense_data_frame(),
        pd.DataFrame(
            {
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([1.0, 0.0, 2.5], dtype=np.float64),
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
        ends.dense_counts_for_group(0),
        np.array([1.0, 0.0, 2.5], dtype=np.float64),
    )
    np.testing.assert_array_equal(
        ends.dense_counts_for_group("long_group"),
        np.array([0.0, 4.25, 0.0], dtype=np.float64),
    )
    group_coo = ends.sparse_coo_for_group("long_group")
    assert group_coo.shape == (1, 3)
    np.testing.assert_array_equal(group_coo.col, np.array([1], dtype=np.int32))
    pd.testing.assert_frame_equal(
        ends.dense_data_frame_for_group("long_group"),
        pd.DataFrame(
            {
                "motif_index": MOTIF_INDEX,
                "motif": MOTIF_NAMES,
                "count": np.array([0.0, 4.25, 0.0], dtype=np.float64),
                "group_idx": np.array([1, 1, 1], dtype=np.int64),
                "group_name": np.array(
                    ["long_group", "long_group", "long_group"],
                    dtype=object,
                ),
                "eligible_windows": np.array([2, 2, 2], dtype=np.int64),
                "blacklisted_fraction": np.array([0.125, 0.125, 0.125]),
            }
        ),
    )


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
    _create_array(root, "row_start_bp", np.array([10, 40], dtype=np.uint64), chunks=(2,))
    _create_array(root, "row_end_bp", np.array([20, 60], dtype=np.uint64), chunks=(2,))
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
    _create_array(root, "row_start_bp", np.array([10, 40], dtype=np.uint64), chunks=(2,))
    _create_array(root, "row_end_bp", np.array([20, 60], dtype=np.uint64), chunks=(2,))
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
    _create_array(root, "motif_index", MOTIF_INDEX, chunks=(len(labels),))
    motif_width = len(labels[0]) if len(labels) else 0
    _create_array(
        root,
        "motif_byte",
        np.arange(motif_width, dtype=np.int32),
        chunks=(max(motif_width, 1),),
    )
    motif_ascii = np.frombuffer("".join(labels.tolist()).encode("ascii"), dtype=np.uint8)
    motif_ascii = motif_ascii.reshape((len(labels), motif_width))
    _create_array(
        root,
        "motif_ascii",
        motif_ascii,
        chunks=(max(len(labels), 1), max(motif_width, 1)),
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
