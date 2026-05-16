from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import zarr

import cfdnalab


GROUP_NAMES = np.array(["alpha", "beta long", "gamma-unicode-aa"], dtype=object)
GROUP_AXIS = np.array([0, 1, 2], dtype=np.int32)
ELIGIBLE_INTERVALS = np.array([4, 7, 11], dtype=np.uint32)
LENGTH_BIN = np.array([0, 1], dtype=np.int32)
LENGTH_START_BP = np.array([30, 60], dtype=np.uint32)
LENGTH_END_BP = np.array([60, 90], dtype=np.uint32)
POSITION = np.array([0, 1, 2, 3], dtype=np.int32)
POSITION_BIN_START_BP = np.array([-4, -2, 0, 2], dtype=np.int32)
POSITION_BIN_END_BP = np.array([-2, 0, 2, 4], dtype=np.int32)

# Count value = group * 100 + length_bin * 10 + position. This makes every
# axis visible in assertions and catches accidental transposes or wrong joins.
COUNTS = np.array(
    [
        [[0.0, 1.0, 2.0, 3.0], [10.0, 11.0, 12.0, 13.0]],
        [[100.0, 101.0, 102.0, 103.0], [110.0, 111.0, 112.0, 113.0]],
        [[200.0, 201.0, 202.0, 203.0], [210.0, 211.0, 212.0, 213.0]],
    ],
    dtype=np.float32,
)


def test_load_midpoints_reads_metadata_and_arrays(tmp_path: Path) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")

    profiles = cfdnalab.load_midpoints(store_path)

    assert profiles.path == store_path
    assert profiles.group_names() == GROUP_NAMES.tolist()
    assert profiles.eligible_intervals() == [4, 7, 11]
    pd.testing.assert_frame_equal(
        profiles.groups(),
        pd.DataFrame(
            {
                "group_idx": GROUP_AXIS,
                "group_name": GROUP_NAMES,
                "eligible_intervals": ELIGIBLE_INTERVALS,
            }
        ),
    )
    pd.testing.assert_frame_equal(
        profiles.length_bins(),
        pd.DataFrame(
            {
                "length_bin": LENGTH_BIN,
                "length_start_bp": LENGTH_START_BP,
                "length_end_bp": LENGTH_END_BP,
            }
        ),
    )
    pd.testing.assert_frame_equal(
        profiles.positions(),
        pd.DataFrame(
            {
                "position": POSITION,
                "position_bin_start_bp": POSITION_BIN_START_BP,
                "position_bin_end_bp": POSITION_BIN_END_BP,
            }
        ),
    )
    np.testing.assert_array_equal(profiles.array(), COUNTS)


def test_profile_dataframe_uses_selected_group_length_and_positions(
    tmp_path: Path,
) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")
    profiles = cfdnalab.load_midpoints(store_path)

    data_frame = profiles.data_frame_for_profile(group_idx=1, length_bin_idx=0)

    pd.testing.assert_frame_equal(
        data_frame,
        pd.DataFrame(
            {
                "group_idx": [1, 1, 1, 1],
                "group_name": ["beta long"] * 4,
                "eligible_intervals": [7, 7, 7, 7],
                "length_bin": [0, 0, 0, 0],
                "length_start_bp": [30, 30, 30, 30],
                "length_end_bp": [60, 60, 60, 60],
                "position": POSITION,
                "position_bin_start_bp": POSITION_BIN_START_BP,
                "position_bin_end_bp": POSITION_BIN_END_BP,
                "count": np.array([100.0, 101.0, 102.0, 103.0], dtype=np.float32),
            }
        ),
    )


def test_group_dataframe_preserves_length_position_grid(tmp_path: Path) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")
    profiles = cfdnalab.load_midpoints(store_path)

    data_frame = profiles.data_frame_from_group("gamma-unicode-aa")

    assert data_frame.shape == (8, 10)
    assert data_frame["group_idx"].tolist() == [2] * 8
    assert data_frame["group_name"].tolist() == ["gamma-unicode-aa"] * 8
    assert data_frame["eligible_intervals"].tolist() == [11] * 8
    assert data_frame["length_bin"].tolist() == [0, 0, 0, 0, 1, 1, 1, 1]
    assert data_frame["position"].tolist() == [0, 1, 2, 3, 0, 1, 2, 3]
    assert data_frame["count"].tolist() == [
        200.0,
        201.0,
        202.0,
        203.0,
        210.0,
        211.0,
        212.0,
        213.0,
    ]


def test_length_dataframe_preserves_group_position_grid(tmp_path: Path) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")
    profiles = cfdnalab.load_midpoints(store_path)

    data_frame = profiles.data_frame_from_length(60)

    assert data_frame.shape == (12, 10)
    assert data_frame["group_idx"].tolist() == [
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        2,
    ]
    assert data_frame["group_name"].tolist() == [
        "alpha",
        "alpha",
        "alpha",
        "alpha",
        "beta long",
        "beta long",
        "beta long",
        "beta long",
        "gamma-unicode-aa",
        "gamma-unicode-aa",
        "gamma-unicode-aa",
        "gamma-unicode-aa",
    ]
    assert data_frame["length_bin"].tolist() == [1] * 12
    assert data_frame["length_start_bp"].tolist() == [60] * 12
    assert data_frame["length_end_bp"].tolist() == [90] * 12
    assert data_frame["position"].tolist() == [0, 1, 2, 3] * 3
    assert data_frame["count"].tolist() == [
        10.0,
        11.0,
        12.0,
        13.0,
        110.0,
        111.0,
        112.0,
        113.0,
        210.0,
        211.0,
        212.0,
        213.0,
    ]


def test_array_helpers_return_expected_slices(tmp_path: Path) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")
    profiles = cfdnalab.load_midpoints(store_path)

    np.testing.assert_array_equal(profiles.array_for_profile(2, 1), COUNTS[2, 1, :])
    np.testing.assert_array_equal(profiles.array_from_group("beta long"), COUNTS[1, :, :])
    np.testing.assert_array_equal(profiles.array_from_group_idx(0), COUNTS[0, :, :])
    np.testing.assert_array_equal(profiles.array_from_length(59), COUNTS[:, 0, :])
    np.testing.assert_array_equal(profiles.array_from_length_bin(1), COUNTS[:, 1, :])


def test_lookup_helpers_use_names_and_half_open_length_bins(tmp_path: Path) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")
    profiles = cfdnalab.load_midpoints(store_path)

    assert profiles.group_idx("alpha") == 0
    assert profiles.group_idx("gamma-unicode-aa") == 2
    assert profiles.length_bin_idx(30) == 0
    assert profiles.length_bin_idx(59) == 0
    assert profiles.length_bin_idx(60) == 1
    assert profiles.length_bin_idx(89) == 1

    with pytest.raises(KeyError, match="Unknown midpoint group name"):
        profiles.group_idx("missing")
    with pytest.raises(KeyError, match="No midpoint length bin contains length 90"):
        profiles.length_bin_idx(90)
    with pytest.raises(ValueError, match="Fragment length must be non-negative"):
        profiles.length_bin_idx(-1)


def test_invalid_indices_raise_helpful_errors(tmp_path: Path) -> None:
    store_path = _write_midpoint_store(tmp_path / "sample.midpoint_profiles.zarr")
    profiles = cfdnalab.load_midpoints(store_path)

    with pytest.raises(TypeError, match="group_idx must be an integer"):
        profiles.array_for_profile("0", 0)  # type: ignore[arg-type]
    with pytest.raises(IndexError, match="group_idx 3 is outside 0..2"):
        profiles.array_for_profile(3, 0)
    with pytest.raises(IndexError, match="length_bin_idx 2 is outside 0..1"):
        profiles.array_for_profile(0, 2)


def test_loader_rejects_invalid_paths(tmp_path: Path) -> None:
    missing_path = tmp_path / "missing.midpoint_profiles.zarr"
    file_path = tmp_path / "not_a_directory.zarr"
    file_path.write_text("not a zarr store")
    wrong_suffix_path = tmp_path / "sample.not_zarr"
    wrong_suffix_path.mkdir()

    with pytest.raises(FileNotFoundError, match="does not exist"):
        cfdnalab.load_midpoints(missing_path)
    with pytest.raises(NotADirectoryError, match="not a directory"):
        cfdnalab.load_midpoints(file_path)
    with pytest.raises(ValueError, match="should end with '.zarr'"):
        cfdnalab.load_midpoints(wrong_suffix_path)


def test_loader_rejects_schema_and_shape_problems(tmp_path: Path) -> None:
    wrong_schema = _write_midpoint_store(
        tmp_path / "wrong_schema.midpoint_profiles.zarr",
        schema="other",
    )
    wrong_version = _write_midpoint_store(
        tmp_path / "wrong_version.midpoint_profiles.zarr",
        schema_version=99,
    )
    missing_array = _write_midpoint_store(
        tmp_path / "missing_array.midpoint_profiles.zarr",
        omit={"eligible_intervals"},
    )
    wrong_counts_dimensions = _write_midpoint_store(
        tmp_path / "wrong_dimensions.midpoint_profiles.zarr",
        counts_dimension_names=("length_bin", "group", "position"),
    )
    noncontiguous_axis = _write_midpoint_store(
        tmp_path / "noncontiguous_axis.midpoint_profiles.zarr",
        group_axis=np.array([0, 2, 3], dtype=np.int32),
    )
    shape_mismatch = _write_midpoint_store(
        tmp_path / "shape_mismatch.midpoint_profiles.zarr",
        counts=COUNTS[:, :, :3],
    )
    missing_group_labels = _write_midpoint_store(
        tmp_path / "missing_group_labels.midpoint_profiles.zarr",
        group_labels=None,
    )

    with pytest.raises(ValueError, match="Expected cfdnalab_schema"):
        cfdnalab.load_midpoints(wrong_schema)
    with pytest.raises(ValueError, match="Unsupported midpoint schema version"):
        cfdnalab.load_midpoints(wrong_version)
    with pytest.raises(ValueError, match="missing arrays: \\['eligible_intervals'\\]"):
        cfdnalab.load_midpoints(missing_array)
    with pytest.raises(ValueError, match="counts dimensions must be"):
        cfdnalab.load_midpoints(wrong_counts_dimensions)
    with pytest.raises(ValueError, match="group axis must be contiguous"):
        cfdnalab.load_midpoints(noncontiguous_axis)
    with pytest.raises(ValueError, match="counts shape does not match"):
        cfdnalab.load_midpoints(shape_mismatch)
    with pytest.raises(ValueError, match="group array is missing group-name labels"):
        cfdnalab.load_midpoints(missing_group_labels)


def _write_midpoint_store(
    path: Path,
    *,
    schema: str = "midpoint_profiles",
    schema_version: int = 1,
    omit: set[str] | None = None,
    group_axis: np.ndarray = GROUP_AXIS,
    group_labels: list[str] | None = GROUP_NAMES.tolist(),
    counts: np.ndarray = COUNTS,
    counts_dimension_names: tuple[str, str, str] = ("group", "length_bin", "position"),
) -> Path:
    omitted_arrays = omit or set()
    root = zarr.open_group(str(path), mode="w", zarr_format=3)
    root.attrs["cfdnalab_schema"] = schema
    root.attrs["cfdnalab_schema_version"] = schema_version

    _create_array(
        root,
        "counts",
        counts,
        chunks=(2, 2, 3),
        dimension_names=counts_dimension_names,
    )
    group = _create_array(root, "group", group_axis, chunks=(len(group_axis),))
    group.attrs["label_field"] = "group_name"
    if group_labels is not None:
        group.attrs["labels"] = group_labels
    _create_array(root, "eligible_intervals", ELIGIBLE_INTERVALS, chunks=(3,))
    _create_array(root, "length_bin", LENGTH_BIN, chunks=(2,))
    _create_array(root, "length_start_bp", LENGTH_START_BP, chunks=(2,))
    _create_array(root, "length_end_bp", LENGTH_END_BP, chunks=(2,))
    _create_array(root, "position", POSITION, chunks=(4,))
    _create_array(root, "position_bin_start_bp", POSITION_BIN_START_BP, chunks=(4,))
    _create_array(root, "position_bin_end_bp", POSITION_BIN_END_BP, chunks=(4,))

    for array_name in omitted_arrays:
        del root[array_name]

    return path


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
