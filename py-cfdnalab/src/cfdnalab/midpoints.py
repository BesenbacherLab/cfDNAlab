"""
Load cfDNAlab midpoint profile Zarr outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
import numbers
import pathlib
from typing import Any, Sequence

import numpy as np
import pandas as pd
import zarr

from ._helpers import (
    normalize_strings,
    normalize_length_bin_selector,
    normalize_zero_based_indices,
    resolve_unique_match,
    validate_fragment_length,
)

MIDPOINT_MIN_SUPPORTED_SCHEMA_VERSION = 1
MIDPOINT_MAX_SUPPORTED_SCHEMA_VERSION = 1
REQUIRED_ARRAYS = {
    "counts",
    "group",
    "eligible_intervals",
    "length_bin",
    "length_start_bp",
    "length_end_bp",
    "position",
    "position_bin_start_bp",
    "position_bin_end_bp",
}


@dataclass
class LoadedProfile:
    """Validated midpoint Zarr handles and axis metadata."""

    store: Any
    counts: zarr.Array
    group_idx: np.ndarray
    group_names: np.ndarray
    eligible_intervals: np.ndarray
    length_bin: np.ndarray
    length_start_bp: np.ndarray
    length_end_bp: np.ndarray
    position: np.ndarray
    position_bin_start_bp: np.ndarray
    position_bin_end_bp: np.ndarray


class MidpointProfiles:
    """
    Helper for loading and slicing midpoint profile Zarr output.

    Midpoint profiles store counts as `(group, length_bin, position)`. The class
    exposes metadata as pandas data frames and count slices as NumPy arrays.
    """

    def __init__(self, path: pathlib.Path | str) -> None:
        """
        Load a midpoint profile Zarr store.

        Parameters
        ----------
        path
            Path to a `<prefix>.midpoint_profiles.zarr` directory.
        """
        self.path = pathlib.Path(path)
        self.profiles = MidpointProfiles._load_zarr(self.path)

    def __repr__(self) -> str:
        """
        Return a compact summary with path, schema version, and count shape.
        """
        shape = tuple(self.profiles.counts.shape)
        schema_version = self.profiles.store.attrs.get("cfdnalab_schema_version")
        return (
            "MidpointProfiles("
            f"path={str(self.path)!r}, "
            f"schema_version={schema_version!r}, "
            f"shape={shape!r}"
            ")"
        )

    @staticmethod
    def _load_zarr(path: pathlib.Path | str) -> LoadedProfile:
        """
        Open and validate a midpoint profile Zarr store.

        Parameters
        ----------
        path
            Path to a `<prefix>.midpoint_profiles.zarr` directory.

        Returns
        -------
        LoadedProfile
            Loaded Zarr handles and axis metadata.
        """
        path = pathlib.Path(path)
        _validate_zarr_store_path(path)

        try:
            store = zarr.open_group(str(path), mode="r", zarr_format=3)
        except Exception as error:
            raise ValueError(f"Could not open midpoint Zarr store at {path}") from error
        _validate_root_metadata(store)
        _validate_required_arrays(store)

        counts = store["counts"]
        _validate_counts_dimensions(counts)

        # Load small coordinate arrays eagerly so later helpers only slice counts
        group_idx = _read_array(store, "group")
        length_bin = _read_array(store, "length_bin")
        position = _read_array(store, "position")
        eligible_intervals = _read_array(store, "eligible_intervals")
        length_start_bp = _read_array(store, "length_start_bp")
        length_end_bp = _read_array(store, "length_end_bp")
        position_bin_start_bp = _read_array(store, "position_bin_start_bp")
        position_bin_end_bp = _read_array(store, "position_bin_end_bp")
        group_names = _read_group_names(store)

        # The count tensor must align with every public coordinate axis
        expected_shape = (len(group_idx), len(length_bin), len(position))
        if tuple(counts.shape) != expected_shape:
            raise ValueError(
                "counts shape does not match coordinate arrays: "
                f"counts={counts.shape}, coordinates={expected_shape}"
            )

        _validate_axis(group_idx, "group")
        _validate_axis(length_bin, "length_bin")
        _validate_axis(position, "position")
        _validate_same_length(group_names, group_idx, "group_name", "group")
        _validate_same_length(
            eligible_intervals, group_idx, "eligible_intervals", "group"
        )
        _validate_same_length(
            length_start_bp, length_bin, "length_start_bp", "length_bin"
        )
        _validate_same_length(length_end_bp, length_bin, "length_end_bp", "length_bin")
        _validate_same_length(
            position_bin_start_bp, position, "position_bin_start_bp", "position"
        )
        _validate_same_length(
            position_bin_end_bp, position, "position_bin_end_bp", "position"
        )

        return LoadedProfile(
            store=store,
            counts=counts,
            group_idx=group_idx,
            group_names=group_names,
            eligible_intervals=eligible_intervals,
            length_bin=length_bin,
            length_start_bp=length_start_bp,
            length_end_bp=length_end_bp,
            position=position,
            position_bin_start_bp=position_bin_start_bp,
            position_bin_end_bp=position_bin_end_bp,
        )

    def group_idx(self, group_name: str) -> int:
        """
        Find the midpoint group index for a group name.

        Parameters
        ----------
        group_name
            Group name to resolve.

        Returns
        -------
        int
            Group index.
        """
        return self._resolve_group_name(group_name)

    def length_bin_idx(self, length: int) -> int:
        """
        Find the length-bin index whose interval contains a fragment length.

        Parameters
        ----------
        length
            Fragment length in bp.

        Returns
        -------
        int
            Length-bin index.
        """
        return self._resolve_length(length)

    def group_metadata(self) -> pd.DataFrame:
        """
        Return midpoint group labels and eligible interval counts.

        Returns
        -------
        pandas.DataFrame
            Columns are `group_idx`, `group_name`, and `eligible_intervals`.
        """
        return pd.DataFrame(
            {
                "group_idx": self.profiles.group_idx,
                "group_name": self.profiles.group_names,
                "eligible_intervals": self.profiles.eligible_intervals,
            }
        )

    def counts_array(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
    ) -> np.ndarray:
        """
        Return midpoint counts as a dense NumPy array.

        The result keeps the midpoint count dimensions in the same order as
        the file: group, length bin, then position. Scalar selectors keep their
        axis as length one, so the shape is always
        `(selected groups, selected length bins, positions)`.

        Parameters
        ----------
        groups
            `None` for all groups, one group name, or a sequence of group names.
            Use either `groups` or `group_idxs`, not both.
        group_idxs
            `None` for all groups, one group index, or a sequence of group
            indices. Use either `groups` or `group_idxs`, not both.
        with_lengths
            Fragment length or lengths in bp. Counts are returned for the
            length bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`. Counts are
            returned for whole length bins that overlap this range.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(group, length_bin, position)`.
        """
        group_indices = self._resolve_group_selector(groups, group_idxs)
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        counts = np.empty(
            (
                len(group_indices),
                len(length_bin_indices),
                len(self.profiles.position),
            ),
            dtype=self.profiles.counts.dtype,
        )
        for output_group_idx, group_idx in enumerate(group_indices):
            for output_length_bin_idx, length_bin_idx in enumerate(length_bin_indices):
                counts[output_group_idx, output_length_bin_idx, :] = np.asarray(
                    self.profiles.counts[int(group_idx), int(length_bin_idx), :]
                )
        return counts

    def length_bins(self) -> pd.DataFrame:
        """
        Get the fragment length bins available in this midpoint-profile output.

        Length bins are half-open intervals. A bin with `length_start_bp=30`
        and `length_end_bp=50` contains fragment lengths `30 <= length < 50`.

        Returns
        -------
        pandas.DataFrame
            Columns are `length_bin`, `length_start_bp`, and `length_end_bp`.
        """
        return pd.DataFrame(
            {
                "length_bin": self.profiles.length_bin,
                "length_start_bp": self.profiles.length_start_bp,
                "length_end_bp": self.profiles.length_end_bp,
            }
        )

    def positions(self) -> pd.DataFrame:
        """
        Get the midpoint position bins available in this output.

        Returns
        -------
        pandas.DataFrame
            Columns are `position`, `position_bin_start_bp`, and
            `position_bin_end_bp`.
        """
        return pd.DataFrame(
            {
                "position": self.profiles.position,
                "position_bin_start_bp": self.profiles.position_bin_start_bp,
                "position_bin_end_bp": self.profiles.position_bin_end_bp,
            }
        )

    def data_frame(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of midpoint profile counts.

        Use this for tabular analysis of the midpoint count array. The result
        expands the selected group and length-bin axes across all midpoint
        position bins, with group, length-bin, and position metadata on each
        row.

        Parameters
        ----------
        groups
            `None` for all groups, one group name, or a sequence of group names.
            Use either `groups` or `group_idxs`, not both.
        group_idxs
            `None` for all groups, one group index, or a sequence of group
            indices. Use either `groups` or `group_idxs`, not both.
        with_lengths
            Fragment length or lengths in bp. The returned rows use the length
            bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`. Returned
            rows use whole length bins that overlap this range.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.

        Returns
        -------
        pandas.DataFrame
            One row per selected group, length bin, and midpoint position bin.
        """
        group_indices = self._resolve_group_selector(groups, group_idxs)
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self._data_frame_for_indices(group_indices, length_bin_indices)

    def _data_frame_for_indices(
        self, group_indices: np.ndarray, length_bin_indices: np.ndarray
    ) -> pd.DataFrame:
        """
        Build a long data frame for selected group and length-bin indices.
        """
        if len(group_indices) == 0 or len(length_bin_indices) == 0:
            return self._empty_data_frame()

        frames: list[pd.DataFrame] = []
        counts = self.counts_array(
            group_idxs=group_indices,
            length_bin_idxs=length_bin_indices,
        )
        for output_group_idx, group_idx in enumerate(group_indices):
            for output_length_bin_idx, length_bin_idx in enumerate(length_bin_indices):
                profile = counts[output_group_idx, output_length_bin_idx, :]
                frames.append(
                    pd.DataFrame(
                        {
                            "group_idx": int(self.profiles.group_idx[group_idx]),
                            "group_name": self.profiles.group_names[group_idx],
                            "eligible_intervals": int(
                                self.profiles.eligible_intervals[group_idx]
                            ),
                            "length_bin": int(
                                self.profiles.length_bin[length_bin_idx]
                            ),
                            "length_start_bp": int(
                                self.profiles.length_start_bp[length_bin_idx]
                            ),
                            "length_end_bp": int(
                                self.profiles.length_end_bp[length_bin_idx]
                            ),
                            "position": self.profiles.position,
                            "position_bin_start_bp": self.profiles.position_bin_start_bp,
                            "position_bin_end_bp": self.profiles.position_bin_end_bp,
                            "count": profile,
                        }
                    )
                )
        return pd.concat(frames, ignore_index=True)

    def _empty_data_frame(self) -> pd.DataFrame:
        """
        Return an empty midpoint data frame with public columns.
        """
        return pd.DataFrame(
            {
                "group_idx": self.profiles.group_idx[:0],
                "group_name": self.profiles.group_names[:0],
                "eligible_intervals": self.profiles.eligible_intervals[:0],
                "length_bin": self.profiles.length_bin[:0],
                "length_start_bp": self.profiles.length_start_bp[:0],
                "length_end_bp": self.profiles.length_end_bp[:0],
                "position": self.profiles.position[:0],
                "position_bin_start_bp": self.profiles.position_bin_start_bp[:0],
                "position_bin_end_bp": self.profiles.position_bin_end_bp[:0],
                "count": np.asarray([], dtype=float),
            }
        )

    def _resolve_group_selector(
        self,
        groups: str | Sequence[str] | None,
        group_idxs: int | Sequence[int] | None,
    ) -> np.ndarray:
        """
        Normalize optional group selectors to group-axis indices.
        """
        if groups is not None and group_idxs is not None:
            raise ValueError("Use either groups or group_idxs, not both")
        if groups is None and group_idxs is None:
            return np.arange(len(self.profiles.group_idx), dtype=np.int64)
        if groups is not None:
            group_names = normalize_strings(groups, name="groups")
            return np.asarray(
                [self._resolve_group_name(group_name) for group_name in group_names],
                dtype=np.int64,
            )
        return normalize_zero_based_indices(
            group_idxs,
            size=len(self.profiles.group_idx),
            name="group_idxs",
            index_name="group_idx",
        )

    def _resolve_length_bin_selector(
        self,
        with_lengths: int | Sequence[int] | None,
        with_length_range: Sequence[int] | None,
        length_bin_idxs: int | Sequence[int] | None,
    ) -> np.ndarray:
        """
        Normalize optional fragment length or length-bin selectors.
        """
        return normalize_length_bin_selector(
            with_lengths=with_lengths,
            with_length_range=with_length_range,
            length_bin_idxs=length_bin_idxs,
            length_start_bp=self.profiles.length_start_bp,
            length_end_bp=self.profiles.length_end_bp,
            selector_context="midpoint",
        )

    def _resolve_group_name(self, group_name: str) -> int:
        """
        Resolve a group name to its group index.

        Parameters
        ----------
        group_name
            Group name to resolve.

        Returns
        -------
        int
            Group index.
        """
        return resolve_unique_match(
            self.profiles.group_names == group_name,
            missing_message=f"Unknown midpoint group name: {group_name!r}",
            duplicate_message=f"Midpoint group name is not unique: {group_name!r}",
        )

    def _resolve_length(self, length: int) -> int:
        """
        Resolve a fragment length to its length-bin index.

        Parameters
        ----------
        length
            Fragment length in bp.

        Returns
        -------
        int
            Length-bin index.
        """
        length = validate_fragment_length(length)
        return resolve_unique_match(
            (self.profiles.length_start_bp <= length)
            & (length < self.profiles.length_end_bp),
            missing_message=f"No midpoint length bin contains length {length}",
            duplicate_message=f"Multiple midpoint length bins contain length {length}",
        )


def read_midpoints(path: pathlib.Path | str) -> MidpointProfiles:
    """
    Open a cfDNAlab midpoint profile Zarr store.

    Parameters
    ----------
    path
        Path to a `.midpoint_profiles.zarr` directory.

    Returns
    -------
    MidpointProfiles
        Loaded midpoint profile helper.
    """
    return MidpointProfiles(path)


def _validate_zarr_store_path(path: pathlib.Path) -> None:
    """
    Validate that a path points to a Zarr V3 midpoint store directory.
    """
    if not path.exists():
        raise FileNotFoundError(f"Midpoint Zarr store does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(
            f"Midpoint Zarr store path exists but is not a directory: {path}"
        )
    if path.suffix != ".zarr":
        raise ValueError(
            f"Midpoint Zarr store path should end with '.zarr', got: {path}"
        )
    if not (path / "zarr.json").is_file():
        raise ValueError(
            f"Midpoint Zarr store is missing Zarr V3 metadata file: {path / 'zarr.json'}"
        )


def _validate_root_metadata(store: Any) -> None:
    """
    Validate root attributes that identify the midpoint schema version.
    """
    schema = store.attrs.get("cfdnalab_schema")
    if schema != "midpoint_profiles":
        raise ValueError(
            f"Expected cfdnalab_schema='midpoint_profiles', found {schema!r}"
        )

    schema_version = store.attrs.get("cfdnalab_schema_version")
    if not isinstance(schema_version, numbers.Integral) or not (
        MIDPOINT_MIN_SUPPORTED_SCHEMA_VERSION
        <= int(schema_version)
        <= MIDPOINT_MAX_SUPPORTED_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unsupported midpoint schema version: "
            f"{schema_version!r}. Supported range: "
            f"{MIDPOINT_MIN_SUPPORTED_SCHEMA_VERSION}..{MIDPOINT_MAX_SUPPORTED_SCHEMA_VERSION}"
        )


def _validate_required_arrays(store: Any) -> None:
    """
    Require every array needed by the public midpoint API.
    """
    missing = sorted(name for name in REQUIRED_ARRAYS if not _has_array(store, name))
    if missing:
        raise ValueError(f"Midpoint Zarr store is missing arrays: {missing}")


def _has_array(store: Any, name: str) -> bool:
    """
    Return whether a Zarr group contains an array path.
    """
    try:
        store[name]
    except Exception:
        return False
    return True


def _validate_counts_dimensions(counts: Any) -> None:
    """
    Require midpoint counts to use the expected named axes.
    """
    dimension_names = tuple(getattr(counts.metadata, "dimension_names", ()) or ())
    expected = ("group", "length_bin", "position")
    if dimension_names != expected:
        raise ValueError(
            f"counts dimensions must be {expected}, found {dimension_names}"
        )


def _read_array(store: Any, name: str) -> np.ndarray:
    """
    Load a Zarr array fully into a NumPy array.
    """
    return np.asarray(store[name][:])


def _read_group_names(store: Any) -> np.ndarray:
    """
    Read group-name labels from the group axis metadata.
    """
    label_field = store["group"].attrs.get("label_field")
    if label_field != "group_name":
        raise ValueError(
            f"group labels must have label_field='group_name', found {label_field!r}"
        )
    labels = store["group"].attrs.get("labels")
    if labels is None:
        raise ValueError("group array is missing group-name labels")
    return np.asarray(labels, dtype=str)


def _validate_axis(values: np.ndarray, name: str) -> None:
    """
    Require an axis coordinate array to be contiguous zero-based indices.
    """
    expected = np.arange(len(values), dtype=values.dtype)
    if not np.array_equal(values, expected):
        raise ValueError(f"{name} axis must be contiguous 0-based indices")


def _validate_same_length(
    values: np.ndarray, axis_values: np.ndarray, values_name: str, axis_name: str
) -> None:
    """
    Require a metadata vector to have the same length as its axis.
    """
    if len(values) != len(axis_values):
        raise ValueError(
            f"{values_name} length ({len(values)}) does not match "
            f"{axis_name} length ({len(axis_values)})"
        )
