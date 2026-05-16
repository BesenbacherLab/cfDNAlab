"""
Classes for loading and interacting with the midpoints .zarr output.
"""

from __future__ import annotations

from dataclasses import dataclass
import numbers
import pathlib
from typing import Any, List

import numpy as np
import pandas as pd
import zarr

MIN_SUPPORTED_SCHEMA_VERSION = 1
MAX_SUPPORTED_SCHEMA_VERSION = 1
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
    def __init__(self, path: pathlib.Path | str) -> None:
        """
        Load a midpoint profile Zarr store.

        Parameters
        ----------
        path
            Path to a `<prefix>.midpoint_profiles.zarr` directory.

        Returns
        -------
        None
        """
        self.path = pathlib.Path(path)
        self.profiles = MidpointProfiles._load_zarr(self.path)

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

        group_idx = _read_array(store, "group")
        length_bin = _read_array(store, "length_bin")
        position = _read_array(store, "position")
        eligible_intervals = _read_array(store, "eligible_intervals")
        length_start_bp = _read_array(store, "length_start_bp")
        length_end_bp = _read_array(store, "length_end_bp")
        position_bin_start_bp = _read_array(store, "position_bin_start_bp")
        position_bin_end_bp = _read_array(store, "position_bin_end_bp")
        group_names = _read_group_names(store)

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

    def group_names(self) -> List[str]:
        """
        Return the group names.

        Returns
        -------
        list[str]
            Group names in count-array row order.
        """
        return self.profiles.group_names.tolist()

    def eligible_intervals(self) -> List[int]:
        """
        Return the number of eligible intervals for each group.

        Returns
        -------
        list[int]
            Eligible interval counts in count-array row order.
        """
        return self.profiles.eligible_intervals.astype(int).tolist()

    def group_idx(self, group_name: str) -> int:
        """
        Return the index for a group name.

        Parameters
        ----------
        group_name
            Group name to resolve.

        Returns
        -------
        int
            Zero-based group index.
        """
        return self._resolve_group_name(group_name)

    def length_bin_idx(self, length: int) -> int:
        """
        Return the length-bin index for a fragment length.

        Parameters
        ----------
        length
            Fragment length in bp.

        Returns
        -------
        int
            Zero-based length-bin index.
        """
        return self._resolve_length(length)

    def groups(self) -> pd.DataFrame:
        """
        Return group metadata.

        Returns
        -------
        pandas.DataFrame
            One row per group.
        """
        return pd.DataFrame(
            {
                "group_idx": self.profiles.group_idx,
                "group_name": self.profiles.group_names,
                "eligible_intervals": self.profiles.eligible_intervals,
            }
        )

    def length_bins(self) -> pd.DataFrame:
        """
        Return fragment length-bin metadata.

        Returns
        -------
        pandas.DataFrame
            One row per length bin.
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
        Return midpoint position-bin metadata.

        Returns
        -------
        pandas.DataFrame
            One row per midpoint position bin.
        """
        return pd.DataFrame(
            {
                "position": self.profiles.position,
                "position_bin_start_bp": self.profiles.position_bin_start_bp,
                "position_bin_end_bp": self.profiles.position_bin_end_bp,
            }
        )

    def data_frame_for_profile(
        self, group_idx: int, length_bin_idx: int
    ) -> pd.DataFrame:
        """
        Build a long data frame for one group and one length bin.

        Use `.group_idx(group_name=)` and `.length_bin_idx(length=)` to get the
        indices for a group name and fragment length.

        Parameters
        ----------
        group_idx
            Zero-based group index to extract.
        length_bin_idx
            Zero-based length-bin index to extract.

        Returns
        -------
        pandas.DataFrame
            Counts and position metadata for one profile.
        """
        group_idx = self._validate_group_idx(group_idx)
        length_bin_idx = self._validate_length_bin_idx(length_bin_idx)
        profile = self.array_for_profile(group_idx, length_bin_idx)

        return pd.DataFrame(
            {
                "group_idx": int(self.profiles.group_idx[group_idx]),
                "group_name": self.profiles.group_names[group_idx],
                "eligible_intervals": int(self.profiles.eligible_intervals[group_idx]),
                "length_bin": int(self.profiles.length_bin[length_bin_idx]),
                "length_start_bp": int(self.profiles.length_start_bp[length_bin_idx]),
                "length_end_bp": int(self.profiles.length_end_bp[length_bin_idx]),
                "position": self.profiles.position,
                "position_bin_start_bp": self.profiles.position_bin_start_bp,
                "position_bin_end_bp": self.profiles.position_bin_end_bp,
                "count": profile,
            }
        )

    def data_frame_from_group(self, group_name: str) -> pd.DataFrame:
        """
        Build a long data frame for one group.

        Parameters
        ----------
        group_name
            Group name to extract.

        Returns
        -------
        pandas.DataFrame
            Counts and length/position metadata for the group.
        """
        group_idx = self._resolve_group_name(group_name)
        return self.data_frame_from_group_idx(group_idx)

    def data_frame_from_group_idx(self, group_idx: int) -> pd.DataFrame:
        """
        Build a long data frame for one group index.

        Parameters
        ----------
        group_idx
            Zero-based group index to extract.

        Returns
        -------
        pandas.DataFrame
            Counts and length/position metadata for the group.
        """
        group_idx = self._validate_group_idx(group_idx)
        profile = self.array_from_group_idx(group_idx)
        length_index, position_index = np.indices(profile.shape)

        return pd.DataFrame(
            {
                "group_idx": int(self.profiles.group_idx[group_idx]),
                "group_name": self.profiles.group_names[group_idx],
                "eligible_intervals": int(self.profiles.eligible_intervals[group_idx]),
                "length_bin": self.profiles.length_bin[length_index.ravel()],
                "length_start_bp": self.profiles.length_start_bp[length_index.ravel()],
                "length_end_bp": self.profiles.length_end_bp[length_index.ravel()],
                "position": self.profiles.position[position_index.ravel()],
                "position_bin_start_bp": self.profiles.position_bin_start_bp[
                    position_index.ravel()
                ],
                "position_bin_end_bp": self.profiles.position_bin_end_bp[
                    position_index.ravel()
                ],
                "count": profile.ravel(),
            }
        )

    def data_frame_from_length(self, length: int) -> pd.DataFrame:
        """
        Build a long data frame for the bin containing a fragment length.

        Parameters
        ----------
        length
            Fragment length in bp.

        Returns
        -------
        pandas.DataFrame
            Counts and group/position metadata for the matching length bin.
        """
        length_bin_idx = self._resolve_length(length)
        return self.data_frame_from_length_bin(length_bin_idx)

    def data_frame_from_length_bin(self, length_bin_idx: int) -> pd.DataFrame:
        """
        Build a long data frame for one length-bin index.

        Parameters
        ----------
        length_bin_idx
            Zero-based length-bin index to extract.

        Returns
        -------
        pandas.DataFrame
            Counts and group/position metadata for the length bin.
        """
        length_bin_idx = self._validate_length_bin_idx(length_bin_idx)
        profile = self.array_from_length_bin(length_bin_idx)
        group_index, position_index = np.indices(profile.shape)

        return pd.DataFrame(
            {
                "group_idx": self.profiles.group_idx[group_index.ravel()],
                "group_name": self.profiles.group_names[group_index.ravel()],
                "eligible_intervals": self.profiles.eligible_intervals[
                    group_index.ravel()
                ],
                "length_bin": int(self.profiles.length_bin[length_bin_idx]),
                "length_start_bp": int(self.profiles.length_start_bp[length_bin_idx]),
                "length_end_bp": int(self.profiles.length_end_bp[length_bin_idx]),
                "position": self.profiles.position[position_index.ravel()],
                "position_bin_start_bp": self.profiles.position_bin_start_bp[
                    position_index.ravel()
                ],
                "position_bin_end_bp": self.profiles.position_bin_end_bp[
                    position_index.ravel()
                ],
                "count": profile.ravel(),
            }
        )

    def array_for_profile(self, group_idx: int, length_bin_idx: int) -> np.ndarray:
        """
        Load counts for one group and one length bin.

        Use `.group_idx(group_name=)` and `.length_bin_idx(length=)` to get the
        indices for a group name and fragment length.

        Parameters
        ----------
        group_idx
            Zero-based group index to extract.
        length_bin_idx
            Zero-based length-bin index to extract.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(position,)`.
        """
        group_idx = self._validate_group_idx(group_idx)
        length_bin_idx = self._validate_length_bin_idx(length_bin_idx)
        return np.asarray(self.profiles.counts[group_idx, length_bin_idx, :])

    def array(self) -> np.ndarray:
        """
        Load the full midpoint count tensor.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(group, length_bin, position)`.
        """
        return np.asarray(self.profiles.counts[:])

    def array_from_group(self, group_name: str) -> np.ndarray:
        """
        Load counts for one group.

        Parameters
        ----------
        group_name
            Group name to extract.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(length_bin, position)`.
        """
        group_idx = self._resolve_group_name(group_name)
        return self.array_from_group_idx(group_idx)

    def array_from_group_idx(self, group_idx: int) -> np.ndarray:
        """
        Load counts for one group index.

        Parameters
        ----------
        group_idx
            Zero-based group index to extract.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(length_bin, position)`.
        """
        group_idx = self._validate_group_idx(group_idx)
        return np.asarray(self.profiles.counts[group_idx, :, :])

    def array_from_length(self, length: int) -> np.ndarray:
        """
        Load counts for the bin containing a fragment length.

        Parameters
        ----------
        length
            Fragment length in bp.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(group, position)`.
        """
        length_bin_idx = self._resolve_length(length)
        return self.array_from_length_bin(length_bin_idx)

    def array_from_length_bin(self, length_bin_idx: int) -> np.ndarray:
        """
        Load counts for one length-bin index.

        Parameters
        ----------
        length_bin_idx
            Zero-based length-bin index to extract.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(group, position)`.
        """
        length_bin_idx = self._validate_length_bin_idx(length_bin_idx)
        return np.asarray(self.profiles.counts[:, length_bin_idx, :])

    def _resolve_group_name(self, group_name: str) -> int:
        """
        Resolve a group name to its zero-based group index.

        Parameters
        ----------
        group_name
            Group name to resolve.

        Returns
        -------
        int
            Zero-based group index.
        """
        matches = np.flatnonzero(self.profiles.group_names == group_name)
        if len(matches) == 0:
            raise KeyError(f"Unknown midpoint group name: {group_name!r}")
        if len(matches) > 1:
            raise ValueError(f"Midpoint group name is not unique: {group_name!r}")
        return int(matches[0])

    def _resolve_length(self, length: int) -> int:
        """
        Resolve a fragment length to its zero-based length-bin index.

        Parameters
        ----------
        length
            Fragment length in bp.

        Returns
        -------
        int
            Zero-based length-bin index.
        """
        if length < 0:
            raise ValueError(f"Fragment length must be non-negative, got {length}")
        matches = np.flatnonzero(
            (self.profiles.length_start_bp <= length)
            & (length < self.profiles.length_end_bp)
        )
        if len(matches) == 0:
            raise KeyError(f"No midpoint length bin contains length {length}")
        if len(matches) > 1:
            raise ValueError(f"Multiple midpoint length bins contain length {length}")
        return int(matches[0])

    def _validate_group_idx(self, group_idx: int) -> int:
        """
        Validate and normalize a group index.

        Parameters
        ----------
        group_idx
            Zero-based group index.

        Returns
        -------
        int
            Validated group index.
        """
        return _validate_index(group_idx, len(self.profiles.group_idx), "group_idx")

    def _validate_length_bin_idx(self, length_bin_idx: int) -> int:
        """
        Validate and normalize a length-bin index.

        Parameters
        ----------
        length_bin_idx
            Zero-based length-bin index.

        Returns
        -------
        int
            Validated length-bin index.
        """
        return _validate_index(
            length_bin_idx, len(self.profiles.length_bin), "length_bin_idx"
        )


def load_midpoints(path: pathlib.Path | str) -> MidpointProfiles:
    """
    Load a midpoint profile Zarr store.

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
    schema = store.attrs.get("cfdnalab_schema")
    if schema != "midpoint_profiles":
        raise ValueError(
            f"Expected cfdnalab_schema='midpoint_profiles', found {schema!r}"
        )

    schema_version = store.attrs.get("cfdnalab_schema_version")
    if not isinstance(schema_version, numbers.Integral) or not (
        MIN_SUPPORTED_SCHEMA_VERSION <= int(schema_version) <= MAX_SUPPORTED_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unsupported midpoint schema version: "
            f"{schema_version!r}. Supported range: "
            f"{MIN_SUPPORTED_SCHEMA_VERSION}..{MAX_SUPPORTED_SCHEMA_VERSION}"
        )


def _validate_required_arrays(store: Any) -> None:
    missing = sorted(name for name in REQUIRED_ARRAYS if not _has_array(store, name))
    if missing:
        raise ValueError(f"Midpoint Zarr store is missing arrays: {missing}")


def _has_array(store: Any, name: str) -> bool:
    try:
        store[name]
    except Exception:
        return False
    return True


def _validate_counts_dimensions(counts: Any) -> None:
    dimension_names = tuple(getattr(counts.metadata, "dimension_names", ()) or ())
    expected = ("group", "length_bin", "position")
    if dimension_names != expected:
        raise ValueError(
            f"counts dimensions must be {expected}, found {dimension_names}"
        )


def _read_array(store: Any, name: str) -> np.ndarray:
    return np.asarray(store[name][:])


def _read_group_names(store: Any) -> np.ndarray:
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
    expected = np.arange(len(values), dtype=values.dtype)
    if not np.array_equal(values, expected):
        raise ValueError(f"{name} axis must be contiguous 0-based indices")


def _validate_same_length(
    values: np.ndarray, axis_values: np.ndarray, values_name: str, axis_name: str
) -> None:
    if len(values) != len(axis_values):
        raise ValueError(
            f"{values_name} length ({len(values)}) does not match "
            f"{axis_name} length ({len(axis_values)})"
        )


def _validate_index(index: int, size: int, name: str) -> int:
    if not isinstance(index, numbers.Integral):
        raise TypeError(f"{name} must be an integer, got {type(index).__name__}")
    index = int(index)
    if index < 0 or index >= size:
        raise IndexError(f"{name} {index} is outside 0..{size - 1}")
    return index
