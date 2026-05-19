"""
Load cfDNAlab length-count TSV outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
import io
import os
import pathlib
import re
import shutil
import subprocess
from typing import Sequence

import numpy as np
import pandas as pd
import zstandard

from ._helpers import (
    filter_blacklisted_fraction,
    normalize_length_bin_selector,
    normalize_strings,
    normalize_zero_based_indices,
    resolve_unique_match,
    validate_fragment_length,
    validate_scalar_bool,
)

COUNT_COLUMN_PATTERN = re.compile(r"^count_([0-9]+)(?:_([0-9]+))?$")
VALID_LENGTH_VALUES = ("count", "fraction", "density")
VALID_DENOMINATORS = {"all_bins", "selected_bins"}


@dataclass
class LoadedLengths:
    """Validated length-count table split into axes, row metadata, and counts."""

    mode: str
    length_bin: np.ndarray
    length_start_bp: np.ndarray
    length_end_bp: np.ndarray
    count_columns: list[str]
    counts: np.ndarray
    row_metadata: pd.DataFrame


class LengthCounts:
    """Common API for global, windowed, and grouped length-count outputs."""

    def __init__(
        self, path: pathlib.Path | str, loaded_lengths: LoadedLengths | None = None
    ) -> None:
        """
        Create a length-count loader for a TSV path.

        Use `read_lengths()` for normal construction so the returned object
        matches the output mode.
        """
        self.path = _normalize_length_counts_path(path)
        self.lengths = loaded_lengths or self._load_tsv(self.path)

    def __repr__(self) -> str:
        """
        Return a compact summary with the source path, mode, and count shape.
        """
        return (
            f"{self.__class__.__name__}("
            f"path={str(self.path)!r}, "
            f"mode={self.lengths.mode!r}, "
            f"shape={self.lengths.counts.shape!r}"
            ")"
        )

    @staticmethod
    def _load_tsv(path: pathlib.Path) -> LoadedLengths:
        """
        Read a length-count TSV and validate its public schema.

        The TSV itself carries the metadata needed by the Python API. Count
        columns define the length bins, and leading metadata columns define
        whether the file is global, windowed, or grouped.
        """
        _validate_length_counts_path(path)
        header_columns = _read_length_count_header(path)
        table = _read_length_counts_table(path)
        if table.shape[1] == 0:
            raise ValueError("Length-count TSV must contain at least one count column")

        count_columns = _length_count_columns(header_columns)
        bin_metadata = _parse_length_count_columns(count_columns)
        mode = _infer_length_count_mode(header_columns, count_columns)
        row_metadata = _length_row_metadata(table, mode, count_columns[0])
        counts = _length_counts_matrix_from_table(table, count_columns)
        if mode == "global" and counts.shape[0] != 1:
            raise ValueError(
                f"Global length-count output must contain exactly one row, found {counts.shape[0]}"
            )
        if len(row_metadata) != counts.shape[0]:
            raise ValueError("Length-count row metadata does not match count row count")

        return LoadedLengths(
            mode=mode,
            length_bin=bin_metadata["length_bin"],
            length_start_bp=bin_metadata["length_start_bp"],
            length_end_bp=bin_metadata["length_end_bp"],
            count_columns=count_columns,
            counts=counts,
            row_metadata=row_metadata,
        )

    def length_bins(self) -> pd.DataFrame:
        """
        Return fragment length bin definitions used by the count columns.

        Length bins are half-open intervals. A bin with `length_start_bp=30`
        and `length_end_bp=50` contains fragment lengths `30 <= length < 50`.

        Returns
        -------
        pandas.DataFrame
            Columns are `length_bin`, `length_start_bp`, `length_end_bp`,
            `length_midpoint_bp`, and `length_width_bp`.
        """
        return pd.DataFrame(
            {
                "length_bin": self.lengths.length_bin,
                "length_start_bp": self.lengths.length_start_bp,
                "length_end_bp": self.lengths.length_end_bp,
                "length_midpoint_bp": (
                    self.lengths.length_start_bp + self.lengths.length_end_bp
                )
                / 2,
                "length_width_bp": self.lengths.length_end_bp
                - self.lengths.length_start_bp,
            }
        )

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

        Raises
        ------
        KeyError
            If no length bin contains `length`.
        """
        length = validate_fragment_length(length)
        return resolve_unique_match(
            (self.lengths.length_start_bp <= length)
            & (length < self.lengths.length_end_bp),
            missing_message=f"No length-count bin contains length {length}",
            duplicate_message=f"Multiple length-count bins contain length {length}",
        )

    def counts_array(
        self,
        *,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
    ) -> np.ndarray:
        """
        Return raw length counts as a dense NumPy array.

        Use `with_lengths`, `with_length_range`, or `length_bin_idxs` to select
        length bins. Range selection uses whole bins overlapping the half-open
        `[start, end)` bp range.

        Parameters
        ----------
        with_lengths
            Fragment length or lengths in bp. Counts are returned for the
            length bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(output row, length_bin)`. Output rows are
            windows for windowed output, groups for grouped output, and the
            single global summary row for global output.
        """
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self.lengths.counts[:, length_bin_indices].copy()

    def _data_frame_for_rows(
        self,
        row_indices: np.ndarray,
        *,
        value: str,
        denominator: str,
        keep_wide: bool,
        max_blacklisted_fraction: float,
        length_bin_indices: np.ndarray,
    ) -> pd.DataFrame:
        """
        Assemble selected output rows as a long or wide data frame.

        Row selection happens before blacklist filtering so explicit selectors
        and `max_blacklisted_fraction` compose the same way for windows and
        groups.
        """
        row_indices = filter_blacklisted_fraction(
            self.lengths.row_metadata,
            max_blacklisted_fraction,
            row_indices=row_indices,
        )
        value = _validate_value(value)
        denominator = _validate_denominator(denominator)
        keep_wide = validate_scalar_bool(keep_wide, "keep_wide")
        values = self._value_matrix(
            row_indices,
            length_bin_indices,
            value,
            denominator,
        )
        if keep_wide:
            return self._wide_data_frame(row_indices, length_bin_indices, values, value)
        return self._long_data_frame(row_indices, length_bin_indices, values, value)

    def _value_matrix(
        self,
        row_indices: np.ndarray,
        length_bin_indices: np.ndarray,
        value: str,
        denominator: str,
    ) -> np.ndarray:
        """
        Convert raw counts to counts, within-row fractions, or densities.

        Density is the selected fraction divided by the length-bin width. Rows
        with zero denominator counts get `NaN` fractions and densities rather
        than silently reporting zero.
        """
        counts = self.lengths.counts[row_indices, :].astype(float, copy=False)
        selected_counts = counts[:, length_bin_indices].copy()
        if value == "count":
            return selected_counts

        row_totals = (
            selected_counts.sum(axis=1)
            if denominator == "selected_bins"
            else counts.sum(axis=1)
        )
        fractions = np.full(selected_counts.shape, np.nan, dtype=float)
        positive_rows = row_totals > 0
        fractions[positive_rows, :] = (
            selected_counts[positive_rows, :] / row_totals[positive_rows, None]
        )
        if value == "fraction":
            return fractions

        widths = (
            self.lengths.length_end_bp[length_bin_indices]
            - self.lengths.length_start_bp[length_bin_indices]
        )
        return fractions / widths

    def _wide_data_frame(
        self,
        row_indices: np.ndarray,
        length_bin_indices: np.ndarray,
        values: np.ndarray,
        value: str,
    ) -> pd.DataFrame:
        """
        Join selected row metadata to one value column per length bin.
        """
        metadata = self.lengths.row_metadata.iloc[row_indices].reset_index(drop=True)
        value_data_frame = pd.DataFrame(
            values,
            columns=_value_column_names(
                [self.lengths.count_columns[index] for index in length_bin_indices],
                value,
            ),
        )
        return pd.concat([metadata, value_data_frame], axis=1)

    def _long_data_frame(
        self,
        row_indices: np.ndarray,
        length_bin_indices: np.ndarray,
        values: np.ndarray,
        value: str,
    ) -> pd.DataFrame:
        """
        Repeat selected row metadata across bins and attach one value column.
        """
        num_rows = len(row_indices)
        num_bins = len(length_bin_indices)
        metadata = self.lengths.row_metadata.iloc[row_indices].reset_index(drop=True)
        metadata = metadata.loc[metadata.index.repeat(num_bins)].reset_index(drop=True)
        # Preserve column order and dtypes after blacklist filtering removes all rows
        selected_bins = self.length_bins().iloc[length_bin_indices].reset_index(
            drop=True
        )
        bins = selected_bins.iloc[[]].copy() if num_rows == 0 else pd.concat(
            [selected_bins] * num_rows, ignore_index=True
        )
        data_frame = pd.concat([metadata, bins], axis=1)
        data_frame[value] = values.ravel()
        return data_frame

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
            length_start_bp=self.lengths.length_start_bp,
            length_end_bp=self.lengths.length_end_bp,
            selector_context="length-count",
        )


class GlobalLengthCounts(LengthCounts):
    """Length counts for global output."""

    def data_frame(
        self,
        *,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
        value: str = "count",
        denominator: str = "all_bins",
        keep_wide: bool = False,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame for the global fragment length distribution.

        Long output has one row per length bin with bin metadata. Wide output
        has one row with one value column per length bin.

        Parameters
        ----------
        with_lengths
            Fragment length or lengths in bp. Returned values use the length
            bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`. Returned
            values use whole length bins that overlap this range.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.
        value
            One of `"count"`, `"fraction"`, or `"density"`. Fractions are
            within the global row. Densities are fractions divided by the
            length-bin width.
        denominator
            For `"fraction"` and `"density"`, `"all_bins"` divides by the row
            total over all length bins, while `"selected_bins"` divides by the
            total over the returned length bins. Ignored for `"count"`.
        keep_wide
            If `False`, return one row per length bin. If `True`, return one
            row with one value column per length bin.

        Returns
        -------
        pandas.DataFrame
            Global length-count values with length-bin metadata for long output
            or value-prefixed columns for wide output.
        """
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self._data_frame_for_rows(
            np.arange(self.lengths.counts.shape[0]),
            value=value,
            denominator=denominator,
            keep_wide=keep_wide,
            max_blacklisted_fraction=1.0,
            length_bin_indices=length_bin_indices,
        )


class WindowedLengthCounts(LengthCounts):
    """Length counts for fixed-size or BED-window output."""

    def window_metadata(self) -> pd.DataFrame:
        """
        Return genomic window metadata for this length-count output.

        Returns
        -------
        pandas.DataFrame
            Columns are `window_idx`, `chrom`, `start`, `end`, and optionally
            `blacklisted_fraction`.
        """
        return self.lengths.row_metadata.copy()

    def counts_array(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
    ) -> np.ndarray:
        """
        Return raw length counts as a dense NumPy array.

        Scalar selectors keep their axis as length one, so the shape is always
        `(selected windows, length_bin)`.

        Parameters
        ----------
        window_idxs
            `None` for all windows, one window index, or a sequence of window
            indices.
        with_lengths
            Fragment length or lengths in bp. Counts are returned for the
            length bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(window, length_bin)`.
        """
        row_indices = normalize_zero_based_indices(
            window_idxs,
            size=self.lengths.counts.shape[0],
            name="window_idxs",
            index_name="window_idx",
        )
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self.lengths.counts[np.ix_(row_indices, length_bin_indices)].copy()

    def data_frame(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
        value: str = "count",
        denominator: str = "all_bins",
        keep_wide: bool = False,
        max_blacklisted_fraction: float = 1.0,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of fragment length distributions for windows.

        Use `window_idxs` to keep only selected genomic windows. Long output has
        one row per selected window and length bin. Wide output has one row per
        selected window with one value column per length bin.

        Parameters
        ----------
        window_idxs
            `None` for all windows, a window index, or a sequence of window indices.
        with_lengths
            Fragment length or lengths in bp. Returned values use the length
            bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`. Returned
            values use whole length bins that overlap this range.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.
        value
            One of `"count"`, `"fraction"`, or `"density"`. Fractions are
            within each selected window. Densities are fractions divided by the
            length-bin width.
        denominator
            For `"fraction"` and `"density"`, `"all_bins"` divides by each
            row's total over all length bins, while `"selected_bins"` divides
            by the total over the returned length bins. Ignored for `"count"`.
        keep_wide
            If `False`, return one row per selected window and length bin. If
            `True`, return one row per selected window with one value column per
            length bin.
        max_blacklisted_fraction
            Maximum `blacklisted_fraction` in 0..1 to keep. The default `1.0`
            keeps all selected windows.

        Returns
        -------
        pandas.DataFrame
            Window metadata and length-count values.
        """
        row_indices = normalize_zero_based_indices(
            window_idxs,
            size=self.lengths.counts.shape[0],
            name="window_idxs",
            index_name="window_idx",
        )
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self._data_frame_for_rows(
            row_indices,
            value=value,
            denominator=denominator,
            keep_wide=keep_wide,
            max_blacklisted_fraction=max_blacklisted_fraction,
            length_bin_indices=length_bin_indices,
        )


class GroupedLengthCounts(LengthCounts):
    """Length counts for grouped BED output."""

    def group_metadata(self) -> pd.DataFrame:
        """
        Return grouped BED metadata for this length-count output.

        Returns
        -------
        pandas.DataFrame
            Columns are `group_idx`, `group_name`, `eligible_windows`, and
            optionally `blacklisted_fraction`.
        """
        return self.lengths.row_metadata.copy()

    def group_idx(self, group_name: str) -> int:
        """
        Find the count-row index for a group name.

        Parameters
        ----------
        group_name
            Group name to resolve.

        Returns
        -------
        int
            Group index.
        """
        group_names = self.lengths.row_metadata["group_name"].to_numpy()
        return resolve_unique_match(
            group_names == group_name,
            missing_message=f"Unknown length-count group name: {group_name!r}",
            duplicate_message=(
                f"Length-count group name is not unique: {group_name!r}"
            ),
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
        Return raw length counts as a dense NumPy array.

        Scalar selectors keep their axis as length one, so the shape is always
        `(selected groups, length_bin)`.

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
            Two bp bounds defining a half-open range `[start, end)`.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(group, length_bin)`.
        """
        row_indices = self._resolve_group_selector(groups, group_idxs)
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self.lengths.counts[np.ix_(row_indices, length_bin_indices)].copy()

    def data_frame(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        with_lengths: int | Sequence[int] | None = None,
        with_length_range: Sequence[int] | None = None,
        length_bin_idxs: int | Sequence[int] | None = None,
        value: str = "count",
        denominator: str = "all_bins",
        keep_wide: bool = False,
        max_blacklisted_fraction: float = 1.0,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of fragment length distributions for groups.

        Use `groups` or `group_idxs` to keep only selected grouped BED rows.
        Long output has one row per selected group and length bin. Wide output
        has one row per selected group with one value column per length bin.

        Parameters
        ----------
        groups
            `None` for all groups, one group name, or a sequence of group names.
            Use either `groups` or `group_idxs`, not both.
        group_idxs
            `None` for all groups, one group index, or a sequence of group
            indices. Use either `groups` or `group_idxs`, not both.
        with_lengths
            Fragment length or lengths in bp. Returned values use the length
            bins containing these lengths. Multiple lengths must select
            distinct length bins.
        with_length_range
            Two bp bounds defining a half-open range `[start, end)`. Returned
            values use whole length bins that overlap this range.
        length_bin_idxs
            `None` for all length bins, one length-bin index, or a sequence of
            length-bin indices. Use only one of `with_lengths`,
            `with_length_range`, or `length_bin_idxs`.
        value
            One of `"count"`, `"fraction"`, or `"density"`. Fractions are
            within each selected group. Densities are fractions divided by the
            length-bin width.
        denominator
            For `"fraction"` and `"density"`, `"all_bins"` divides by each
            row's total over all length bins, while `"selected_bins"` divides
            by the total over the returned length bins. Ignored for `"count"`.
        keep_wide
            If `False`, return one row per selected group and length bin. If
            `True`, return one row per selected group with one value column per
            length bin.
        max_blacklisted_fraction
            Maximum `blacklisted_fraction` in 0..1 to keep. The default `1.0`
            keeps all selected groups.

        Returns
        -------
        pandas.DataFrame
            Group metadata and length-count values.
        """
        row_indices = self._resolve_group_selector(groups, group_idxs)
        length_bin_indices = self._resolve_length_bin_selector(
            with_lengths, with_length_range, length_bin_idxs
        )
        return self._data_frame_for_rows(
            row_indices,
            value=value,
            denominator=denominator,
            keep_wide=keep_wide,
            max_blacklisted_fraction=max_blacklisted_fraction,
            length_bin_indices=length_bin_indices,
        )

    def _resolve_group_selector(
        self,
        groups: str | Sequence[str] | None,
        group_idxs: int | Sequence[int] | None,
    ) -> np.ndarray:
        """
        Normalize group-name or group-index selectors to row indices.
        """
        if groups is not None and group_idxs is not None:
            raise ValueError("Use either groups or group_idxs, not both")
        if groups is None and group_idxs is None:
            return np.arange(self.lengths.counts.shape[0])
        if groups is not None:
            group_names = normalize_strings(groups, name="groups")
            return np.asarray(
                [self.group_idx(group_name) for group_name in group_names],
                dtype=np.int64,
            )
        return normalize_zero_based_indices(
            group_idxs,
            size=self.lengths.counts.shape[0],
            name="group_idxs",
            index_name="group_idx",
        )


def read_lengths(
    path: pathlib.Path | str,
) -> GlobalLengthCounts | WindowedLengthCounts | GroupedLengthCounts:
    """
    Read a cfDNAlab length-count TSV and return the matching loader class.

    Parameters
    ----------
    path
        Path to a `.length_counts.tsv` or `.length_counts.tsv.zst` file.

    Returns
    -------
    LengthCounts
        `GlobalLengthCounts`, `WindowedLengthCounts`, or
        `GroupedLengthCounts`, depending on the TSV metadata columns.
    """
    path = _normalize_length_counts_path(path)
    loaded = LengthCounts._load_tsv(path)
    constructors = {
        "global": GlobalLengthCounts,
        "windowed": WindowedLengthCounts,
        "grouped": GroupedLengthCounts,
    }
    return constructors[loaded.mode](path, loaded)


def _normalize_length_counts_path(path: pathlib.Path | str) -> pathlib.Path:
    """
    Normalize user-supplied length-count paths before file validation.
    """
    if not isinstance(path, (str, os.PathLike)):
        raise TypeError(
            f"Length-count path must be a path-like string, got {type(path).__name__}"
        )
    return pathlib.Path(path)


def _validate_length_counts_path(path: pathlib.Path) -> None:
    """
    Validate that a path points to a supported length-count TSV file.
    """
    if not path.exists():
        raise FileNotFoundError(f"Length-count TSV does not exist: {path}")
    if path.is_dir():
        raise IsADirectoryError(f"Length-count path exists but is a directory: {path}")
    path_text = str(path).lower()
    if not (path_text.endswith(".tsv") or path_text.endswith(".tsv.zst")):
        raise ValueError(
            f"Length-count path must end in '.tsv' or '.tsv.zst', got: {path}"
        )


def _read_length_counts_table(path: pathlib.Path) -> pd.DataFrame:
    """
    Read a plain or Zstandard-compressed length-count TSV into pandas.

    pandas is the primary reader. The command-line fallback is only used when
    pandas cannot read a `.zst` file in the active environment.
    """
    try:
        # Pandas uses zstandard for .zst input when compression inference works.
        table = pd.read_csv(path, sep="\t", compression="infer", keep_default_na=False)
        return _normalize_length_count_table(table)
    except Exception as error:
        if not str(path).lower().endswith(".zst"):
            raise ValueError(f"Could not read length-count TSV at {path}") from error

        # Keep a CLI fallback for older or externally constrained pandas setups
        zstd = shutil.which("zstd")
        if not zstd:
            raise ValueError(
                "Could not read .tsv.zst length-count file with pandas, and the "
                "zstd command-line tool was not found"
            ) from error

        try:
            result = subprocess.run(
                [zstd, "-dc", str(path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            table = pd.read_csv(
                io.BytesIO(result.stdout),
                sep="\t",
                compression=None,
                keep_default_na=False,
            )
            return _normalize_length_count_table(table)
        except subprocess.CalledProcessError as subprocess_error:
            stderr = subprocess_error.stderr.decode("utf-8", errors="replace").strip()
            raise ValueError(
                f"Could not decompress length-count TSV with zstd: {stderr}"
            ) from subprocess_error
        except Exception as fallback_error:
            raise ValueError(
                f"Could not read length-count TSV at {path}"
            ) from fallback_error


def _read_length_count_header(path: pathlib.Path) -> list[str]:
    """
    Read TSV header columns before pandas can mangle duplicate names.
    """
    try:
        if str(path).lower().endswith(".zst"):
            with zstandard.open(path, mode="rt", encoding="utf-8", newline="") as file:
                header = file.readline()
        else:
            with path.open("r", encoding="utf-8", newline="") as file:
                header = file.readline()
    except Exception as error:
        raise ValueError(f"Could not read length-count TSV header at {path}") from error

    if header == "":
        raise ValueError("Length-count TSV must contain a header row")
    return header.rstrip("\r\n").split("\t")


def _normalize_length_count_table(table: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize inferred pandas column types that are public metadata.
    """
    if "chrom" in table.columns:
        table["chrom"] = table["chrom"].astype(str)
    return table


def _length_count_columns(column_names: list[str]) -> list[str]:
    """
    Return count columns and require them to be the final contiguous block.

    This keeps the parser unambiguous. All leading columns are metadata, and all
    trailing columns are length-bin counts.
    """
    count_columns = [name for name in column_names if COUNT_COLUMN_PATTERN.match(name)]
    if not count_columns:
        raise ValueError(
            "Length-count TSV must contain count columns named "
            "count_<length> or count_<start>_<end>"
        )
    if len(set(column_names)) != len(column_names):
        raise ValueError("Length-count TSV column names must be unique")
    first_count_column = column_names.index(count_columns[0])
    expected_count_columns = column_names[first_count_column:]
    if count_columns != expected_count_columns:
        raise ValueError(
            "Length-count TSV count columns must be contiguous and follow metadata columns"
        )
    return count_columns


def _infer_length_count_mode(column_names: list[str], count_columns: list[str]) -> str:
    """
    Infer global, windowed, or grouped output from leading metadata columns.
    """
    first_count_column = column_names.index(count_columns[0])
    metadata_columns = column_names[:first_count_column]
    if metadata_columns == []:
        return "global"
    if metadata_columns in (
        ["chrom", "start", "end"],
        ["chrom", "start", "end", "blacklisted_fraction"],
    ):
        return "windowed"
    if metadata_columns in (
        ["group_name", "eligible_windows"],
        ["group_name", "eligible_windows", "blacklisted_fraction"],
    ):
        return "grouped"
    raise ValueError(
        "Could not infer length-count output mode from metadata columns: "
        + ", ".join(metadata_columns)
    )


def _parse_length_count_columns(count_columns: list[str]) -> dict[str, np.ndarray]:
    """
    Parse `count_*` column names into half-open length-bin metadata.

    `count_167` means `[167, 168)`. `count_150_200` means `[150, 200)`.
    Duplicate intervals are rejected because lookups by fragment length would
    become ambiguous.
    """
    starts: list[int] = []
    ends: list[int] = []
    for count_column in count_columns:
        match = COUNT_COLUMN_PATTERN.match(count_column)
        if match is None:
            raise ValueError(f"Invalid length-count column name: {count_column}")
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) is not None else start + 1
        starts.append(start)
        ends.append(end)

    starts_array = np.asarray(starts, dtype=np.int64)
    ends_array = np.asarray(ends, dtype=np.int64)
    _validate_half_open_intervals(starts_array, ends_array, "length bin")
    intervals = list(zip(starts, ends))
    if len(set(intervals)) != len(intervals):
        raise ValueError("Length-count TSV contains duplicate length bins")

    return {
        "length_bin": np.arange(len(count_columns), dtype=np.int32),
        "length_start_bp": starts_array,
        "length_end_bp": ends_array,
    }


def _length_row_metadata(
    table: pd.DataFrame, mode: str, first_count_column: str
) -> pd.DataFrame:
    """
    Build public row metadata for the inferred length-count mode.

    The Rust TSV does not include explicit row indices. Python exposes
    row indices derived from file order.
    """
    metadata_columns = table.columns[: table.columns.get_loc(first_count_column)]
    metadata = table.loc[:, metadata_columns]
    if mode == "global":
        return pd.DataFrame(index=np.arange(len(table)))

    if mode == "windowed":
        chrom = _required_string_array(metadata["chrom"], "chrom")
        start = _required_integer_array(metadata["start"], "start", nonnegative=True)
        end = _required_integer_array(metadata["end"], "end", nonnegative=True)
        _validate_half_open_intervals(start, end, "window")
        out = pd.DataFrame(
            {
                "window_idx": np.arange(len(table), dtype=np.int32),
                "chrom": chrom,
                "start": start,
                "end": end,
            }
        )
    else:
        group_name = _required_string_array(metadata["group_name"], "group_name")
        eligible_windows = _required_integer_array(
            metadata["eligible_windows"], "eligible_windows", nonnegative=True
        )
        out = pd.DataFrame(
            {
                "group_idx": np.arange(len(table), dtype=np.int32),
                "group_name": group_name,
                "eligible_windows": eligible_windows,
            }
        )

    if "blacklisted_fraction" in metadata.columns:
        out["blacklisted_fraction"] = _required_fraction_array(
            metadata["blacklisted_fraction"], "blacklisted_fraction"
        )
    return out


def _length_counts_matrix_from_table(
    table: pd.DataFrame, count_columns: list[str]
) -> np.ndarray:
    """
    Extract non-negative length counts from count columns as a numeric matrix.
    """
    values = []
    for count_column in count_columns:
        column = _required_numeric_array(table[count_column], count_column)
        if np.any(column < 0):
            raise ValueError(f"{count_column} must contain non-negative values")
        values.append(column)
    return np.column_stack(values).astype(float, copy=False)


def _required_numeric_array(values: pd.Series, name: str) -> np.ndarray:
    """
    Convert a TSV column to finite numeric values or fail with column context.
    """
    try:
        array = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    except Exception as error:
        raise ValueError(f"{name} must contain numeric values") from error
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array


def _required_integer_array(
    values: pd.Series, name: str, *, nonnegative: bool = False
) -> np.ndarray:
    """
    Convert a TSV column to integer values and optionally require non-negativity.
    """
    array = _required_numeric_array(values, name)
    if np.any(array != np.floor(array)):
        raise ValueError(f"{name} must contain integer values")
    if nonnegative and np.any(array < 0):
        raise ValueError(f"{name} must contain non-negative integer values")
    return array.astype(np.int64)


def _required_fraction_array(values: pd.Series, name: str) -> np.ndarray:
    """
    Convert a TSV column to finite fractions in the closed interval 0..1.
    """
    array = _required_numeric_array(values, name)
    if np.any((array < 0) | (array > 1)):
        raise ValueError(f"{name} must contain finite fractions in 0..1")
    return array


def _required_string_array(values: pd.Series, name: str) -> np.ndarray:
    """
    Convert a TSV column to non-empty strings.
    """
    array = values.to_numpy(dtype=object)
    if any(not isinstance(value, str) or value == "" for value in array):
        raise ValueError(f"{name} must contain non-empty character strings")
    return array


def _validate_half_open_intervals(
    starts: np.ndarray, ends: np.ndarray, label: str
) -> None:
    """
    Require each half-open interval to have positive width.
    """
    if np.any(starts >= ends):
        raise ValueError(f"{label} start must be smaller than {label} end")


def _validate_value(value: str) -> str:
    """
    Validate the requested length-count value scale.
    """
    if value not in VALID_LENGTH_VALUES:
        raise ValueError(
            "value must be one of "
            + ", ".join(repr(valid_value) for valid_value in VALID_LENGTH_VALUES)
        )
    return value


def _validate_denominator(denominator: str) -> str:
    """
    Validate the row-total basis used for fractions and densities.
    """
    if denominator not in VALID_DENOMINATORS:
        valid = ", ".join(repr(value) for value in sorted(VALID_DENOMINATORS))
        raise ValueError(f"denominator must be one of {valid}")
    return denominator


def _value_column_names(count_columns: list[str], value: str) -> list[str]:
    """
    Build wide-output column names for counts, fractions, or densities.
    """
    if value == "count":
        return count_columns
    return [column.replace("count_", f"{value}_", 1) for column in count_columns]
