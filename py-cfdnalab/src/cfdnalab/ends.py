"""
Load cfDNAlab end-motif Zarr outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
import numbers
import pathlib
from typing import Any, List, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import zarr

from ._helpers import (
    filter_blacklisted_fraction,
    normalize_strings,
    normalize_zero_based_indices,
    resolve_unique_match,
    validate_scalar_bool,
    validate_zero_based_index,
)

END_MOTIF_MIN_SUPPORTED_SCHEMA_VERSION = 1
END_MOTIF_MAX_SUPPORTED_SCHEMA_VERSION = 1
VALID_STORAGE_MODES = {"dense", "sparse_coo"}
VALID_ROW_MODES = {"global", "size", "bed", "grouped_bed"}


@dataclass
class LoadedEndMotifs:
    """Validated end-motif Zarr handles and row or motif metadata."""

    store: Any
    storage_mode: str
    row_mode: str
    motif_index: np.ndarray
    motif_names: np.ndarray
    row: np.ndarray
    counts: zarr.Array | None
    sparse_row: np.ndarray | None
    sparse_motif: np.ndarray | None
    sparse_count: np.ndarray | None
    sparse_shape: np.ndarray | None
    row_labels: np.ndarray | None
    chromosome: np.ndarray | None
    chromosome_names: np.ndarray | None
    row_chromosome: np.ndarray | None
    row_start_bp: np.ndarray | None
    row_end_bp: np.ndarray | None
    blacklisted_fraction: np.ndarray | None
    group_idx: np.ndarray | None
    group_names: np.ndarray | None
    eligible_windows: np.ndarray | None


class EndMotifCounts:
    """Common API for global, windowed, and grouped end-motif outputs."""

    def __init__(
        self,
        path: pathlib.Path | str,
        loaded_end_motifs: LoadedEndMotifs | None = None,
    ) -> None:
        """
        Load an end-motif count Zarr store.

        Parameters
        ----------
        path
            Path to a `<prefix>.end_motifs.zarr` directory.
        loaded_end_motifs
            Preloaded store data used by `read_end_motifs`.
        """
        self.path = pathlib.Path(path)
        if loaded_end_motifs is None:
            loaded_end_motifs = EndMotifCounts._load_zarr(self.path)
        self.end_motifs = loaded_end_motifs

    def __repr__(self) -> str:
        """
        Return a compact summary with path, schema version, modes, and shape.
        """
        schema_version = self.end_motifs.store.attrs.get("cfdnalab_schema_version")
        if self.end_motifs.storage_mode == "dense":
            shape = tuple(_required(self.end_motifs.counts, "counts").shape)
        else:
            shape = tuple(
                _required(self.end_motifs.sparse_shape, "sparse/shape").astype(int)
            )
        return (
            f"{self.__class__.__name__}("
            f"path={str(self.path)!r}, "
            f"schema_version={schema_version!r}, "
            f"storage_mode={self.end_motifs.storage_mode!r}, "
            f"row_mode={self.end_motifs.row_mode!r}, "
            f"shape={shape!r}"
            ")"
        )

    @staticmethod
    def _load_zarr(path: pathlib.Path | str) -> LoadedEndMotifs:
        """
        Open and validate an end-motif count Zarr store.

        Parameters
        ----------
        path
            Path to a `<prefix>.end_motifs.zarr` directory.

        Returns
        -------
        LoadedEndMotifs
            Loaded Zarr handles and axis metadata.
        """
        path = pathlib.Path(path)
        _validate_zarr_store_path(path)

        try:
            store = zarr.open_group(str(path), mode="r", zarr_format=3)
        except Exception as error:
            raise ValueError(
                f"Could not open end-motif Zarr store at {path}"
            ) from error

        storage_mode, row_mode = _validate_root_metadata(store)
        _validate_required_arrays(store, storage_mode, row_mode)

        # Motif and row axes are small coordinate arrays, so keep them in memory
        motif_index = _read_array(store, "motif_index")
        motif_names = _read_motif_ascii_labels(store, len(motif_index))
        row = _read_array(store, "row")
        _validate_axis(motif_index, "motif_index")
        _validate_axis(row, "row")

        counts = None
        sparse_row = None
        sparse_motif = None
        sparse_count = None
        sparse_shape = None
        if storage_mode == "dense":
            # Dense counts stay as a Zarr array handle until a method asks for data
            counts = store["counts"]
            _validate_dense_counts(counts, len(row), len(motif_index))
        else:
            # Sparse coordinates are the count data, so validate them eagerly
            sparse_row = _read_array(store, "sparse/row")
            sparse_motif = _read_array(store, "sparse/motif")
            sparse_count = _read_array(store, "sparse/count")
            sparse_shape = _read_array(store, "sparse/shape")
            _validate_sparse_arrays(
                store,
                sparse_row,
                sparse_motif,
                sparse_count,
                sparse_shape,
                len(row),
                len(motif_index),
            )

        row_labels = None
        chromosome = None
        chromosome_names = None
        row_chromosome = None
        row_start_bp = None
        row_end_bp = None
        blacklisted_fraction = None
        group_idx = None
        group_names = None
        eligible_windows = None

        # Row-mode metadata determines which subclass read_end_motifs returns
        if row_mode == "global":
            row_labels = _read_labels(store["row"], "row_label", len(row), "row")
        elif row_mode in {"size", "bed"}:
            chromosome = _read_array(store, "chromosome")
            chromosome_names = _read_labels(
                store["chromosome"],
                "chromosome_name",
                len(chromosome),
                "chromosome",
            )
            row_chromosome = _read_array(store, "row_chromosome")
            row_start_bp = _read_array(store, "row_start_bp")
            row_end_bp = _read_array(store, "row_end_bp")
            blacklisted_fraction = _read_array(store, "blacklisted_fraction")
            _validate_axis(chromosome, "chromosome")
            _validate_same_length(row_chromosome, row, "row_chromosome", "row")
            _validate_same_length(row_start_bp, row, "row_start_bp", "row")
            _validate_same_length(row_end_bp, row, "row_end_bp", "row")
            _validate_same_length(
                blacklisted_fraction, row, "blacklisted_fraction", "row"
            )
            if np.any(row_chromosome < 0) or np.any(row_chromosome >= len(chromosome)):
                raise ValueError(
                    "row_chromosome contains an index outside the chromosome axis"
                )
        elif row_mode == "grouped_bed":
            group_idx = _read_array(store, "group")
            group_names = _read_labels(
                store["group"], "group_name", len(group_idx), "group"
            )
            eligible_windows = _read_array(store, "eligible_windows")
            blacklisted_fraction = _read_array(store, "blacklisted_fraction")
            _validate_same_length(group_idx, row, "group", "row")
            _validate_same_length(group_names, row, "group_name labels", "row")
            _validate_same_length(eligible_windows, row, "eligible_windows", "row")
            _validate_same_length(
                blacklisted_fraction, row, "blacklisted_fraction", "row"
            )
            _validate_axis(group_idx, "group")

        return LoadedEndMotifs(
            store=store,
            storage_mode=storage_mode,
            row_mode=row_mode,
            motif_index=motif_index,
            motif_names=motif_names,
            row=row,
            counts=counts,
            sparse_row=sparse_row,
            sparse_motif=sparse_motif,
            sparse_count=sparse_count,
            sparse_shape=sparse_shape,
            row_labels=row_labels,
            chromosome=chromosome,
            chromosome_names=chromosome_names,
            row_chromosome=row_chromosome,
            row_start_bp=row_start_bp,
            row_end_bp=row_end_bp,
            blacklisted_fraction=blacklisted_fraction,
            group_idx=group_idx,
            group_names=group_names,
            eligible_windows=eligible_windows,
        )

    def storage_mode(self) -> str:
        """
        Return how end-motif counts are stored on disk.

        Returns
        -------
        str
            Either `"dense"` or `"sparse_coo"`.
        """
        return self.end_motifs.storage_mode

    def row_mode(self) -> str:
        """
        Return what each end-motif count row represents.

        Returns
        -------
        str
            One of `"global"`, `"size"`, `"bed"`, or `"grouped_bed"`.
        """
        return self.end_motifs.row_mode

    def motifs(self) -> List[str]:
        """
        Return end-motif labels in motif-axis order.

        Returns
        -------
        list[str]
            Motif labels in motif-axis order.
        """
        return self.end_motifs.motif_names.tolist()

    def motif_idx(self, motif: str) -> int:
        """
        Find the motif-axis index for a motif label.

        Parameters
        ----------
        motif
            Motif label to resolve.

        Returns
        -------
        int
            Motif index.
        """
        return self._resolve_motif(motif)

    def has_motif(self, motif: str) -> bool:
        """
        Return whether a motif label exists in this output.

        Sparse output only stores observed motifs, so an unobserved motif will
        return `False` even if it is part of the theoretical motif universe.

        Parameters
        ----------
        motif
            Motif label to check.

        Returns
        -------
        bool
            Whether the motif can be resolved in this output.
        """
        return bool(np.any(self.end_motifs.motif_names == motif))

    def motif_metadata(self) -> pd.DataFrame:
        """
        Get the motif labels and motif indices available in this output.

        Returns
        -------
        pandas.DataFrame
            Columns are `motif_index` and `motif`.
        """
        return pd.DataFrame(
            {
                "motif_index": self.end_motifs.motif_index,
                "motif": self.end_motifs.motif_names,
            }
        )

    def sparse_coo(self) -> sparse.coo_matrix:
        """
        Return end-motif counts as a SciPy COO matrix.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(output row, motif)`. Output rows are
            windows for windowed output, groups for grouped output, and the
            single global summary row for global output.
        """
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return sparse.coo_matrix(np.asarray(counts[:]))

        row_index = _required(self.end_motifs.sparse_row, "sparse/row")
        motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
        count = _required(self.end_motifs.sparse_count, "sparse/count")
        shape = tuple(
            _required(self.end_motifs.sparse_shape, "sparse/shape").astype(int)
        )
        return sparse.coo_matrix(
            (
                count,
                (
                    row_index.astype(np.int64, copy=False),
                    motif_index.astype(np.int64, copy=False),
                ),
            ),
            shape=shape,
        )

    def _data_frame(
        self,
        *,
        densify: bool = False,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
    ) -> pd.DataFrame:
        """
        Shared implementation behind mode-specific public data_frame() methods.
        """
        densify = validate_scalar_bool(densify, "densify")
        row_indices = self._resolve_row_selector(window_idxs, groups, group_idxs)
        motif_indices = self._resolve_motif_selector(motifs, motif_idxs)
        row_indices = filter_blacklisted_fraction(
            self._row_metadata_data_frame(),
            max_blacklisted_fraction,
            row_indices=row_indices,
        )
        if self.end_motifs.storage_mode == "sparse_coo" and not densify:
            return self._stored_data_frame_for_indices(row_indices, motif_indices)
        return self._complete_data_frame_for_indices(row_indices, motif_indices, densify)

    def sparse_coo_for_motif(self, motif: str) -> sparse.coo_matrix:
        """
        Return sparse counts for one motif across output rows.

        Parameters
        ----------
        motif
            Motif label to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(output row, 1)`.
        """
        return self.sparse_coo_for_motif_idx(self._resolve_motif(motif))

    def sparse_coo_for_motif_idx(self, motif_idx: int) -> sparse.coo_matrix:
        """
        Return sparse counts for one motif index across output rows.

        Parameters
        ----------
        motif_idx
            Motif index to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(output row, 1)`.
        """
        motif_idx = self._validate_motif_idx(motif_idx)
        if self.end_motifs.storage_mode == "sparse_coo":
            sparse_row_index = _required(self.end_motifs.sparse_row, "sparse/row")
            sparse_motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
            sparse_count = _required(self.end_motifs.sparse_count, "sparse/count")
            matches = sparse_motif_index == motif_idx
            row_index = sparse_row_index[matches].astype(np.int64, copy=False)
            col_index = np.zeros(len(row_index), dtype=np.int64)
            count = sparse_count[matches]
            return sparse.coo_matrix(
                (count, (row_index, col_index)),
                shape=(len(self.end_motifs.row), 1),
            )
        return self.sparse_coo().tocsr()[:, motif_idx].tocoo()

    def dense_counts_zarr_array(self) -> zarr.Array:
        """
        Return the lazy Zarr counts array for dense output.

        This returns the on-disk Zarr array handle without loading the full
        dense matrix into memory. Sparse output has no dense `counts` array.

        Returns
        -------
        zarr.Array
            Dense count array with shape `(output row, motif)`.
        """
        if self.end_motifs.storage_mode != "dense":
            raise ValueError(
                "dense_counts_zarr_array() is only available for dense end-motif output"
            )
        return _required(self.end_motifs.counts, "counts")

    def dense_counts_matrix(self, allow_densify: bool = False) -> np.ndarray:
        """
        Return end-motif counts as a dense NumPy matrix.

        Sparse stores are only densified when `allow_densify=True`.

        Parameters
        ----------
        allow_densify
            If `True`, allow sparse stores to be converted to a dense in-memory
            matrix.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(output row, motif)`.
        """
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[:])
        _require_densify(allow_densify, "dense_counts_matrix")
        return self.sparse_coo().toarray()

    def dense_counts_for_motif(
        self, motif: str, allow_densify: bool = False
    ) -> np.ndarray:
        """
        Return dense counts for one motif across output rows.

        Parameters
        ----------
        motif
            Motif label to extract.
        allow_densify
            If `True`, allow sparse stores to be converted to a dense vector.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(output row,)`.
        """
        return self.dense_counts_for_motif_idx(
            self._resolve_motif(motif), allow_densify=allow_densify
        )

    def dense_counts_for_motif_idx(
        self, motif_idx: int, allow_densify: bool = False
    ) -> np.ndarray:
        """
        Return dense counts for one motif index across output rows.

        Sparse stores are only densified when `allow_densify=True`.

        Parameters
        ----------
        motif_idx
            Motif index to extract.
        allow_densify
            If `True`, allow sparse stores to be converted to a dense vector.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(output row,)`.
        """
        motif_idx = self._validate_motif_idx(motif_idx)
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[:, motif_idx])

        _require_densify(allow_densify, "dense_counts_for_motif_idx")
        values = np.zeros(len(self.end_motifs.row), dtype=float)
        sparse_row_index = _required(self.end_motifs.sparse_row, "sparse/row")
        sparse_motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
        sparse_count = _required(self.end_motifs.sparse_count, "sparse/count")
        matches = sparse_motif_index == motif_idx
        values[sparse_row_index[matches].astype(int)] = sparse_count[matches]
        return values

    def _complete_data_frame_for_indices(
        self, row_indices: np.ndarray, motif_indices: np.ndarray, densify: bool
    ) -> pd.DataFrame:
        """
        Build rows for every selected output row and motif.
        """
        row_metadata = self._row_metadata_data_frame().iloc[row_indices].reset_index(
            drop=True
        )
        motif_metadata = self.motif_metadata().iloc[motif_indices].reset_index(
            drop=True
        )
        row_count = len(row_indices)
        motif_count = len(motif_indices)
        if row_count == 0 or motif_count == 0:
            return self._empty_data_frame()

        counts = self.dense_counts_matrix(allow_densify=densify)
        selected_counts = counts[np.ix_(row_indices, motif_indices)]
        repeated_rows = row_metadata.loc[
            row_metadata.index.repeat(motif_count)
        ].reset_index(drop=True)
        repeated_motifs = pd.concat(
            [motif_metadata] * row_count, ignore_index=True
        )
        data_frame = pd.concat([repeated_rows, repeated_motifs], axis=1)
        data_frame["count"] = selected_counts.ravel()
        return data_frame

    def _stored_data_frame_for_indices(
        self, row_indices: np.ndarray, motif_indices: np.ndarray
    ) -> pd.DataFrame:
        """
        Build a selected data frame from stored COO rows.
        """
        if len(row_indices) == 0 or len(motif_indices) == 0:
            return self._empty_data_frame()

        sparse_rows = _required(self.end_motifs.sparse_row, "sparse/row").astype(int)
        sparse_motifs = _required(self.end_motifs.sparse_motif, "sparse/motif").astype(
            int
        )
        sparse_counts = _required(self.end_motifs.sparse_count, "sparse/count")
        row_positions = {int(row_idx): order for order, row_idx in enumerate(row_indices)}
        motif_positions = {
            int(motif_idx): order for order, motif_idx in enumerate(motif_indices)
        }
        matches = np.asarray(
            [
                row_idx in row_positions and motif_idx in motif_positions
                for row_idx, motif_idx in zip(sparse_rows, sparse_motifs)
            ],
            dtype=bool,
        )
        if not np.any(matches):
            return self._empty_data_frame()

        matched_rows = sparse_rows[matches]
        matched_motifs = sparse_motifs[matches]
        sort_order = np.argsort(
            [
                row_positions[int(row_idx)] * len(motif_indices)
                + motif_positions[int(motif_idx)]
                for row_idx, motif_idx in zip(matched_rows, matched_motifs)
            ],
            kind="stable",
        )
        matched_rows = matched_rows[sort_order]
        matched_motifs = matched_motifs[sort_order]
        matched_counts = sparse_counts[matches][sort_order]

        row_metadata = self._row_metadata_data_frame().iloc[matched_rows].reset_index(
            drop=True
        )
        motif_metadata = self.motif_metadata().iloc[matched_motifs].reset_index(
            drop=True
        )
        data_frame = pd.concat([row_metadata, motif_metadata], axis=1)
        data_frame["count"] = matched_counts
        return data_frame

    def _empty_data_frame(self) -> pd.DataFrame:
        """
        Return an empty end-motif data frame with public columns.
        """
        row_metadata = self._row_metadata_data_frame().iloc[[]].reset_index(drop=True)
        motif_metadata = self.motif_metadata().iloc[[]].reset_index(drop=True)
        data_frame = pd.concat([row_metadata, motif_metadata], axis=1)
        data_frame["count"] = np.asarray([], dtype=float)
        return data_frame

    def _row_metadata_data_frame(self) -> pd.DataFrame:
        """
        Build mode-specific row metadata for the active end-motif row mode.

        Global output has one labeled row. Windowed output exposes genomic
        coordinates, and grouped output exposes group names and eligible-window
        counts.
        """
        if self.end_motifs.row_mode == "global":
            data_frame = pd.DataFrame()
            data_frame["row_label"] = self.end_motifs.row_labels
        elif self.end_motifs.row_mode in {"size", "bed"}:
            data_frame = pd.DataFrame()
            row_chromosome = _required(self.end_motifs.row_chromosome, "row_chromosome")
            chromosome_names = _required(
                self.end_motifs.chromosome_names, "chromosome_names"
            )
            data_frame["window_idx"] = self.end_motifs.row
            data_frame["chrom"] = chromosome_names[row_chromosome.astype(int)]
            data_frame["start"] = _required(self.end_motifs.row_start_bp, "row_start_bp")
            data_frame["end"] = _required(self.end_motifs.row_end_bp, "row_end_bp")
            data_frame["blacklisted_fraction"] = _required(
                self.end_motifs.blacklisted_fraction, "blacklisted_fraction"
            )
        elif self.end_motifs.row_mode == "grouped_bed":
            data_frame = pd.DataFrame()
            data_frame["group_idx"] = _required(self.end_motifs.group_idx, "group")
            data_frame["group_name"] = _required(self.end_motifs.group_names, "group_name")
            data_frame["eligible_windows"] = _required(
                self.end_motifs.eligible_windows, "eligible_windows"
            )
            data_frame["blacklisted_fraction"] = _required(
                self.end_motifs.blacklisted_fraction, "blacklisted_fraction"
            )
        else:
            data_frame = pd.DataFrame()
        return data_frame

    def _resolve_row_selector(
        self,
        window_idxs: int | Sequence[int] | None,
        groups: str | Sequence[str] | None,
        group_idxs: int | Sequence[int] | None,
    ) -> np.ndarray:
        """
        Normalize row selectors for the active row mode.
        """
        if self.end_motifs.row_mode in {"size", "bed"}:
            if groups is not None or group_idxs is not None:
                raise ValueError("Grouped selectors can only be used with grouped output")
            return normalize_zero_based_indices(
                window_idxs,
                size=len(self.end_motifs.row),
                name="window_idxs",
                index_name="window_idx",
            )
        if self.end_motifs.row_mode == "grouped_bed":
            if window_idxs is not None:
                raise ValueError("window_idxs can only be used with windowed output")
            if groups is not None and group_idxs is not None:
                raise ValueError("Use either groups or group_idxs, not both")
            if groups is not None:
                group_names = normalize_strings(groups, name="groups")
                return np.asarray(
                    [self.group_idx(group_name) for group_name in group_names],
                    dtype=np.int64,
                )
            return normalize_zero_based_indices(
                group_idxs,
                size=len(self.end_motifs.row),
                name="group_idxs",
                index_name="group_idx",
            )
        if window_idxs is not None or groups is not None or group_idxs is not None:
            raise ValueError("Row selectors cannot be used with global output")
        return np.arange(len(self.end_motifs.row), dtype=np.int64)

    def _resolve_motif_selector(
        self,
        motifs: str | Sequence[str] | None,
        motif_idxs: int | Sequence[int] | None,
    ) -> np.ndarray:
        """
        Normalize motif-name or motif-index selectors.
        """
        if motifs is not None and motif_idxs is not None:
            raise ValueError("Use either motifs or motif_idxs, not both")
        if motifs is None and motif_idxs is None:
            return np.arange(len(self.end_motifs.motif_index), dtype=np.int64)
        if motifs is not None:
            motif_names = normalize_strings(motifs, name="motifs")
            return np.asarray(
                [self._resolve_motif(motif) for motif in motif_names],
                dtype=np.int64,
            )
        return normalize_zero_based_indices(
            motif_idxs,
            size=len(self.end_motifs.motif_index),
            name="motif_idxs",
            index_name="motif_idx",
        )

    def _resolve_motif(self, motif: str) -> int:
        """
        Resolve a motif label to its motif index.
        """
        return resolve_unique_match(
            self.end_motifs.motif_names == motif,
            missing_message=f"Unknown end-motif label: {motif!r}",
            duplicate_message=f"End-motif label is not unique: {motif!r}",
        )

    def _validate_row(self, row: int) -> int:
        """
        Validate and normalize one row index.
        """
        return validate_zero_based_index(row, len(self.end_motifs.row), "row")

    def _validate_motif_idx(self, motif_idx: int) -> int:
        """
        Validate and normalize one motif index.
        """
        return validate_zero_based_index(
            motif_idx, len(self.end_motifs.motif_index), "motif_idx"
        )


class GlobalEndMotifCounts(EndMotifCounts):
    """End-motif counts for global output."""

    def data_frame(
        self,
        *,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame for global end-motif counts.

        Sparse outputs return stored non-zero motif counts unless
        `densify=True`. Densifying adds explicit zero-count rows for selected
        observed motifs. Dense outputs always include zero counts.

        Parameters
        ----------
        densify
            If `True`, sparse outputs add explicit zero-count rows for selected
            observed motifs. Dense outputs ignore this option.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.

        Returns
        -------
        pandas.DataFrame
            Global row metadata, motif metadata, and `count`.
        """
        return self._data_frame(
            densify=densify,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def dense_counts_vec(self, allow_densify: bool = False) -> np.ndarray:
        """
        Return global end-motif counts as a dense vector.

        Sparse stores are only densified when `allow_densify=True`.

        Parameters
        ----------
        allow_densify
            If `True`, allow sparse stores to be converted to dense counts.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(motif,)`.
        """
        return self.dense_counts_matrix(allow_densify=allow_densify)[0, :]


class WindowedEndMotifCounts(EndMotifCounts):
    """End-motif counts for fixed-size or BED-window output."""

    def data_frame(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of end-motif counts for genomic windows.

        Use `window_idxs` to keep only selected windows and `motifs` or
        `motif_idxs` to keep only selected motifs. Sparse outputs return stored
        non-zero rows unless `densify=True`. Densifying adds explicit
        zero-count rows for selected observed motifs. Dense outputs always
        include zero counts.

        Parameters
        ----------
        window_idxs
            `None` for all windows, one window index, or a sequence of window
            indices.
        densify
            If `True`, sparse outputs add explicit zero-count rows for selected
            observed motifs. Dense outputs ignore this option.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.
        max_blacklisted_fraction
            Maximum row `blacklisted_fraction` in 0..1 to retain before counts
            are returned. The default `1.0` keeps all selected windows.

        Returns
        -------
        pandas.DataFrame
            Window metadata, motif metadata, and `count`.
        """
        return self._data_frame(
            densify=densify,
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
        )

    def windows(self) -> pd.DataFrame:
        """
        Get the genomic windows available in this end-motif output.

        Public genomic window metadata uses `window_idx`, `chrom`, `start`,
        and `end` columns.

        Returns
        -------
        pandas.DataFrame
            Columns are `window_idx`, `chrom`, `start`, `end`, and
            `blacklisted_fraction`.
        """
        return self._row_metadata_data_frame()

    def sparse_coo_for_window(self, window_idx: int) -> sparse.coo_matrix:
        """
        Return sparse motif counts for one genomic window.

        Parameters
        ----------
        window_idx
            Window index to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(1, motif)`.
        """
        window_idx = self._validate_row(window_idx)
        if self.end_motifs.storage_mode == "sparse_coo":
            sparse_row_index = _required(self.end_motifs.sparse_row, "sparse/row")
            sparse_motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
            sparse_count = _required(self.end_motifs.sparse_count, "sparse/count")
            matches = sparse_row_index == window_idx
            row_index = np.zeros(int(matches.sum()), dtype=np.int64)
            motif_index = sparse_motif_index[matches].astype(np.int64, copy=False)
            return sparse.coo_matrix(
                (sparse_count[matches], (row_index, motif_index)),
                shape=(1, len(self.end_motifs.motif_index)),
            )
        return self.sparse_coo().tocsr()[window_idx, :].tocoo()

    def dense_counts_for_window(
        self, window_idx: int, allow_densify: bool = False
    ) -> np.ndarray:
        """
        Return motif counts for one genomic window as a dense vector.

        Sparse stores are only densified when `allow_densify=True`.

        Parameters
        ----------
        window_idx
            Window index to extract.
        allow_densify
            If `True`, allow sparse stores to be converted to dense counts.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(motif,)`.
        """
        window_idx = self._validate_row(window_idx)
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[window_idx, :])
        _require_densify(allow_densify, "dense_counts_for_window")
        return self.sparse_coo_for_window(window_idx).toarray()[0, :]


class GroupedEndMotifCounts(EndMotifCounts):
    """End-motif counts for grouped BED output."""

    def data_frame(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of end-motif counts for grouped BED rows.

        Use `groups` or `group_idxs` to keep only selected groups and `motifs`
        or `motif_idxs` to keep only selected motifs. Sparse outputs return
        stored non-zero rows unless `densify=True`. Densifying adds explicit
        zero-count rows for selected observed motifs. Dense outputs always
        include zero counts.

        Parameters
        ----------
        groups
            `None` for all groups, one group name, or a sequence of group names.
            Use either `groups` or `group_idxs`, not both.
        group_idxs
            `None` for all groups, one group index, or a sequence of group
            indices. Use either `groups` or `group_idxs`, not both.
        densify
            If `True`, sparse outputs add explicit zero-count rows for selected
            observed motifs. Dense outputs ignore this option.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.
        max_blacklisted_fraction
            Maximum row `blacklisted_fraction` in 0..1 to retain before counts
            are returned. The default `1.0` keeps all selected groups.

        Returns
        -------
        pandas.DataFrame
            Group metadata, motif metadata, and `count`.
        """
        return self._data_frame(
            densify=densify,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
        )

    def groups(self) -> pd.DataFrame:
        """
        Get the BED groups available in this end-motif output.

        Returns
        -------
        pandas.DataFrame
            Columns are `group_idx`, `group_name`, `eligible_windows`, and
            `blacklisted_fraction`.
        """
        return self._row_metadata_data_frame()

    def group_idx(self, group_name: str) -> int:
        """
        Find the end-motif row index for a group name.

        Parameters
        ----------
        group_name
            Group name to resolve.

        Returns
        -------
        int
            Group index.
        """
        group_names = _required(self.end_motifs.group_names, "group_name")
        return resolve_unique_match(
            group_names == group_name,
            missing_message=f"Unknown end-motif group name: {group_name!r}",
            duplicate_message=f"End-motif group name is not unique: {group_name!r}",
        )

    def sparse_coo_for_group(self, group: int | str) -> sparse.coo_matrix:
        """
        Return sparse motif counts for one grouped BED row.

        Parameters
        ----------
        group
            Group index or group name to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(1, motif)`.
        """
        group_idx = self._resolve_group(group)
        if self.end_motifs.storage_mode == "sparse_coo":
            sparse_row_index = _required(self.end_motifs.sparse_row, "sparse/row")
            sparse_motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
            sparse_count = _required(self.end_motifs.sparse_count, "sparse/count")
            matches = sparse_row_index == group_idx
            row_index = np.zeros(int(matches.sum()), dtype=np.int64)
            motif_index = sparse_motif_index[matches].astype(np.int64, copy=False)
            return sparse.coo_matrix(
                (sparse_count[matches], (row_index, motif_index)),
                shape=(1, len(self.end_motifs.motif_index)),
            )
        return self.sparse_coo().tocsr()[group_idx, :].tocoo()

    def dense_counts_for_group(
        self, group: int | str, allow_densify: bool = False
    ) -> np.ndarray:
        """
        Return motif counts for one grouped BED row as a dense vector.

        Sparse stores are only densified when `allow_densify=True`.

        Parameters
        ----------
        group
            Group index or group name to extract.
        allow_densify
            If `True`, allow sparse stores to be converted to dense counts.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(motif,)`.
        """
        group_idx = self._resolve_group(group)
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[group_idx, :])
        _require_densify(allow_densify, "dense_counts_for_group")
        return self.sparse_coo_for_group(group_idx).toarray()[0, :]

    def _resolve_group(self, group: int | str) -> int:
        """
        Resolve a group name or group index to a row index.
        """
        if isinstance(group, str):
            return self.group_idx(group)
        return validate_zero_based_index(group, len(self.end_motifs.row), "group")


def read_end_motifs(path: pathlib.Path | str) -> EndMotifCounts:
    """
    Open a cfDNAlab end-motif count Zarr store.

    Parameters
    ----------
    path
        Path to an `.end_motifs.zarr` directory.

    Returns
    -------
    EndMotifCounts
        Mode-specific end-motif count helper.
    """
    path = pathlib.Path(path)
    loaded = EndMotifCounts._load_zarr(path)
    if loaded.row_mode == "global":
        return GlobalEndMotifCounts(path, loaded)
    if loaded.row_mode in {"size", "bed"}:
        return WindowedEndMotifCounts(path, loaded)
    if loaded.row_mode == "grouped_bed":
        return GroupedEndMotifCounts(path, loaded)
    raise ValueError(f"Unsupported end-motif row mode: {loaded.row_mode!r}")


def _validate_zarr_store_path(path: pathlib.Path) -> None:
    """
    Validate that a path points to a Zarr V3 end-motif store directory.
    """
    if not path.exists():
        raise FileNotFoundError(f"End-motif Zarr store does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(
            f"End-motif Zarr store path exists but is not a directory: {path}"
        )
    if path.suffix != ".zarr":
        raise ValueError(
            f"End-motif Zarr store path should end with '.zarr', got: {path}"
        )
    if not (path / "zarr.json").is_file():
        raise ValueError(
            f"End-motif Zarr store is missing Zarr V3 metadata file: {path / 'zarr.json'}"
        )


def _validate_root_metadata(store: Any) -> tuple[str, str]:
    """
    Validate root attributes and return storage and row modes.
    """
    schema = store.attrs.get("cfdnalab_schema")
    if schema != "end_motif_counts":
        raise ValueError(
            f"Expected cfdnalab_schema='end_motif_counts', found {schema!r}"
        )

    schema_version = store.attrs.get("cfdnalab_schema_version")
    if not isinstance(schema_version, numbers.Integral) or not (
        END_MOTIF_MIN_SUPPORTED_SCHEMA_VERSION
        <= int(schema_version)
        <= END_MOTIF_MAX_SUPPORTED_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unsupported end-motif schema version: "
            f"{schema_version!r}. Supported range: "
            f"{END_MOTIF_MIN_SUPPORTED_SCHEMA_VERSION}..{END_MOTIF_MAX_SUPPORTED_SCHEMA_VERSION}"
        )

    storage_mode = store.attrs.get("storage_mode")
    if storage_mode not in VALID_STORAGE_MODES:
        raise ValueError(f"Unsupported end-motif storage mode: {storage_mode!r}")

    row_mode = store.attrs.get("row_mode")
    if row_mode not in VALID_ROW_MODES:
        raise ValueError(f"Unsupported end-motif row mode: {row_mode!r}")

    return storage_mode, row_mode


def _validate_required_arrays(store: Any, storage_mode: str, row_mode: str) -> None:
    """
    Require all arrays needed by the active storage mode and row mode.
    """
    required = {"motif_index", "motif_byte", "motif_ascii", "row"}
    if storage_mode == "dense":
        required.add("counts")
    else:
        required.update(
            {
                "sparse/row",
                "sparse/motif",
                "sparse/count",
                "sparse/shape",
                "sparse/sparse_dimension",
            }
        )

    if row_mode in {"size", "bed"}:
        required.update(
            {
                "chromosome",
                "row_chromosome",
                "row_start_bp",
                "row_end_bp",
                "blacklisted_fraction",
            }
        )
    elif row_mode == "grouped_bed":
        required.update({"group", "eligible_windows", "blacklisted_fraction"})

    missing = sorted(name for name in required if not _has_array(store, name))
    if missing:
        raise ValueError(f"End-motif Zarr store is missing arrays: {missing}")


def _has_array(store: Any, name: str) -> bool:
    """
    Return whether a Zarr group contains an array path.
    """
    try:
        store[name]
    except Exception:
        return False
    return True


def _validate_dense_counts(counts: Any, n_rows: int, n_motifs: int) -> None:
    """
    Validate dense end-motif count array dimensions and shape.
    """
    _validate_array_dimensions(counts, ("row", "motif"), "dense counts")
    if tuple(counts.shape) != (n_rows, n_motifs):
        raise ValueError(
            "dense counts shape does not match row and motif axes: "
            f"counts={counts.shape}, coordinates={(n_rows, n_motifs)}"
        )


def _validate_sparse_arrays(
    store: Any,
    row: np.ndarray,
    motif: np.ndarray,
    count: np.ndarray,
    shape: np.ndarray,
    n_rows: int,
    n_motifs: int,
) -> None:
    """
    Validate sparse COO arrays against row and motif axes.

    Sparse output stores row indices, motif indices, and counts as parallel
    arrays. Shape and dimension labels are validated separately so malformed COO
    stores fail before users convert them to SciPy matrices.
    """
    _validate_array_dimensions(store["sparse/row"], ("nnz",), "sparse/row")
    _validate_array_dimensions(store["sparse/motif"], ("nnz",), "sparse/motif")
    _validate_array_dimensions(store["sparse/count"], ("nnz",), "sparse/count")
    _validate_array_dimensions(
        store["sparse/shape"], ("sparse_dimension",), "sparse/shape"
    )
    _validate_array_dimensions(
        store["sparse/sparse_dimension"],
        ("sparse_dimension",),
        "sparse/sparse_dimension",
    )
    _validate_same_length(row, motif, "sparse/row", "sparse/motif")
    _validate_same_length(row, count, "sparse/row", "sparse/count")
    if len(shape) != 2:
        raise ValueError(f"sparse/shape must have length 2, found {len(shape)}")
    if tuple(shape.astype(int)) != (n_rows, n_motifs):
        raise ValueError(
            "sparse/shape does not match row and motif axes: "
            f"shape={tuple(shape)}, coordinates={(n_rows, n_motifs)}"
        )
    if len(row) > 0:
        if np.any(row < 0) or np.any(row >= n_rows):
            raise ValueError("sparse/row contains an index outside the row axis")
        if np.any(motif < 0) or np.any(motif >= n_motifs):
            raise ValueError("sparse/motif contains an index outside the motif axis")

    sparse_dimension = _read_array(store, "sparse/sparse_dimension")
    _validate_axis(sparse_dimension, "sparse_dimension")
    # Label order matters because sparse COO columns are interpreted by position
    sparse_dimension_labels = _read_labels(
        store["sparse/sparse_dimension"],
        "sparse_dimension_name",
        len(sparse_dimension),
        "sparse_dimension",
    )
    if sparse_dimension_labels.tolist() != ["row", "motif"]:
        raise ValueError(
            "sparse_dimension labels must be ['row', 'motif'], "
            f"found {sparse_dimension_labels.tolist()!r}"
        )


def _validate_array_dimensions(
    array: Any, expected: tuple[str, ...], array_name: str
) -> None:
    """
    Require a Zarr array to expose expected named dimensions.
    """
    dimension_names = tuple(getattr(array.metadata, "dimension_names", ()) or ())
    if dimension_names != expected:
        raise ValueError(
            f"{array_name} dimensions must be {expected}, found {dimension_names}"
        )


def _read_array(store: Any, name: str) -> np.ndarray:
    """
    Load a Zarr array fully into a NumPy array.
    """
    return np.asarray(store[name][:])


def _read_motif_ascii_labels(store: Any, expected_len: int) -> np.ndarray:
    """
    Decode fixed-width ASCII motif labels into Python strings.

    The schema stores motif labels as bytes so labels are portable across
    languages and do not depend on object-string array support.
    """
    motif_byte = _read_array(store, "motif_byte")
    _validate_axis(motif_byte, "motif_byte")

    motif_ascii = _read_array(store, "motif_ascii")
    if motif_ascii.ndim != 2:
        raise ValueError(f"motif_ascii must have rank 2, found rank {motif_ascii.ndim}")
    expected_shape = (expected_len, len(motif_byte))
    if tuple(motif_ascii.shape) != expected_shape:
        raise ValueError(
            "motif_ascii shape does not match motif axes: "
            f"motif_ascii={motif_ascii.shape}, expected={expected_shape}"
        )
    if motif_ascii.dtype != np.uint8:
        raise ValueError(
            f"motif_ascii must have dtype uint8, found {motif_ascii.dtype}"
        )

    # Each row is one fixed-width ASCII motif label
    try:
        labels = [bytes(row).decode("ascii") for row in motif_ascii]
    except UnicodeDecodeError as error:
        raise ValueError("motif_ascii contains non-ASCII bytes") from error
    return np.asarray(labels, dtype=str)


def _read_labels(
    array: Any, field_name: str, expected_len: int, array_name: str
) -> np.ndarray:
    """
    Read string labels from a coordinate array's metadata.
    """
    label_field = array.attrs.get("label_field")
    if label_field != field_name:
        raise ValueError(
            f"{array_name} labels must have label_field={field_name!r}, "
            f"found {label_field!r}"
        )
    labels = array.attrs.get("labels")
    if labels is None:
        raise ValueError(f"{array_name} array is missing labels")
    labels = np.asarray(labels, dtype=str)
    if len(labels) != expected_len:
        raise ValueError(
            f"{array_name} labels length ({len(labels)}) does not match "
            f"axis length ({expected_len})"
        )
    return labels


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


def _require_densify(allow_densify: bool, method_name: str) -> None:
    """
    Require explicit opt-in before expanding sparse output to a dense array.
    """
    if not allow_densify:
        raise ValueError(
            f"{method_name}() would densify a sparse end-motif store. "
            "Use sparse_coo() or pass allow_densify=True."
        )


def _required(
    value: np.ndarray | zarr.Array | None, name: str
) -> np.ndarray | zarr.Array:
    """
    Return a loaded value or fail with schema context when it is unavailable.
    """
    if value is None:
        raise ValueError(f"End-motif Zarr store does not contain {name}")
    return value
