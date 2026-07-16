"""
Load cfDNAlab end-motif Zarr outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import numbers
import pathlib
from typing import Any, Sequence

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
)
from .ref_kmers import RefKmerFrequencies

END_MOTIF_MIN_SUPPORTED_SCHEMA_VERSION = 1
END_MOTIF_MAX_SUPPORTED_SCHEMA_VERSION = 2
VALID_STORAGE_MODES = {"dense", "sparse_coo"}
VALID_ROW_MODES = {"global", "size", "bed", "grouped_bed"}
VALID_MOTIF_AXIS_KINDS = {"motif", "motif_group"}


@dataclass
class LoadedEndMotifs:
    """Validated end-motif Zarr handles and row or motif metadata."""

    store: Any
    storage_mode: str
    row_mode: str
    motif_axis_kind: str
    motif_index: np.ndarray
    motif_names: np.ndarray | None
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
        _reject_empty_end_motif_counts_from_metadata(path)

        try:
            store = zarr.open_group(str(path), mode="r", zarr_format=3)
        except Exception as error:
            raise ValueError(
                f"Could not open end-motif Zarr store at {path}"
            ) from error

        storage_mode, row_mode, motif_axis_kind = _validate_root_metadata(store)
        _validate_required_arrays(store, storage_mode, row_mode, motif_axis_kind)

        # Motif and row axes are small coordinate arrays, so keep them in memory
        motif_index = _read_array(store, "motif_index")
        if motif_axis_kind == "motif":
            motif_names = _read_motif_ascii_labels(store, len(motif_index))
        else:
            motif_names = _read_labels(
                store["motif_index"],
                "motif_group",
                len(motif_index),
                "motif_index",
            )
        _validate_unique_labels(motif_names, "end-motif")
        row = _read_array(store, "row")
        _validate_axis(motif_index, "motif_index")
        _validate_axis(row, "row")
        if len(motif_index) == 0:
            _raise_no_end_motif_counts_available()

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
            if len(sparse_count) == 0:
                _raise_no_end_motif_counts_available()
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
            if len(row) != 1 or row_labels.tolist() != ["global"]:
                raise ValueError(
                    "global end-motif output must contain exactly one row "
                    "labeled 'global'"
                )
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
            motif_axis_kind=motif_axis_kind,
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

    def motifs_metadata(self) -> pd.DataFrame:
        """
        Return motif-axis labels and motif indices available in this output.

        For grouped motifs-file output, the `motif` labels are the group names
        used during counting.

        Returns
        -------
        pandas.DataFrame
            Columns are `motif_index` and `motif`.
        """
        return pd.DataFrame(
            {
                "motif_index": self.end_motifs.motif_index,
                "motif": _required(self.end_motifs.motif_names, "motif_names"),
            }
        )

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
        motif_names = _required(self.end_motifs.motif_names, "motif_names")
        return bool(np.any(motif_names == motif))

    def _sparse_counts_matrix_all(self) -> sparse.coo_matrix:
        """
        Return all end-motif counts as a SciPy COO matrix.
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

    def _sparse_counts_matrix_for_indices(
        self, row_indices: np.ndarray, motif_indices: np.ndarray
    ) -> sparse.coo_matrix:
        """
        Return selected rows and motifs as a SciPy COO matrix.

        SciPy sparse indexing preserves the requested row and motif order, so
        callers get a matrix whose axes match the selector order rather than
        the original file order.
        """
        if len(row_indices) == 0 or len(motif_indices) == 0:
            return sparse.coo_matrix((len(row_indices), len(motif_indices)))
        return (
            self._sparse_counts_matrix_all()
            .tocsr()[row_indices, :][:, motif_indices]
            .tocoo()
        )

    def _sparse_counts_matrix(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Shared implementation behind mode-specific sparse count methods.
        """
        row_indices = self._resolve_row_selector(window_idxs, groups, group_idxs)
        motif_indices = self._resolve_motif_selector(motifs, motif_idxs)
        return self._sparse_counts_matrix_for_indices(row_indices, motif_indices)

    def _dense_counts_array_for_indices(
        self,
        row_indices: np.ndarray,
        motif_indices: np.ndarray,
        *,
        allow_densify: bool,
    ) -> np.ndarray:
        """
        Return selected rows and motifs as a dense NumPy array.
        """
        allow_densify = validate_scalar_bool(allow_densify, "allow_densify")
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(
                counts.get_orthogonal_selection((row_indices, motif_indices))
            )

        _require_densify(allow_densify, "dense_counts_array")
        return self._sparse_counts_matrix_for_indices(
            row_indices, motif_indices
        ).toarray()

    def _dense_counts_array(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Shared implementation behind mode-specific dense count methods.
        """
        row_indices = self._resolve_row_selector(window_idxs, groups, group_idxs)
        motif_indices = self._resolve_motif_selector(motifs, motif_idxs)
        return self._dense_counts_array_for_indices(
            row_indices,
            motif_indices,
            allow_densify=allow_densify,
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
        return self._complete_data_frame_for_indices(
            row_indices, motif_indices, densify
        )

    def _data_frame_with_optional_reference_correction(
        self,
        *,
        ref_kmers: RefKmerFrequencies | None = None,
        densify: bool = False,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> pd.DataFrame:
        """
        Dispatch row-mode data frame methods to raw or corrected loading.

        Without `ref_kmers`, this preserves the ordinary selectors and returns
        raw rows. A two-sided correction choice is rejected on that path because
        there is no reference to correct against. With `ref_kmers`, the same row,
        motif, densification, and blacklist selectors are forwarded to the
        reference-correction implementation.
        """
        if ref_kmers is None:
            if two_sided_correction is not None:
                raise ValueError("two_sided_correction requires ref_kmers")
            return self._data_frame(
                densify=densify,
                window_idxs=window_idxs,
                groups=groups,
                group_idxs=group_idxs,
                motifs=motifs,
                motif_idxs=motif_idxs,
                max_blacklisted_fraction=max_blacklisted_fraction,
            )
        from .reference_correction import _reference_corrected_data_frame

        return _reference_corrected_data_frame(
            self,
            ref_kmers,
            window_idxs=window_idxs,
            groups=groups,
            group_idxs=group_idxs,
            densify=densify,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def _corrected_counts_array(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> np.ndarray:
        """
        Forward dense corrected-count loading through the shared implementation.

        Keeping this forwarding path shared across row modes ensures that row
        and motif selectors, fixed-shape validation, densification rules, and
        two-sided correction choices have the same behavior for every caller.
        """
        from .reference_correction import _reference_corrected_counts_array

        return _reference_corrected_counts_array(
            self,
            ref_kmers,
            window_idxs=window_idxs,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def _sparse_corrected_counts_matrix(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> sparse.coo_matrix:
        """
        Forward sparse corrected-count loading through the shared implementation.

        Keeping this forwarding path shared across row modes preserves selector
        validation, the fixed motif axis, unsupported-reference handling, and
        two-sided correction behavior without constructing a dense matrix.
        """
        from .reference_correction import _sparse_reference_corrected_counts_matrix

        return _sparse_reference_corrected_counts_matrix(
            self,
            ref_kmers,
            window_idxs=window_idxs,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def corrected_motifs_metadata(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        use_global_bias: bool = False,
        two_sided_correction: str | None = None,
    ) -> pd.DataFrame:
        """
        Return the motif axis used by reference-corrected matrices.

        Use this method to interpret the columns returned by
        `corrected_counts_array()` and `sparse_corrected_counts_matrix()`. Rows
        are in matrix-column order. `matrix_column` is the zero-based
        column in the returned matrix. `motif_index` is zero-based and refers to
        the full correction-mode axis described by `motif`, not necessarily the
        motif axis stored in the end-motif file.

        For `"joint"`, `"split"`, and one-sided correction, this is the
        selected stored motif axis. For `"outside"` and `"inside"`, repeated
        side labels are deduplicated in their first stored-motif occurrence
        order. Label selection returns labels in the requested order. Motif
        index selection is not available for `"outside"` or `"inside"` because
        those modes create a new axis.

        Parameters
        ----------
        ref_kmers
            Loaded reference k-mer output used for correction.
        motifs
            Motif label or labels on the correction-mode axis. Use either
            `motifs` or `motif_idxs`, not both.
        motif_idxs
            Stored motif index or indices. This is only available for
            `"joint"`, `"split"`, and one-sided correction.
        use_global_bias
            Whether a global reference k-mer output may be applied to a
            non-global end-motif output.
        two_sided_correction
            Required for two-sided motif labels. Use `"joint"`, `"split"`,
            `"outside"`, or `"inside"`. Leave as `None` for one-sided motifs
            or motif groups.

        Returns
        -------
        pandas.DataFrame
            Matrix-column metadata with `matrix_column`, `motif_index`, and
            `motif` columns.
        """
        from .reference_correction import (
            _resolve_correction_mode,
            _selected_mode_axis,
            _validate_reference_correction_inputs,
            _validate_two_sided_correction,
        )

        two_sided_correction = _validate_two_sided_correction(two_sided_correction)
        use_global_bias = validate_scalar_bool(use_global_bias, "use_global_bias")
        _validate_reference_correction_inputs(self, ref_kmers, use_global_bias)
        correction_mode = _resolve_correction_mode(
            self,
            ref_kmers,
            two_sided_correction,
        )
        motif_labels = _selected_mode_axis(
            self,
            correction_mode,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

        # Map selected labels back to their position on the complete mode axis
        # because matrix-column order can differ from mode-axis order
        full_mode_labels = (
            self.motifs_metadata()["motif"].astype(str).tolist()
            if correction_mode.mode in {"exact", "split"}
            else list(correction_mode.side_labels)
        )
        mode_index_by_label = {
            label: index for index, label in enumerate(full_mode_labels)
        }
        motif_indices = np.asarray(
            [mode_index_by_label[label] for label in motif_labels],
            dtype=np.int64,
        )
        return pd.DataFrame(
            {
                "matrix_column": np.arange(len(motif_labels), dtype=np.int64),
                "motif_index": motif_indices,
                "motif": motif_labels,
            }
        )

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

    def _complete_data_frame_for_indices(
        self, row_indices: np.ndarray, motif_indices: np.ndarray, densify: bool
    ) -> pd.DataFrame:
        """
        Build rows for every selected output row and motif.
        """
        row_metadata = (
            self._row_metadata_data_frame().iloc[row_indices].reset_index(drop=True)
        )
        motif_metadata = (
            self._motif_axis_metadata_data_frame()
            .iloc[motif_indices]
            .reset_index(drop=True)
        )
        row_count = len(row_indices)
        motif_count = len(motif_indices)
        if row_count == 0 or motif_count == 0:
            return self._empty_data_frame()

        selected_counts = self._dense_counts_array_for_indices(
            row_indices,
            motif_indices,
            allow_densify=densify,
        )
        repeated_rows = row_metadata.loc[
            row_metadata.index.repeat(motif_count)
        ].reset_index(drop=True)
        repeated_motifs = pd.concat([motif_metadata] * row_count, ignore_index=True)
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
        row_positions = {
            int(row_idx): order for order, row_idx in enumerate(row_indices)
        }
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

        row_metadata = (
            self._row_metadata_data_frame().iloc[matched_rows].reset_index(drop=True)
        )
        motif_metadata = (
            self._motif_axis_metadata_data_frame()
            .iloc[matched_motifs]
            .reset_index(drop=True)
        )
        data_frame = pd.concat([row_metadata, motif_metadata], axis=1)
        data_frame["count"] = matched_counts
        return data_frame

    def _empty_data_frame(self) -> pd.DataFrame:
        """
        Return an empty end-motif data frame with public columns.
        """
        row_metadata = self._row_metadata_data_frame().iloc[[]].reset_index(drop=True)
        motif_metadata = (
            self._motif_axis_metadata_data_frame().iloc[[]].reset_index(drop=True)
        )
        data_frame = pd.concat([row_metadata, motif_metadata], axis=1)
        data_frame["count"] = np.asarray([], dtype=float)
        return data_frame

    def _motif_axis_metadata_data_frame(self) -> pd.DataFrame:
        """
        Build metadata for the active count-column axis.
        """
        return self.motifs_metadata()

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
            data_frame["start"] = _required(
                self.end_motifs.row_start_bp, "row_start_bp"
            )
            data_frame["end"] = _required(self.end_motifs.row_end_bp, "row_end_bp")
            data_frame["blacklisted_fraction"] = _required(
                self.end_motifs.blacklisted_fraction, "blacklisted_fraction"
            )
        elif self.end_motifs.row_mode == "grouped_bed":
            data_frame = pd.DataFrame()
            data_frame["group_idx"] = _required(self.end_motifs.group_idx, "group")
            data_frame["group_name"] = _required(
                self.end_motifs.group_names, "group_name"
            )
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
                raise ValueError(
                    "Grouped selectors can only be used with grouped output"
                )
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
        Normalize column selectors for the active motif axis.
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
            _required(self.end_motifs.motif_names, "motif_names") == motif,
            missing_message=f"Unknown end-motif label: {motif!r}",
            duplicate_message=f"End-motif label is not unique: {motif!r}",
        )


class GlobalEndMotifCounts(EndMotifCounts):
    """End-motif counts for global output."""

    def data_frame(
        self,
        *,
        ref_kmers: RefKmerFrequencies | None = None,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame for global end-motif counts.

        Sparse outputs return stored non-zero motif counts unless
        `densify=True`. Densifying adds explicit zero-count rows for selected
        observed motifs. Dense outputs always include zero counts. Pass
        `ref_kmers` to add reference-corrected counts.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        `corrected_frequency` is normalized from `corrected_count` over the full
        correction-mode motif axis for each output row. Motif selection filters
        those frequencies afterward and does not renormalize them. A selected
        subset can therefore sum to less than 1. If the corrected total is
        zero, finite frequencies are zero. With `unsupported_motifs="keep_na"`,
        one undefined positive corrected count makes all frequencies in that
        output row `NaN`.

        Parameters
        ----------
        ref_kmers
            Optional loaded reference k-mer output used for correction.
        densify
            If `True`, sparse outputs add explicit zero-count rows for selected
            observed motifs. Dense outputs ignore this option.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.
        use_global_bias
            Whether a global reference k-mer output may be applied to every row.
        unsupported_motifs
            What to do when an observed sample motif has no positive correction
            factor under the selected mode. Use `"error"`, `"drop"`, or `"keep_na"`.
        two_sided_correction
            Required for two-sided motif labels such as `"AC_GT"` when
            `ref_kmers` is passed. Use `"joint"`, `"split"`, `"outside"`, or
            `"inside"`. Leave as `None` for one-sided motifs or motif groups.

        Returns
        -------
        pandas.DataFrame
            Global row metadata, motif metadata, and `count`. If `ref_kmers`
            is passed, also includes `corrected_count` and
            `corrected_frequency`.
        """
        return self._data_frame_with_optional_reference_correction(
            ref_kmers=ref_kmers,
            densify=densify,
            motifs=motifs,
            motif_idxs=motif_idxs,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def dense_counts_array(
        self,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return global end-motif counts as a dense NumPy array.

        Sparse stores are only densified when `allow_densify=True`. Scalar
        motif selectors keep their axis as length one, so the shape is always
        `(1, selected motifs)`.

        Parameters
        ----------
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.
        allow_densify
            If `True`, allow sparse stores to be converted to dense counts.

        Returns
        -------
        numpy.ndarray
            Dense count array with shape `(global row, motif)`.
        """
        return self._dense_counts_array(
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def sparse_counts_matrix(
        self,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return global end-motif counts as a SciPy sparse matrix.

        Scalar motif selectors keep their axis as length one, so the shape is
        always `(1, selected motifs)`.

        Parameters
        ----------
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse count matrix with shape `(global row, motif)`.
        """
        return self._sparse_counts_matrix(
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def corrected_counts_array(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> np.ndarray:
        """
        Return global reference-corrected end-motif counts as a dense array.

        The result has one column per motif on the selected correction-mode
        axis. `"joint"`, `"split"`, and one-sided correction retain the selected
        stored motif axis. `"outside"` and `"inside"` replace it with a
        deduplicated side axis, so the number of columns can differ from
        `dense_counts_array()`. Use `corrected_motifs_metadata()` to inspect the
        exact column labels and order.

        Sparse end-motif stores are not densified unless
        `allow_densify=True`. Use `sparse_corrected_counts_matrix()` to keep a
        sparse result.

        `unsupported_motifs="drop"` is not allowed because arrays have a fixed
        row and motif shape. Use
        `data_frame(ref_kmers=..., unsupported_motifs="drop")` when
        unsupported motifs should be omitted.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        An observed sample motif with a positive count is unsupported when it
        has no positive correction factor under the selected mode.
        `unsupported_motifs="keep_na"` keeps
        that matrix cell as `NaN`. The `"drop"` policy is unavailable for
        matrices because it would change a fixed result axis.

        """
        return self._corrected_counts_array(
            ref_kmers,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def sparse_corrected_counts_matrix(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> sparse.coo_matrix:
        """
        Return global reference-corrected end-motif counts as a sparse matrix.

        The result has one column per motif on the selected correction-mode
        axis. `"outside"` and `"inside"` can therefore have fewer columns than
        `sparse_counts_matrix()`. Use `corrected_motifs_metadata()` to inspect
        the exact column labels and order. Corrected zeroes are not stored.
        Corrected `NaN` values from
        `unsupported_motifs="keep_na"` are stored so they remain visible.

        `unsupported_motifs="drop"` is not allowed because sparse matrices
        still have a fixed row and motif shape. Use
        `data_frame(ref_kmers=..., unsupported_motifs="drop")` when
        unsupported motifs should be omitted.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        An observed sample motif with a positive count is unsupported when it
        has no positive correction factor under the selected mode.
        `unsupported_motifs="keep_na"` keeps
        that matrix cell as `NaN`. The `"drop"` policy is unavailable for
        matrices because it would change a fixed result axis.
        """
        return self._sparse_corrected_counts_matrix(
            ref_kmers,
            motifs=motifs,
            motif_idxs=motif_idxs,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )


class WindowedEndMotifCounts(EndMotifCounts):
    """End-motif counts for fixed-size or BED-window output."""

    def data_frame(
        self,
        *,
        ref_kmers: RefKmerFrequencies | None = None,
        window_idxs: int | Sequence[int] | None = None,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of end-motif counts for genomic windows.

        Use `window_idxs` to keep only selected windows and `motifs` or
        `motif_idxs` to keep only selected motifs. Sparse outputs return stored
        non-zero rows unless `densify=True`. Densifying adds explicit
        zero-count rows for selected observed motifs. Dense outputs always
        include zero counts. Pass `ref_kmers` to add reference-corrected counts.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        `corrected_frequency` is normalized from `corrected_count` over the full
        correction-mode motif axis for each output row. Motif selection filters
        those frequencies afterward and does not renormalize them. A selected
        subset can therefore sum to less than 1. If the corrected total is
        zero, finite frequencies are zero. With `unsupported_motifs="keep_na"`,
        one undefined positive corrected count makes all frequencies in that
        output row `NaN`.

        Parameters
        ----------
        ref_kmers
            Optional loaded reference k-mer output used for correction.
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
        use_global_bias
            Whether a global reference k-mer output may be applied to every row.
        unsupported_motifs
            What to do when an observed sample motif has no positive correction
            factor under the selected mode. Use `"error"`, `"drop"`, or `"keep_na"`.
        two_sided_correction
            Required for two-sided motif labels such as `"AC_GT"` when
            `ref_kmers` is passed. Use `"joint"`, `"split"`, `"outside"`, or
            `"inside"`. Leave as `None` for one-sided motifs or motif groups.

        Returns
        -------
        pandas.DataFrame
            Window metadata, motif metadata, and `count`. If `ref_kmers` is
            passed, also includes `corrected_count` and
            `corrected_frequency`.
        """
        return self._data_frame_with_optional_reference_correction(
            ref_kmers=ref_kmers,
            densify=densify,
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def window_metadata(self) -> pd.DataFrame:
        """
        Return genomic window metadata for this end-motif output.

        Public genomic window metadata uses `window_idx`, `chrom`, `start`,
        and `end` columns.

        Returns
        -------
        pandas.DataFrame
            Columns are `window_idx`, `chrom`, `start`, `end`, and
            `blacklisted_fraction`.
        """
        return self._row_metadata_data_frame()

    def dense_counts_array(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return windowed end-motif counts as a dense NumPy array.

        Sparse stores are only densified when `allow_densify=True`. Scalar
        selectors keep their axes as length one, so the shape is always
        `(selected windows, selected motifs)`.

        Parameters
        ----------
        window_idxs
            `None` for all windows, one window index, or a sequence of window
            indices.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.
        allow_densify
            If `True`, allow sparse stores to be converted to dense counts.

        Returns
        -------
        numpy.ndarray
            Dense count array with shape `(window, motif)`.
        """
        return self._dense_counts_array(
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def sparse_counts_matrix(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return windowed end-motif counts as a SciPy sparse matrix.

        Scalar selectors keep their axes as length one, so the shape is always
        `(selected windows, selected motifs)`.

        Parameters
        ----------
        window_idxs
            `None` for all windows, one window index, or a sequence of window
            indices.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse count matrix with shape `(window, motif)`.
        """
        return self._sparse_counts_matrix(
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def corrected_counts_array(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> np.ndarray:
        """
        Return windowed reference-corrected end-motif counts as a dense array.

        The result has one row per selected window and one column per motif on
        the selected correction-mode axis. `"outside"` and `"inside"` replace
        the stored joint axis with a deduplicated side axis. Use
        `corrected_motifs_metadata()` to inspect the exact column labels and
        order. Sparse end-motif stores are not densified unless
        `allow_densify=True`. Use
        `sparse_corrected_counts_matrix()` to keep a sparse result.

        `unsupported_motifs="drop"` is not allowed because arrays have a fixed
        row and motif shape. Use
        `data_frame(ref_kmers=..., unsupported_motifs="drop")` when
        unsupported motifs should be omitted.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        An observed sample motif with a positive count is unsupported when it
        has no positive correction factor under the selected mode.
        `unsupported_motifs="keep_na"` keeps
        that matrix cell as `NaN`. The `"drop"` policy is unavailable for
        matrices because it would change a fixed result axis.

        """
        return self._corrected_counts_array(
            ref_kmers,
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def sparse_corrected_counts_matrix(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> sparse.coo_matrix:
        """
        Return windowed reference-corrected end-motif counts as a sparse matrix.

        The result has one row per selected window and one column per motif on
        the selected correction-mode axis. `"outside"` and `"inside"` can have
        fewer columns than `sparse_counts_matrix()`. Use
        `corrected_motifs_metadata()` to inspect the exact column labels and
        order. Corrected zeroes are not stored. Corrected `NaN` values from
        `unsupported_motifs="keep_na"` are stored so they remain visible.

        `unsupported_motifs="drop"` is not allowed because sparse matrices
        still have a fixed row and motif shape. Use
        `data_frame(ref_kmers=..., unsupported_motifs="drop")` when
        unsupported motifs should be omitted.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        An observed sample motif with a positive count is unsupported when it
        has no positive correction factor under the selected mode.
        `unsupported_motifs="keep_na"` keeps
        that matrix cell as `NaN`. The `"drop"` policy is unavailable for
        matrices because it would change a fixed result axis.
        """
        return self._sparse_corrected_counts_matrix(
            ref_kmers,
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )


class GroupedEndMotifCounts(EndMotifCounts):
    """End-motif counts for grouped BED output."""

    def data_frame(
        self,
        *,
        ref_kmers: RefKmerFrequencies | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame of end-motif counts for grouped BED rows.

        Use `groups` or `group_idxs` to keep only selected groups and `motifs`
        or `motif_idxs` to keep only selected motifs. Sparse outputs return
        stored non-zero rows unless `densify=True`. Densifying adds explicit
        zero-count rows for selected observed motifs. Dense outputs always
        include zero counts. Pass `ref_kmers` to add reference-corrected counts.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        `corrected_frequency` is normalized from `corrected_count` over the full
        correction-mode motif axis for each output row. Motif selection filters
        those frequencies afterward and does not renormalize them. A selected
        subset can therefore sum to less than 1. If the corrected total is
        zero, finite frequencies are zero. With `unsupported_motifs="keep_na"`,
        one undefined positive corrected count makes all frequencies in that
        output row `NaN`.

        Parameters
        ----------
        ref_kmers
            Optional loaded reference k-mer output used for correction.
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
        use_global_bias
            Whether a global reference k-mer output may be applied to every row.
        unsupported_motifs
            What to do when an observed sample motif has no positive correction
            factor under the selected mode. Use `"error"`, `"drop"`, or `"keep_na"`.
        two_sided_correction
            Required for two-sided motif labels such as `"AC_GT"` when
            `ref_kmers` is passed. Use `"joint"`, `"split"`, `"outside"`, or
            `"inside"`. Leave as `None` for one-sided motifs or motif groups.

        Returns
        -------
        pandas.DataFrame
            Group metadata, motif metadata, and `count`. If `ref_kmers` is
            passed, also includes `corrected_count` and
            `corrected_frequency`.
        """
        return self._data_frame_with_optional_reference_correction(
            ref_kmers=ref_kmers,
            densify=densify,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def group_metadata(self) -> pd.DataFrame:
        """
        Return grouped BED metadata for this end-motif output.

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

    def dense_counts_array(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return grouped end-motif counts as a dense NumPy array.

        Sparse stores are only densified when `allow_densify=True`. Scalar
        selectors keep their axes as length one, so the shape is always
        `(selected groups, selected motifs)`.

        Parameters
        ----------
        groups
            `None` for all groups, one group name, or a sequence of group names.
            Use either `groups` or `group_idxs`, not both.
        group_idxs
            `None` for all groups, one group index, or a sequence of group
            indices. Use either `groups` or `group_idxs`, not both.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.
        allow_densify
            If `True`, allow sparse stores to be converted to dense counts.

        Returns
        -------
        numpy.ndarray
            Dense count array with shape `(group, motif)`.
        """
        return self._dense_counts_array(
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def sparse_counts_matrix(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return grouped end-motif counts as a SciPy sparse matrix.

        Scalar selectors keep their axes as length one, so the shape is always
        `(selected groups, selected motifs)`.

        Parameters
        ----------
        groups
            `None` for all groups, one group name, or a sequence of group names.
            Use either `groups` or `group_idxs`, not both.
        group_idxs
            `None` for all groups, one group index, or a sequence of group
            indices. Use either `groups` or `group_idxs`, not both.
        motifs
            Motif label or labels. Use either `motifs` or `motif_idxs`, not both.
        motif_idxs
            Motif index or indices. Use either `motifs` or `motif_idxs`, not both.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse count matrix with shape `(group, motif)`.
        """
        return self._sparse_counts_matrix(
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def corrected_counts_array(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> np.ndarray:
        """
        Return grouped reference-corrected end-motif counts as a dense array.

        The result has one row per selected group and one column per motif on
        the selected correction-mode axis. `"outside"` and `"inside"` replace
        the stored joint axis with a deduplicated side axis. Use
        `corrected_motifs_metadata()` to inspect the exact column labels and
        order. Sparse end-motif stores are not densified unless
        `allow_densify=True`. Use
        `sparse_corrected_counts_matrix()` to keep a sparse result.

        `unsupported_motifs="drop"` is not allowed because arrays have a fixed
        row and motif shape. Use
        `data_frame(ref_kmers=..., unsupported_motifs="drop")` when
        unsupported motifs should be omitted.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        An observed sample motif with a positive count is unsupported when it
        has no positive correction factor under the selected mode.
        `unsupported_motifs="keep_na"` keeps
        that matrix cell as `NaN`. The `"drop"` policy is unavailable for
        matrices because it would change a fixed result axis.
        """
        return self._corrected_counts_array(
            ref_kmers,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )

    def sparse_corrected_counts_matrix(
        self,
        ref_kmers: RefKmerFrequencies,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        max_blacklisted_fraction: float = 1.0,
        use_global_bias: bool = False,
        unsupported_motifs: str = "error",
        two_sided_correction: str | None = None,
    ) -> sparse.coo_matrix:
        """
        Return grouped reference-corrected end-motif counts as a sparse matrix.

        The result has one row per selected group and one column per motif on
        the selected correction-mode axis. `"outside"` and `"inside"` can have
        fewer columns than `sparse_counts_matrix()`. Use
        `corrected_motifs_metadata()` to inspect the exact column labels and
        order. Corrected zeroes are not stored. Corrected `NaN` values from
        `unsupported_motifs="keep_na"` are stored so they remain visible.

        `unsupported_motifs="drop"` is not allowed because sparse matrices
        still have a fixed row and motif shape. Use
        `data_frame(ref_kmers=..., unsupported_motifs="drop")` when
        unsupported motifs should be omitted.

        Reference correction
        --------------------
        Reference correction divides each observed end-motif count by a
        reference-based correction factor for the matched row. This factor is
        computed from the motif frequencies in the reference k-mer output and
        normalized so a uniform reference composition leaves counts unchanged.
        Motifs that are common in the reference row are scaled down. Motifs
        that are rare in the reference row are scaled up. Only motifs with a
        positive reference frequency contribute to the row's correction
        support.

        Two-sided correction modes
        --------------------------
        When motif labels contain both outside and inside bases, such as
        `"AC_GT"`, `two_sided_correction` chooses both the motif labels in the
        result and the correction factor used for each returned count.

        - `"joint"` keeps full labels such as `"AC_GT"` and corrects each
          count using the exact reference k-mer `"ACGT"`.

        - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
          correction factor from the two sides separately. For `"AC_GT"`,
          separate correction factors are calculated for outside label `"AC"`
          and inside label `"GT"`. Those two correction factors are multiplied
          and applied to the observed `"AC_GT"` count. Use this when you want
          full two-sided motif labels in the result, but the exact full
          reference k-mers are too sparse or you want the reference correction
          to treat outside and inside sequence composition separately.

        - `"outside"` returns outside labels such as `"AC_"`. For each outside
          label, all full motif counts with that outside label are summed first.
          For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
          count. That summed count is corrected using the outside label `"AC"`.

        - `"inside"` returns inside labels such as `"_GT"`. For each inside
          label, all full motif counts with that inside label are summed first.
          For example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"`
          count. That summed count is corrected using the inside label `"GT"`.

        For `"split"`, `"outside"`, and `"inside"`, side-specific reference
        frequencies are calculated from the loaded full-length reference
        k-mers. For example, the outside frequency for `"AC"` is the sum of
        frequencies for loaded k-mers with prefix `"AC"`, such as `"ACTG"` and
        `"ACAA"`. The inside frequency for `"TG"` is the corresponding sum over
        loaded k-mers with suffix `"TG"`. Separate shorter reference k-mer runs
        are not required.

        A motifs file used for the reference output restricts these sums to the
        k-mers in that file. Without a motifs file, all k-mers in the reference
        output can contribute, including k-mers absent from the sample
        end-motif output.

        An observed sample motif with a positive count is unsupported when it
        has no positive correction factor under the selected mode.
        `unsupported_motifs="keep_na"` keeps
        that matrix cell as `NaN`. The `"drop"` policy is unavailable for
        matrices because it would change a fixed result axis.
        """
        return self._sparse_corrected_counts_matrix(
            ref_kmers,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
            use_global_bias=use_global_bias,
            unsupported_motifs=unsupported_motifs,
            two_sided_correction=two_sided_correction,
        )


def read_end_motifs(
    path: pathlib.Path | str,
) -> GlobalEndMotifCounts | WindowedEndMotifCounts | GroupedEndMotifCounts:
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


def _raise_no_end_motif_counts_available() -> None:
    """
    Raise the public error for sparse stores with no stored counts.
    """
    raise ValueError(
        "No end-motif counts are available in this store. "
        "If you expected motifs or groups from `--motifs-file` with zero counts to remain in the output, "
        "rerun `cfdna ends` with `--all-motifs`."
    )


def _reject_empty_end_motif_counts_from_metadata(path: pathlib.Path) -> None:
    """
    Reject sparse no-count stores before the Zarr library opens every array node.
    """
    root_attributes = _raw_zarr_attributes(path)
    if root_attributes.get("cfdnalab_schema") != "end_motif_counts":
        return

    motif_shape = _raw_zarr_array_shape(path, "motif_index")
    if motif_shape is not None and len(motif_shape) >= 1 and motif_shape[0] == 0:
        _raise_no_end_motif_counts_available()

    if root_attributes.get("storage_mode") == "sparse_coo":
        sparse_count_shape = _raw_zarr_array_shape(path, "sparse/count")
        if (
            sparse_count_shape is not None
            and len(sparse_count_shape) >= 1
            and sparse_count_shape[0] == 0
        ):
            _raise_no_end_motif_counts_available()


def _raw_zarr_attributes(path: pathlib.Path) -> dict[str, Any]:
    """
    Read root attributes directly from `zarr.json`.
    """
    try:
        metadata = json.loads((path / "zarr.json").read_text())
    except OSError:
        return {}
    attributes = metadata.get("attributes")
    return attributes if isinstance(attributes, dict) else {}


def _raw_zarr_array_shape(
    path: pathlib.Path, array_name: str
) -> tuple[int, ...] | None:
    """
    Read one array shape directly from its metadata file.
    """
    metadata_path = path.joinpath(*array_name.split("/"), "zarr.json")
    if not metadata_path.is_file():
        return None
    metadata = json.loads(metadata_path.read_text())
    shape = metadata.get("shape")
    if not isinstance(shape, list):
        return None
    parsed_shape = []
    for dimension in shape:
        if not isinstance(dimension, numbers.Integral) or int(dimension) < 0:
            raise ValueError(
                f"{array_name} metadata shape must be a non-negative integer vector"
            )
        parsed_shape.append(int(dimension))
    return tuple(parsed_shape)


def _validate_root_metadata(store: Any) -> tuple[str, str, str]:
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

    motif_axis_kind = store.attrs.get("motif_axis_kind")
    if motif_axis_kind is None:
        if int(schema_version) != 1:
            raise ValueError("end-motif schema v2 stores must declare motif_axis_kind")
        motif_axis_kind = "motif"
    if motif_axis_kind not in VALID_MOTIF_AXIS_KINDS:
        raise ValueError(f"Unsupported end-motif motif axis kind: {motif_axis_kind!r}")

    return storage_mode, row_mode, motif_axis_kind


def _validate_required_arrays(
    store: Any, storage_mode: str, row_mode: str, motif_axis_kind: str
) -> None:
    """
    Require all arrays needed by the active storage mode and row mode.
    """
    required = {"motif_index", "row"}
    if motif_axis_kind == "motif":
        required.update({"motif_byte", "motif_ascii"})
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
    if expected_len > 0 and len(motif_byte) == 0:
        raise ValueError(
            "motif_ascii cannot decode non-empty motif axis with zero motif_byte width"
        )

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
    for motif_index, label in enumerate(labels):
        if _contains_control_character(label):
            raise ValueError(
                f"motif_ascii row {motif_index} contains a control character"
            )
    return np.asarray(labels, dtype=str)


def _validate_unique_labels(labels: np.ndarray, label_name: str) -> None:
    """
    Reject duplicate labels on a public selector axis.
    """
    if len(set(labels)) != len(labels):
        raise ValueError(f"duplicate {label_name} label")


def _contains_control_character(label: str) -> bool:
    """
    Return whether a label contains a Unicode control character.
    """
    return any(
        ord(character) < 32 or 0x7F <= ord(character) <= 0x9F for character in label
    )


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
    if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence):
        raise ValueError(f"{array_name} labels must be a list of character strings")

    validated_labels: list[str] = []
    for label in labels:
        if not isinstance(label, str):
            raise ValueError(f"{array_name} labels must be character strings")
        if _contains_control_character(label):
            raise ValueError(f"{array_name} labels must not contain control characters")
        validated_labels.append(label)

    label_array = np.asarray(validated_labels, dtype=str)
    if len(label_array) != expected_len:
        raise ValueError(
            f"{array_name} labels length ({len(label_array)}) does not match "
            f"axis length ({expected_len})"
        )
    return label_array


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
            "Use sparse_counts_matrix() or pass allow_densify=True."
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
