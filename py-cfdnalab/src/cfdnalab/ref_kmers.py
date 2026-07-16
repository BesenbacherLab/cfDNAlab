"""
Load cfDNAlab reference k-mer Zarr outputs.
"""

from __future__ import annotations

from collections.abc import Sequence
import copy
from dataclasses import dataclass
import json
import numbers
import pathlib
import unicodedata
from typing import Any

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

REF_KMER_SUPPORTED_SCHEMA_VERSION = 1
VALID_STORAGE_MODES = {"dense", "sparse_coo"}
VALID_ROW_MODES = {"global", "size", "bed", "grouped_bed"}
VALID_MOTIF_AXIS_KINDS = {"motif", "motif_group"}
EXPECTED_VALUE_UNITS = "reference_kmer_frequency"
EXPECTED_COUNT_UNITS = "reference_kmer_count"
EXPECTED_ROW_SCALING_ARRAY = "row_scaling_factor"
EXPECTED_COUNT_RECONSTRUCTION = (
    "reference_kmer_count = frequency * row_scaling_factor[row]"
)


@dataclass
class LoadedRefKmers:
    """Validated reference k-mer Zarr handles and row or motif metadata."""

    store: Any
    storage_mode: str
    row_mode: str
    motif_axis_kind: str
    kmer_size: int
    canonical: bool
    all_motifs: bool
    assign_by: str
    motif_index: np.ndarray
    motif_names: np.ndarray
    row: np.ndarray
    row_scaling_factor: np.ndarray
    reference_contig_footprint: Any
    frequencies: zarr.Array | None
    sparse_row: np.ndarray | None
    sparse_motif: np.ndarray | None
    sparse_frequency: np.ndarray | None
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


class RefKmerFrequencies:
    """
    Common API for global, windowed, and grouped reference k-mer outputs.

    A reference k-mer output describes the expected background k-mer
    composition for one or more rows. A row can be the whole reference, a
    genomic window, a BED interval, or a grouped BED entry.

    The file stores frequencies. Count helpers reconstruct comparable counts as
    `frequency * row_scaling_factor[row]`, where the scaling factor is the
    number of reference k-mer positions represented by that row.

    Public Python selectors use zero-based indices. For example, `motif_idxs=0`
    selects the first motif in `motifs_metadata()`.
    """

    def __init__(
        self,
        path: pathlib.Path | str,
        loaded_ref_kmers: LoadedRefKmers | None = None,
    ) -> None:
        """
        Load a reference k-mer frequency output directory.

        The directory is a Zarr store on disk, but most users can work through
        the data frame, NumPy array, and SciPy sparse matrix methods without
        using Zarr directly.

        Parameters
        ----------
        path
            Path to a `<prefix>.ref_kmers.zarr` directory.
        loaded_ref_kmers
            Preloaded store data used by `read_ref_kmers`.
        """
        self.path = pathlib.Path(path)
        if loaded_ref_kmers is None:
            loaded_ref_kmers = RefKmerFrequencies._load_zarr(self.path)
        self.ref_kmers = loaded_ref_kmers

    def __repr__(self) -> str:
        """
        Return a compact summary with path, schema version, modes, and shape.
        """
        schema_version = self.ref_kmers.store.attrs.get("cfdnalab_schema_version")
        if self.ref_kmers.storage_mode == "dense":
            shape = tuple(_required(self.ref_kmers.frequencies, "frequencies").shape)
        else:
            shape = tuple(
                _required(self.ref_kmers.sparse_shape, "sparse/shape").astype(int)
            )
        return (
            f"{self.__class__.__name__}("
            f"path={str(self.path)!r}, "
            f"schema_version={schema_version!r}, "
            f"storage_mode={self.ref_kmers.storage_mode!r}, "
            f"row_mode={self.ref_kmers.row_mode!r}, "
            f"motif_axis_kind={self.ref_kmers.motif_axis_kind!r}, "
            f"shape={shape!r}"
            ")"
        )

    @staticmethod
    def _load_zarr(path: pathlib.Path | str) -> LoadedRefKmers:
        """
        Open and validate a reference k-mer frequency Zarr store.
        """
        path = pathlib.Path(path)
        _validate_zarr_store_path(path)

        try:
            store = zarr.open_group(str(path), mode="r", zarr_format=3)
        except Exception as error:
            raise ValueError(
                f"Could not open reference k-mer Zarr store at {path}"
            ) from error

        metadata = _validate_root_metadata(store)
        storage_mode = metadata["storage_mode"]
        row_mode = metadata["row_mode"]
        motif_axis_kind = metadata["motif_axis_kind"]
        _validate_required_arrays(store, storage_mode, row_mode, motif_axis_kind)

        _validate_array_dimensions(store["motif_index"], ("motif",), "motif_index")
        motif_index = _read_array(store, "motif_index")
        _validate_axis(motif_index, "motif_index")
        if motif_axis_kind == "motif":
            motif_names = _read_motif_ascii_labels(store, len(motif_index))
            _validate_concrete_motif_labels(
                motif_names,
                metadata["kmer_size"],
                metadata["canonical"],
            )
        else:
            motif_names = _read_labels(
                store["motif_index"],
                "motif_group",
                len(motif_index),
                "motif_index",
            )
            _validate_unique_labels(motif_names, "reference k-mer motif-group")

        _validate_array_dimensions(store["row"], ("row",), "row")
        row = _read_array(store, "row")
        _validate_axis(row, "row")

        row_scaling_factor = _read_array(store, "row_scaling_factor")
        _validate_array_dimensions(
            store["row_scaling_factor"],
            ("row",),
            "row_scaling_factor",
        )
        _validate_same_length(row_scaling_factor, row, "row_scaling_factor", "row")
        _validate_row_scaling_factors(row_scaling_factor)

        reference_contig_footprint = _read_reference_contig_footprint(store)

        frequencies = None
        sparse_row = None
        sparse_motif = None
        sparse_frequency = None
        sparse_shape = None
        if storage_mode == "dense":
            frequencies = store["frequencies"]
            _validate_dense_frequencies(frequencies, len(row), len(motif_index))
        else:
            sparse_row = _read_array(store, "sparse/row")
            sparse_motif = _read_array(store, "sparse/motif")
            sparse_frequency = _read_array(store, "sparse/frequency")
            sparse_shape = _read_array(store, "sparse/shape")
            _validate_sparse_arrays(
                store,
                sparse_row,
                sparse_motif,
                sparse_frequency,
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

        if row_mode == "global":
            row_labels = _read_labels(store["row"], "row_label", len(row), "row")
            if len(row) != 1 or row_labels.tolist() != ["global"]:
                raise ValueError(
                    "global reference k-mer output must contain exactly one row "
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
            _validate_array_dimensions(
                store["chromosome"], ("chromosome",), "chromosome"
            )
            _validate_array_dimensions(
                store["row_chromosome"], ("row",), "row_chromosome"
            )
            _validate_array_dimensions(store["row_start_bp"], ("row",), "row_start_bp")
            _validate_array_dimensions(store["row_end_bp"], ("row",), "row_end_bp")
            _validate_array_dimensions(
                store["blacklisted_fraction"], ("row",), "blacklisted_fraction"
            )
            _validate_axis(chromosome, "chromosome")
            _validate_same_length(row_chromosome, row, "row_chromosome", "row")
            _validate_same_length(row_start_bp, row, "row_start_bp", "row")
            _validate_same_length(row_end_bp, row, "row_end_bp", "row")
            _validate_same_length(
                blacklisted_fraction, row, "blacklisted_fraction", "row"
            )
            _validate_index_values(row_chromosome, len(chromosome), "row_chromosome")
            _validate_half_open_intervals(
                row_start_bp, row_end_bp, "row_start_bp", "row_end_bp"
            )
            _validate_fraction_values(blacklisted_fraction, "blacklisted_fraction")
        elif row_mode == "grouped_bed":
            group_idx = _read_array(store, "group")
            group_names = _read_labels(
                store["group"], "group_name", len(group_idx), "group"
            )
            eligible_windows = _read_array(store, "eligible_windows")
            blacklisted_fraction = _read_array(store, "blacklisted_fraction")
            _validate_array_dimensions(store["group"], ("row",), "group")
            _validate_array_dimensions(
                store["eligible_windows"], ("row",), "eligible_windows"
            )
            _validate_array_dimensions(
                store["blacklisted_fraction"], ("row",), "blacklisted_fraction"
            )
            _validate_same_length(group_idx, row, "group", "row")
            _validate_same_length(group_names, row, "group_name labels", "row")
            _validate_same_length(eligible_windows, row, "eligible_windows", "row")
            _validate_same_length(
                blacklisted_fraction, row, "blacklisted_fraction", "row"
            )
            _validate_axis(group_idx, "group")
            _validate_nonnegative_integer_values(eligible_windows, "eligible_windows")
            _validate_fraction_values(blacklisted_fraction, "blacklisted_fraction")

        return LoadedRefKmers(
            store=store,
            storage_mode=storage_mode,
            row_mode=row_mode,
            motif_axis_kind=motif_axis_kind,
            kmer_size=metadata["kmer_size"],
            canonical=metadata["canonical"],
            all_motifs=metadata["all_motifs"],
            assign_by=metadata["assign_by"],
            motif_index=motif_index,
            motif_names=motif_names,
            row=row,
            row_scaling_factor=row_scaling_factor,
            reference_contig_footprint=reference_contig_footprint,
            frequencies=frequencies,
            sparse_row=sparse_row,
            sparse_motif=sparse_motif,
            sparse_frequency=sparse_frequency,
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
        Return whether the output is saved as dense values or sparse values.

        Dense output has a value for every row and motif. Sparse output stores
        only non-zero values and is usually better for large outputs.
        """
        return self.ref_kmers.storage_mode

    def row_mode(self) -> str:
        """
        Return what each reference k-mer frequency row represents.

        `global` has one row for the whole reference. `size` and `bed` rows are
        genomic windows. `grouped_bed` rows are BED groups.
        """
        return self.ref_kmers.row_mode

    def motif_axis_kind(self) -> str:
        """
        Return whether columns represent concrete k-mers or motif groups.
        """
        return self.ref_kmers.motif_axis_kind

    def kmer_size(self) -> int:
        """
        Return the k-mer size used by the reference k-mer command.
        """
        return self.ref_kmers.kmer_size

    def canonical(self) -> bool:
        """
        Return whether reverse-complement k-mers were collapsed.
        """
        return self.ref_kmers.canonical

    def all_motifs(self) -> bool:
        """
        Return whether the command kept every requested motif target.

        For full k-mer output, this means every A/C/G/T k-mer for the
        requested k. For motifs-file output, this means every target from the
        motifs file.
        """
        return self.ref_kmers.all_motifs

    def assign_by(self) -> str:
        """
        Return the window assignment rule used by the command.
        """
        return self.ref_kmers.assign_by

    def reference_contig_footprint(self) -> Any:
        """
        Return which reference contigs contributed to the output.

        The result is decoded metadata from the output file. It is useful for
        checking that the reference used for loading matches the expected
        genome or contig subset.
        """
        return copy.deepcopy(self.ref_kmers.reference_contig_footprint)

    def motifs_metadata(self) -> pd.DataFrame:
        """
        Return the k-mer labels and motif indices available in this output.

        For grouped motifs-file output, the `motif` labels are the motif-group
        names used during counting. The returned `motif_index` values are
        zero-based and can be passed to `motif_idxs`.

        If `all_motifs()` is false, the motif axis is the combined set of
        motifs or motifs-file targets observed anywhere in the output.
        Densifying sparse output fills zeroes only across this listed axis. It
        does not add every possible k-mer.
        """
        return pd.DataFrame(
            {
                "motif_index": self.ref_kmers.motif_index,
                "motif": self.ref_kmers.motif_names,
            }
        )

    def motif_idx(self, motif: str) -> int:
        """
        Find the zero-based motif-axis index for a motif label.
        """
        return self._resolve_motif(motif)

    def has_motif(self, motif: str) -> bool:
        """
        Return whether a motif label exists in this output.

        Observed-only output can omit motifs that were not observed anywhere in
        the output, so this can return `False` even for a valid A/C/G/T k-mer.
        """
        return bool(np.any(self.ref_kmers.motif_names == motif))

    def row_scaling_factors(self) -> pd.DataFrame:
        """
        Return row metadata and factors used to reconstruct k-mer counts.

        Frequencies are fractions. Multiplying a row's frequency by its
        `row_scaling_factor` gives the reconstructed count for that row.
        """
        data_frame = self._row_metadata_data_frame().copy()
        data_frame["row_scaling_factor"] = self.ref_kmers.row_scaling_factor
        return data_frame

    def dense_frequencies_zarr_array(self) -> zarr.Array:
        """
        Return the on-disk frequency array for advanced dense-output workflows.

        Most users should use `dense_frequencies_array()`, which returns a
        NumPy array, or `data_frame()`, which returns a pandas data frame. This
        method returns the underlying Zarr array so advanced users can slice it
        without first reading the whole array into memory.

        This is only available for dense output. Sparse output does not have an
        on-disk dense frequency array. Use `sparse_frequencies_matrix()` for
        sparse output, or call a dense helper with `allow_densify=True` when the
        full selected result is small enough to hold in memory.
        """
        if self.ref_kmers.storage_mode != "dense":
            raise ValueError(
                "dense_frequencies_zarr_array() is only available for dense "
                "reference k-mer output"
            )
        return _required(self.ref_kmers.frequencies, "frequencies")

    def _sparse_frequencies_matrix_all(self) -> sparse.coo_matrix:
        """
        Return all reference k-mer frequencies as a SciPy COO matrix.
        """
        if self.ref_kmers.storage_mode == "dense":
            frequencies = np.asarray(
                _required(self.ref_kmers.frequencies, "frequencies")[:]
            )
            _validate_frequency_values(frequencies, "frequencies")
            return sparse.coo_matrix(frequencies)

        row_index = _required(self.ref_kmers.sparse_row, "sparse/row")
        motif_index = _required(self.ref_kmers.sparse_motif, "sparse/motif")
        frequency = _required(self.ref_kmers.sparse_frequency, "sparse/frequency")
        shape = tuple(
            _required(self.ref_kmers.sparse_shape, "sparse/shape").astype(int)
        )
        return sparse.coo_matrix(
            (
                frequency,
                (
                    row_index.astype(np.int64, copy=False),
                    motif_index.astype(np.int64, copy=False),
                ),
            ),
            shape=shape,
        )

    def _sparse_counts_matrix_all(self) -> sparse.coo_matrix:
        """
        Return all reconstructed reference k-mer counts as a SciPy COO matrix.
        """
        if self.ref_kmers.storage_mode == "dense":
            frequencies = np.asarray(
                _required(self.ref_kmers.frequencies, "frequencies")[:]
            )
            _validate_frequency_values(frequencies, "frequencies")
            counts = frequencies * self.ref_kmers.row_scaling_factor[:, np.newaxis]
            return sparse.coo_matrix(counts)

        row_index = _required(self.ref_kmers.sparse_row, "sparse/row")
        motif_index = _required(self.ref_kmers.sparse_motif, "sparse/motif")
        frequency = _required(self.ref_kmers.sparse_frequency, "sparse/frequency")
        shape = tuple(
            _required(self.ref_kmers.sparse_shape, "sparse/shape").astype(int)
        )
        count = frequency * self.ref_kmers.row_scaling_factor[row_index.astype(int)]
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

    def _sparse_frequencies_matrix_for_indices(
        self, row_indices: np.ndarray, motif_indices: np.ndarray
    ) -> sparse.coo_matrix:
        """
        Return selected rows and motifs as a frequency COO matrix.
        """
        if len(row_indices) == 0 or len(motif_indices) == 0:
            return sparse.coo_matrix((len(row_indices), len(motif_indices)))
        return (
            self._sparse_frequencies_matrix_all()
            .tocsr()[row_indices, :][:, motif_indices]
            .tocoo()
        )

    def _sparse_counts_matrix_for_indices(
        self, row_indices: np.ndarray, motif_indices: np.ndarray
    ) -> sparse.coo_matrix:
        """
        Return selected rows and motifs as a reconstructed count COO matrix.
        """
        if len(row_indices) == 0 or len(motif_indices) == 0:
            return sparse.coo_matrix((len(row_indices), len(motif_indices)))
        return (
            self._sparse_counts_matrix_all()
            .tocsr()[row_indices, :][:, motif_indices]
            .tocoo()
        )

    def _dense_frequencies_array_for_indices(
        self,
        row_indices: np.ndarray,
        motif_indices: np.ndarray,
        *,
        allow_densify: bool,
    ) -> np.ndarray:
        """
        Return selected rows and motifs as a dense frequency array.
        """
        allow_densify = validate_scalar_bool(allow_densify, "allow_densify")
        if self.ref_kmers.storage_mode == "dense":
            frequencies = _required(self.ref_kmers.frequencies, "frequencies")
            selected = np.asarray(
                frequencies.get_orthogonal_selection((row_indices, motif_indices))
            )
            _validate_frequency_values(selected, "selected frequencies")
            return selected

        _require_densify(allow_densify, "dense_frequencies_array")
        return self._sparse_frequencies_matrix_for_indices(
            row_indices, motif_indices
        ).toarray()

    def _dense_counts_array_for_indices(
        self,
        row_indices: np.ndarray,
        motif_indices: np.ndarray,
        *,
        allow_densify: bool,
    ) -> np.ndarray:
        """
        Return selected rows and motifs as a dense reconstructed count array.
        """
        frequencies = self._dense_frequencies_array_for_indices(
            row_indices,
            motif_indices,
            allow_densify=allow_densify,
        )
        row_scaling_factor = self.ref_kmers.row_scaling_factor[row_indices]
        return frequencies * row_scaling_factor[:, np.newaxis]

    def _sparse_frequencies_matrix(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Shared implementation behind mode-specific sparse frequency methods.
        """
        row_indices = self._resolve_row_selector(window_idxs, groups, group_idxs)
        motif_indices = self._resolve_motif_selector(motifs, motif_idxs)
        return self._sparse_frequencies_matrix_for_indices(row_indices, motif_indices)

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

    def _dense_frequencies_array(
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
        Shared implementation behind mode-specific dense frequency methods.
        """
        row_indices = self._resolve_row_selector(window_idxs, groups, group_idxs)
        motif_indices = self._resolve_motif_selector(motifs, motif_idxs)
        return self._dense_frequencies_array_for_indices(
            row_indices,
            motif_indices,
            allow_densify=allow_densify,
        )

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
        if self.ref_kmers.storage_mode == "sparse_coo" and not densify:
            return self._stored_data_frame_for_indices(row_indices, motif_indices)
        return self._complete_data_frame_for_indices(
            row_indices, motif_indices, densify
        )

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

        selected_frequencies = self._dense_frequencies_array_for_indices(
            row_indices,
            motif_indices,
            allow_densify=densify,
        )
        selected_counts = (
            selected_frequencies
            * self.ref_kmers.row_scaling_factor[row_indices][:, np.newaxis]
        )
        repeated_rows = row_metadata.loc[
            row_metadata.index.repeat(motif_count)
        ].reset_index(drop=True)
        repeated_motifs = pd.concat([motif_metadata] * row_count, ignore_index=True)
        data_frame = pd.concat([repeated_rows, repeated_motifs], axis=1)
        data_frame["frequency"] = selected_frequencies.ravel()
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

        sparse_rows = _required(self.ref_kmers.sparse_row, "sparse/row").astype(int)
        sparse_motifs = _required(self.ref_kmers.sparse_motif, "sparse/motif").astype(
            int
        )
        sparse_frequencies = _required(
            self.ref_kmers.sparse_frequency, "sparse/frequency"
        )
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
        matched_frequencies = sparse_frequencies[matches][sort_order]
        matched_counts = (
            matched_frequencies * self.ref_kmers.row_scaling_factor[matched_rows]
        )

        row_metadata = (
            self._row_metadata_data_frame().iloc[matched_rows].reset_index(drop=True)
        )
        motif_metadata = (
            self._motif_axis_metadata_data_frame()
            .iloc[matched_motifs]
            .reset_index(drop=True)
        )
        data_frame = pd.concat([row_metadata, motif_metadata], axis=1)
        data_frame["frequency"] = matched_frequencies
        data_frame["count"] = matched_counts
        return data_frame

    def _empty_data_frame(self) -> pd.DataFrame:
        """
        Return an empty reference k-mer data frame with public columns.
        """
        row_metadata = self._row_metadata_data_frame().iloc[[]].reset_index(drop=True)
        motif_metadata = (
            self._motif_axis_metadata_data_frame().iloc[[]].reset_index(drop=True)
        )
        data_frame = pd.concat([row_metadata, motif_metadata], axis=1)
        data_frame["frequency"] = np.asarray([], dtype=float)
        data_frame["count"] = np.asarray([], dtype=float)
        return data_frame

    def _motif_axis_metadata_data_frame(self) -> pd.DataFrame:
        """
        Build metadata for the active motif axis.
        """
        return self.motifs_metadata()

    def _row_metadata_data_frame(self) -> pd.DataFrame:
        """
        Build mode-specific row metadata for the active reference k-mer row mode.
        """
        if self.ref_kmers.row_mode == "global":
            data_frame = pd.DataFrame()
            data_frame["row_label"] = self.ref_kmers.row_labels
        elif self.ref_kmers.row_mode in {"size", "bed"}:
            data_frame = pd.DataFrame()
            row_chromosome = _required(self.ref_kmers.row_chromosome, "row_chromosome")
            chromosome_names = _required(
                self.ref_kmers.chromosome_names, "chromosome_names"
            )
            data_frame["window_idx"] = self.ref_kmers.row
            data_frame["chrom"] = chromosome_names[row_chromosome.astype(int)]
            data_frame["start"] = _required(self.ref_kmers.row_start_bp, "row_start_bp")
            data_frame["end"] = _required(self.ref_kmers.row_end_bp, "row_end_bp")
            data_frame["blacklisted_fraction"] = _required(
                self.ref_kmers.blacklisted_fraction, "blacklisted_fraction"
            )
        elif self.ref_kmers.row_mode == "grouped_bed":
            data_frame = pd.DataFrame()
            data_frame["group_idx"] = _required(self.ref_kmers.group_idx, "group")
            data_frame["group_name"] = _required(
                self.ref_kmers.group_names, "group_name"
            )
            data_frame["eligible_windows"] = _required(
                self.ref_kmers.eligible_windows, "eligible_windows"
            )
            data_frame["blacklisted_fraction"] = _required(
                self.ref_kmers.blacklisted_fraction, "blacklisted_fraction"
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
        if self.ref_kmers.row_mode in {"size", "bed"}:
            if groups is not None or group_idxs is not None:
                raise ValueError(
                    "Grouped selectors can only be used with grouped output"
                )
            return normalize_zero_based_indices(
                window_idxs,
                size=len(self.ref_kmers.row),
                name="window_idxs",
                index_name="window_idx",
            )
        if self.ref_kmers.row_mode == "grouped_bed":
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
                size=len(self.ref_kmers.row),
                name="group_idxs",
                index_name="group_idx",
            )
        if window_idxs is not None or groups is not None or group_idxs is not None:
            raise ValueError("Row selectors cannot be used with global output")
        return np.arange(len(self.ref_kmers.row), dtype=np.int64)

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
            return np.arange(len(self.ref_kmers.motif_index), dtype=np.int64)
        if motifs is not None:
            motif_names = normalize_strings(motifs, name="motifs")
            return np.asarray(
                [self._resolve_motif(motif) for motif in motif_names],
                dtype=np.int64,
            )
        return normalize_zero_based_indices(
            motif_idxs,
            size=len(self.ref_kmers.motif_index),
            name="motif_idxs",
            index_name="motif_idx",
        )

    def _resolve_motif(self, motif: str) -> int:
        """
        Resolve a motif label to its motif index.
        """
        return resolve_unique_match(
            self.ref_kmers.motif_names == motif,
            missing_message=f"Unknown reference k-mer motif label: {motif!r}",
            duplicate_message=f"Reference k-mer motif label is not unique: {motif!r}",
        )


class GlobalRefKmerFrequencies(RefKmerFrequencies):
    """Reference k-mer frequencies for global output."""

    def data_frame(
        self,
        *,
        densify: bool = False,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """
        Create a pandas data frame for global reference k-mer frequencies.

        The data frame contains one row per selected motif with `frequency` and
        reconstructed `count`. Sparse output returns only non-zero stored values
        unless `densify=True`, which adds zero-frequency rows for the selected
        motifs in `motifs_metadata()`. For observed-only output, those selected
        labels are the combined set observed anywhere in the output. Densifying
        does not add every possible k-mer unless `all_motifs()` is true.
        """
        return self._data_frame(
            densify=densify,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def dense_frequencies_array(
        self,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return global reference k-mer frequencies as a dense NumPy array.

        Sparse output requires `allow_densify=True` because the method creates
        an in-memory array with one value for every selected motif in
        `motifs_metadata()`, including zeroes that were not stored on disk.
        """
        return self._dense_frequencies_array(
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def dense_counts_array(
        self,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return reconstructed global reference k-mer counts as a dense array.

        Counts are reconstructed from frequencies with the global row scaling
        factor. Sparse output requires `allow_densify=True` because the method
        creates an in-memory array with explicit zeroes for selected motifs in
        `motifs_metadata()`.
        """
        return self._dense_counts_array(
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def sparse_frequencies_matrix(
        self,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return global reference k-mer frequencies as a SciPy sparse matrix.
        """
        return self._sparse_frequencies_matrix(
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def sparse_counts_matrix(
        self,
        *,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return reconstructed global reference k-mer counts as a SciPy sparse matrix.
        """
        return self._sparse_counts_matrix(
            motifs=motifs,
            motif_idxs=motif_idxs,
        )


class WindowedRefKmerFrequencies(RefKmerFrequencies):
    """Reference k-mer frequencies for fixed-size or BED-window output."""

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
        Create a pandas data frame for windowed reference k-mer frequencies.

        `window_idxs` and `motif_idxs` are zero-based indices. The data frame
        includes window metadata, motif metadata, `frequency`, and reconstructed
        `count`.

        Sparse output returns only non-zero stored values unless `densify=True`,
        which adds zero-frequency rows for the selected windows and the
        selected motifs in `motifs_metadata()`. For observed-only output, those
        selected labels are the combined set observed anywhere in the output.
        Densifying does not add every possible k-mer unless `all_motifs()` is
        true.
        """
        return self._data_frame(
            densify=densify,
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
        )

    def window_metadata(self) -> pd.DataFrame:
        """
        Return genomic window metadata for this reference k-mer output.

        `window_idx` is the zero-based row index accepted by window selectors.
        `start` and `end` are half-open genomic coordinates.
        """
        return self._row_metadata_data_frame()

    def dense_frequencies_array(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return windowed reference k-mer frequencies as a dense NumPy array.

        Sparse output requires `allow_densify=True` because the method creates
        an in-memory array with one value for every selected window and every
        selected motif in `motifs_metadata()`, including zeroes that were not
        stored on disk.
        """
        return self._dense_frequencies_array(
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def dense_counts_array(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return reconstructed windowed reference k-mer counts as a dense array.

        Counts are reconstructed row-wise from frequencies and
        `row_scaling_factor`. Sparse output requires `allow_densify=True`
        because the method creates an in-memory array with explicit zeroes for
        selected motifs in `motifs_metadata()`.
        """
        return self._dense_counts_array(
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def sparse_frequencies_matrix(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return windowed reference k-mer frequencies as a SciPy sparse matrix.
        """
        return self._sparse_frequencies_matrix(
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )

    def sparse_counts_matrix(
        self,
        *,
        window_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return reconstructed windowed reference k-mer counts as a SciPy sparse matrix.
        """
        return self._sparse_counts_matrix(
            window_idxs=window_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )


class GroupedRefKmerFrequencies(RefKmerFrequencies):
    """Reference k-mer frequencies for grouped BED output."""

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
        Create a pandas data frame for grouped reference k-mer frequencies.

        `group_idxs` and `motif_idxs` are zero-based indices. `groups` selects
        rows by group name. The data frame includes group metadata, motif
        metadata, `frequency`, and reconstructed `count`.

        Sparse output returns only non-zero stored values unless `densify=True`,
        which adds zero-frequency rows for the selected groups and the selected
        motifs in `motifs_metadata()`. For observed-only output, those selected
        labels are the combined set observed anywhere in the output. Densifying
        does not add every possible k-mer unless `all_motifs()` is true.
        """
        return self._data_frame(
            densify=densify,
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            max_blacklisted_fraction=max_blacklisted_fraction,
        )

    def group_metadata(self) -> pd.DataFrame:
        """
        Return grouped BED metadata for this reference k-mer output.

        `group_idx` is the zero-based row index accepted by group selectors.
        """
        return self._row_metadata_data_frame()

    def group_idx(self, group_name: str) -> int:
        """
        Find the zero-based reference k-mer row index for a group name.
        """
        group_names = _required(self.ref_kmers.group_names, "group_name")
        return resolve_unique_match(
            group_names == group_name,
            missing_message=f"Unknown reference k-mer group name: {group_name!r}",
            duplicate_message=(
                f"Reference k-mer group name is not unique: {group_name!r}"
            ),
        )

    def dense_frequencies_array(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
        allow_densify: bool = False,
    ) -> np.ndarray:
        """
        Return grouped reference k-mer frequencies as a dense NumPy array.

        Sparse output requires `allow_densify=True` because the method creates
        an in-memory array with one value for every selected group and every
        selected motif in `motifs_metadata()`, including zeroes that were not
        stored on disk.
        """
        return self._dense_frequencies_array(
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
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
        Return reconstructed grouped reference k-mer counts as a dense array.

        Counts are reconstructed row-wise from frequencies and
        `row_scaling_factor`. Sparse output requires `allow_densify=True`
        because the method creates an in-memory array with explicit zeroes for
        selected motifs in `motifs_metadata()`.
        """
        return self._dense_counts_array(
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
            allow_densify=allow_densify,
        )

    def sparse_frequencies_matrix(
        self,
        *,
        groups: str | Sequence[str] | None = None,
        group_idxs: int | Sequence[int] | None = None,
        motifs: str | Sequence[str] | None = None,
        motif_idxs: int | Sequence[int] | None = None,
    ) -> sparse.coo_matrix:
        """
        Return grouped reference k-mer frequencies as a SciPy sparse matrix.
        """
        return self._sparse_frequencies_matrix(
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
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
        Return reconstructed grouped reference k-mer counts as a SciPy sparse matrix.
        """
        return self._sparse_counts_matrix(
            groups=groups,
            group_idxs=group_idxs,
            motifs=motifs,
            motif_idxs=motif_idxs,
        )


def read_ref_kmers(
    path: pathlib.Path | str,
) -> GlobalRefKmerFrequencies | WindowedRefKmerFrequencies | GroupedRefKmerFrequencies:
    """
    Open a cfDNAlab reference k-mer frequency output directory.
    """
    path = pathlib.Path(path)
    loaded = RefKmerFrequencies._load_zarr(path)
    if loaded.row_mode == "global":
        return GlobalRefKmerFrequencies(path, loaded)
    if loaded.row_mode in {"size", "bed"}:
        return WindowedRefKmerFrequencies(path, loaded)
    if loaded.row_mode == "grouped_bed":
        return GroupedRefKmerFrequencies(path, loaded)
    raise ValueError(f"Unsupported reference k-mer row mode: {loaded.row_mode!r}")


def _validate_zarr_store_path(path: pathlib.Path) -> None:
    """
    Validate that a path points to a Zarr V3 reference k-mer store directory.
    """
    if not path.exists():
        raise FileNotFoundError(f"Reference k-mer Zarr store does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(
            f"Reference k-mer Zarr store path exists but is not a directory: {path}"
        )
    if path.suffix != ".zarr":
        raise ValueError(
            f"Reference k-mer Zarr store path should end with '.zarr', got: {path}"
        )
    if not (path / "zarr.json").is_file():
        raise ValueError(
            "Reference k-mer Zarr store is missing Zarr V3 metadata file: "
            f"{path / 'zarr.json'}"
        )


def _validate_root_metadata(store: Any) -> dict[str, Any]:
    """
    Validate root attributes and return metadata used by the API.
    """
    schema = store.attrs.get("cfdnalab_schema")
    if schema != "ref_kmer_frequencies":
        raise ValueError(
            f"Expected cfdnalab_schema='ref_kmer_frequencies', found {schema!r}"
        )

    schema_version = store.attrs.get("cfdnalab_schema_version")
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, numbers.Integral)
        or int(schema_version) != REF_KMER_SUPPORTED_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unsupported reference k-mer schema version: "
            f"{schema_version!r}. Supported version: "
            f"{REF_KMER_SUPPORTED_SCHEMA_VERSION}"
        )

    storage_mode = store.attrs.get("storage_mode")
    if storage_mode not in VALID_STORAGE_MODES:
        raise ValueError(f"Unsupported reference k-mer storage mode: {storage_mode!r}")

    row_mode = store.attrs.get("row_mode")
    if row_mode not in VALID_ROW_MODES:
        raise ValueError(f"Unsupported reference k-mer row mode: {row_mode!r}")

    motif_axis_kind = store.attrs.get("motif_axis_kind")
    if motif_axis_kind not in VALID_MOTIF_AXIS_KINDS:
        raise ValueError(
            f"Unsupported reference k-mer motif axis kind: {motif_axis_kind!r}"
        )

    _validate_root_string_attr(store, "value_units", EXPECTED_VALUE_UNITS)
    _validate_root_string_attr(store, "count_units", EXPECTED_COUNT_UNITS)
    _validate_root_string_attr(
        store,
        "row_scaling_factor_array",
        EXPECTED_ROW_SCALING_ARRAY,
    )
    _validate_root_string_attr(
        store,
        "count_reconstruction",
        EXPECTED_COUNT_RECONSTRUCTION,
    )
    if storage_mode == "dense":
        _validate_root_string_attr(store, "primary_array", "frequencies")
    else:
        _validate_root_string_attr(store, "primary_group", "sparse")
        _validate_root_string_attr(store, "sparse_format", "coo")
        sparse_indices_base = store.attrs.get("sparse_indices_base")
        if (
            isinstance(sparse_indices_base, bool)
            or not isinstance(sparse_indices_base, numbers.Integral)
            or int(sparse_indices_base) != 0
        ):
            raise ValueError(
                "sparse_indices_base must be 0 for reference k-mer sparse COO output"
            )

    kmer_size = store.attrs.get("kmer_size")
    if (
        isinstance(kmer_size, bool)
        or not isinstance(kmer_size, numbers.Integral)
        or int(kmer_size) <= 0
    ):
        raise ValueError(f"kmer_size must be a positive integer, found {kmer_size!r}")

    canonical = store.attrs.get("canonical")
    if not isinstance(canonical, bool):
        raise ValueError(f"canonical must be a boolean, found {canonical!r}")

    all_motifs = store.attrs.get("all_motifs")
    if not isinstance(all_motifs, bool):
        raise ValueError(f"all_motifs must be a boolean, found {all_motifs!r}")

    assign_by = store.attrs.get("assign_by")
    if not isinstance(assign_by, str) or not assign_by:
        raise ValueError(f"assign_by must be a non-empty string, found {assign_by!r}")

    return {
        "storage_mode": storage_mode,
        "row_mode": row_mode,
        "motif_axis_kind": motif_axis_kind,
        "kmer_size": int(kmer_size),
        "canonical": canonical,
        "all_motifs": all_motifs,
        "assign_by": assign_by,
    }


def _validate_root_string_attr(store: Any, name: str, expected: str) -> None:
    """
    Require a root string attribute to match the schema value.
    """
    value = store.attrs.get(name)
    if value != expected:
        raise ValueError(
            f"Reference k-mer root attribute {name!r} must be {expected!r}, "
            f"found {value!r}"
        )


def _validate_required_arrays(
    store: Any, storage_mode: str, row_mode: str, motif_axis_kind: str
) -> None:
    """
    Require all arrays needed by the active storage mode and row mode.
    """
    required = {
        "motif_index",
        "row",
        "row_scaling_factor",
        "reference_contig_footprint_json",
    }
    if motif_axis_kind == "motif":
        required.update({"motif_byte", "motif_ascii"})
    if storage_mode == "dense":
        required.add("frequencies")
    else:
        required.update(
            {
                "sparse/row",
                "sparse/motif",
                "sparse/frequency",
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
        raise ValueError(f"Reference k-mer Zarr store is missing arrays: {missing}")


def _has_array(store: Any, name: str) -> bool:
    """
    Return whether a Zarr group contains an array path.
    """
    try:
        store[name]
    except Exception:
        return False
    return True


def _validate_dense_frequencies(frequencies: Any, n_rows: int, n_motifs: int) -> None:
    """
    Validate dense reference k-mer frequency array dimensions and shape.
    """
    _validate_array_dimensions(frequencies, ("row", "motif"), "dense frequencies")
    if tuple(frequencies.shape) != (n_rows, n_motifs):
        raise ValueError(
            "dense frequencies shape does not match row and motif axes: "
            f"frequencies={frequencies.shape}, coordinates={(n_rows, n_motifs)}"
        )


def _validate_sparse_arrays(
    store: Any,
    row: np.ndarray,
    motif: np.ndarray,
    frequency: np.ndarray,
    shape: np.ndarray,
    n_rows: int,
    n_motifs: int,
) -> None:
    """
    Validate sparse COO arrays against row and motif axes.
    """
    _validate_array_dimensions(store["sparse/row"], ("nnz",), "sparse/row")
    _validate_array_dimensions(store["sparse/motif"], ("nnz",), "sparse/motif")
    _validate_array_dimensions(store["sparse/frequency"], ("nnz",), "sparse/frequency")
    _validate_array_dimensions(
        store["sparse/shape"], ("sparse_dimension",), "sparse/shape"
    )
    _validate_array_dimensions(
        store["sparse/sparse_dimension"],
        ("sparse_dimension",),
        "sparse/sparse_dimension",
    )
    _validate_same_length(row, motif, "sparse/row", "sparse/motif")
    _validate_same_length(row, frequency, "sparse/row", "sparse/frequency")
    _validate_integer_values(row, "sparse/row")
    _validate_integer_values(motif, "sparse/motif")
    _validate_integer_values(shape, "sparse/shape")
    _validate_frequency_values(frequency, "sparse/frequency")
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
        previous_coordinate: tuple[int, int] | None = None
        for row_index, motif_index in zip(row.astype(int), motif.astype(int)):
            coordinate = (int(row_index), int(motif_index))
            if previous_coordinate is not None and coordinate <= previous_coordinate:
                raise ValueError(
                    "sparse COO entries must be sorted and unique by row, motif"
                )
            previous_coordinate = coordinate

    sparse_dimension = _read_array(store, "sparse/sparse_dimension")
    _validate_axis(sparse_dimension, "sparse_dimension")
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


def _read_reference_contig_footprint(store: Any) -> Any:
    """
    Decode the JSON reference contig footprint.
    """
    array = store["reference_contig_footprint_json"]
    _validate_array_dimensions(
        array,
        ("json_byte",),
        "reference_contig_footprint_json",
    )
    raw = _read_array(store, "reference_contig_footprint_json")
    if raw.dtype != np.uint8:
        raise ValueError(
            f"reference_contig_footprint_json must have dtype uint8, found {raw.dtype}"
        )
    try:
        return json.loads(bytes(raw).decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(
            "reference_contig_footprint_json must contain UTF-8 JSON"
        ) from error


def _read_motif_ascii_labels(store: Any, expected_len: int) -> np.ndarray:
    """
    Decode fixed-width ASCII motif labels into Python strings.
    """
    _validate_array_dimensions(store["motif_byte"], ("motif_byte",), "motif_byte")
    motif_byte = _read_array(store, "motif_byte")
    _validate_axis(motif_byte, "motif_byte")

    _validate_array_dimensions(
        store["motif_ascii"], ("motif", "motif_byte"), "motif_ascii"
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
    _validate_integer_values(values, name)
    expected = np.arange(len(values), dtype=values.dtype)
    if not np.array_equal(values, expected):
        raise ValueError(f"{name} axis must be contiguous 0-based indices")


def _validate_integer_values(values: np.ndarray, name: str) -> None:
    """
    Require an array to contain integer values.
    """
    if not np.issubdtype(values.dtype, np.integer):
        raise ValueError(f"{name} must contain integer values")


def _validate_nonnegative_integer_values(values: np.ndarray, name: str) -> None:
    """
    Require an array to contain non-negative integer values.
    """
    _validate_integer_values(values, name)
    if np.any(values < 0):
        raise ValueError(f"{name} must contain non-negative integer values")


def _validate_index_values(values: np.ndarray, axis_length: int, name: str) -> None:
    """
    Require an array to contain valid zero-based indices into an axis.
    """
    _validate_nonnegative_integer_values(values, name)
    if np.any(values >= axis_length):
        raise ValueError(f"{name} contains an index outside the referenced axis")


def _validate_half_open_intervals(
    starts: np.ndarray, ends: np.ndarray, start_name: str, end_name: str
) -> None:
    """
    Require paired start and end arrays to contain non-empty half-open intervals.
    """
    _validate_nonnegative_integer_values(starts, start_name)
    _validate_nonnegative_integer_values(ends, end_name)
    _validate_same_length(starts, ends, start_name, end_name)
    if np.any(starts >= ends):
        raise ValueError(f"{start_name} must be smaller than {end_name}")


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


def _validate_frequency_values(values: np.ndarray, name: str) -> None:
    """
    Require frequency values to be finite fractions in 0..1.
    """
    if values.size == 0:
        return
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError(f"{name} must contain finite frequency values in 0..1")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values > 1.0):
        raise ValueError(f"{name} must contain finite frequency values in 0..1")


def _validate_fraction_values(values: np.ndarray, name: str) -> None:
    """
    Require values to be finite fractions in 0..1.
    """
    if values.size == 0:
        return
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError(f"{name} must contain finite fractions in 0..1")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values > 1.0):
        raise ValueError(f"{name} must contain finite fractions in 0..1")


def _validate_row_scaling_factors(values: np.ndarray) -> None:
    """
    Require row scaling factors to be finite non-negative values.
    """
    if values.size == 0:
        return
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError("row_scaling_factor must contain finite non-negative values")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("row_scaling_factor must contain finite non-negative values")


def _validate_concrete_motif_labels(
    motif_names: np.ndarray, kmer_size: int, canonical: bool
) -> None:
    """
    Validate concrete A/C/G/T motif labels against root metadata.
    """
    _validate_unique_labels(motif_names, "reference k-mer motif")
    for motif_index, motif in enumerate(motif_names):
        if len(motif) != kmer_size:
            raise ValueError(
                "reference k-mer motif label "
                f"{motif!r} at motif_index {motif_index} has length {len(motif)}, "
                f"expected {kmer_size}"
            )
        invalid_bases = sorted(set(motif) - {"A", "C", "G", "T"})
        if invalid_bases:
            raise ValueError(
                "reference k-mer motif label "
                f"{motif!r} at motif_index {motif_index} contains invalid "
                f"bases {invalid_bases!r}"
            )
        if canonical:
            canonical_motif = _canonical_ref_kmer(motif)
            if motif != canonical_motif:
                raise ValueError(
                    "canonical reference k-mer motif label "
                    f"{motif!r} at motif_index {motif_index} should be "
                    f"{canonical_motif!r}"
                )


def _validate_unique_labels(labels: np.ndarray, label_name: str) -> None:
    """
    Reject duplicate labels on a public axis.
    """
    if len(set(labels)) != len(labels):
        raise ValueError(f"duplicate {label_name} label")


def _contains_control_character(label: str) -> bool:
    """
    Return whether a label contains a Unicode control character.
    """
    return any(unicodedata.category(character) == "Cc" for character in label)


def _canonical_ref_kmer(motif: str) -> str:
    """
    Return the canonical label used by cfDNAlab reference k-mer output.
    """
    reverse_complement = _reverse_complement(motif)
    if len(motif) % 2 == 0:
        return min(motif, reverse_complement)
    middle_base = motif[len(motif) // 2]
    if middle_base in {"A", "C"}:
        return motif
    return reverse_complement


def _reverse_complement(motif: str) -> str:
    """
    Return the reverse complement of a concrete A/C/G/T motif.
    """
    return motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _require_densify(allow_densify: bool, method_name: str) -> None:
    """
    Require explicit opt-in before expanding sparse output to a dense array.
    """
    if not allow_densify:
        raise ValueError(
            f"{method_name}() would turn sparse reference k-mer output into "
            "a dense in-memory array. "
            "Use sparse_frequencies_matrix(), sparse_counts_matrix(), "
            "or pass allow_densify=True."
        )


def _required(
    value: np.ndarray | zarr.Array | None, name: str
) -> np.ndarray | zarr.Array:
    """
    Return a loaded value or fail with schema context when it is unavailable.
    """
    if value is None:
        raise ValueError(f"Reference k-mer Zarr store does not contain {name}")
    return value
