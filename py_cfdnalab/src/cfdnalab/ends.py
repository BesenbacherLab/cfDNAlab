"""
Classes for loading and interacting with the ends .zarr output.
"""

from __future__ import annotations

from dataclasses import dataclass
import numbers
import pathlib
from typing import Any, List

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import zarr

SUPPORTED_SCHEMA_VERSION = 1
VALID_STORAGE_MODES = {"dense", "sparse_coo"}
VALID_ROW_MODES = {"global", "size", "bed", "grouped_bed"}


@dataclass
class LoadedEndMotifs:
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
    """Base class for end-motif count helpers."""

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
            Preloaded store data used by `load_end_motifs`.

        Returns
        -------
        None
        """
        self.path = pathlib.Path(path)
        if loaded_end_motifs is None:
            loaded_end_motifs = EndMotifCounts._load_zarr(self.path)
        self.end_motifs = loaded_end_motifs

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
            raise ValueError(f"Could not open end-motif Zarr store at {path}") from error

        storage_mode, row_mode = _validate_root_metadata(store)
        _validate_required_arrays(store, storage_mode, row_mode)

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
            counts = store["counts"]
            _validate_dense_counts(counts, len(row), len(motif_index))
        else:
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
                raise ValueError("row_chromosome contains an index outside the chromosome axis")
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
        Return the count storage mode.

        Returns
        -------
        str
            Either `"dense"` or `"sparse_coo"`.
        """
        return self.end_motifs.storage_mode

    def row_mode(self) -> str:
        """
        Return the row metadata mode.

        Returns
        -------
        str
            One of `"global"`, `"size"`, `"bed"`, or `"grouped_bed"`.
        """
        return self.end_motifs.row_mode

    def motifs(self) -> List[str]:
        """
        Return motif labels.

        Returns
        -------
        list[str]
            Motif labels in count-column order.
        """
        return self.end_motifs.motif_names.tolist()

    def motif_idx(self, motif: str) -> int:
        """
        Return the index for a motif label.

        Parameters
        ----------
        motif
            Motif label to resolve.

        Returns
        -------
        int
            Zero-based motif index.
        """
        return self._resolve_motif(motif)

    def motif_metadata(self) -> pd.DataFrame:
        """
        Return motif metadata.

        Returns
        -------
        pandas.DataFrame
            One row per motif.
        """
        return pd.DataFrame(
            {
                "motif_index": self.end_motifs.motif_index,
                "motif": self.end_motifs.motif_names,
            }
        )

    def sparse_coo(self) -> sparse.coo_matrix:
        """
        Return counts as a SciPy COO sparse matrix.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(row, motif)`.
        """
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return sparse.coo_matrix(np.asarray(counts[:]))

        row_index = _required(self.end_motifs.sparse_row, "sparse/row")
        motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
        count = _required(self.end_motifs.sparse_count, "sparse/count")
        shape = tuple(_required(self.end_motifs.sparse_shape, "sparse/shape").astype(int))
        return sparse.coo_matrix((count, (row_index, motif_index)), shape=shape)

    def sparse_coo_data_frame(self) -> pd.DataFrame:
        """
        Return non-zero COO entries as a data frame.

        Returns
        -------
        pandas.DataFrame
            One row per stored non-zero count.
        """
        if self.end_motifs.storage_mode == "sparse_coo":
            row_index = _required(self.end_motifs.sparse_row, "sparse/row")
            motif_index = _required(self.end_motifs.sparse_motif, "sparse/motif")
            count = _required(self.end_motifs.sparse_count, "sparse/count")
        else:
            coo = self.sparse_coo()
            row_index = coo.row
            motif_index = coo.col
            count = coo.data
        motif_lookup_index = motif_index.astype(int)
        return pd.DataFrame(
            {
                "row": row_index,
                "motif_index": motif_index,
                "motif": self.end_motifs.motif_names[motif_lookup_index],
                "count": count,
            }
        )

    def sparse_coo_for_motif(self, motif: str) -> sparse.coo_matrix:
        """
        Return a sparse column matrix for one motif.

        Parameters
        ----------
        motif
            Motif label to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(row, 1)`.
        """
        return self.sparse_coo_for_motif_idx(self._resolve_motif(motif))

    def sparse_coo_for_motif_idx(self, motif_idx: int) -> sparse.coo_matrix:
        """
        Return a sparse column matrix for one motif index.

        Parameters
        ----------
        motif_idx
            Zero-based motif index to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(row, 1)`.
        """
        motif_idx = self._validate_motif_idx(motif_idx)
        return self.sparse_coo().tocsr()[:, motif_idx].tocoo()

    def dense_array(self) -> np.ndarray:
        """
        Load or reconstruct the full dense count matrix.

        Returns
        -------
        numpy.ndarray
            Count array with shape `(row, motif)`.
        """
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[:])
        return self.sparse_coo().toarray()

    def dense_array_for_motif(self, motif: str) -> np.ndarray:
        """
        Return dense counts for one motif.

        Parameters
        ----------
        motif
            Motif label to extract.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(row,)`.
        """
        return self.dense_array_for_motif_idx(self._resolve_motif(motif))

    def dense_array_for_motif_idx(self, motif_idx: int) -> np.ndarray:
        """
        Return dense counts for one motif index.

        Parameters
        ----------
        motif_idx
            Zero-based motif index to extract.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(row,)`.
        """
        motif_idx = self._validate_motif_idx(motif_idx)
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[:, motif_idx])

        values = np.zeros(len(self.end_motifs.row), dtype=float)
        sparse_row = _required(self.end_motifs.sparse_row, "sparse/row")
        sparse_motif = _required(self.end_motifs.sparse_motif, "sparse/motif")
        sparse_count = _required(self.end_motifs.sparse_count, "sparse/count")
        matches = sparse_motif == motif_idx
        values[sparse_row[matches].astype(int)] = sparse_count[matches]
        return values

    def dense_data_frame_for_motif(self, motif: str) -> pd.DataFrame:
        """
        Build a dense data frame for one motif across all output rows.

        Parameters
        ----------
        motif
            Motif label to extract.

        Returns
        -------
        pandas.DataFrame
            Row metadata and counts for one motif.
        """
        return self.dense_data_frame_for_motif_idx(self._resolve_motif(motif))

    def dense_data_frame_for_motif_idx(self, motif_idx: int) -> pd.DataFrame:
        """
        Build a dense data frame for one motif index across all output rows.

        Parameters
        ----------
        motif_idx
            Zero-based motif index to extract.

        Returns
        -------
        pandas.DataFrame
            Row metadata and counts for one motif.
        """
        motif_idx = self._validate_motif_idx(motif_idx)
        frame = self._row_metadata_frame()
        frame["motif_index"] = int(self.end_motifs.motif_index[motif_idx])
        frame["motif"] = self.end_motifs.motif_names[motif_idx]
        frame["count"] = self.dense_array_for_motif_idx(motif_idx)
        return frame

    def _row_metadata_frame(self) -> pd.DataFrame:
        if self.end_motifs.row_mode == "global":
            frame = pd.DataFrame()
            frame["row_label"] = self.end_motifs.row_labels
        elif self.end_motifs.row_mode in {"size", "bed"}:
            frame = pd.DataFrame()
            row_chromosome = _required(self.end_motifs.row_chromosome, "row_chromosome")
            chromosome_names = _required(self.end_motifs.chromosome_names, "chromosome_names")
            frame["window_idx"] = self.end_motifs.row
            frame["chromosome"] = row_chromosome
            frame["chromosome_name"] = chromosome_names[row_chromosome.astype(int)]
            frame["window_start_bp"] = _required(
                self.end_motifs.row_start_bp, "row_start_bp"
            )
            frame["window_end_bp"] = _required(self.end_motifs.row_end_bp, "row_end_bp")
            frame["blacklisted_fraction"] = _required(
                self.end_motifs.blacklisted_fraction, "blacklisted_fraction"
            )
        elif self.end_motifs.row_mode == "grouped_bed":
            frame = pd.DataFrame()
            frame["group_idx"] = _required(self.end_motifs.group_idx, "group")
            frame["group_name"] = _required(self.end_motifs.group_names, "group_name")
            frame["eligible_windows"] = _required(
                self.end_motifs.eligible_windows, "eligible_windows"
            )
            frame["blacklisted_fraction"] = _required(
                self.end_motifs.blacklisted_fraction, "blacklisted_fraction"
            )
        else:
            frame = pd.DataFrame()
        return frame

    def _resolve_motif(self, motif: str) -> int:
        matches = np.flatnonzero(self.end_motifs.motif_names == motif)
        if len(matches) == 0:
            raise KeyError(f"Unknown end-motif label: {motif!r}")
        if len(matches) > 1:
            raise ValueError(f"End-motif label is not unique: {motif!r}")
        return int(matches[0])

    def _validate_row(self, row: int) -> int:
        return _validate_index(row, len(self.end_motifs.row), "row")

    def _validate_motif_idx(self, motif_idx: int) -> int:
        return _validate_index(motif_idx, len(self.end_motifs.motif_index), "motif_idx")


class GlobalEndMotifCounts(EndMotifCounts):
    """End-motif counts for global output."""

    def counts(self) -> np.ndarray:
        """
        Return the global dense count vector.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(motif,)`.
        """
        return self.dense_array()[0, :]

    def data_frame(self) -> pd.DataFrame:
        """
        Build a dense motif count data frame for global output.

        Returns
        -------
        pandas.DataFrame
            One row per motif.
        """
        frame = self.motif_metadata()
        frame["count"] = self.counts()
        return frame


class WindowedEndMotifCounts(EndMotifCounts):
    """End-motif counts for fixed-size or BED-window output."""

    def windows(self) -> pd.DataFrame:
        """
        Return window metadata.

        Returns
        -------
        pandas.DataFrame
            One row per output window.
        """
        return self._row_metadata_frame()

    def sparse_coo_for_window(self, window_idx: int) -> sparse.coo_matrix:
        """
        Return a sparse row matrix for one window.

        Parameters
        ----------
        window_idx
            Zero-based window index to extract.

        Returns
        -------
        scipy.sparse.coo_matrix
            Sparse matrix with shape `(1, motif)`.
        """
        window_idx = self._validate_row(window_idx)
        return self.sparse_coo().tocsr()[window_idx, :].tocoo()

    def dense_array_for_window(self, window_idx: int) -> np.ndarray:
        """
        Return dense counts for one window.

        Parameters
        ----------
        window_idx
            Zero-based window index to extract.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(motif,)`.
        """
        window_idx = self._validate_row(window_idx)
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[window_idx, :])
        return self.sparse_coo_for_window(window_idx).toarray()[0, :]

    def dense_data_frame_for_window(self, window_idx: int) -> pd.DataFrame:
        """
        Build a dense motif count data frame for one window.

        Parameters
        ----------
        window_idx
            Zero-based window index to extract.

        Returns
        -------
        pandas.DataFrame
            One row per motif.
        """
        window_idx = self._validate_row(window_idx)
        frame = self.motif_metadata()
        frame["count"] = self.dense_array_for_window(window_idx)
        window_metadata = self.windows().iloc[window_idx].to_dict()
        for name, value in window_metadata.items():
            frame[name] = value
        return frame


class GroupedEndMotifCounts(EndMotifCounts):
    """End-motif counts for grouped BED output."""

    def groups(self) -> pd.DataFrame:
        """
        Return group metadata.

        Returns
        -------
        pandas.DataFrame
            One row per group.
        """
        return self._row_metadata_frame()

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
        group_names = _required(self.end_motifs.group_names, "group_name")
        matches = np.flatnonzero(group_names == group_name)
        if len(matches) == 0:
            raise KeyError(f"Unknown end-motif group name: {group_name!r}")
        if len(matches) > 1:
            raise ValueError(f"End-motif group name is not unique: {group_name!r}")
        return int(matches[0])

    def sparse_coo_for_group(self, group: int | str) -> sparse.coo_matrix:
        """
        Return a sparse row matrix for one group.

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
        return self.sparse_coo().tocsr()[group_idx, :].tocoo()

    def dense_array_for_group(self, group: int | str) -> np.ndarray:
        """
        Return dense counts for one group.

        Parameters
        ----------
        group
            Group index or group name to extract.

        Returns
        -------
        numpy.ndarray
            Dense count vector with shape `(motif,)`.
        """
        group_idx = self._resolve_group(group)
        if self.end_motifs.storage_mode == "dense":
            counts = _required(self.end_motifs.counts, "counts")
            return np.asarray(counts[group_idx, :])
        return self.sparse_coo_for_group(group_idx).toarray()[0, :]

    def dense_data_frame_for_group(self, group: int | str) -> pd.DataFrame:
        """
        Build a dense motif count data frame for one group.

        Parameters
        ----------
        group
            Group index or group name to extract.

        Returns
        -------
        pandas.DataFrame
            One row per motif.
        """
        group_idx = self._resolve_group(group)
        frame = self.motif_metadata()
        frame["count"] = self.dense_array_for_group(group_idx)
        group_metadata = self.groups().iloc[group_idx].to_dict()
        for name, value in group_metadata.items():
            frame[name] = value
        return frame

    def _resolve_group(self, group: int | str) -> int:
        if isinstance(group, str):
            return self.group_idx(group)
        return _validate_index(group, len(self.end_motifs.row), "group")


def load_end_motifs(path: pathlib.Path | str) -> EndMotifCounts:
    """
    Load an end-motif count Zarr store.

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
    schema = store.attrs.get("cfdnalab_schema")
    if schema != "end_motif_counts":
        raise ValueError(
            f"Expected cfdnalab_schema='end_motif_counts', found {schema!r}"
        )

    schema_version = store.attrs.get("cfdnalab_schema_version")
    if schema_version != SUPPORTED_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported end-motif schema version: "
            f"{schema_version!r}. Supported version: {SUPPORTED_SCHEMA_VERSION}"
        )

    storage_mode = store.attrs.get("storage_mode")
    if storage_mode not in VALID_STORAGE_MODES:
        raise ValueError(f"Unsupported end-motif storage mode: {storage_mode!r}")

    row_mode = store.attrs.get("row_mode")
    if row_mode not in VALID_ROW_MODES:
        raise ValueError(f"Unsupported end-motif row mode: {row_mode!r}")

    return storage_mode, row_mode


def _validate_required_arrays(store: Any, storage_mode: str, row_mode: str) -> None:
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
    try:
        store[name]
    except Exception:
        return False
    return True


def _validate_dense_counts(counts: Any, n_rows: int, n_motifs: int) -> None:
    dimension_names = tuple(getattr(counts.metadata, "dimension_names", ()) or ())
    expected = ("row", "motif")
    if dimension_names != expected:
        raise ValueError(
            f"dense counts dimensions must be {expected}, found {dimension_names}"
        )
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
    _read_labels(
        store["sparse/sparse_dimension"],
        "sparse_dimension_name",
        len(sparse_dimension),
        "sparse_dimension",
    )


def _read_array(store: Any, name: str) -> np.ndarray:
    return np.asarray(store[name][:])


def _read_motif_ascii_labels(store: Any, expected_len: int) -> np.ndarray:
    motif_byte = _read_array(store, "motif_byte")
    _validate_axis(motif_byte, "motif_byte")

    motif_ascii = _read_array(store, "motif_ascii")
    if motif_ascii.ndim != 2:
        raise ValueError(
            f"motif_ascii must have rank 2, found rank {motif_ascii.ndim}"
        )
    expected_shape = (expected_len, len(motif_byte))
    if tuple(motif_ascii.shape) != expected_shape:
        raise ValueError(
            "motif_ascii shape does not match motif axes: "
            f"motif_ascii={motif_ascii.shape}, expected={expected_shape}"
        )
    if motif_ascii.dtype != np.uint8:
        raise ValueError(f"motif_ascii must have dtype uint8, found {motif_ascii.dtype}")

    try:
        labels = [bytes(row).decode("ascii") for row in motif_ascii]
    except UnicodeDecodeError as error:
        raise ValueError("motif_ascii contains non-ASCII bytes") from error
    return np.asarray(labels, dtype=str)


def _read_labels(
    array: Any, field_name: str, expected_len: int, array_name: str
) -> np.ndarray:
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


def _required(value: np.ndarray | zarr.Array | None, name: str) -> np.ndarray | zarr.Array:
    if value is None:
        raise ValueError(f"End-motif Zarr store does not contain {name}")
    return value
