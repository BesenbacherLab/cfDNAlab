from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sparse

from ._helpers import filter_blacklisted_fraction, validate_scalar_bool
from .ends import EndMotifCounts
from .ref_kmers import RefKmerFrequencies

UNSUPPORTED_MOTIF_POLICIES = {"error", "drop", "keep_na"}


def _reference_corrected_data_frame(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    *,
    window_idxs: int | Sequence[int] | None = None,
    groups: str | Sequence[str] | None = None,
    group_idxs: int | Sequence[int] | None = None,
    densify: bool = False,
    motifs: str | Sequence[str] | None = None,
    motif_idxs: int | Sequence[int] | None = None,
    max_blacklisted_fraction: float = 1.0,
    use_global_bias: bool = False,
    unsupported_motifs: str = "error",
) -> pd.DataFrame:
    """
    Correct end-motif counts for matching reference k-mer frequencies.

    The returned data frame starts from `ends.data_frame()`. Each count is
    divided by `reference_frequency * correction_motif_count`.
    `correction_motif_count` is computed separately for each reference row from
    motifs with positive reference frequency. A uniform reference within that
    row leaves counts unchanged.

    Concrete end-motif labels are matched to reference k-mers by removing the
    `_` separator, for example `AT_CG -> ATCG`. Motif-group outputs are matched
    directly by group label.

    Reference k-mer output is read without densifying. For sparse reference
    output, omitted row/motif pairs are treated as zero frequency.

    Sample motifs can be absent from the reference genome, or can have zero
    reference frequency in the matched row. Positive end-motif counts for those
    motifs cannot be divided by a reference bias. By default this is an error.
    Set `unsupported_motifs="drop"` to omit rows whose reference motif has no
    positive reference frequency, or `unsupported_motifs="keep_na"` to keep
    them with `NaN` corrected counts.

    By default, end-motif and reference k-mer rows must match exactly. If
    `ref_kmers` is global and `ends` is windowed or grouped, pass
    `use_global_bias=True` to apply the global reference composition to every
    end-motif row.

    `cfdna ends` and `cfdna ref-kmers` both write forward-oriented motif labels.
    Right-end motifs have already been reverse-complemented by `cfdna ends`.

    Window, group, and motif selectors follow the same rules as
    `ends.data_frame()`. Motif selectors choose the returned end-motif rows.
    They do not change the reference support used for scaling, so selecting a
    motif gives the same corrected value as filtering the full corrected data
    frame afterward.

    Parameters
    ----------
    ends
        Loaded end-motif output.
    ref_kmers
        Loaded reference k-mer output built with matching reference settings.
    window_idxs
        Window index or indices for windowed outputs.
    groups
        Group name or names for grouped outputs. Use either `groups` or
        `group_idxs`, not both.
    group_idxs
        Group index or indices for grouped outputs. Use either `groups` or
        `group_idxs`, not both.
    densify
        Whether to include explicit zero-count end-motif rows for sparse
        end-motif output. Reference k-mer output is not densified.
    motifs
        End-motif label or labels to return. Use either `motifs` or
        `motif_idxs`, not both.
    motif_idxs
        End-motif index or indices to return. Use either `motifs` or
        `motif_idxs`, not both.
    max_blacklisted_fraction
        Maximum row `blacklisted_fraction` in 0..1 to keep before correction.
        For matched windowed or grouped references, this filters end-motif rows
        first. Reference rows are then matched to the remaining end-motif rows,
        not filtered independently by the reference file's blacklist fractions.
    use_global_bias
        Whether a global reference k-mer output may be applied to every
        non-global end-motif row.
    unsupported_motifs
        What to do when an observed sample motif has no positive reference
        frequency. Use `"error"`, `"drop"`, or `"keep_na"`.

    Returns
    -------
    pandas.DataFrame
        End-motif rows with `reference_motif`, `reference_frequency`,
        `correction_motif_count`, `reference_scale`, and
        `reference_corrected_count`.
    """
    use_global_bias = validate_scalar_bool(use_global_bias, "use_global_bias")
    unsupported_motifs = _validate_unsupported_motifs(unsupported_motifs)
    _validate_reference_correction_inputs(ends, ref_kmers, use_global_bias)
    row_columns = _reference_correction_row_columns(ends.row_mode())
    reference_row_columns = _reference_correction_reference_row_columns(
        ends,
        ref_kmers,
        use_global_bias,
    )
    _validate_matching_rows(ends, ref_kmers, row_columns, reference_row_columns)

    end_rows = ends._data_frame(
        densify=densify,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        motifs=motifs,
        motif_idxs=motif_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
    )
    if end_rows.empty:
        return _add_empty_reference_correction_columns(end_rows)

    ref_row_indices = _reference_row_indices_for_end_rows(
        ref_kmers,
        end_rows,
        reference_row_columns,
    )
    ref_rows = _reference_rows_for_indices(
        ref_kmers,
        ref_row_indices,
    )

    end_rows = end_rows.copy()
    ref_rows = ref_rows.copy()
    if ends.end_motifs.motif_axis_kind == "motif_group":
        end_rows["reference_motif"] = end_rows["motif"]
    else:
        end_rows["reference_motif"] = end_rows["motif"].str.replace(
            "_", "", regex=False
        )

    ref_rows = ref_rows.rename(
        columns={"motif": "reference_motif", "frequency": "reference_frequency"}
    )
    ref_columns = reference_row_columns + ["reference_motif", "reference_frequency"]
    ref_rows = ref_rows[ref_columns]
    ref_rows = _filter_reference_rows_to_end_rows(
        ref_rows,
        end_rows,
        reference_row_columns,
    )
    merge_columns = reference_row_columns + ["reference_motif"]
    duplicated_reference = ref_rows.duplicated(merge_columns)
    if duplicated_reference.any():
        raise ValueError("Reference k-mer rows are not unique for row and motif labels")
    reference_support_counts = _reference_support_counts(
        ref_rows,
        reference_row_columns,
    )

    corrected = end_rows.merge(
        ref_rows,
        on=merge_columns,
        how="left",
        sort=False,
    )
    missing_reference = corrected["reference_frequency"].isna()
    if missing_reference.any():
        corrected.loc[missing_reference, "reference_frequency"] = 0.0

    corrected = _add_correction_motif_count(
        corrected,
        reference_support_counts,
        reference_row_columns,
    )

    unsupported_reference = corrected["reference_frequency"].le(0.0)
    positive_unsupported_reference = unsupported_reference & corrected["count"].gt(0.0)
    if positive_unsupported_reference.any() and unsupported_motifs == "error":
        unsupported_labels = sorted(
            corrected.loc[positive_unsupported_reference, "reference_motif"].unique()
        )
        raise ValueError(
            "Positive-count end motifs have no positive reference frequency: "
            f"{unsupported_labels!r}. Pass unsupported_motifs='drop' to omit "
            "those rows, or unsupported_motifs='keep_na' to keep them with "
            "NaN corrected counts."
        )
    if unsupported_motifs == "drop":
        corrected = corrected.loc[~unsupported_reference].copy()

    corrected["reference_scale"] = (
        corrected["reference_frequency"] * corrected["correction_motif_count"]
    )
    corrected["reference_corrected_count"] = 0.0
    supported_reference = corrected["reference_scale"].gt(0.0)
    corrected.loc[supported_reference, "reference_corrected_count"] = (
        corrected.loc[supported_reference, "count"]
        / corrected.loc[supported_reference, "reference_scale"]
    )
    if unsupported_motifs == "keep_na":
        corrected.loc[
            positive_unsupported_reference.reindex(corrected.index, fill_value=False),
            "reference_corrected_count",
        ] = np.nan
    return corrected


def _reference_corrected_counts_array(
    ends: EndMotifCounts,
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
) -> np.ndarray:
    """
    Return reference-corrected end-motif counts as a dense NumPy array.
    """
    allow_densify = validate_scalar_bool(allow_densify, "allow_densify")
    _validate_fixed_shape_unsupported_policy(
        unsupported_motifs,
        "corrected_counts_array",
    )
    if ends.storage_mode() == "sparse_coo" and not allow_densify:
        raise ValueError(
            "corrected_counts_array() would turn sparse end-motif output into "
            "a dense in-memory array. Use sparse_corrected_counts_matrix() or "
            "pass allow_densify=True."
        )
    row_indices, motif_indices = _selected_end_axes(
        ends,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        motifs=motifs,
        motif_idxs=motif_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
    )
    corrected = _reference_corrected_data_frame(
        ends,
        ref_kmers,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        densify=True,
        motifs=motifs,
        motif_idxs=motif_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
        use_global_bias=use_global_bias,
        unsupported_motifs=unsupported_motifs,
    )
    return corrected["reference_corrected_count"].to_numpy(dtype=float).reshape(
        (len(row_indices), len(motif_indices))
    )


def _sparse_reference_corrected_counts_matrix(
    ends: EndMotifCounts,
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
) -> sparse.coo_matrix:
    """
    Return reference-corrected end-motif counts as a SciPy COO matrix.
    """
    _validate_fixed_shape_unsupported_policy(
        unsupported_motifs,
        "sparse_corrected_counts_matrix",
    )
    row_indices, motif_indices = _selected_end_axes(
        ends,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        motifs=motifs,
        motif_idxs=motif_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
    )
    if len(row_indices) == 0 or len(motif_indices) == 0:
        return sparse.coo_matrix((len(row_indices), len(motif_indices)))

    corrected = _reference_corrected_data_frame(
        ends,
        ref_kmers,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        densify=False,
        motifs=motifs,
        motif_idxs=motif_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
        use_global_bias=use_global_bias,
        unsupported_motifs=unsupported_motifs,
    )
    if corrected.empty:
        return sparse.coo_matrix((len(row_indices), len(motif_indices)))

    row_positions = _selected_row_positions(ends, row_indices)
    row_columns = _reference_correction_row_columns(ends.row_mode())
    motif_positions = _selected_motif_positions(ends, motif_indices)
    corrected_values = corrected["reference_corrected_count"].to_numpy(dtype=float)
    stored = (corrected_values != 0.0) | np.isnan(corrected_values)
    if not np.any(stored):
        return sparse.coo_matrix((len(row_indices), len(motif_indices)))

    matrix_rows = np.asarray(
        [
            _row_position_for_corrected_row(row_positions, row, row_columns)
            for _, row in corrected.iterrows()
        ],
        dtype=np.int64,
    )[stored]
    matrix_columns = np.asarray(
        [motif_positions[row["motif"]] for _, row in corrected.iterrows()],
        dtype=np.int64,
    )[stored]
    return sparse.coo_matrix(
        (
            corrected_values[stored],
            (matrix_rows, matrix_columns),
        ),
        shape=(len(row_indices), len(motif_indices)),
    )


def _validate_unsupported_motifs(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError(
            f"unsupported_motifs must be a string, got {type(value).__name__}"
        )
    if value not in UNSUPPORTED_MOTIF_POLICIES:
        choices = "', '".join(sorted(UNSUPPORTED_MOTIF_POLICIES))
        raise ValueError(f"unsupported_motifs must be one of '{choices}'")
    return value


def _validate_fixed_shape_unsupported_policy(value: str, method_name: str) -> None:
    unsupported_motifs = _validate_unsupported_motifs(value)
    if unsupported_motifs == "drop":
        raise ValueError(
            f"unsupported_motifs='drop' cannot be represented in a fixed-shape "
            f"{method_name}() result. Use data_frame(ref_kmers=..., "
            "unsupported_motifs='drop') or unsupported_motifs='keep_na'."
        )


def _selected_end_axes(
    ends: EndMotifCounts,
    *,
    window_idxs: int | Sequence[int] | None,
    groups: str | Sequence[str] | None,
    group_idxs: int | Sequence[int] | None,
    motifs: str | Sequence[str] | None,
    motif_idxs: int | Sequence[int] | None,
    max_blacklisted_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    row_indices = ends._resolve_row_selector(window_idxs, groups, group_idxs)
    row_indices = filter_blacklisted_fraction(
        ends._row_metadata_data_frame(),
        max_blacklisted_fraction,
        row_indices=row_indices,
    )
    motif_indices = ends._resolve_motif_selector(motifs, motif_idxs)
    return row_indices, motif_indices


def _selected_row_positions(
    ends: EndMotifCounts,
    row_indices: np.ndarray,
) -> dict[tuple[object, ...], int]:
    row_columns = _reference_correction_row_columns(ends.row_mode())
    row_metadata = (
        ends._row_metadata_data_frame().iloc[row_indices].reset_index(drop=True)
    )
    return {
        tuple(row[column] for column in row_columns): position
        for position, (_, row) in enumerate(row_metadata.iterrows())
    }


def _selected_motif_positions(
    ends: EndMotifCounts,
    motif_indices: np.ndarray,
) -> dict[object, int]:
    motif_metadata = ends.motifs_metadata().iloc[motif_indices].reset_index(drop=True)
    return {
        row["motif"]: position
        for position, (_, row) in enumerate(motif_metadata.iterrows())
    }


def _row_position_for_corrected_row(
    row_positions: dict[tuple[object, ...], int],
    row: pd.Series,
    row_columns: list[str],
) -> int:
    return row_positions[tuple(row[column] for column in row_columns)]


def _validate_reference_correction_inputs(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    use_global_bias: bool,
) -> None:
    if not isinstance(ends, EndMotifCounts):
        raise TypeError("ends must be an EndMotifCounts object")
    if not isinstance(ref_kmers, RefKmerFrequencies):
        raise TypeError("ref_kmers must be a RefKmerFrequencies object")
    if use_global_bias and ref_kmers.row_mode() != "global":
        raise ValueError(
            "use_global_bias=True requires a global reference k-mer output"
        )
    if ends.row_mode() != ref_kmers.row_mode():
        if ref_kmers.row_mode() == "global" and ends.row_mode() != "global":
            if not use_global_bias:
                raise ValueError(
                    "Reference k-mer output is global but end-motif output is "
                    f"{ends.row_mode()!r}. Pass use_global_bias=True to apply the "
                    "global reference bias to every end-motif row."
                )
        else:
            raise ValueError(
                "End-motif and reference k-mer row modes must match: "
                f"{ends.row_mode()!r} != {ref_kmers.row_mode()!r}"
            )
    if ends.end_motifs.motif_axis_kind == "motif_group":
        if ref_kmers.motif_axis_kind() != "motif_group":
            raise ValueError(
                "Grouped end-motif output requires grouped reference k-mer output"
            )
        return

    if ref_kmers.motif_axis_kind() != "motif":
        raise ValueError(
            "Concrete end-motif output requires concrete reference k-mer output"
        )
    if ref_kmers.canonical():
        raise ValueError(
            "Reference correction requires non-canonical reference k-mer output"
        )

    reference_motif_lengths = (
        ends.motifs_metadata()["motif"].str.replace("_", "", regex=False).str.len()
    )
    if not (reference_motif_lengths == ref_kmers.kmer_size()).all():
        raise ValueError(
            "End-motif width must match reference k-mer size "
            f"({ref_kmers.kmer_size()})"
        )


def _validate_matching_rows(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    row_columns: list[str],
    reference_row_columns: list[str],
) -> None:
    if not reference_row_columns:
        return

    end_row_keys = (
        ends._row_metadata_data_frame()[row_columns]
        .drop_duplicates()
        .sort_values(row_columns)
        .reset_index(drop=True)
    )
    ref_row_keys = (
        ref_kmers._row_metadata_data_frame()[reference_row_columns]
        .drop_duplicates()
        .sort_values(reference_row_columns)
        .reset_index(drop=True)
    )
    if len(end_row_keys) != len(ends._row_metadata_data_frame()):
        raise ValueError("End-motif row labels are not unique enough for correction")
    if len(ref_row_keys) != len(ref_kmers._row_metadata_data_frame()):
        raise ValueError(
            "Reference k-mer row labels are not unique enough for correction"
        )
    if not end_row_keys.equals(ref_row_keys):
        raise ValueError(
            "End-motif and reference k-mer rows do not match. "
            "Run ref-kmers with the same windowing or grouping."
        )


def _reference_correction_row_columns(row_mode: str) -> list[str]:
    if row_mode == "global":
        return ["row_label"]
    if row_mode in {"size", "bed"}:
        return ["window_idx", "chrom", "start", "end"]
    if row_mode == "grouped_bed":
        return ["group_name"]
    raise ValueError(f"Unsupported end-motif row mode for correction: {row_mode!r}")


def _reference_correction_reference_row_columns(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    use_global_bias: bool,
) -> list[str]:
    if (
        use_global_bias
        and ref_kmers.row_mode() == "global"
        and ends.row_mode() != "global"
    ):
        return []
    return _reference_correction_row_columns(ref_kmers.row_mode())


def _reference_row_indices_for_end_rows(
    ref_kmers: RefKmerFrequencies,
    end_rows: pd.DataFrame,
    reference_row_columns: list[str],
) -> np.ndarray:
    if not reference_row_columns:
        return np.arange(len(ref_kmers._row_metadata_data_frame()), dtype=np.int64)

    reference_metadata = ref_kmers._row_metadata_data_frame()
    reference_indices_by_key = {
        row_key: row_index
        for row_index, row_key in enumerate(
            _row_key_tuples(reference_metadata, reference_row_columns)
        )
    }
    selected_row_keys = dict.fromkeys(
        _row_key_tuples(end_rows, reference_row_columns)
    )
    try:
        return np.asarray(
            [reference_indices_by_key[row_key] for row_key in selected_row_keys],
            dtype=np.int64,
        )
    except KeyError as error:
        raise ValueError(
            "Selected end-motif row has no matching reference k-mer row"
        ) from error


def _reference_rows_for_indices(
    ref_kmers: RefKmerFrequencies,
    row_indices: np.ndarray,
) -> pd.DataFrame:
    motif_indices = np.arange(len(ref_kmers.motifs_metadata()), dtype=np.int64)
    if ref_kmers.storage_mode() == "sparse_coo":
        return ref_kmers._stored_data_frame_for_indices(row_indices, motif_indices)
    return ref_kmers._complete_data_frame_for_indices(
        row_indices,
        motif_indices,
        densify=False,
    )


def _reference_support_counts(
    ref_rows: pd.DataFrame,
    reference_row_columns: list[str],
) -> pd.DataFrame | int:
    positive_ref_rows = ref_rows.loc[ref_rows["reference_frequency"].gt(0.0)]
    if not reference_row_columns:
        return len(positive_ref_rows)
    if positive_ref_rows.empty:
        return pd.DataFrame(
            columns=reference_row_columns + ["correction_motif_count"]
        )
    return (
        positive_ref_rows.groupby(reference_row_columns, sort=False)
        .size()
        .reset_index(name="correction_motif_count")
    )


def _filter_reference_rows_to_end_rows(
    ref_rows: pd.DataFrame,
    end_rows: pd.DataFrame,
    reference_row_columns: list[str],
) -> pd.DataFrame:
    if not reference_row_columns:
        return ref_rows
    selected_row_keys = set(_row_key_tuples(end_rows, reference_row_columns))
    if not selected_row_keys:
        return ref_rows.iloc[[]].reset_index(drop=True)
    ref_row_keys = _row_key_tuples(ref_rows, reference_row_columns)
    keep = np.asarray(
        [row_key in selected_row_keys for row_key in ref_row_keys],
        dtype=bool,
    )
    return ref_rows.loc[keep].reset_index(drop=True)


def _row_key_tuples(
    data_frame: pd.DataFrame,
    row_columns: list[str],
) -> list[tuple[object, ...]]:
    return list(data_frame[row_columns].itertuples(index=False, name=None))


def _add_correction_motif_count(
    corrected: pd.DataFrame,
    reference_support_counts: pd.DataFrame | int,
    reference_row_columns: list[str],
) -> pd.DataFrame:
    if not reference_row_columns:
        corrected = corrected.copy()
        corrected["correction_motif_count"] = int(reference_support_counts)
        return corrected

    corrected = corrected.merge(
        reference_support_counts,
        on=reference_row_columns,
        how="left",
        sort=False,
    )
    corrected["correction_motif_count"] = (
        corrected["correction_motif_count"].fillna(0).astype(np.int64)
    )
    return corrected


def _add_empty_reference_correction_columns(data_frame: pd.DataFrame) -> pd.DataFrame:
    data_frame = data_frame.copy()
    data_frame["reference_motif"] = np.asarray([], dtype=object)
    data_frame["reference_frequency"] = np.asarray([], dtype=float)
    data_frame["correction_motif_count"] = np.asarray([], dtype=np.int64)
    data_frame["reference_scale"] = np.asarray([], dtype=float)
    data_frame["reference_corrected_count"] = np.asarray([], dtype=float)
    return data_frame
