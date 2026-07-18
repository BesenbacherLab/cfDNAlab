from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sparse

from ._helpers import (
    filter_blacklisted_fraction,
    normalize_strings,
    validate_scalar_bool,
)
from .ends import EndMotifCounts
from .ref_kmers import RefKmerFrequencies

UNSUPPORTED_MOTIF_POLICIES = {"error", "drop", "keep_na"}
TWO_SIDED_CORRECTION_MODES = {"joint", "split", "outside", "inside"}


@dataclass(frozen=True)
class _ReferenceCorrectionMode:
    """Resolved correction strategy and its side-specific axis information."""

    mode: str
    outside_width: int = 0
    inside_width: int = 0
    side_labels: tuple[str, ...] = ()


@dataclass
class _ReferenceCorrectionContext:
    """
    Validated row and motif axes shared by all corrected output forms.

    `end_row_key_columns` identifies rows in the returned end-count data.
    `reference_row_key_columns` identifies the reference composition joined to
    those rows and is empty when a global composition is broadcast.
    """

    correction_mode: _ReferenceCorrectionMode
    end_row_key_columns: list[str]
    reference_row_key_columns: list[str]
    selected_mode_labels: list[str]


def _prepare_reference_correction(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    *,
    motifs: str | Sequence[str] | None,
    motif_idxs: int | Sequence[int] | None,
    use_global_bias: bool,
    two_sided_correction: str | None,
) -> _ReferenceCorrectionContext:
    """
    Validate shared correction inputs and resolve reusable row and motif axes.

    The returned context records the correction mode, sample and reference row
    keys, and selected output labels. Data frame, dense-array, and sparse-matrix
    paths reuse it so validation and axis selection cannot diverge.
    """
    # Validate options and the relationship between the two stored outputs first
    use_global_bias = validate_scalar_bool(use_global_bias, "use_global_bias")
    two_sided_correction = _validate_two_sided_correction(two_sided_correction)
    _validate_reference_correction_inputs(ends, ref_kmers, use_global_bias)
    # Resolve the output motif axis and the metadata keys used on each side
    correction_mode = _resolve_correction_mode(ends, ref_kmers, two_sided_correction)
    end_row_key_columns = _reference_correction_row_columns(ends.row_mode())
    reference_row_key_columns = _reference_correction_reference_row_columns(
        ends,
        ref_kmers,
        use_global_bias,
    )
    _validate_matching_rows(
        ends,
        ref_kmers,
        end_row_key_columns,
        reference_row_key_columns,
    )
    # Interpret the motif selector on the axis produced by the correction mode
    selected_mode_labels = _selected_mode_axis(
        ends,
        correction_mode,
        motifs=motifs,
        motif_idxs=motif_idxs,
    )
    return _ReferenceCorrectionContext(
        correction_mode=correction_mode,
        end_row_key_columns=end_row_key_columns,
        reference_row_key_columns=reference_row_key_columns,
        selected_mode_labels=selected_mode_labels,
    )


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
    two_sided_correction: str | None = None,
) -> pd.DataFrame:
    """
    Correct end-motif counts for matching reference k-mer frequencies.

    The returned data frame starts from `ends.data_frame()` and adds
    `corrected_count` and `corrected_frequency`. Formula internals such as
    reference frequencies and correction factors are not returned.

    Motif labels are matched to reference k-mers by removing the `_` separator,
    for example `AT_CG -> ATCG`. Motif-group outputs are matched directly by
    group label.

    Reference k-mer output is read without densifying. For sparse reference
    output, omitted row/motif pairs are treated as zero frequency.

    Reference correction
    --------------------
    Reference correction divides each observed end-motif count by a
    reference-based correction factor for the matched row. This factor is
    computed from the motif frequencies in the reference k-mer output and
    normalized so a uniform reference composition leaves counts unchanged.
    Motifs that are common in the reference row are scaled down. Motifs that
    are rare in the reference row are scaled up. Only motifs with a positive
    reference frequency contribute to the row's correction support.

    Two-sided correction modes
    --------------------------
    When motif labels contain both outside and inside bases, such as `"AC_GT"`,
    `two_sided_correction` chooses both the motif labels in the result and the
    correction factor used for each returned count.

    - `"joint"` keeps full labels such as `"AC_GT"` and corrects each count
      using the exact reference k-mer `"ACGT"`.

    - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
      correction factor from the two sides separately. For `"AC_GT"`, separate
      correction factors are calculated for outside label `"AC"` and inside
      label `"GT"`. Those two correction factors are multiplied and applied to
      the observed `"AC_GT"` count. Use this when you want full two-sided motif
      labels in the result, but the exact full reference k-mers are too sparse
      or you want the reference correction to treat outside and inside sequence
      composition separately.

    - `"outside"` returns outside labels such as `"AC_"`. For each outside
      label, all full motif counts with that outside label are summed first.
      For example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"`
      count. That summed count is corrected using the outside label `"AC"`.

    - `"inside"` returns inside labels such as `"_GT"`. For each inside label,
      all full motif counts with that inside label are summed first. For
      example, `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"` count.
      That summed count is corrected using the inside label `"GT"`.

    One-sided outputs do not accept an explicit mode.

    Reference motif set
    ------------------
    For `"split"`, `"outside"`, and `"inside"`, side-specific reference
    frequencies are calculated from the loaded full-length reference k-mers.
    For example, the outside frequency for `"AC"` is the sum of frequencies for
    loaded k-mers with prefix `"AC"`, such as `"ACTG"` and `"ACAA"`. The inside
    frequency for `"TG"` is the corresponding sum over loaded k-mers with suffix
    `"TG"`. Separate shorter reference k-mer runs are not required.

    A motifs file used for the reference output restricts these sums to the
    k-mers in that file. Without a motifs file, all k-mers in the reference
    output can contribute, including k-mers absent from the sample end-motif
    output.

    Corrected frequencies
    ---------------------
    `corrected_frequency` is `corrected_count` divided by the sum of corrected
    counts for the same output row. The sum uses the full motif axis produced by
    the selected correction mode. Motif selection is applied afterward, so a
    selected subset keeps its full-axis frequencies and may sum to less than 1.

    If the finite corrected-count total is zero, finite frequencies are zero.
    Correction fails if dividing a finite count by its correction factor would
    produce a non-finite corrected count.
    With `unsupported_motifs="keep_na"`, one positive count without a positive
    correction factor makes every corrected frequency in that output row
    `NaN`, because the full corrected total is unknown. With
    `unsupported_motifs="drop"`, unsupported units are removed before the total
    is calculated.

    An observed sample motif with a positive count is unsupported when it has no
    positive correction factor under the selected mode. By default this is an
    error. Set
    `unsupported_motifs="drop"` to omit unsupported rows, or
    `unsupported_motifs="keep_na"` to keep them with `NaN` corrected counts.

    By default, end-motif and reference k-mer rows must match exactly. If
    `ref_kmers` is global and `ends` is windowed or grouped, pass
    `use_global_bias=True` to apply the global reference composition to every
    end-motif row.

    `cfdna ends` writes motifs from each fragment end inward. Right-end motifs
    are therefore reverse-complemented relative to the stored reference.
    Reference correction requires `cfdna ref-kmers --orientation both`, which
    averages the reference-forward sequence and its reverse complement.

    Window, group, and motif selectors follow the same rules as
    `ends.data_frame()`. Motif selectors choose the returned end-motif rows.
    They do not change the reference support used for correction, so selecting
    a motif gives the same corrected value as filtering the full corrected data
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
        What to do when an observed sample motif has no positive correction
        factor under the selected mode. Use `"error"`, `"drop"`, or `"keep_na"`.
    two_sided_correction
        Required for two-sided motif labels such as `"AC_GT"`. Use `"joint"`,
        `"split"`, `"outside"`, or `"inside"`. Leave as `None` for one-sided
        motifs or motif groups.

    Returns
    -------
    pandas.DataFrame
        End-motif rows with `corrected_count` and `corrected_frequency`.
    """
    unsupported_motifs = _validate_unsupported_motifs(unsupported_motifs)
    context = _prepare_reference_correction(
        ends,
        ref_kmers,
        motifs=motifs,
        motif_idxs=motif_idxs,
        use_global_bias=use_global_bias,
        two_sided_correction=two_sided_correction,
    )
    return _reference_corrected_data_frame_from_context(
        ends,
        ref_kmers,
        context,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        densify=densify,
        max_blacklisted_fraction=max_blacklisted_fraction,
        unsupported_motifs=unsupported_motifs,
    )


def _reference_corrected_data_frame_from_context(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    context: _ReferenceCorrectionContext,
    *,
    window_idxs: int | Sequence[int] | None,
    groups: str | Sequence[str] | None,
    group_idxs: int | Sequence[int] | None,
    densify: bool,
    max_blacklisted_fraction: float,
    unsupported_motifs: str,
) -> pd.DataFrame:
    """
    Build a corrected data frame from previously validated correction state.

    End counts are loaded without motif filtering because corrected frequencies
    need every represented motif in the output axis. Side-only modes also need
    every full motif that contributes to a selected side. Reference support is
    calculated separately from the full reference motif axis. Motif selection
    is applied only after counts and frequencies are calculated.
    """

    # Load selected sample rows without motif filtering. Exact and
    # split corrected counts could be calculated from selected motifs alone, but
    # their corrected frequencies still need the corrected total over all motifs
    end_rows = ends._data_frame(
        densify=densify,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        motifs=None,
        motif_idxs=None,
        max_blacklisted_fraction=max_blacklisted_fraction,
    )
    if end_rows.empty:
        return _add_empty_reference_correction_columns(end_rows)

    # Preserve the public sample-row and column order across joins and aggregation
    output_columns = end_rows.columns.tolist()
    end_rows = end_rows.copy()
    end_rows["_cfdnalab_row_order"] = _row_order_indices(
        end_rows,
        context.end_row_key_columns,
    )

    # Select only reference rows matching the chosen sample rows, while keeping
    # every reference motif needed to define support and side marginals
    ref_row_indices = _reference_row_indices_for_end_rows(
        ref_kmers,
        end_rows,
        context.reference_row_key_columns,
    )
    ref_rows = _reference_rows_for_indices(
        ref_kmers,
        ref_row_indices,
    )

    ref_rows = _prepared_reference_rows(
        ref_rows,
        context.reference_row_key_columns,
    )

    # Convert the reference frequencies into correction denominators according
    # to whether full motifs or separate motif sides define the correction
    if context.correction_mode.mode == "exact":
        corrected = _correct_exact_label_data_frame(
            ends,
            end_rows,
            ref_rows,
            context.reference_row_key_columns,
            unsupported_motifs,
        )
    elif context.correction_mode.mode == "split":
        corrected = _correct_split_data_frame(
            end_rows,
            ref_rows,
            context.reference_row_key_columns,
            context.correction_mode,
            unsupported_motifs,
        )
    else:
        corrected = _correct_side_data_frame(
            end_rows,
            ref_rows,
            context.reference_row_key_columns,
            context.correction_mode,
            output_columns,
            unsupported_motifs,
        )

    # Normalize over the complete correction-mode motif axis before applying the
    # requested motif selector, so selection cannot change reported frequencies
    corrected = _add_corrected_frequency(
        corrected,
        context.end_row_key_columns,
    )
    corrected = corrected.loc[
        corrected["motif"].isin(context.selected_mode_labels)
    ].copy()

    # Restore the selector's motif order within the original sample-row order
    if context.selected_mode_labels:
        motif_order = {
            motif: index for index, motif in enumerate(context.selected_mode_labels)
        }
        corrected["_cfdnalab_motif_order"] = corrected["motif"].map(motif_order)
        corrected = corrected.sort_values(
            ["_cfdnalab_row_order", "_cfdnalab_motif_order"],
            kind="stable",
        )
        corrected = corrected.drop(columns=["_cfdnalab_motif_order"])
    corrected_columns = output_columns + ["corrected_count", "corrected_frequency"]
    return corrected[corrected_columns].reset_index(drop=True)


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
    two_sided_correction: str | None = None,
) -> np.ndarray:
    """
    Return reference-corrected end-motif counts as a dense NumPy array.

    The method validates that the unsupported-motif policy has a fixed shape,
    resolves selected rows and correction-mode labels, and calculates a complete
    corrected data frame before reshaping counts into row-by-motif order. Sparse
    source output requires explicit permission to densify.
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
    context = _prepare_reference_correction(
        ends,
        ref_kmers,
        motifs=motifs,
        motif_idxs=motif_idxs,
        use_global_bias=use_global_bias,
        two_sided_correction=two_sided_correction,
    )
    row_indices = _selected_rows(
        ends,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
    )
    corrected = _reference_corrected_data_frame_from_context(
        ends,
        ref_kmers,
        context,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        densify=True,
        max_blacklisted_fraction=max_blacklisted_fraction,
        unsupported_motifs=unsupported_motifs,
    )
    return (
        corrected["corrected_count"]
        .to_numpy(dtype=float)
        .reshape((len(row_indices), len(context.selected_mode_labels)))
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
    two_sided_correction: str | None = None,
) -> sparse.coo_matrix:
    """
    Return reference-corrected end-motif counts as a SciPy COO matrix.

    Correction uses the same full-axis support and frequency rules as the data
    frame path, but only nonzero and `NaN` corrected values are stored. Selected
    correction-row keys and mode labels determine matrix coordinates and shape,
    including when no rows or motifs are selected.
    """
    _validate_fixed_shape_unsupported_policy(
        unsupported_motifs,
        "sparse_corrected_counts_matrix",
    )
    context = _prepare_reference_correction(
        ends,
        ref_kmers,
        motifs=motifs,
        motif_idxs=motif_idxs,
        use_global_bias=use_global_bias,
        two_sided_correction=two_sided_correction,
    )
    row_indices = _selected_rows(
        ends,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        max_blacklisted_fraction=max_blacklisted_fraction,
    )
    if len(row_indices) == 0 or len(context.selected_mode_labels) == 0:
        return sparse.coo_matrix((len(row_indices), len(context.selected_mode_labels)))

    corrected = _reference_corrected_data_frame_from_context(
        ends,
        ref_kmers,
        context,
        window_idxs=window_idxs,
        groups=groups,
        group_idxs=group_idxs,
        densify=False,
        max_blacklisted_fraction=max_blacklisted_fraction,
        unsupported_motifs=unsupported_motifs,
    )
    if corrected.empty:
        return sparse.coo_matrix((len(row_indices), len(context.selected_mode_labels)))

    row_positions = _selected_row_positions(ends, row_indices)
    motif_positions = {
        motif_label: position
        for position, motif_label in enumerate(context.selected_mode_labels)
    }
    corrected_values = corrected["corrected_count"].to_numpy(dtype=float)
    stored = (corrected_values != 0.0) | np.isnan(corrected_values)
    if not np.any(stored):
        return sparse.coo_matrix((len(row_indices), len(context.selected_mode_labels)))

    matrix_rows = np.asarray(
        [
            row_positions[row_key]
            for row_key in _row_key_tuples(corrected, context.end_row_key_columns)
        ],
        dtype=np.int64,
    )[stored]
    matrix_columns = (
        corrected["motif"].map(motif_positions).to_numpy(dtype=np.int64)[stored]
    )
    return sparse.coo_matrix(
        (
            corrected_values[stored],
            (matrix_rows, matrix_columns),
        ),
        shape=(len(row_indices), len(context.selected_mode_labels)),
    )


def _validate_unsupported_motifs(value: str) -> str:
    """
    Validate the policy for observed motifs without reference support.

    The value must be exactly `"error"`, `"drop"`, or `"keep_na"`. Returning
    the validated string lets callers use a single validation path.
    """
    if not isinstance(value, str):
        raise TypeError(
            f"unsupported_motifs must be a string, got {type(value).__name__}"
        )
    if value not in UNSUPPORTED_MOTIF_POLICIES:
        choices = "', '".join(sorted(UNSUPPORTED_MOTIF_POLICIES))
        raise ValueError(f"unsupported_motifs must be one of '{choices}'")
    return value


def _validate_fixed_shape_unsupported_policy(value: str, method_name: str) -> None:
    """
    Validate unsupported-motif handling for a fixed-shape output.

    Dropping motifs would give different rows or columns different shapes, so
    matrix and array outputs accept only the error and missing-value policies.
    """
    unsupported_motifs = _validate_unsupported_motifs(value)
    if unsupported_motifs == "drop":
        raise ValueError(
            f"unsupported_motifs='drop' cannot be represented in a fixed-shape "
            f"{method_name}() result. Use data_frame(ref_kmers=..., "
            "unsupported_motifs='drop') or unsupported_motifs='keep_na'."
        )


def _validate_two_sided_correction(value: str | None) -> str | None:
    """
    Validate the requested interpretation of two-sided motif labels.

    `None` is preserved for one-sided or grouped axes. Explicit values must be
    one of the four supported two-sided correction modes.
    """
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(
            f"two_sided_correction must be a string or None, got {type(value).__name__}"
        )
    if value not in TWO_SIDED_CORRECTION_MODES:
        choices = "', '".join(sorted(TWO_SIDED_CORRECTION_MODES))
        raise ValueError(f"two_sided_correction must be one of '{choices}'")
    return value


def _resolve_correction_mode(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    two_sided_correction: str | None,
) -> _ReferenceCorrectionMode:
    """
    Resolve how sample motif labels map to reference correction factors.

    Grouped, empty, one-sided, and `"joint"` axes use exact-label correction.
    Other two-sided modes retain the inferred side widths, and side-only modes
    also store the derived output labels used for aggregation.
    """
    # Motif groups already define the complete correction axis. They have no
    # outside/inside split to reinterpret, so only exact-label correction applies
    if ends.end_motifs.motif_axis_kind == "motif_group":
        if two_sided_correction is not None:
            raise ValueError(
                "Motif-group end-motif outputs do not accept two_sided_correction"
            )
        return _ReferenceCorrectionMode(mode="exact")

    motif_labels = ends.motifs_metadata()["motif"].astype(str).tolist()
    # An empty stored axis has nothing from which to infer side widths. Treat it
    # as exact so empty selections can still return their expected shape
    if not motif_labels:
        return _ReferenceCorrectionMode(mode="exact")

    # Resolve the split from the sample labels, then verify that the reference
    # labels can use exactly that boundary before choosing a correction formula
    outside_width, inside_width = _infer_end_motif_side_widths(
        motif_labels,
        ref_kmers.kmer_size(),
    )
    _validate_reference_labels_split_cleanly(
        ref_kmers.motifs_metadata()["motif"].astype(str),
        outside_width,
        inside_width,
    )
    if outside_width == 0 or inside_width == 0:
        if two_sided_correction is not None:
            raise ValueError(
                "One-sided end-motif outputs do not accept two_sided_correction"
            )
        return _ReferenceCorrectionMode(mode="exact")
    # A true two-sided axis is ambiguous without an explicit choice because the
    # choice controls both the correction formula and, for side modes, output shape
    if two_sided_correction is None:
        raise ValueError(
            "two-sided end-motif labels with both outside and inside bases require two_sided_correction"
        )
    if two_sided_correction == "joint":
        return _ReferenceCorrectionMode(mode="exact")
    # Split mode keeps the stored full-motif axis. Outside and inside modes
    # instead need a stable, deduplicated derived axis for selection and ordering
    side_labels = (
        tuple(_side_axis_labels(motif_labels, two_sided_correction))
        if two_sided_correction in {"outside", "inside"}
        else ()
    )
    return _ReferenceCorrectionMode(
        mode=two_sided_correction,
        outside_width=outside_width,
        inside_width=inside_width,
        side_labels=side_labels,
    )


def _infer_end_motif_side_widths(
    motif_labels: Sequence[str],
    reference_kmer_size: int,
) -> tuple[int, int]:
    """
    Infer the outside and inside widths shared by all end-motif labels.

    Every label must contain the separator, have the reference k-mer's total
    width, and use the same split. An empty motif axis cannot define widths.
    """
    inferred_widths: tuple[int, int] | None = None
    for motif_label in motif_labels:
        outside, inside = _split_end_motif_label(motif_label)
        widths = (len(outside), len(inside))
        if sum(widths) != reference_kmer_size:
            raise ValueError(
                "End-motif width must match reference k-mer size "
                f"({reference_kmer_size}): {motif_label}"
            )
        if inferred_widths is None:
            inferred_widths = widths
        elif inferred_widths != widths:
            raise ValueError(
                "All end-motif labels must use the same outside and inside widths"
            )
    if inferred_widths is None:
        raise ValueError("Cannot infer side widths from an empty motif axis")
    return inferred_widths


def _split_end_motif_label(motif_label: str) -> tuple[str, str]:
    """
    Split an end-motif label into its outside and inside bases.

    Exactly one underscore must separate the two sides. Either side may be
    empty for a valid one-sided motif label.
    """
    parts = motif_label.split("_")
    if len(parts) != 2:
        raise ValueError(
            "End-motif label must contain exactly one '_' to separate "
            f"outside and inside bases: {motif_label}"
        )
    return parts[0], parts[1]


def _validate_reference_labels_split_cleanly(
    reference_labels: Sequence[str],
    outside_width: int,
    inside_width: int,
) -> None:
    """
    Ensure reference labels can be split using the end-motif side widths.

    Reference labels have no underscore, so their length must equal the sum of
    the inferred outside and inside widths.
    """
    expected_width = outside_width + inside_width
    for reference_label in reference_labels:
        if len(reference_label) != expected_width:
            raise ValueError(
                "Reference motif label must split into outside width "
                f"{outside_width} and inside width {inside_width}: {reference_label}"
            )


def _selected_mode_axis(
    ends: EndMotifCounts,
    correction_mode: _ReferenceCorrectionMode,
    *,
    motifs: str | Sequence[str] | None,
    motif_idxs: int | Sequence[int] | None,
) -> list[str]:
    """
    Resolve labels on the axis produced by the correction mode.

    Exact and split modes select the stored motif axis. Side-only modes select
    the derived side-label axis by label and reject stored motif indices because
    those indices do not describe the derived columns.
    """
    mode = correction_mode.mode
    # Exact and split results retain a one-to-one relationship with stored motif
    # indices, so the standard end-motif selector defines their result axis
    if mode in {"exact", "split"}:
        motif_indices = ends._resolve_motif_selector(motifs, motif_idxs)
        return ends.motifs_metadata().iloc[motif_indices]["motif"].astype(str).tolist()

    # Side-only labels were derived by collapsing multiple stored motifs. A
    # stored motif index therefore cannot identify a derived result column
    if motif_idxs is not None:
        raise ValueError(
            "motif index selectors are not supported for outside or inside reference correction"
        )
    full_side_labels = list(correction_mode.side_labels)
    if motifs is None:
        return full_side_labels
    # Validate label selectors against the derived axis while preserving the
    # caller's requested order for later data frame and matrix construction
    requested_labels = normalize_strings(motifs, name="motifs")
    side_label_set = set(full_side_labels)
    for requested_label in requested_labels:
        if requested_label not in side_label_set:
            raise ValueError(f"Side-mode motif axis has no label {requested_label!r}")
    return requested_labels


def _side_axis_labels(motif_labels: Sequence[str], side_mode: str) -> list[str]:
    """
    Build the unique labels for an outside- or inside-only motif axis.

    Outside labels retain a trailing underscore and inside labels retain a
    leading underscore. First occurrence determines output order.
    """
    side_labels: list[str] = []
    seen_labels: set[str] = set()
    for motif_label in motif_labels:
        outside, inside = _split_end_motif_label(motif_label)
        side_label = f"{outside}_" if side_mode == "outside" else f"_{inside}"
        if side_label not in seen_labels:
            seen_labels.add(side_label)
            side_labels.append(side_label)
    return side_labels


def _prepared_reference_rows(
    ref_rows: pd.DataFrame,
    reference_row_columns: list[str],
) -> pd.DataFrame:
    """
    Prepare reference rows for joining to selected end-motif rows.

    The motif and frequency columns receive unambiguous reference names, only
    correction columns are retained. The caller has already selected the exact
    reference rows needed for the chosen end-motif rows.
    """
    ref_rows = ref_rows.copy()
    ref_rows = ref_rows.rename(
        columns={"motif": "reference_motif", "frequency": "reference_frequency"}
    )
    ref_columns = reference_row_columns + ["reference_motif", "reference_frequency"]
    return ref_rows[ref_columns]


def _positive_support_counts(
    frequencies: pd.Series,
    data_frame: pd.DataFrame,
    row_columns: list[str],
) -> pd.Series:
    """
    Count positive frequencies within each reference row.

    Global reference composition receives the same scalar count on every row.
    Keyed output is grouped by its row-identifying metadata. Empty input returns
    an aligned empty integer series.
    """
    positive_frequency = frequencies.gt(0.0).astype(np.int64)
    if not row_columns:
        return pd.Series(
            int(positive_frequency.sum()),
            index=data_frame.index,
            dtype=np.int64,
        )
    row_groups = [data_frame[column] for column in row_columns]
    return positive_frequency.groupby(
        row_groups,
        sort=False,
        dropna=False,
    ).transform("sum")


def _correct_exact_label_data_frame(
    ends: EndMotifCounts,
    end_rows: pd.DataFrame,
    ref_rows: pd.DataFrame,
    reference_row_columns: list[str],
    unsupported_motifs: str,
) -> pd.DataFrame:
    """
    Correct each count with the matching full reference-motif frequency.

    Concrete end-motif labels are matched after removing their underscore,
    while motif-group labels match directly. Within each reference row, the
    correction denominator is the motif frequency multiplied by the number of
    positive-frequency motifs, making a uniform reference denominator equal to
    one. The selected unsupported-motif policy is applied after this join.
    """
    # Translate the sample label to the corresponding reference label before
    # joining. Concrete motif labels differ only by the end separator
    end_rows = end_rows.copy()
    if ends.end_motifs.motif_axis_kind == "motif_group":
        end_rows["reference_motif"] = end_rows["motif"]
    else:
        end_rows["reference_motif"] = end_rows["motif"].str.replace(
            "_", "", regex=False
        )
    merge_columns = reference_row_columns + ["reference_motif"]
    if ref_rows.duplicated(merge_columns).any():
        raise ValueError("Reference k-mer rows are not unique for row and motif labels")

    # Attach support before joining so frequency and support arrive together
    ref_rows = ref_rows.copy()
    ref_rows["number_of_supported_motifs"] = _positive_support_counts(
        ref_rows["reference_frequency"],
        ref_rows,
        reference_row_columns,
    )

    # Attach each motif's frequency and its row's support count in one join
    corrected = end_rows.merge(
        ref_rows,
        on=merge_columns,
        how="left",
        sort=False,
    )
    corrected["reference_frequency"] = corrected["reference_frequency"].fillna(0.0)
    corrected["number_of_supported_motifs"] = (
        corrected["number_of_supported_motifs"].fillna(0).astype(np.int64)
    )
    # Normalize against uniform composition. With N supported motifs, a uniform
    # reference frequency of 1/N gives a denominator of one
    corrected["reference_denominator"] = (
        corrected["reference_frequency"] * corrected["number_of_supported_motifs"]
    )
    return _apply_reference_denominator_policy(
        corrected,
        unsupported_motifs,
    )


def _correct_split_data_frame(
    end_rows: pd.DataFrame,
    ref_rows: pd.DataFrame,
    reference_row_columns: list[str],
    correction_mode: _ReferenceCorrectionMode,
    unsupported_motifs: str,
) -> pd.DataFrame:
    """
    Correct full motifs from separate outside and inside reference marginals.

    Reference frequencies are summed by prefix and suffix, each side is
    normalized against its own positive support, and the two denominators are
    multiplied. Counts and the original full-motif axis are otherwise retained.
    """
    # Expose each full sample motif's two join labels
    outside_width = correction_mode.outside_width
    inside_width = correction_mode.inside_width
    end_rows = _add_end_sides(end_rows, outside_width, inside_width)
    # Collapse the complete reference motif set independently by prefix and suffix
    outside_reference = _side_reference_denominator(
        ref_rows,
        reference_row_columns,
        "outside",
        outside_width,
        inside_width,
    )
    inside_reference = _side_reference_denominator(
        ref_rows,
        reference_row_columns,
        "inside",
        outside_width,
        inside_width,
    )

    # Attach both normalized side denominators, then combine their effects
    corrected = _merge_side_denominator(
        end_rows,
        outside_reference,
        reference_row_columns,
        "outside",
    )
    corrected = _merge_side_denominator(
        corrected,
        inside_reference,
        reference_row_columns,
        "inside",
    )
    corrected["reference_denominator"] = (
        corrected["outside_reference_denominator"]
        * corrected["inside_reference_denominator"]
    )
    return _apply_reference_denominator_policy(
        corrected,
        unsupported_motifs,
    )


def _correct_side_data_frame(
    end_rows: pd.DataFrame,
    ref_rows: pd.DataFrame,
    reference_row_columns: list[str],
    correction_mode: _ReferenceCorrectionMode,
    output_columns: list[str],
    unsupported_motifs: str,
) -> pd.DataFrame:
    """
    Collapse full motifs onto one side before applying reference correction.

    Counts sharing the selected outside or inside label are summed per output
    row. The function assigns indices from the derived side axis, restores row
    metadata, and divides each aggregate by the matching marginal denominator.
    """
    # Replace each full sample motif with its selected side label
    outside_width = correction_mode.outside_width
    inside_width = correction_mode.inside_width
    side_mode = correction_mode.mode
    end_rows = _add_end_sides(end_rows, outside_width, inside_width)
    side_column = "outside" if side_mode == "outside" else "inside"
    end_rows["motif"] = (
        end_rows["outside"] + "_"
        if side_mode == "outside"
        else "_" + end_rows["inside"]
    )
    side_axis = correction_mode.side_labels
    side_index_by_label = {label: index for index, label in enumerate(side_axis)}
    end_rows["motif_index"] = (
        end_rows["motif"].map(side_index_by_label).astype(np.int64)
    )
    # Save row metadata once because the following aggregation keeps only keys,
    # the derived motif label, and the summed count
    row_metadata_columns = [
        column
        for column in output_columns
        if column not in {"motif_index", "motif", "count"}
    ]
    row_metadata = end_rows[
        ["_cfdnalab_row_order"] + row_metadata_columns
    ].drop_duplicates("_cfdnalab_row_order")
    # Sum all full sample motifs that contribute to the same side label
    aggregated = (
        end_rows.groupby(
            ["_cfdnalab_row_order", "motif_index", "motif", side_column],
            sort=False,
        )["count"]
        .sum()
        .reset_index()
    )
    aggregated = aggregated.merge(
        row_metadata,
        on="_cfdnalab_row_order",
        how="left",
        sort=False,
        validate="many_to_one",
    )
    # Derive side frequencies from the complete reference motif set and attach the
    # denominator for the side retained in the result
    side_reference = _side_reference_denominator(
        ref_rows,
        reference_row_columns,
        side_column,
        outside_width,
        inside_width,
    )
    corrected = _merge_side_denominator(
        aggregated,
        side_reference,
        reference_row_columns,
        side_column,
    )
    denominator_column = f"{side_column}_reference_denominator"
    corrected["reference_denominator"] = corrected[denominator_column]
    return _apply_reference_denominator_policy(
        corrected,
        unsupported_motifs,
    )


def _add_end_sides(
    end_rows: pd.DataFrame,
    outside_width: int,
    inside_width: int,
) -> pd.DataFrame:
    """
    Parse motif labels into validated outside and inside columns.

    A copy is returned so callers can reshape it without changing their input.
    Every parsed label must match the widths resolved for the correction mode.
    """
    end_rows = end_rows.copy()
    outside_inside = end_rows["motif"].map(_split_end_motif_label)
    end_rows["outside"] = [outside for outside, _ in outside_inside]
    end_rows["inside"] = [inside for _, inside in outside_inside]
    invalid_width = end_rows["outside"].str.len().ne(outside_width) | end_rows[
        "inside"
    ].str.len().ne(inside_width)
    if invalid_width.any():
        raise ValueError("End-motif label does not match inferred side widths")
    return end_rows


def _side_reference_denominator(
    ref_rows: pd.DataFrame,
    reference_row_columns: list[str],
    side_column: str,
    outside_width: int,
    inside_width: int,
) -> pd.DataFrame:
    """
    Build the correction-denominator table for one motif side.

    Full reference labels are reduced to the requested prefix or suffix. Their
    frequencies are then aggregated within each reference row. Each marginal
    frequency is multiplied by the number of positive-frequency side labels in
    that row, so uniform side composition gives denominator one.
    """
    # Turn each full reference motif into the selected sample-side label
    ref_rows = ref_rows.copy()
    if side_column == "outside":
        ref_rows[side_column] = ref_rows["reference_motif"].str.slice(
            0,
            outside_width,
        )
    else:
        ref_rows[side_column] = ref_rows["reference_motif"].str.slice(
            outside_width,
            outside_width + inside_width,
        )

    # Sum full reference-motif frequencies into a marginal frequency per side
    group_columns = reference_row_columns + [side_column]
    side_frequencies = (
        ref_rows.groupby(group_columns, sort=False)["reference_frequency"]
        .sum()
        .reset_index(name="side_reference_frequency")
    )

    side_frequencies["side_support_count"] = _positive_support_counts(
        side_frequencies["side_reference_frequency"],
        side_frequencies,
        reference_row_columns,
    )

    # Scale relative to a uniform distribution over the supported side labels
    denominator_column = f"{side_column}_reference_denominator"
    side_frequencies[denominator_column] = (
        side_frequencies["side_reference_frequency"]
        * side_frequencies["side_support_count"]
    )
    return side_frequencies[reference_row_columns + [side_column, denominator_column]]


def _merge_side_denominator(
    end_rows: pd.DataFrame,
    side_reference: pd.DataFrame,
    reference_row_columns: list[str],
    side_column: str,
) -> pd.DataFrame:
    """
    Join one side's correction denominator to sample rows.

    Reference row keys and the outside or inside label identify the denominator.
    Missing matches are deliberately preserved for policy handling.
    """
    denominator_column = f"{side_column}_reference_denominator"
    side_columns = reference_row_columns + [side_column, denominator_column]
    return end_rows.merge(
        side_reference[side_columns],
        on=reference_row_columns + [side_column],
        how="left",
        sort=False,
    )


def _apply_reference_denominator_policy(
    corrected: pd.DataFrame,
    unsupported_motifs: str,
) -> pd.DataFrame:
    """
    Divide counts by valid denominators and apply unsupported-motif handling.

    Missing denominators become zero. Positive counts without positive support
    either raise, are dropped, or become `NaN`, while unsupported zero counts
    remain zero. The function also rejects any division producing a non-finite
    corrected count.
    """
    corrected = corrected.copy()
    denominator_column = "reference_denominator"

    # A missing denominator is the result of a left join with no matching
    # reference support. It is semantically unsupported, just like an explicit zero
    corrected[denominator_column] = corrected[denominator_column].fillna(0.0)

    # Keep two masks because unsupported zero counts are harmless under error and
    # keep_na, while drop intentionally removes the complete unsupported axis
    unsupported_reference = corrected[denominator_column].le(0.0)
    positive_unsupported_reference = unsupported_reference & corrected["count"].gt(0.0)

    # Fail before modifying the table so the error reports every affected motif
    if positive_unsupported_reference.any() and unsupported_motifs == "error":
        unsupported_labels = sorted(
            corrected.loc[positive_unsupported_reference, "motif"].unique()
        )
        raise ValueError(
            "Positive-count end motifs have no positive reference-based correction factor: "
            f"{unsupported_labels!r}. Pass unsupported_motifs='drop' to omit "
            "those rows, or unsupported_motifs='keep_na' to keep them with "
            "NaN corrected counts."
        )

    # Drop is based on reference support, not observed count. This removes both
    # positive and zero-count rows and therefore permits a variable-shaped result
    if unsupported_motifs == "drop":
        corrected = corrected.loc[~unsupported_reference].copy()

    # Starting at zero gives unsupported zero-count rows their defined result and
    # lets the division below operate exclusively on positive denominators
    corrected["corrected_count"] = 0.0
    supported_reference = corrected[denominator_column].gt(0.0)
    with np.errstate(over="ignore", invalid="ignore"):
        corrected_values = corrected.loc[supported_reference, "count"].to_numpy(
            dtype=float
        ) / corrected.loc[supported_reference, denominator_column].to_numpy(dtype=float)

    # A denominator can be positive but too small to produce a finite floating-
    # point result. Reject that row instead of allowing infinities into normalization
    if not np.isfinite(corrected_values).all():
        non_finite_labels = sorted(
            corrected.loc[
                supported_reference,
                "motif",
            ]
            .iloc[np.flatnonzero(~np.isfinite(corrected_values))]
            .unique()
        )
        raise ValueError(
            "Reference correction produced non-finite corrected counts for "
            f"motifs {non_finite_labels!r}. Reference correction factors may "
            "be too small for the observed counts."
        )
    corrected.loc[supported_reference, "corrected_count"] = corrected_values

    # Only positive unsupported counts are unknown. A zero observation corrected
    # by an unknown factor is still exactly zero. Frequency handling later expands
    # any resulting NaN to the complete output row because its total is unknown
    if unsupported_motifs == "keep_na":
        corrected.loc[positive_unsupported_reference, "corrected_count"] = np.nan
    return corrected


def _add_corrected_frequency(
    corrected: pd.DataFrame,
    end_row_key_columns: list[str],
) -> pd.DataFrame:
    """
    Calculate corrected motif frequencies within each output row.

    Counts are scaled by the row maximum before summing to avoid overflow. A
    finite zero-total row receives zero frequencies, while `keep_na` makes the
    entire row's frequencies unknown if any positive corrected count is unknown.
    """
    corrected = corrected.copy()
    corrected["corrected_frequency"] = 0.0
    if corrected.empty:
        return corrected

    # Group by the metadata columns that identify an output row, then reuse the
    # resulting index for the maximum, total, and missing-value checks
    row_group_indices = corrected.groupby(
        end_row_key_columns,
        sort=False,
        dropna=False,
    ).ngroup()
    corrected_counts = corrected["corrected_count"]

    # Scaling before summing avoids overflow. Replacing undefined divisions with
    # zero also gives zero-total rows their required zero frequencies
    row_maximum = corrected_counts.groupby(row_group_indices).transform("max")
    scaled_counts = (
        corrected_counts.div(row_maximum).where(row_maximum.gt(0.0), 0.0).fillna(0.0)
    )
    scaled_totals = scaled_counts.groupby(row_group_indices).transform("sum")
    corrected["corrected_frequency"] = scaled_counts.div(scaled_totals).where(
        row_maximum.gt(0.0),
        0.0,
    )

    # One unknown corrected count makes the row total unknown, so no motif in
    # that row has a defensible corrected frequency
    row_has_unknown_count = (
        corrected_counts.isna().groupby(row_group_indices).transform("any")
    )
    corrected.loc[row_has_unknown_count, "corrected_frequency"] = np.nan

    return corrected


def _row_order_indices(
    data_frame: pd.DataFrame,
    end_row_key_columns: list[str],
) -> np.ndarray:
    """
    Assign a stable integer position to each distinct correction row.

    Equal row-key tuples receive the same position, and first occurrence defines
    the order used to restore rows after merges and aggregation.
    """
    row_keys = _row_key_tuples(data_frame, end_row_key_columns)
    row_order_by_key = {
        row_key: row_order for row_order, row_key in enumerate(dict.fromkeys(row_keys))
    }
    return np.asarray(
        [row_order_by_key[row_key] for row_key in row_keys], dtype=np.int64
    )


def _selected_rows(
    ends: EndMotifCounts,
    *,
    window_idxs: int | Sequence[int] | None,
    groups: str | Sequence[str] | None,
    group_idxs: int | Sequence[int] | None,
    max_blacklisted_fraction: float,
) -> np.ndarray:
    """
    Resolve requested end-motif rows and apply the blacklist filter.

    Selection follows the end-motif object's row-mode rules. Filtering happens
    afterward so returned indices describe the rows represented in the result.
    """
    row_indices = ends._resolve_row_selector(window_idxs, groups, group_idxs)
    return filter_blacklisted_fraction(
        ends._row_metadata_data_frame(),
        max_blacklisted_fraction,
        row_indices=row_indices,
    )


def _selected_row_positions(
    ends: EndMotifCounts,
    row_indices: np.ndarray,
) -> dict[tuple[object, ...], int]:
    """
    Map selected correction-row keys to output matrix positions.

    Metadata is taken in selected-index order, so the mapping agrees with the
    dense and sparse matrix row order.
    """
    end_row_key_columns = _reference_correction_row_columns(ends.row_mode())
    row_metadata = (
        ends._row_metadata_data_frame().iloc[row_indices].reset_index(drop=True)
    )
    return {
        row_key: position
        for position, row_key in enumerate(
            _row_key_tuples(row_metadata, end_row_key_columns)
        )
    }


def _validate_reference_correction_inputs(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    use_global_bias: bool,
) -> None:
    """
    Validate objects and axes before any reference correction is calculated.

    Row modes must match unless a global reference is explicitly broadcast.
    Reference output must be non-canonical, motif-group axes must match, and
    concrete end-motif widths must equal the reference k-mer size.
    """
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
    if ref_kmers.canonical():
        raise ValueError(
            "Reference correction requires non-canonical reference k-mer output"
        )
    if ref_kmers.orientation() != "both":
        raise ValueError(
            "Reference correction requires reference k-mer output generated "
            "with `--orientation both`"
        )
    if ends.end_motifs.motif_axis_kind == "motif_group":
        if ref_kmers.motif_axis_kind() != "motif_group":
            raise ValueError(
                "Grouped end-motif output requires grouped reference k-mer output"
            )
        return

    if ref_kmers.motif_axis_kind() != "motif":
        raise ValueError(
            "End-motif output with motif labels requires reference k-mer "
            "output with motif labels"
        )
    reference_motif_lengths = (
        ends.motifs_metadata()["motif"].str.replace("_", "", regex=False).str.len()
    )
    if not (reference_motif_lengths == ref_kmers.kmer_size()).all():
        raise ValueError(
            f"End-motif width must match reference k-mer size ({ref_kmers.kmer_size()})"
        )


def _validate_matching_rows(
    ends: EndMotifCounts,
    ref_kmers: RefKmerFrequencies,
    end_row_key_columns: list[str],
    reference_row_key_columns: list[str],
) -> None:
    """
    Ensure end and reference outputs describe the same correction rows.

    Both sets of row keys must be unique and identical after sorting. No keyed
    comparison is needed when a global reference row is being broadcast.
    """
    if not reference_row_key_columns:
        return

    end_row_keys = (
        ends._row_metadata_data_frame()[end_row_key_columns]
        .drop_duplicates()
        .sort_values(end_row_key_columns)
        .reset_index(drop=True)
    )
    ref_row_keys = (
        ref_kmers._row_metadata_data_frame()[reference_row_key_columns]
        .drop_duplicates()
        .sort_values(reference_row_key_columns)
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
    """
    Return the metadata columns that identify a correction row.

    Global, windowed, and grouped outputs use different keys. Unknown row modes
    fail here instead of producing an ambiguous join.
    """
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
    """
    Choose the reference columns used to join correction rows.

    Broadcasting a global reference to non-global ends uses no join columns,
    which applies that single reference composition to every selected end row.
    Otherwise the reference row mode determines the keys.
    """
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
    """
    Resolve the reference rows needed for the selected end-motif rows.

    Keyed modes preserve the first-seen order of distinct selected row keys and
    require every key to exist. Global broadcasting selects the available global
    reference row directly.
    """
    if not reference_row_columns:
        return np.arange(len(ref_kmers._row_metadata_data_frame()), dtype=np.int64)

    reference_metadata = ref_kmers._row_metadata_data_frame()
    reference_indices_by_key = {
        row_key: row_index
        for row_index, row_key in enumerate(
            _row_key_tuples(reference_metadata, reference_row_columns)
        )
    }
    selected_row_keys = dict.fromkeys(_row_key_tuples(end_rows, reference_row_columns))
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
    """
    Load the complete reference motif axis for selected reference rows.

    All motifs are required because correction support and side marginals use
    the complete reference motif set, not only sample-selected motifs. Sparse stores read
    stored rows, while dense stores return their complete selected coordinates.
    """
    motif_indices = np.arange(len(ref_kmers.motifs_metadata()), dtype=np.int64)
    if ref_kmers.storage_mode() == "sparse_coo":
        return ref_kmers._stored_data_frame_for_indices(row_indices, motif_indices)
    return ref_kmers._complete_data_frame_for_indices(
        row_indices,
        motif_indices,
        densify=False,
    )


def _row_key_tuples(
    data_frame: pd.DataFrame,
    key_columns: list[str],
) -> list[tuple[object, ...]]:
    """
    Convert correction-row columns into ordered, hashable keys.

    The tuples are used consistently for membership checks, row matching, and
    matrix-position maps.
    """
    return list(data_frame[key_columns].itertuples(index=False, name=None))


def _add_empty_reference_correction_columns(data_frame: pd.DataFrame) -> pd.DataFrame:
    """
    Give an empty result the same correction columns as a non-empty result.

    Both columns are explicitly floating point so downstream array and matrix
    construction sees stable dtypes even when no rows were selected.
    """
    data_frame = data_frame.copy()
    data_frame["corrected_count"] = np.asarray([], dtype=float)
    data_frame["corrected_frequency"] = np.asarray([], dtype=float)
    return data_frame
