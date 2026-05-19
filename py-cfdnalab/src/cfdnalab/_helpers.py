"""Internal helpers shared by cfDNAlab loaders."""

from __future__ import annotations

import numbers
from typing import Sequence

import numpy as np
import pandas as pd


def resolve_unique_match(
    condition: np.ndarray, *, missing_message: str, duplicate_message: str
) -> int:
    """
    Return the only index where a boolean condition is true.

    Use this for public lookup helpers where zero matches and multiple matches
    are both errors with different meanings. The caller owns the wording so the
    final error message can stay specific to groups, motifs, length bins, or
    any other lookup target.

    Parameters
    ----------
    condition
        Boolean array marking candidate matches.
    missing_message
        Message for the `KeyError` raised when there are no matches.
    duplicate_message
        Message for the `ValueError` raised when there are multiple matches.

    Returns
    -------
    int
        Index of the single matching element.
    """
    matches = np.flatnonzero(condition)
    if len(matches) == 0:
        raise KeyError(missing_message)
    if len(matches) > 1:
        raise ValueError(duplicate_message)
    return int(matches[0])


def validate_fragment_length(length: int) -> int:
    """
    Validate a fragment length used for length-bin lookup.

    Python booleans are integer subclasses, but treating `True` as length 1
    would silently select the wrong bin.

    Parameters
    ----------
    length
        Fragment length in bp.

    Returns
    -------
    int
        Non-negative integer fragment length.
    """
    if isinstance(length, bool) or not isinstance(length, numbers.Integral):
        raise TypeError(
            f"Fragment length must be an integer, got {type(length).__name__}"
        )
    length = int(length)
    if length < 0:
        raise ValueError(f"Fragment length must be non-negative, got {length}")
    return length


def validate_zero_based_index(index: int, size: int, name: str) -> int:
    """
    Validate a public Python index.

    Parameters
    ----------
    index
        User-supplied index.
    size
        Axis size.
    name
        Index name used in error messages.

    Returns
    -------
    int
        Validated index.
    """
    if isinstance(index, bool) or not isinstance(index, numbers.Integral):
        raise TypeError(f"{name} must be an integer, got {type(index).__name__}")
    index = int(index)
    if index < 0 or index >= size:
        raise IndexError(f"{name} {index} is outside 0..{size - 1}")
    return index


def validate_scalar_bool(value: bool, name: str) -> bool:
    """
    Validate a boolean option.
    """
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a bool, got {type(value).__name__}")
    return value


def validate_unique_values(values: Sequence[object], name: str) -> None:
    """
    Reject duplicate user-requested selector values.

    Duplicated selectors usually indicate a mistake and would otherwise return
    duplicated rows that are hard to distinguish from true repeated metadata.
    """
    seen: set[object] = set()
    for value in values:
        if value in seen:
            raise ValueError(f"{name} contains duplicate values")
        seen.add(value)


def normalize_zero_based_indices(
    indices: int | Sequence[int] | None,
    *,
    size: int,
    name: str,
    index_name: str,
) -> np.ndarray:
    """
    Normalize an optional scalar-or-vector index selector.
    """
    if indices is None:
        return np.arange(size, dtype=np.int64)
    if isinstance(indices, numbers.Integral):
        values = [validate_zero_based_index(indices, size, index_name)]
    else:
        try:
            values = [
                validate_zero_based_index(index, size, index_name)
                for index in indices
            ]
        except TypeError as error:
            raise TypeError(
                f"{name} must be an integer or sequence of integers"
            ) from error
    validate_unique_values(values, name)
    return np.asarray(values, dtype=np.int64)


def normalize_strings(
    values: str | Sequence[str],
    *,
    name: str,
) -> list[str]:
    """
    Normalize a scalar-or-vector string selector.
    """
    if isinstance(values, str):
        items = [values]
    else:
        try:
            items = list(values)
        except TypeError as error:
            raise TypeError(
                f"{name} must be a string or sequence of strings"
            ) from error
    if any(not isinstance(value, str) for value in items):
        raise TypeError(f"{name} must contain strings")
    validate_unique_values(items, name)
    return items


def validate_fraction(value: float, name: str) -> float:
    """
    Validate one scalar fraction in the closed interval 0..1.
    """
    if (
        isinstance(value, bool)
        or not isinstance(value, numbers.Real)
        or not np.isfinite(value)
        or value < 0
        or value > 1
    ):
        raise ValueError(f"{name} must be a single finite fraction in 0..1")
    return float(value)


def filter_blacklisted_fraction(
    data_frame: pd.DataFrame,
    max_blacklisted_fraction: float,
    *,
    row_indices: np.ndarray | None = None,
) -> pd.DataFrame | np.ndarray:
    """
    Filter rows or selected row indices by `blacklisted_fraction`.

    When `row_indices` is supplied, the returned value is a filtered copy of
    those indices. Otherwise the returned value is a filtered data frame.
    """
    max_blacklisted_fraction = validate_fraction(
        max_blacklisted_fraction, "max_blacklisted_fraction"
    )
    if "blacklisted_fraction" not in data_frame.columns:
        if max_blacklisted_fraction == 1:
            return data_frame if row_indices is None else row_indices
        raise ValueError(
            "Cannot filter by max_blacklisted_fraction because this output has no "
            "blacklisted_fraction column"
        )

    blacklist_values = data_frame["blacklisted_fraction"].to_numpy()
    if row_indices is None:
        keep = blacklist_values <= max_blacklisted_fraction
        return data_frame.loc[keep].reset_index(drop=True)

    keep = blacklist_values[row_indices] <= max_blacklisted_fraction
    return row_indices[keep]
