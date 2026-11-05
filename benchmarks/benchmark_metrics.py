"""Numerical checks shared by the benchmark output comparisons.

Bounds allow an absolute overshoot of 1e-5 before clipping. Undefined metrics
are handled by their calculators, so unexpected NaN and infinity always fail here.
Only the small-array Pearson calculation needs NumPy. Coverage can use the
scalar checks without installing it or expanding coverage into genomic arrays.
"""

import math


def checked_metric(value, name, minimum=0.0, maximum=None):
    """Require a finite metric and clip only small errors at the supplied bounds.

    The default requires a nonnegative result. Signed metrics use minimum=None.
    Fractions and other normalized metrics additionally supply maximum=1.0.
    """
    if not math.isfinite(value):
        raise RuntimeError(f"Invalid {name}: {value} (expected a finite value)")
    if minimum is not None:
        if value < minimum - 1e-5:
            raise RuntimeError(f"Invalid {name}: {value} (minimum {minimum})")
        value = max(minimum, value)
    if maximum is not None:
        if value > maximum + 1e-5:
            raise RuntimeError(f"Invalid {name}: {value} (maximum {maximum})")
        value = min(maximum, value)
    return value


def correlation_value(covariance, left_variance, right_variance):
    """Check covariance and variances, then calculate and validate Pearson r.

    Zero variance gives an undefined correlation. Nonfinite accumulators and
    negative variances indicate invalid arithmetic and must raise an error.
    """
    if (
        not all(
            math.isfinite(value)
            for value in (covariance, left_variance, right_variance)
        )
        or left_variance < 0
        or right_variance < 0
    ):
        raise RuntimeError(
            "Invalid correlation value: nonfinite covariance/variance or negative variance"
        )
    if left_variance == 0 or right_variance == 0:
        return "NA"
    # Dividing by each square root avoids overflowing the product of variances
    value = covariance / math.sqrt(left_variance) / math.sqrt(right_variance)
    return checked_metric(value, "correlation value", minimum=-1.0, maximum=1.0)


def pearson_correlation(left, right):
    """Calculate Pearson for small profiles and check r before clipping.

    NumPy's corrcoef clips internally, which would hide excessive overshoots.
    Covariance leaves the raw calculation available for the shared bounds check.
    """
    import numpy as np

    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.ndim != 1 or right.ndim != 1 or left.shape != right.shape or left.size == 0:
        raise ValueError("Expected two nonempty, equally sized profiles for correlation")
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("Nonfinite profile values in correlation")
    # Exact equality identifies constant profiles even when their mean rounds
    if (left == left[0]).all() or (right == right[0]).all():
        return "NA"
    covariance = np.cov(left, right, ddof=0)
    if covariance[0, 0] == 0 or covariance[1, 1] == 0:
        raise RuntimeError("Invalid correlation value: nonconstant profile has zero variance")
    return correlation_value(covariance[0, 1], covariance[0, 0], covariance[1, 1])
