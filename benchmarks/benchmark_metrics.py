"""Numerical checks shared by the benchmark output comparisons.

Bounds allow an absolute overshoot of 1e-5 before clipping. Undefined metrics
are handled by their calculators, so unexpected NaN and infinity always fail here.
Profile arrays use SciPy's Pearson calculation. Coverage uses its accumulated
weighted covariance and variances without expanding genomic runs into arrays.
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
    """Finish coverage's weighted correlation from its streaming accumulators.

    Zero variance gives an undefined correlation. Nonfinite accumulators and
    negative variances indicate invalid arithmetic and must raise an error.
    These checks apply to the coverage accumulator, not SciPy's array calculation.
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
