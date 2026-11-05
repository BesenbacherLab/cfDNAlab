"""Numerical boundary checks with manually derived expected values."""

import importlib.util
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, correlation_value, pearson_correlation


class MetricValueTests(unittest.TestCase):
    """Reject invalid metrics while allowing small boundary rounding errors."""

    def test_accepts_finite_nonnegative_and_signed_metrics(self):
        self.assertEqual(checked_metric(2.5, "MAE"), 2.5)
        self.assertEqual(checked_metric(0.0, "RMSE"), 0.0)
        self.assertEqual(checked_metric(-2.5, "signed mean", minimum=None), -2.5)

    def test_clips_only_small_boundary_overshoots(self):
        # Fractions and normalized diversity range from 0 to 1
        for value, expected in ((-1e-5, 0.0), (0.25, 0.25), (1.00001, 1.0)):
            with self.subTest(value=value):
                self.assertEqual(
                    checked_metric(value, "fraction", maximum=1.0), expected
                )

    def test_rejects_excessive_overshoots_and_nonfinite_metrics(self):
        outside_lower = math.nextafter(-1e-5, -math.inf)
        outside_upper = math.nextafter(1.00001, math.inf)
        for value in (outside_lower, outside_upper, math.nan, math.inf, -math.inf):
            with self.subTest(value=value):
                with self.assertRaisesRegex(RuntimeError, "Invalid fraction"):
                    checked_metric(value, "fraction", maximum=1.0)


class CorrelationValueTests(unittest.TestCase):
    """Check raw coefficients before they can be clipped into the valid range."""

    def test_clips_small_overshoots_on_both_sides(self):
        # Unit variances give denominator 1, so covariance equals the raw r
        for covariance, expected in ((-1.00001, -1.0), (0.5, 0.5), (1.00001, 1.0)):
            with self.subTest(covariance=covariance):
                self.assertEqual(correlation_value(covariance, 1.0, 1.0), expected)

    def test_rejects_large_overshoots_and_invalid_accumulators(self):
        outside_limit = math.nextafter(1.00001, math.inf)
        cases = [(-outside_limit, 1.0, 1.0), (outside_limit, 1.0, 1.0)]
        for invalid in (math.nan, math.inf, -math.inf):
            cases.extend(
                ((invalid, 1.0, 1.0), (0.0, invalid, 1.0), (0.0, 1.0, invalid))
            )
        cases.extend(((0.0, -1.0, 1.0), (0.0, 1.0, -1.0)))
        for covariance, left_variance, right_variance in cases:
            with self.subTest(values=(covariance, left_variance, right_variance)):
                with self.assertRaisesRegex(RuntimeError, "Invalid correlation value"):
                    correlation_value(covariance, left_variance, right_variance)

    def test_returns_na_for_constant_profiles(self):
        for left_variance, right_variance in ((0.0, 1.0), (1.0, 0.0), (0.0, 0.0)):
            with self.subTest(variances=(left_variance, right_variance)):
                self.assertEqual(
                    correlation_value(0.0, left_variance, right_variance), "NA"
                )

    def test_avoids_overflow_from_multiplying_large_variances(self):
        # Equal covariance and variances give r=1. Their product exceeds float range
        self.assertAlmostEqual(correlation_value(1e308, 1e308, 1e308), 1.0)


@unittest.skipUnless(importlib.util.find_spec("numpy") is not None, "Requires NumPy")
class ArrayCorrelationTests(unittest.TestCase):
    """Check the Pearson calculation used by lengths, motifs and midpoints."""

    def test_calculates_known_positive_and_negative_correlations(self):
        # For (0,1,2) and (0,1,1), centered sums are 1, 2 and 2/3
        # The coefficient is 1/sqrt(2 * 2/3) = sqrt(3)/2
        self.assertAlmostEqual(
            pearson_correlation([0, 1, 2], [0, 1, 1]), math.sqrt(3) / 2
        )
        self.assertAlmostEqual(pearson_correlation([1, 2, 3], [3, 2, 1]), -1.0)
        self.assertAlmostEqual(
            pearson_correlation([1e12 + 1, 1e12 + 2, 1e12 + 3], [2, 4, 6]), 1.0
        )

    def test_returns_na_for_exactly_constant_profiles(self):
        # Repeated 0.1 values are constant even if their floating point mean differs
        self.assertEqual(pearson_correlation([0.1] * 11, list(range(11))), "NA")
        self.assertEqual(pearson_correlation([1, 2, 3], [2, 2, 2]), "NA")
        self.assertEqual(pearson_correlation([1], [2]), "NA")

    def test_rejects_invalid_inputs_and_invalid_covariance(self):
        for left, right in (
            ([], []),
            ([1, 2], [1]),
            ([1, math.nan], [1, 2]),
            ([1, 2], [1, math.inf]),
            ([0, 1e308], [0, 1]),
            ([0, 1e-200], [0, 1]),
        ):
            with self.subTest(left=left, right=right):
                with self.assertRaisesRegex((ValueError, RuntimeError), "correlation"):
                    pearson_correlation(left, right)


if __name__ == "__main__":
    unittest.main()
