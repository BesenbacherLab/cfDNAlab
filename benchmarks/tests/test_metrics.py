"""Scalar metric checks with manually derived expected values.

Profile correlations are covered by the length, motif and midpoint comparator tests.
"""

import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, correlation_value


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


if __name__ == "__main__":
    unittest.main()
