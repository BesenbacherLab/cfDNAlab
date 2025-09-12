
// TODO: Check manually - generated but not validated!

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use cfdnalab::utils::{coverage::CoveragePrefix, fragment::Fragment};

    // Simple approx helpers
    fn feq(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }
    fn deq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn doc_example_pipeline() -> Result<()> {
        let length: u32 = 300;
        let mut cp = CoveragePrefix::initialize_coverage_prefix(length);

        // Unweighted and weighted fragments
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 100,
            end: 200,
        })?;
        cp.add_fragment_to_prefix_weighted(
            Fragment {
                tid: 0,
                start: 150,
                end: 250,
            },
            0.87,
        )?;

        // Optional blacklist
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(120, 140)?;
        cp.finalize_blacklist_prefix();

        // Build per-base coverage and indexes
        cp.finalize_coverage();
        cp.build_query_index()?;

        // Coverage length matches sequence length
        let cov = cp.coverage().unwrap();
        assert_eq!(cov.len() as u32, length);

        // Spot check coverage values around boundaries
        assert!(feq(cov[99], 0.0, 1e-6));
        assert!(feq(cov[100], 1.0, 1e-6));
        assert!(feq(cov[119], 1.0, 1e-6));
        assert!(feq(cov[120], 1.0, 1e-6)); // masked later but raw coverage is 1.0
        assert!(feq(cov[149], 1.0, 1e-6));
        assert!(feq(cov[150], 1.87, 1e-6));
        assert!(feq(cov[199], 1.87, 1e-6));
        assert!(feq(cov[200], 0.87, 1e-6));
        assert!(feq(cov[249], 0.87, 1e-6));
        assert!(feq(cov[250], 0.0, 1e-6));

        // Sum and averages over [100, 300)
        // Expected sums:
        // [100,120): 20 * 1.0 = 20.0
        // [120,140): 20 * 1.0 = 20.0 (masked segment)
        // [140,150): 10 * 1.0 = 10.0
        // [150,200): 50 * 1.87 = 93.5
        // [200,250): 50 * 0.87 = 43.5
        // [250,300): 50 * 0.0 = 0.0
        // Total including masked = 187.0; excluding masked = 167.0
        let sum_all = cp.sum_coverage(100, 300, false)?;
        let sum_ok = cp.sum_coverage(100, 300, true)?;
        assert!(deq(sum_all, 187.0, 1e-9));
        assert!(deq(sum_ok, 167.0, 1e-9));

        let avg_all = cp.avg_coverage(100, 300, false)?;
        let avg_ok = cp.avg_coverage(100, 300, true)?;
        assert!(feq(avg_all, 187.0 / 200.0, 1e-6));
        assert!(feq(avg_ok, 167.0 / 180.0, 1e-6)); // 20 masked bases removed from denominator

        // Position queries with NaN for masked
        let ys = cp.coverage_at_positions_nan(&[119, 120, 139, 140, 150])?;
        assert!(!ys[0].is_nan());
        assert!(ys[1].is_nan()); // 120 masked
        assert!(ys[2].is_nan()); // 139 masked
        assert!(!ys[3].is_nan()); // 140 not masked
        assert!(!ys[4].is_nan());

        // Mask at positions
        let ms = cp.mask_at_positions(&[119, 120, 139, 140])?;
        assert_eq!(ms, vec![0, 1, 1, 0]);

        Ok(())
    }

    #[test]
    fn errors_before_finalize_and_mask_requirements() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(1000);

        // Using averages before finalize_coverage should error
        let err = cp.avg_coverage(0, 10, false).unwrap_err();
        assert!(format!("{err}").contains("coverage not finalized"));

        // Add coverage and finalize
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 10,
        })?;
        cp.finalize_coverage();

        // Excluding blacklisted requires finalized mask if a blacklist exists
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(2, 5)?;
        let err = cp.avg_coverage(0, 10, true).unwrap_err();
        assert!(format!("{err}").contains("blacklist present but not finalized"));

        // Finalize mask and now it should work
        cp.finalize_blacklist_prefix();
        let _ = cp.avg_coverage(0, 10, true)?; // no panic

        Ok(())
    }

    #[test]
    fn add_fragment_after_finalize_requires_refinalize() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(100);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        cp.finalize_coverage();
        cp.build_query_index()?;

        // Add another fragment; coverage should be invalidated
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 20,
            end: 30,
        })?;
        // Now any query should complain coverage not finalized
        let err = cp.sum_coverage(0, 40, false).unwrap_err();
        assert!(format!("{err}").contains("coverage not finalized"));

        // Re-finalize and query again
        cp.finalize_coverage();
        cp.build_query_index()?;
        let sum = cp.sum_coverage(0, 40, false)?;
        // Expected sum = 10 + 10 = 20
        assert!(deq(sum, 20.0, 1e-9));
        Ok(())
    }

    #[test]
    fn drop_prefix_blocks_additions() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(50);
        cp.drop_prefix();
        let err = cp
            .add_fragment_to_prefix(Fragment {
                tid: 0,
                start: 0,
                end: 1,
            })
            .unwrap_err();
        assert!(format!("{err}").contains("prefix was dropped"));
        Ok(())
    }

    #[test]
    fn bulk_queries_parallel_and_serial_match() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(1000);
        // Create a few simple fragments
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 110,
        })?;
        cp.add_fragment_to_prefix_weighted(
            Fragment {
                tid: 0,
                start: 200,
                end: 400,
            },
            0.5,
        )?;
        cp.finalize_coverage();
        // No blacklist
        cp.build_query_index()?;

        let intervals = vec![(0, 10), (10, 110), (100, 300), (350, 450), (900, 1000)];
        let sums_ser = cp.bulk_sum_coverage(&intervals, false, false)?;
        let sums_par = cp.bulk_sum_coverage(&intervals, false, true)?;
        assert_eq!(sums_ser.len(), intervals.len());
        assert_eq!(sums_par.len(), intervals.len());
        for i in 0..intervals.len() {
            assert!(deq(sums_ser[i], sums_par[i], 1e-9));
        }

        let avgs_ser = cp.bulk_avg_coverage(&intervals, false, false)?;
        let avgs_par = cp.bulk_avg_coverage(&intervals, false, true)?;
        assert_eq!(avgs_ser.len(), intervals.len());
        assert_eq!(avgs_par.len(), intervals.len());
        for i in 0..intervals.len() {
            assert!(feq(avgs_ser[i], avgs_par[i], 1e-6));
        }

        Ok(())
    }

    #[test]
    fn blacklist_finalize_is_non_destructive_and_affects_queries() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(100);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        cp.finalize_coverage();

        // Build a blacklist delta and clone it
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(12, 18)?;
        let bl_before = cp._get_bl_delta().clone(); // Access is allowed in submodule tests

        // Finalize mask; delta should remain unchanged
        cp.finalize_blacklist_prefix();
        let bl_after = cp._get_bl_delta().clone();
        assert_eq!(bl_before, bl_after);

        // Build indexes and check effect on queries
        cp.build_query_index()?;
        let sum_all = cp.sum_coverage(10, 20, false)?;
        let sum_ok = cp.sum_coverage(10, 20, true)?;
        // Sum without excluding is 10 * 1.0 = 10. Excluding removes [12,18) length 6
        assert!(deq(sum_all, 10.0, 1e-9));
        assert!(deq(sum_ok, 4.0, 1e-9));
        Ok(())
    }

    #[test]
    fn coverage_at_positions_and_mask() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(60);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        cp.finalize_coverage();

        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(12, 15)?;
        cp.finalize_blacklist_prefix();

        let vals = cp.coverage_at_positions(&[9, 10, 12, 14, 15, 19, 20])?;
        assert_eq!(vals, vec![0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]);

        let vals_nan = cp.coverage_at_positions_nan(&[12, 13, 14, 15])?;
        assert!(vals_nan[0].is_nan());
        assert!(vals_nan[1].is_nan());
        assert!(vals_nan[2].is_nan());
        assert!(!vals_nan[3].is_nan());

        let mask = cp.mask_at_positions(&[11, 12, 14, 15])?;
        assert_eq!(mask, vec![0, 1, 1, 0]);
        Ok(())
    }

    #[test]
    fn invalid_inputs_are_rejected() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(50);

        // Negative weight
        let err = cp
            .add_fragment_to_prefix_weighted(
                Fragment {
                    tid: 0,
                    start: 0,
                    end: 10,
                },
                -0.1,
            )
            .unwrap_err();
        assert!(format!("{err}").contains("invalid weight"));

        // Start >= end for fragment
        let err = cp
            .add_fragment_to_prefix(Fragment {
                tid: 0,
                start: 10,
                end: 10,
            })
            .unwrap_err();
        assert!(format!("{err}").contains("start 10 >= end 10"));

        // Out-of-bounds blacklist
        cp.initialize_blacklist_prefix();
        let err = cp.add_blacklist_to_prefix(45, 60).unwrap_err();
        assert!(format!("{err}").contains("out of bounds"));

        // Bounds check in queries
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 10,
        })?;
        cp.finalize_coverage();
        let err = cp.sum_coverage(10, 60, false).unwrap_err();
        assert!(format!("{err}").contains("exceeds sequence length"));

        Ok(())
    }

    // Append these to your existing `#[cfg(test)] mod tests`
    //
    // Note: These tests assume the current API semantics in your snippet,
    // including that `finalize_coverage()` is infallible and will panic if the prefix was dropped.

    use std::panic::{AssertUnwindSafe, catch_unwind};

    #[test]
    fn empty_sequence_finalize_and_query() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(0);
        // Finalize and build indexes on empty sequence
        let cov = cp.finalize_coverage();
        assert_eq!(cov.len(), 0);
        cp.build_query_index()?;
        // Bulk queries on empty set of intervals
        let sums = cp.bulk_sum_coverage(&[], false, false)?;
        let avgs = cp.bulk_avg_coverage(&[], false, false)?;
        assert!(sums.is_empty());
        assert!(avgs.is_empty());
        Ok(())
    }

    #[test]
    fn single_base_coverage_and_queries() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(1);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 1,
        })?;
        let cov = cp.finalize_coverage();
        assert_eq!(cov, &[1.0]);
        cp.build_query_index()?;
        assert!(deq(cp.sum_coverage(0, 1, false)?, 1.0, 1e-12));
        assert!(feq(cp.avg_coverage(0, 1, false)?, 1.0, 1e-6));
        // Zero-width interval average
        assert!(feq(cp.avg_coverage(1, 1, false)?, 0.0, 1e-6));
        Ok(())
    }

    #[test]
    fn exclude_without_blacklist_equals_include() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(200);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        cp.finalize_coverage();
        cp.build_query_index()?;
        // No blacklist present, excluding should equal including
        let a = cp.sum_coverage(0, 200, false)?;
        let b = cp.sum_coverage(0, 200, true)?;
        assert!(deq(a, b, 1e-12));
        let a = cp.avg_coverage(0, 200, false)?;
        let b = cp.avg_coverage(0, 200, true)?;
        assert!(feq(a, b, 1e-6));
        Ok(())
    }

    #[test]
    fn finalize_blacklist_without_intervals_is_noop() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(100);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 10,
        })?;
        cp.finalize_coverage();
        cp.initialize_blacklist_prefix();
        cp.finalize_blacklist_prefix(); // no intervals added
        cp.build_query_index()?;
        let inc = cp.sum_coverage(0, 100, false)?;
        let exc = cp.sum_coverage(0, 100, true)?;
        assert!(deq(inc, exc, 1e-12));
        Ok(())
    }

    #[test]
    fn idempotent_build_query_index() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(50);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 5,
            end: 15,
        })?;
        cp.finalize_coverage();
        cp.build_query_index()?;
        let s1 = cp.sum_coverage(0, 50, false)?;
        let a1 = cp.avg_coverage(0, 50, false)?;
        // Rebuild indexes again
        cp.build_query_index()?;
        let s2 = cp.sum_coverage(0, 50, false)?;
        let a2 = cp.avg_coverage(0, 50, false)?;
        assert!(deq(s1, s2, 1e-12));
        assert!(feq(a1, a2, 1e-6));
        Ok(())
    }

    #[test]
    fn bulk_empty_intervals() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(10);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 2,
            end: 5,
        })?;
        cp.finalize_coverage();
        cp.build_query_index()?;
        let sums = cp.bulk_sum_coverage(&[], false, false)?;
        let avgs = cp.bulk_avg_coverage(&[], false, false)?;
        assert!(sums.is_empty());
        assert!(avgs.is_empty());
        Ok(())
    }

    #[test]
    fn position_bounds_and_errors() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(5);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 1,
            end: 4,
        })?;
        cp.finalize_coverage();
        // In-bounds positions
        let vals = cp.coverage_at_positions(&[0, 1, 3, 4])?;
        assert_eq!(vals, vec![0.0, 1.0, 1.0, 0.0]);
        // Out-of-bounds should error
        let err = cp.coverage_at_positions(&[5]).unwrap_err();
        assert!(format!("{err}").contains("out of bounds"));
        Ok(())
    }

    #[test]
    fn blacklist_affects_queries_as_expected() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(100);
        // Coverage segments: [10,30)=1.0 and [40,90)=0.5
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 30,
        })?;
        cp.add_fragment_to_prefix_weighted(
            Fragment {
                tid: 0,
                start: 40,
                end: 90,
            },
            0.5,
        )?;
        cp.finalize_coverage();
        // Blacklist [20,25) and [80,100)
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(20, 25)?;
        cp.add_blacklist_to_prefix(80, 100)?;
        cp.finalize_blacklist_prefix();
        cp.build_query_index()?;

        // Sum including mask
        let s_all = cp.sum_coverage(0, 100, false)?;
        // Manual sum: 20*1 + 50*0.5 = 45.0
        assert!(deq(s_all, 45.0, 1e-12));

        // Excluding masked removes [20,25) and [80,100) from numerator
        let s_exc = cp.sum_coverage(0, 100, true)?;
        // Manual excluding removes 5*1 + 10*0.5 = 10
        assert!(deq(s_exc, 35.0, 1e-12));

        // Averages
        let a_all = cp.avg_coverage(0, 100, false)?;
        assert!(feq(a_all, 45.0 / 100.0, 1e-6));
        let a_exc = cp.avg_coverage(0, 100, true)?;
        // Denominator excludes 5 + 20 = 25 masked bases
        assert!(feq(a_exc, 35.0 / 75.0, 1e-6));

        Ok(())
    }

    #[test]
    fn updating_blacklist_requires_refinalize_to_affect_results() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(100);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 100,
        })?;
        cp.finalize_coverage();

        // First mask [10,20)
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(10, 20)?;
        cp.finalize_blacklist_prefix();
        cp.build_query_index()?;

        let s_exc1 = cp.sum_coverage(0, 100, true)?;
        assert!(deq(s_exc1, 90.0, 1e-12));

        // Edit blacklist to also include [30,40) but do not finalize yet
        cp.add_blacklist_to_prefix(30, 40)?;
        // Rebuild indexes now uses stale mask since we didn't finalize again
        // We call build_query_index to simulate callers that rebuild for safety
        cp.build_query_index()?;

        // Panics as we did not finalize blacklist again after adding
        let err = cp.sum_coverage(0, 100, true).unwrap_err();
        assert!(format!("{err}").contains("blacklist present but not finalized"));

        // Finalize mask and rebuild indexes to apply the change
        cp.finalize_blacklist_prefix();
        cp.build_query_index()?;
        let s_exc2 = cp.sum_coverage(0, 100, true)?;
        assert!(deq(s_exc2, 80.0, 1e-12));
        Ok(())
    }

    #[test]
    fn finalize_twice_is_stable() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(30);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        let c1 = cp.finalize_coverage().to_vec();
        let c2 = cp.finalize_coverage().to_vec();
        assert_eq!(c1, c2);
        Ok(())
    }

    #[test]
    fn bulk_parallel_vs_serial_equivalence_with_mask() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(1000);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 500,
        })?;
        cp.add_fragment_to_prefix_weighted(
            Fragment {
                tid: 0,
                start: 250,
                end: 750,
            },
            0.5,
        )?;
        cp.finalize_coverage();

        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(400, 450)?;
        cp.add_blacklist_to_prefix(700, 900)?;
        cp.finalize_blacklist_prefix();

        cp.build_query_index()?;

        let intervals: Vec<(u32, u32)> = vec![
            (0, 0),
            (0, 1000),
            (250, 260),
            (390, 410),  // spans into masked region
            (440, 460),  // mostly masked
            (700, 900),  // fully masked
            (900, 1000), // masked then unmasked
        ];

        let s_ser = cp.bulk_sum_coverage(&intervals, true, false)?;
        let s_par = cp.bulk_sum_coverage(&intervals, true, true)?;
        for i in 0..intervals.len() {
            assert!(deq(s_ser[i], s_par[i], 1e-12));
        }

        let a_ser = cp.bulk_avg_coverage(&intervals, true, false)?;
        let a_par = cp.bulk_avg_coverage(&intervals, true, true)?;
        for i in 0..intervals.len() {
            assert!(feq(a_ser[i], a_par[i], 1e-6));
        }
        Ok(())
    }

    #[test]
    fn finalize_after_drop_prefix_panics() {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(10);
        cp.drop_prefix();
        // finalize_coverage currently assumes the prefix exists and will panic
        let panicked = catch_unwind(AssertUnwindSafe(|| {
            let _ = cp.finalize_coverage();
        }))
        .is_err();
        assert!(panicked);
    }

    #[test]
    fn queries_cover_edges_exactly() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(10);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 0,
            end: 10,
        })?;
        cp.finalize_coverage();
        cp.build_query_index()?;
        // Full range
        assert!(deq(cp.sum_coverage(0, 10, false)?, 10.0, 1e-12));
        // Left edge zero-width
        assert!(feq(cp.avg_coverage(0, 0, false)?, 0.0, 1e-6));
        // Right edge zero-width
        assert!(feq(cp.avg_coverage(10, 10, false)?, 0.0, 1e-6));
        // Single base at last position
        assert!(feq(cp.avg_coverage(9, 10, false)?, 1.0, 1e-6));
        Ok(())
    }

    #[test]
    fn manual_vs_indexed_sum_consistency() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(200);
        // Coverage 1.0 on [20,60) and 0.5 on [100,150)
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 20,
            end: 60,
        })?;
        cp.add_fragment_to_prefix_weighted(
            Fragment {
                tid: 0,
                start: 100,
                end: 150,
            },
            0.5,
        )?;
        let cov = cp.finalize_coverage().to_vec();
        cp.build_query_index()?;

        let intervals = vec![(0, 200), (0, 20), (20, 60), (50, 120), (140, 160)];
        for &(a, b) in &intervals {
            let sum_manual: f64 = cov[a as usize..b as usize].iter().map(|&x| x as f64).sum();
            let sum_indexed = cp.sum_coverage(a, b, false)?;
            assert!(deq(sum_manual, sum_indexed, 1e-9));
        }
        Ok(())
    }

    #[test]
    fn mask_positions_nan_semantics() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(30);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        cp.finalize_coverage();

        // No blacklist yields no NaNs
        let v = cp.coverage_at_positions_nan(&[9, 10, 19, 20])?;
        for x in &v {
            assert!(!x.is_nan());
        }

        // With blacklist, NaNs appear inside masked region
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_to_prefix(12, 15)?;
        cp.finalize_blacklist_prefix();
        let v = cp.coverage_at_positions_nan(&[11, 12, 13, 14, 15])?;
        assert!(!v[0].is_nan());
        assert!(v[1].is_nan());
        assert!(v[2].is_nan());
        assert!(v[3].is_nan());
        assert!(!v[4].is_nan());
        Ok(())
    }

    #[test]
    fn add_fragment_after_indexes_invalidates_and_requires_refinalize() -> Result<()> {
        let mut cp = CoveragePrefix::initialize_coverage_prefix(100);
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 10,
            end: 20,
        })?;
        cp.finalize_coverage();
        cp.build_query_index()?;
        let s1 = cp.sum_coverage(0, 100, false)?;
        assert!(deq(s1, 10.0, 1e-12));

        // Add more fragments, which should invalidate coverage and indexes
        cp.add_fragment_to_prefix(Fragment {
            tid: 0,
            start: 20,
            end: 30,
        })?;
        // Now querying should fail because coverage not finalized
        let err = cp.avg_coverage(0, 100, false).unwrap_err();
        assert!(format!("{err}").contains("coverage not finalized"));

        // Finalize again and rebuild
        cp.finalize_coverage();
        cp.build_query_index()?;
        let s2 = cp.sum_coverage(0, 100, false)?;
        assert!(deq(s2, 20.0, 1e-12));
        Ok(())
    }
}
