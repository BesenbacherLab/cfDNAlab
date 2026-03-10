#![cfg(feature = "cmd_fcoverage")]

// TODO: Check manually - generated but not validated!

#[cfg(test)]
mod tests_coverage_prefix {
    use anyhow::Result;
    use cfdnalab::shared::{
        coverage::Coverage,
        fragment::{
            minimal_fragment::Fragment,
            segment_fragment::{SegmentedReadInfo, collect_fragment_with_segments},
        },
    };
    use std::panic::{AssertUnwindSafe, catch_unwind};

    // Simple approx helpers
    fn feq(a: f32, b: f32, eps: f32) -> bool {
        (a - b).abs() <= eps
    }
    fn deq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    // SegmentedReadInfo creator
    fn sri(
        tid: i32,
        pos: u32,
        end: u32,
        is_reverse: bool,
        has_ref_gap: bool,
        max_ref_gap: u32,
        segs: &[(u32, u32)],
    ) -> SegmentedReadInfo {
        SegmentedReadInfo {
            tid,
            pos,
            end,
            is_reverse,
            has_ref_gap,
            max_ref_gap,
            ref_mapped_segments: segs.to_vec(),
            gc_tag: Default::default(),
        }
    }

    fn new_cp(len: u32) -> Coverage {
        Coverage::new(len)
    }

    #[test]
    fn doc_example_pipeline() -> Result<()> {
        let length: u32 = 300;
        let mut cp = Coverage::new(length);

        // Unweighted and weighted fragments
        cp.add_fragment(Fragment {
            tid: 0,
            start: 100,
            end: 200,
            gc_tag: Default::default(),
        })?;
        cp.add_fragment_weighted(
            Fragment {
                tid: 0,
                start: 150,
                end: 250,
                gc_tag: Default::default(),
            },
            0.87,
        )?;

        // Optional blacklist
        cp.set_blacklist_mask(&vec![(120, 140)])?;

        // Build per-base coverage and indexes
        cp.finalize_coverage(true);
        cp.build_indexes(false)?;

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

        let cov_vec = cp.coverage().unwrap().to_vec(); // clones, borrow ends here
        let manual = cov_vec[100..300].iter().map(|&x| x as f64).sum::<f64>();
        assert!(deq(sum_all, manual, 1e-9)); // should match 187.0 here

        println!("{:?}", sum_all);
        assert!(deq(sum_all, 187.0, 1e-6));
        assert!(deq(sum_ok, 167.0, 1e-6));

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
    fn add_fragment_after_finalize_requires_refinalize() -> Result<()> {
        let mut cp = Coverage::new(100);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(false);
        cp.build_indexes(false)?;

        // Add another fragment; coverage should be invalidated
        cp.add_fragment(Fragment {
            tid: 0,
            start: 20,
            end: 30,
            gc_tag: Default::default(),
        })?;
        // Now any query should complain coverage not finalized
        let err = cp.sum_coverage(0, 40, false).unwrap_err();
        assert!(format!("{err}").contains("coverage not finalized"));

        // Re-finalize and query again
        cp.finalize_coverage(false);
        cp.build_indexes(false)?;
        let sum = cp.sum_coverage(0, 40, false)?;
        // Expected sum = 10 + 10 = 20
        assert!(deq(sum, 20.0, 1e-9));
        Ok(())
    }

    #[test]
    fn drop_deltas_blocks_additions() -> Result<()> {
        let mut cp = Coverage::new(50);
        cp.drop_deltas();
        let err = cp
            .add_fragment(Fragment {
                tid: 0,
                start: 0,
                end: 1,
                gc_tag: Default::default(),
            })
            .unwrap_err();
        assert!(format!("{err}").contains("prefix was dropped"));
        Ok(())
    }

    #[test]
    fn bulk_queries_parallel_and_serial_match() -> Result<()> {
        let mut cp = Coverage::new(1000);
        // Create a few simple fragments
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 110,
            gc_tag: Default::default(),
        })?;
        cp.add_fragment_weighted(
            Fragment {
                tid: 0,
                start: 200,
                end: 400,
                gc_tag: Default::default(),
            },
            0.5,
        )?;
        cp.finalize_coverage(true);
        // No blacklist
        cp.build_indexes(true)?;

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
    fn coverage_at_positions_and_mask() -> Result<()> {
        let mut cp = Coverage::new(60);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);

        cp.set_blacklist_mask(&vec![(12, 15)])?;

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
        let mut cp = Coverage::new(50);

        // Negative weight
        let err = cp
            .add_fragment_weighted(
                Fragment {
                    tid: 0,
                    start: 0,
                    end: 10,
                    gc_tag: Default::default(),
                },
                -0.1,
            )
            .unwrap_err();
        assert!(format!("{err}").contains("invalid weight"));

        // Start >= end for fragment
        let err = cp
            .add_fragment(Fragment {
                tid: 0,
                start: 10,
                end: 10,
                gc_tag: Default::default(),
            })
            .unwrap_err();
        assert!(format!("{err}").contains("start 10 >= end 10"));

        // Out-of-bounds blacklist
        let err = cp.set_blacklist_mask(&vec![(45, 60)]).unwrap_err();
        assert!(format!("{err}").contains("out of bounds"));

        // Bounds check in queries
        cp.add_fragment(Fragment {
            tid: 0,
            start: 0,
            end: 10,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);
        let err = cp.sum_coverage(10, 60, false).unwrap_err();
        assert!(format!("{err}").contains("exceeds sequence length"));

        Ok(())
    }

    #[test]
    fn empty_sequence_finalize_and_query() -> Result<()> {
        let mut cp = Coverage::new(0);
        // Finalize and build indexes on empty sequence
        let cov = cp.finalize_coverage(true);
        assert_eq!(cov.len(), 0);
        cp.build_indexes(true)?;
        // Bulk queries on empty set of intervals
        let sums = cp.bulk_sum_coverage(&[], false, false)?;
        let avgs = cp.bulk_avg_coverage(&[], false, false)?;
        assert!(sums.is_empty());
        assert!(avgs.is_empty());
        Ok(())
    }

    #[test]
    fn single_base_coverage_and_queries() -> Result<()> {
        let mut cp = Coverage::new(1);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 0,
            end: 1,
            gc_tag: Default::default(),
        })?;
        let cov = cp.finalize_coverage(true);
        assert_eq!(cov, &[1.0]);
        cp.build_indexes(true)?;
        assert!(deq(cp.sum_coverage(0, 1, false)?, 1.0, 1e-12));
        assert!(feq(cp.avg_coverage(0, 1, false)?, 1.0, 1e-6));
        // Zero-width interval average
        assert!(feq(cp.avg_coverage(1, 1, false)?, 0.0, 1e-6));
        Ok(())
    }

    #[test]
    fn exclude_without_blacklist_equals_include() -> Result<()> {
        let mut cp = Coverage::new(200);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;
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
    fn idempotent_build_indexes() -> Result<()> {
        let mut cp = Coverage::new(50);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 5,
            end: 15,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);
        cp.build_indexes(false)?;
        let s1 = cp.sum_coverage(0, 50, false)?;
        let a1 = cp.avg_coverage(0, 50, false)?;
        // Rebuild indexes again
        cp.build_indexes(true)?;
        let s2 = cp.sum_coverage(0, 50, false)?;
        let a2 = cp.avg_coverage(0, 50, false)?;
        assert!(deq(s1, s2, 1e-12));
        assert!(feq(a1, a2, 1e-6));
        Ok(())
    }

    #[test]
    fn bulk_empty_intervals() -> Result<()> {
        let mut cp = Coverage::new(10);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 2,
            end: 5,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;
        let sums = cp.bulk_sum_coverage(&[], false, false)?;
        let avgs = cp.bulk_avg_coverage(&[], false, false)?;
        assert!(sums.is_empty());
        assert!(avgs.is_empty());
        Ok(())
    }

    #[test]
    fn position_bounds_and_errors() -> Result<()> {
        let mut cp = Coverage::new(5);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 1,
            end: 4,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);
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
        let mut cp = Coverage::new(100);
        // Coverage segments: [10,30)=1.0 and [40,90)=0.5
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 30,
            gc_tag: Default::default(),
        })?;
        cp.add_fragment_weighted(
            Fragment {
                tid: 0,
                start: 40,
                end: 90,
                gc_tag: Default::default(),
            },
            0.5,
        )?;
        cp.finalize_coverage(true);
        // Blacklist [20,25) and [80,100)
        cp.set_blacklist_mask(&vec![(20, 25), (80, 100)])?;
        cp.build_indexes(true)?;

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
    fn finalize_twice_is_stable() -> Result<()> {
        let mut cp = Coverage::new(30);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        // Assumes delta is NOT dropped!
        let c1 = cp.finalize_coverage(false).to_vec();
        let c2 = cp.finalize_coverage(false).to_vec();
        assert_eq!(c1, c2);
        Ok(())
    }

    #[test]
    fn bulk_parallel_vs_serial_equivalence_with_mask() -> Result<()> {
        let mut cp = Coverage::new(1000);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 0,
            end: 500,
            gc_tag: Default::default(),
        })?;
        cp.add_fragment_weighted(
            Fragment {
                tid: 0,
                start: 250,
                end: 750,
                gc_tag: Default::default(),
            },
            0.5,
        )?;
        cp.finalize_coverage(true);
        cp.set_blacklist_mask(&vec![(400, 450), (700, 900)])?;

        cp.build_indexes(true)?;

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
    fn finalize_after_drop_deltas_panics() {
        let mut cp = Coverage::new(10);
        cp.drop_deltas();
        // finalize_coverage currently assumes the prefix exists and will panic
        let panicked = catch_unwind(AssertUnwindSafe(|| {
            let _ = cp.finalize_coverage(false);
        }))
        .is_err();
        assert!(panicked);
    }

    #[test]
    fn queries_cover_edges_exactly() -> Result<()> {
        let mut cp = Coverage::new(10);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 0,
            end: 10,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;
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
        let mut cp = Coverage::new(200);
        // Coverage 1.0 on [20,60) and 0.5 on [100,150)
        cp.add_fragment(Fragment {
            tid: 0,
            start: 20,
            end: 60,
            gc_tag: Default::default(),
        })?;
        cp.add_fragment_weighted(
            Fragment {
                tid: 0,
                start: 100,
                end: 150,
                gc_tag: Default::default(),
            },
            0.5,
        )?;
        let cov = cp.finalize_coverage(true).to_vec();
        cp.build_indexes(true)?;

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
        let mut cp = Coverage::new(30);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);

        // No blacklist yields no NaNs
        let v = cp.coverage_at_positions_nan(&[9, 10, 19, 20])?;
        for x in &v {
            assert!(!x.is_nan());
        }

        // With blacklist, NaNs appear inside masked region
        cp.set_blacklist_mask(&vec![(12, 15)])?;

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
        let mut cp = Coverage::new(100);
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(false); // Cannot refinalize if delta is dropped
        cp.build_indexes(false)?;
        let s1 = cp.sum_coverage(0, 100, false)?;
        assert!(deq(s1, 10.0, 1e-12));

        // Add more fragments, which should invalidate coverage and indexes
        cp.add_fragment(Fragment {
            tid: 0,
            start: 20,
            end: 30,
            gc_tag: Default::default(),
        })?;
        // Now querying should fail because coverage not finalized
        let err = cp.avg_coverage(0, 100, false).unwrap_err();
        assert!(format!("{err}").contains("coverage not finalized"));

        // Finalize again and rebuild
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;
        let s2 = cp.sum_coverage(0, 100, false)?;
        assert!(deq(s2, 20.0, 1e-12));
        Ok(())
    }

    // Segmented fragments (handles deletions and gaps)

    #[test]
    fn coverage_no_gaps_exclude_inter_mate_gap() -> Result<()> {
        // Two non-overlapping mates; exclude inter-mate gap
        let fwd = sri(0, 10, 20, false, false, 0, &[]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, false).unwrap();

        let mut cp = new_cp(100);
        cp.add_fragment_with_segments(fws, 1.0)?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;

        // Only read spans counted: 10..20 (10 bp) and 40..50 (10 bp) => 20
        let s = cp.sum_coverage(10, 50, false)?;
        assert!(deq(s, 20.0, 1e-9));

        // Gap 20..40 should be zero
        let gap = cp.sum_coverage(20, 40, false)?;
        assert!(deq(gap, 0.0, 1e-9));
        Ok(())
    }

    #[test]
    fn coverage_no_gaps_include_inter_mate_gap() -> Result<()> {
        // Include the inter-mate gap -> full fragment 10..50 (40 bp)
        let fwd = sri(0, 10, 20, false, false, 0, &[]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, true).unwrap();

        let mut cp = new_cp(100);
        cp.add_fragment_with_segments(fws, 1.0)?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;

        let s = cp.sum_coverage(10, 50, false)?;
        assert!(deq(s, 40.0, 1e-9));
        Ok(())
    }

    #[test]
    fn coverage_with_ref_gap_include_inter_mate_gap() -> Result<()> {
        // forward with internal deletion: [10..20], [25..30]
        // reverse [40..50]; include inter-mate gap -> becomes [10..20], [25..50]
        let fwd = sri(0, 10, 30, false, true, 5, &[(0, 10), (15, 5)]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, true).unwrap();

        let mut cp = new_cp(100);
        cp.add_fragment_with_segments(fws, 1.0)?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;

        // 10..20 (10) + 25..50 (25) = 35
        let s = cp.sum_coverage(10, 50, false)?;
        assert!(deq(s, 35.0, 1e-9));

        // The deletion hole 20..25 is zero
        let hole = cp.sum_coverage(20, 25, false)?;
        assert!(deq(hole, 0.0, 1e-9));
        Ok(())
    }

    #[test]
    fn coverage_with_ref_gap_exclude_inter_mate_gap() -> Result<()> {
        let fwd = sri(0, 10, 30, false, true, 5, &[(0, 10), (15, 5)]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, false).unwrap();

        let mut cp = new_cp(100);
        cp.add_fragment_with_segments(fws, 1.0)?;
        cp.finalize_coverage(true);
        cp.build_indexes(true)?;

        // 10..20 (10) + 25..30 (5) + 40..50 (10) = 25
        let s = cp.sum_coverage(10, 50, false)?;
        assert!(deq(s, 25.0, 1e-9));
        Ok(())
    }

    #[test]
    fn has_blacklist_false_when_no_blacklist() -> Result<()> {
        let mut cp = Coverage::new(100);

        // Add a simple fragment so we can finalize coverage
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);

        // No blacklist configured at all
        assert!(!cp.has_blacklist());
        Ok(())
    }
}

#[cfg(test)]
mod tests_window_results {
    use anyhow::{Result, anyhow};
    use cfdnalab::{
        commands::fcoverage::window_results::{
            CoverageOutput, CoverageWindowAction, WindowValue, compute_window_outputs,
        },
        shared::{coverage::Coverage, fragment::minimal_fragment::Fragment},
    };

    fn deq_f32(a: f32, b: f32, tol: f32) -> bool {
        (a - b).abs() <= tol
    }
    fn deq_f64(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn make_cp_with_simple_fragments(len: u32) -> Result<Coverage> {
        let mut cp = Coverage::new(len);
        // Two 10-bp blocks: [10,20) and [30,40)
        cp.add_fragment(Fragment {
            tid: 0,
            start: 10,
            end: 20,
            gc_tag: Default::default(),
        })?;
        cp.add_fragment(Fragment {
            tid: 0,
            start: 30,
            end: 40,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(false);
        Ok(cp)
    }

    #[test]
    fn compute_windows_average() -> Result<()> {
        let mut cp = make_cp_with_simple_fragments(100)?;

        // Average across [10,20) and [30,40) should both be 1.0
        let windows = vec![(10_u64, 20_u64, 0_u64), (30, 40, 1)];
        let out = compute_window_outputs(
            &mut cp,
            Some(&windows),
            CoverageWindowAction::Average,
            false,
        )?;

        match out {
            CoverageOutput::PerWindow { action, results } => {
                assert_eq!(action, CoverageWindowAction::Average);
                assert_eq!(results.len(), 2);
                assert_eq!(results[0].start, 10);
                assert_eq!(results[0].end, 20);
                assert_eq!(results[0].original_idx, 0);
                match results[0].value {
                    WindowValue::Average(v) => assert!(deq_f32(v, 1.0, 1e-6)),
                    _ => return Err(anyhow!("unexpected payload for Average")),
                }
                match results[1].value {
                    WindowValue::Average(v) => assert!(deq_f32(v, 1.0, 1e-6)),
                    _ => return Err(anyhow!("unexpected payload for Average")),
                }
            }
            _ => return Err(anyhow!("expected PerWindow output")),
        }

        Ok(())
    }

    #[test]
    fn compute_windows_total() -> Result<()> {
        let mut cp = make_cp_with_simple_fragments(100)?;

        // Totals should be window length since coverage is 1.0 in-block
        let windows = vec![(10_u64, 20_u64, 0_u64), (30, 40, 1)];
        let out =
            compute_window_outputs(&mut cp, Some(&windows), CoverageWindowAction::Total, false)?;

        match out {
            CoverageOutput::PerWindow { action, results } => {
                assert_eq!(action, CoverageWindowAction::Total);
                assert_eq!(results.len(), 2);
                match results[0].value {
                    WindowValue::Total(v) => assert!(deq_f64(v, 10.0, 1e-9)),
                    _ => return Err(anyhow!("unexpected payload for Total")),
                }
                match results[1].value {
                    WindowValue::Total(v) => assert!(deq_f64(v, 10.0, 1e-9)),
                    _ => return Err(anyhow!("unexpected payload for Total")),
                }
            }
            _ => return Err(anyhow!("expected PerWindow output")),
        }

        Ok(())
    }

    #[test]
    fn compute_windows_positions_with_nan_blacklist() -> Result<()> {
        let mut cp = Coverage::new(60);
        // One fragment spanning [5, 15)
        cp.add_fragment(Fragment {
            tid: 0,
            start: 5,
            end: 15,
            gc_tag: Default::default(),
        })?;
        cp.finalize_coverage(true);

        // Blacklist [9, 12) so indices 9,10,11 are masked
        cp.set_blacklist_mask(&vec![(9, 12)])?;

        // Window [8, 13) -> positions 8,9,10,11,12
        let windows = vec![(8_u64, 13_u64, 0_u64)];
        let out = compute_window_outputs(
            &mut cp,
            Some(&windows),
            CoverageWindowAction::OnlyIncludeThesePositionsIndexed,
            true,
        )?;

        match out {
            CoverageOutput::PerWindow { action, results } => {
                assert_eq!(
                    action,
                    CoverageWindowAction::OnlyIncludeThesePositionsIndexed
                );
                assert_eq!(results.len(), 1);
                let vals = match &results[0].value {
                    WindowValue::Positions(v) => v.clone(),
                    _ => return Err(anyhow!("unexpected payload for Positions")),
                };
                assert_eq!(vals.len(), 5);

                // Expected: [1.0, NaN, NaN, NaN, 1.0]
                assert!(deq_f32(vals[0], 1.0, 1e-6));
                assert!(vals[1].is_nan());
                assert!(vals[2].is_nan());
                assert!(vals[3].is_nan());
                assert!(deq_f32(vals[4], 1.0, 1e-6));
            }
            _ => return Err(anyhow!("expected PerWindow output")),
        }

        Ok(())
    }

    #[test]
    fn compute_windows_whole_positional_when_none() -> Result<()> {
        let mut cp = make_cp_with_simple_fragments(50)?;
        let out = compute_window_outputs(
            &mut cp,
            None,
            CoverageWindowAction::OnlyIncludeThesePositionsIndexed,
            false,
        )?;

        match out {
            CoverageOutput::WholePositional { start, end, values } => {
                assert_eq!(start, 0);
                assert_eq!(end, 50);
                assert_eq!(values.len(), 50);
                // Spot check
                assert!(deq_f32(values[9], 0.0, 1e-6)); // Just before first fragment
                assert!(deq_f32(values[10], 1.0, 1e-6)); // Start of first fragment
                assert!(deq_f32(values[20], 0.0, 1e-6)); // Just after first fragment
                assert!(deq_f32(values[30], 1.0, 1e-6)); // Start of second fragment
            }
            _ => return Err(anyhow!("expected WholePositional output")),
        }

        Ok(())
    }

    #[test]
    fn compute_windows_errors_if_coverage_not_finalized() -> Result<()> {
        let mut cp = Coverage::new(30);
        // Do not finalize_coverage here
        let windows = vec![(0_u64, 10_u64, 0_u64)];
        let res =
            compute_window_outputs(&mut cp, Some(&windows), CoverageWindowAction::Total, false);
        assert!(res.is_err());
        Ok(())
    }

    #[test]
    fn compute_windows_errors_on_out_of_bounds() -> Result<()> {
        let mut cp = make_cp_with_simple_fragments(50)?;
        let windows = vec![(0_u64, 51_u64, 0_u64)]; // end > length
        let res =
            compute_window_outputs(&mut cp, Some(&windows), CoverageWindowAction::Total, false);
        assert!(res.is_err());
        Ok(())
    }
}
