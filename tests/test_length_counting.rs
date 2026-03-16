#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use anyhow::{Context, Result};
    use cfdnalab::shared::{
        fragment::indel_counting_fragment::{IndelReadInfo, collect_fragment_with_indel_counts},
        interval::IndexedInterval,
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        scale_genome::{compute_window_scaling_over_fragment, compute_window_scaling_over_overlap},
    };

    // ------- helpers -------

    fn approx_eq(a: f64, b: f64, eps: f64) {
        assert!(
            (a - b).abs() <= eps,
            "expected ~{b}, got {a} (|Δ|={})",
            (a - b).abs()
        );
    }

    fn bed_windows(windows: &[(u64, u64)]) -> Vec<IndexedInterval<u64>> {
        // Ignore original_idx here. These tests only need valid half-open windows.
        windows
            .iter()
            .enumerate()
            .map(|(window_idx, &(start, end))| {
                IndexedInterval::new(start, end, window_idx as u64)
                    .expect("test windows should be valid non-empty intervals")
            })
            .collect()
    }

    fn scaling_indices_for_fragment(
        scaling_chr: &[(u64, u64, f32)],
        frag_start: u64,
        frag_end: u64,
    ) -> Result<Vec<usize>> {
        let chrom_len = 1_000u64;
        let mut sf_ptr = 0usize;

        // Build a BED-like view of scaling bins (we only need starts/ends for overlap finding).
        let scaling_bed: Vec<IndexedInterval<u64>> = scaling_chr
            .iter()
            .enumerate()
            .map(|(window_idx, &(start, end, _))| {
                IndexedInterval::new(start, end, window_idx as u64)
                    .expect("scaling bins in tests should be valid non-empty intervals")
            })
            .collect();

        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut sf_ptr,
            Some(&scaling_bed),
            None, // by_size
            frag_start,
            frag_end,
            0.0,   // accept any overlap
            1_000, // look_back (large enough)
        )?
        .context("expected >=1 overlapping scaling bin")?;

        Ok(overlaps.windows.iter().map(|w| w.idx as usize).collect())
    }

    // ------- tests -------

    #[test]
    fn single_window_single_bin_full_overlap_yields_bin_weight() -> Result<()> {
        // Window covers the region, single scaling bin with weight=1.25 (already inverted).
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 100)]);
        let mut wd_ptr = 0usize;

        let frag_start = 20;
        let frag_end = 80; // fragment_len=60
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr: Vec<(u64, u64, f32)> = vec![(0, 100, 1.25)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        // (idx, avg_scaling_over_overlap, overlap_fraction)
        let w = compute_window_scaling_over_overlap(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 1);
        approx_eq(w[0].1, 1.25, 1e-6); // full-window overlap -> avg scaling is the bin's weight
        approx_eq(w[0].2, 1.0, 1e-6); // overlap_fraction relative to fragment length
        Ok(())
    }

    #[test]
    fn full_window_two_bins_averages_by_length() -> Result<()> {
        // Window [0,10), scaling bins: [0,5)->2.0, [5,10)->1.0, fragment [0,10)
        // Expected average over overlap: (5*2 + 5*1) / 10 = 1.5
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 10)]);
        let mut wd_ptr = 0usize;

        let frag_start = 0;
        let frag_end = 10;
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr = vec![(0, 5, 2.0), (5, 10, 1.0)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let w = compute_window_scaling_over_overlap(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 1);
        approx_eq(w[0].1, 1.5, 1e-6); // average scaling
        approx_eq(w[0].2, 1.0, 1e-6); // full overlap of fragment with window
        Ok(())
    }

    #[test]
    fn two_count_windows_partition_additivity() -> Result<()> {
        // Count windows partition [0,10) into [0,5), [5,10)
        // Scaling bins: [0,5)->2.0, [5,10)->1.0
        // Fragment [3,9) crosses both windows (fragment len = 6).
        //
        // For each overlapped window we now get:
        //   - avg_scaling_over_overlap (over that window∩fragment span)
        //   - overlap_fraction = overlap_len / fragment_len
        //
        // Combined weight used by the caller for CountOverlap is:
        //   combined = avg_scaling_over_overlap * overlap_fraction
        //
        // Expected: [3,5) => (2/6)*2.0 = 0.666..., [5,9) => (4/6)*1.0 = 0.666..., sum=1.333...
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 5), (5, 10)]);
        let mut wd_ptr = 0usize;

        let frag_start = 3;
        let frag_end = 9; // len=6
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr = vec![(0, 5, 2.0), (5, 10, 1.0)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        // (idx, avg_scaling_over_overlap, overlap_fraction)
        let w = compute_window_scaling_over_overlap(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 2);

        // Map by window idx for stable assertions, and multiply avg_scaling by overlap_fraction
        let mut by_idx = std::collections::BTreeMap::new();
        for (idx, avg_scaling, overlap_fraction) in w {
            by_idx.insert(idx, avg_scaling * overlap_fraction);
        }

        approx_eq(*by_idx.get(&0).unwrap(), (2.0 / 6.0) * 2.0, 1e-6); // [3,5) in first window
        approx_eq(*by_idx.get(&1).unwrap(), (4.0 / 6.0) * 1.0, 1e-6); // [5,9) in second window
        approx_eq(by_idx.values().copied().sum::<f64>(), 1.3333333333, 1e-6);
        Ok(())
    }

    #[test]
    fn multi_bin_average_over_partial_window() -> Result<()> {
        // Scaling bins: [0,3)->1.0, [3,6)->2.0, [6,9)->0.5
        // Window [0,9), Fragment [2,8) => segments: [2,3):1*1, [3,6):3*2, [6,8):2*0.5
        // Weighted sum = 1 + 6 + 1 = 8; avg = 8/6 = 1.333...
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 9)]);
        let mut wd_ptr = 0usize;

        let frag_start = 2;
        let frag_end = 8; // len=6
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr = vec![(0, 3, 1.0), (3, 6, 2.0), (6, 9, 0.5)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let w = compute_window_scaling_over_overlap(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 1);
        approx_eq(w[0].1, 8.0 / 6.0, 1e-6); // avg scaling over the overlapped span
        // Overlap_fraction is 1.0 here because the single window completely contains [2,8)
        approx_eq(w[0].2, 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn error_on_empty_scaling_indices() -> Result<()> {
        // Build a minimal valid count overlap, but pass empty sf indices
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 10)]);
        let mut wd_ptr = 0usize;

        let frag_start = 2;
        let frag_end = 8;
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .expect("count overlaps must be Some");

        let scaling_chr = vec![(0, 10, 1.0)];
        let sf_idx: Vec<usize> = Vec::new();

        let err =
            compute_window_scaling_over_overlap(&overlaps, &sf_idx, &scaling_chr).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("scaling_bin_indices is empty"),
            "unexpected error: {msg}"
        );
        Ok(())
    }

    #[test]
    fn no_scaling_single_window_uses_overlap_fraction_full_overlap() -> Result<()> {
        // One count window fully covering the fragment -> overlap_fraction == 1.0
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 100)]);
        let mut wd_ptr = 0usize;

        let frag_start = 20;
        let frag_end = 80; // len = 60
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        assert_eq!(overlaps.windows.len(), 1);
        approx_eq(overlaps.windows[0].overlap_fraction as f64, 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn no_scaling_two_windows_uses_overlap_fraction_partition() -> Result<()> {
        // Count windows partition [0,10) into [0,5), [5,10)
        // Fragment [3,9) -> overlaps are lengths 2 and 4 over total len 6
        // Expected fractions: 2/6 and 4/6, summing to 1.0
        let chrom_len = 100;
        let count_wins = bed_windows(&[(0, 5), (5, 10)]);
        let mut wd_ptr = 0usize;

        let frag_start = 3;
        let frag_end = 9; // len = 6
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        assert_eq!(overlaps.windows.len(), 2);

        let mut by_idx = std::collections::BTreeMap::new();
        for ow in overlaps.windows {
            by_idx.insert(ow.idx, ow.overlap_fraction as f64);
        }

        approx_eq(*by_idx.get(&0).unwrap(), 2.0 / 6.0, 1e-6); // [3,5) in first window
        approx_eq(*by_idx.get(&1).unwrap(), 4.0 / 6.0, 1e-6); // [5,9) in second window
        approx_eq(by_idx.values().copied().sum::<f64>(), 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn compute_window_scaling_over_fragment_returns_same_avg_for_all_overlaps() -> Result<()> {
        // Two count windows, one fragment crossing both.
        // Average scaling is computed over the FULL fragment, so both rows share the same avg.
        let count_wins = bed_windows(&[(0, 5), (5, 10)]);
        let chrom_len = 100;
        let mut wd_ptr = 0usize;

        let frag_start = 2;
        let frag_end = 9; // len = 7
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        // Three bins over [0,10): [0,4)->2.0, [4,7)->1.0, [7,10)->3.0.
        // Avg over fragment [2,9): ([2,4)=2*2 + [4,7)=3*1 + [7,9)=2*3) / 7 = (4 + 3 + 6)/7 = 13/7
        let scaling_chr = vec![(0, 4, 2.0), (4, 7, 1.0), (7, 10, 3.0)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let rows = compute_window_scaling_over_fragment(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(rows.len(), 2);
        for (_idx, avg, overlap_fraction) in rows {
            approx_eq(avg, 13.0 / 7.0, 1e-6);
            approx_eq(overlap_fraction, 1.0, 1e-6);
        }
        Ok(())
    }

    #[test]
    fn find_overlapping_windows_respects_min_overlap_threshold() -> Result<()> {
        // One window [0,100).
        // Note: min_overlap_fraction is the fraction of the *interval/fragment* overlapped by the window.
        let chrom_len = 1_000;
        let wins = bed_windows(&[(0, 100)]);
        let mut wd_ptr;

        // 1) Threshold 0.0 -> both intervals should be considered overlapping
        //    - [90,100): interval length 10, overlapped 10 -> 100%
        //    - [99,199): interval length 100, overlapped 1  -> 1%
        wd_ptr = 0;
        {
            let f1 = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                90,
                100,
                0.0,
                500,
            )?
            .context("should find overlap [90,100)")?;
            assert_eq!(f1.windows.len(), 1);

            let f2 = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                99,
                199,
                0.0,
                500,
            )?
            .context("should find tiny (1%) interval overlap [99,199)")?;
            assert_eq!(f2.windows.len(), 1);
        }

        // 2) Threshold 0.05 (5%) -> [99,199) (1%) should be rejected; [90,100) (100%) passes
        wd_ptr = 0;
        {
            let f1 = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                90,
                100,
                0.05,
                500,
            )?
            .context("should find overlap [90,100) at 5%")?;
            assert_eq!(f1.windows.len(), 1);

            let f2 = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                99,
                199,
                0.05,
                500,
            )?;
            assert!(
                f2.is_none(),
                "1% of the interval is <5%, so it should be rejected"
            );
        }

        // 3) Threshold 1.0 (100%) -> only exact full-coverage intervals pass.
        //    - [0,100) is fully covered by [0,100) -> pass
        //    - [50,150) is not fully covered (only 50%) -> reject
        wd_ptr = 0;
        {
            let f_full = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                0,
                100,
                1.0,
                500,
            )?;
            assert!(
                f_full.is_some(),
                "exact full overlap [0,100) should pass 100%"
            );

            let f_full_smaller = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                20,
                80,
                1.0,
                500,
            )?;
            assert!(
                f_full_smaller.is_some(),
                "smaller full overlap [20,80) should pass 100%"
            );

            let f_partial = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                50,
                150,
                1.0,
                500,
            )?;
            assert!(f_partial.is_none(), "partial overlap should fail at 100%");
        }

        Ok(())
    }

    #[test]
    fn midpoint_random_even_with_thread_rng_returns_one_of_two_centers() {
        // Even-length fragment [start, start+len), len=6 -> centers are start+2 or start+3.
        let start = 1000u32;
        let len = 6u32;

        // Sample a few times; each must be one of the two positions.
        for _ in 0..100 {
            let m = midpoint_random_even_with_thread_rng(start, len);
            assert!(
                m == start + len / 2 || m == start + len / 2 - 1,
                "midpoint {m} not one of the two centers"
            );
        }

        // Midpoint from odd-length fragment
        let len = len + 1; // 7
        let m = midpoint_random_even_with_thread_rng(start, len);
        assert!(m == start + len / 2, "midpoint {m} not the center element");
    }

    #[test]
    fn indel_fragment_builder_fast_paths_and_counts() {
        // Forward and reverse on same tid, inward.
        // Overlap will be [120,150) so ref pos 125 lies in the overlap.
        let fwd = IndelReadInfo {
            tid: 0,
            pos: 100,
            end: 150,
            is_reverse: false,
            // One deletion in non-overlap (inside fragment), one within overlap
            deletions: vec![(110, 115), (120, 122)],
            // One insertion in non-overlap (inside fragment) and one in overlap at ref 125
            insertions: vec![(105, 2), (125, 3)],
        };
        let rev = IndelReadInfo {
            tid: 0,
            pos: 120,
            end: 200,
            is_reverse: true,
            // Deletion overlaps [120,122) by 2 bp (121..122), rest are discarded as non-consensus
            deletions: vec![(121, 124)],
            // Insertion at same overlap ref pos 125 but length 1 (min rule will pick 1)
            insertions: vec![(125, 1)],
        };

        let frag = collect_fragment_with_indel_counts(&fwd, &rev, false, true).unwrap();

        // Non-overlap deletions: (110..115)=5 from fwd and none for rev
        assert_eq!(frag.deletions_nonoverlap, 5);

        // Overlap deletions supported: intersection of (fwd 120..122) with (rev 121..124) is (121..122) => 1.
        assert_eq!(frag.deletions_overlap_supported, 1);

        // Non-overlap insertions: (105,+2) is non-overlap (before 120) => 2
        // The (125,*) insertions are in the overlap, so not counted here
        assert_eq!(frag.insertions_nonoverlap, 2);

        // Overlap insertions supported: both at ref 125 => min(3,1)=1
        assert_eq!(frag.insertions_overlap_supported, 1);

        // Length adjusted = ref_len + inserts_total - dels_total = 100 + (5+1) - (2+1) = 97
        assert_eq!(frag.len_ref(), 100);
        assert_eq!(frag.len_indel_adjusted(), 97);

        // Skip and no-counts

        // Skip_indels = true -> None
        let r = collect_fragment_with_indel_counts(&fwd, &rev, true, true);
        assert!(r.is_none());

        // Count_indels = false -> Zero adjustments returned
        let r = collect_fragment_with_indel_counts(&fwd, &rev, false, false).unwrap();
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert_eq!(r.deletions_nonoverlap, 0);
        assert_eq!(r.insertions_nonoverlap, 0);
        assert_eq!(r.deletions_overlap_supported, 0);
        assert_eq!(r.insertions_overlap_supported, 0);
        assert_eq!(r.len_indel_adjusted(), r.len_ref());
    }

    #[test]
    fn scaling_over_fragment_is_constant_per_fragment_even_with_many_windows() -> Result<()> {
        // Three windows tile [0,30); fragment [5,25) overlaps all three.
        let chrom_len = 1000;
        let count_wins = bed_windows(&[(0, 10), (10, 20), (20, 30)]);
        let mut wd_ptr = 0usize;

        let frag_start = 5;
        let frag_end = 25; // len=20
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&count_wins),
            None,
            frag_start,
            frag_end,
            0.0,
            1000,
        )?
        .context("count overlaps")?;

        // Scaling bins across [0,30): [0,10)->1.0, [10,20)->2.0, [20,30)->3.0
        // Average over [5,25) = (5*1 + 10*2 + 5*3)/20 = (5 + 20 + 15)/20 = 40/20 = 2.0
        let scaling_chr = vec![(0, 10, 1.0), (10, 20, 2.0), (20, 30, 3.0)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let rows = compute_window_scaling_over_fragment(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(rows.len(), 3);
        for (_idx, avg, overlap_fraction) in rows {
            approx_eq(avg, 2.0, 1e-6);
            approx_eq(overlap_fraction, 1.0, 1e-6);
        }
        Ok(())
    }

    #[test]
    fn synthetic_no_scaling_midpoint_assignment_counts_full() -> Result<()> {
        // Midpoint assignment ignores overlap fraction and counts as 1.0 per overlapped window,
        // but here we simulate it by building the OverlappingWindows with a 1 bp interval.
        let chrom_len = 100;
        let wins = bed_windows(&[(10, 20), (20, 30)]);
        let mut wd_ptr = 0usize;

        // Fragment midpoint hits exactly 19 -> falls in first window only.
        let midpoint_start = 19;
        let midpoint_end = 20;
        let overlaps = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&wins),
            None,
            midpoint_start,
            midpoint_end,
            0.99, // Require nearly full (for 1 bp, this is fine)
            500,
        )?
        .context("midpoint overlaps")?;

        assert_eq!(overlaps.windows.len(), 1);
        assert_eq!(overlaps.windows[0].idx, 0); // First window
        Ok(())
    }

    #[test]
    fn size_mode_basic_three_bins() -> Result<()> {
        // Bins [0,10), [10,20), [20,30) via size-mode.
        // Fragment [3,27) overlaps 3 bins with lengths 7,10,7 -> fractions over fragment_len=24.
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;
        let by_size = Some(10);

        let frag_start = 3;
        let frag_end = 27; // len = 24

        let res = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            None,    // BED-mode windows not used
            by_size, // fixed-size bins
            frag_start,
            frag_end,
            0.0, // no min fraction threshold
            1_000,
        )?
        .context("expect overlaps in size-mode")?;

        // Expect bins with indices 0,1,2 (i.e., [0,10),[10,20),[20,30))
        assert_eq!(res.windows.len(), 3);

        // Collect overlap fractions in a stable order by idx.
        let mut by_idx = BTreeMap::new();
        for w in res.windows {
            by_idx.insert(w.idx, w.overlap_fraction as f64);
        }

        // Fractions are overlap_len / fragment_len = {7/24, 10/24, 7/24}
        approx_eq(*by_idx.get(&0).unwrap(), 7.0 / 24.0, 1e-6);
        approx_eq(*by_idx.get(&1).unwrap(), 10.0 / 24.0, 1e-6);
        approx_eq(*by_idx.get(&2).unwrap(), 7.0 / 24.0, 1e-6);
        Ok(())
    }

    #[test]
    fn size_mode_min_fraction_filters_bins() -> Result<()> {
        // Same setup as above but require >= 0.3 fraction of the fragment per bin.
        // 7/24 ≈ 0.2917 (filtered), 10/24 ≈ 0.4167 (kept), 7/24 filtered.
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;

        let res = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            None,
            Some(10),
            3,
            27,
            0.30, // threshold
            1_000,
        )?
        .context("size-mode with threshold")?;

        // Only the middle bin remains.
        assert_eq!(res.windows.len(), 1);
        assert_eq!(res.windows[0].idx, 1);
        approx_eq(res.windows[0].overlap_fraction as f64, 10.0 / 24.0, 1e-6);
        Ok(())
    }

    #[test]
    fn bed_mode_edge_touching_and_tiny_overlap() -> Result<()> {
        // Window [20,30). A fragment [10,20) just touches -> no overlap.
        // A fragment [19,20) 1 bp long, [19,20)∩[20,30)=∅ -> no overlap.
        // A fragment [19,21) overlaps 1 bp -> fraction = 1/2.
        let chrom_len = 1_000;
        let wins = bed_windows(&[(20, 30)]);

        // Touching at 20: no overlap
        let mut wd_ptr = 0usize;
        let none1 =
            find_overlapping_windows(chrom_len, &mut wd_ptr, Some(&wins), None, 10, 20, 0.0, 500)?;
        assert!(none1.is_none(), "pure edge touch should not overlap");

        // Still no overlap for [19,20)
        wd_ptr = 0;
        let none2 =
            find_overlapping_windows(chrom_len, &mut wd_ptr, Some(&wins), None, 19, 20, 0.0, 500)?;
        assert!(
            none2.is_none(),
            "1bp ending at 20 should not overlap [20,30)"
        );

        // [19,21): 1 bp overlap with window -> fraction = 1/2 of the fragment
        wd_ptr = 0;
        let some =
            find_overlapping_windows(chrom_len, &mut wd_ptr, Some(&wins), None, 19, 21, 0.0, 500)?
                .context("expected a small overlap")?;
        assert_eq!(some.windows.len(), 1);
        approx_eq(some.windows[0].overlap_fraction as f64, 1.0 / 2.0, 1e-6);
        Ok(())
    }

    #[test]
    fn midpoint_mode_bed() -> Result<()> {
        // “Midpoint” usage simulated by interval [m, m+1) and high min-overlap threshold.
        // Windows: [0,10), [10,20), [20,30). Midpoint m=17 hits only [10,20).
        let chrom_len = 10_000;
        let wins = bed_windows(&[(0, 10), (10, 20), (20, 30)]);
        let mut wd_ptr = 0usize;

        let midpoint = 17u64;
        let res = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&wins),
            None,
            midpoint,
            midpoint + 1,
            0.99, // effectively require 1/1 bp overlap
            1_000,
        )?
        .context("midpoint overlap")?;

        assert_eq!(res.windows.len(), 1);
        assert_eq!(res.windows[0].idx, 1); // [10,20)
        approx_eq(res.windows[0].overlap_fraction as f64, 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn midpoint_mode_size_bins() -> Result<()> {
        // Same as above but with size-mode bins of 10 bp.
        let chrom_len = 10_000;
        let mut wd_ptr = 0usize;

        let m = 17u64;
        let res = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            None,
            Some(10),
            m,
            m + 1,
            0.99,
            1_000,
        )?
        .context("midpoint in size-mode")?;

        assert_eq!(res.windows.len(), 1);
        // Bin index = floor(m/10) = 1
        assert_eq!(res.windows[0].idx, 1);
        approx_eq(res.windows[0].overlap_fraction as f64, 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn streaming_pointer_advances_past_lookback() -> Result<()> {
        // Two windows: [0,100), [100,200). We’ll feed two fragments in ascending order and
        // ensure wd_ptr advances so the second call doesn’t re-check the first window.
        let chrom_len = 10_000;
        let wins = bed_windows(&[(0, 100), (100, 200)]);
        let mut wd_ptr = 0usize;
        let look_back = 0u64; // strict advancement

        // Fragment overlaps only the first window
        let _ = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&wins),
            None,
            10,
            20,
            0.0,
            look_back,
        )?
        .context("first call should overlap [0,100)")?;
        // wd_ptr may still be 0 here (we didn’t push beyond), so drive it with a distant fragment

        // Now a fragment that only hits the second window; with look_back=0 we
        // should skip the first quickly (pointer should not regress)
        let res2 = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&wins),
            None,
            150,
            160,
            0.0,
            look_back,
        )?
        .context("second call should overlap [100,200)")?;
        assert_eq!(res2.windows.len(), 1);
        assert_eq!(res2.windows[0].idx, 1);
        Ok(())
    }

    #[test]
    fn size_mode_one_bp_on_left_boundary_hits_left_bin() -> Result<()> {
        // Bins [0,10) [10,20) ...
        // Fragment [9,10) (len=1) overlaps bin 0 only by 1 bp.
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;

        let res =
            find_overlapping_windows(chrom_len, &mut wd_ptr, None, Some(10), 9, 10, 0.0, 1_000)?
                .context("expect overlap at left boundary")?;

        assert_eq!(res.windows.len(), 1);
        assert_eq!(res.windows[0].idx, 0);
        approx_eq(res.windows[0].overlap_fraction as f64, 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn size_mode_one_bp_on_right_boundary_hits_right_bin() -> Result<()> {
        // Bins [0,10) [10,20) ...
        // Fragment [10,11) (len=1) overlaps bin 1 only by 1 bp.
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;

        let res =
            find_overlapping_windows(chrom_len, &mut wd_ptr, None, Some(10), 10, 11, 0.0, 1_000)?
                .context("expect overlap at right boundary")?;

        assert_eq!(res.windows.len(), 1);
        assert_eq!(res.windows[0].idx, 1);
        approx_eq(res.windows[0].overlap_fraction as f64, 1.0, 1e-6);
        Ok(())
    }

    #[test]
    fn size_mode_one_bp_boundary_with_high_threshold() -> Result<()> {
        // With min_overlap_fraction ~1.0 the 1 bp fragment still passes (1/1).
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;

        let res =
            find_overlapping_windows(chrom_len, &mut wd_ptr, None, Some(10), 10, 11, 0.99, 1_000)?
                .context("expect overlap at boundary with strict threshold")?;

        assert_eq!(res.windows.len(), 1);
        assert_eq!(res.windows[0].idx, 1);
        approx_eq(res.windows[0].overlap_fraction as f64, 1.0, 1e-6);
        Ok(())
    }

    // ---------- Look-back permitting slight backtracking ----------

    #[test]
    fn look_back_allows_small_backtracking_without_losing_windows() -> Result<()> {
        // BED windows [100,200) and [200,300).
        // First fragment [210,220) (hits bin 1).
        // Second fragment [150,160) (goes "backwards" but within look_back=100).
        // wd_ptr should not skip [100,200) based on the look_back logic.
        let chrom_len = 10_000;
        let wins = bed_windows(&[(100, 200), (200, 300)]);
        let look_back = 100u64;
        let mut wd_ptr = 0usize;

        // Forward-ish fragment (later)
        let later = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&wins),
            None,
            210,
            220,
            0.0,
            look_back,
        )?
        .context("later fragment should overlap [200,300)")?;
        assert_eq!(later.windows.len(), 1);
        assert_eq!(later.windows[0].idx, 1);

        // Slight backtrack, but within look_back
        let earlier = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&wins),
            None,
            150,
            160,
            0.0,
            look_back,
        )?
        .context("earlier fragment should still see [100,200)")?;
        assert_eq!(earlier.windows.len(), 1);
        assert_eq!(earlier.windows[0].idx, 0);
        Ok(())
    }

    // ---------- Large fragment spanning many size-bins + threshold ----------

    #[test]
    fn big_fragment_spanning_many_bins_with_threshold_filters_edges() -> Result<()> {
        // Size-mode bins of 10 bp. Fragment [5,95) has len=90 and touches bins 0..9.
        // Overlap lengths: edge bins have 5 bp each; middle bins have 10 bp.
        // With threshold 0.11:
        //   - Middle bins: 10/90 ≈ 0.111... >= 0.11 -> keep
        //   - Edge bins:   5/90  ≈ 0.0555... < 0.11  -> drop
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;

        let res =
            find_overlapping_windows(chrom_len, &mut wd_ptr, None, Some(10), 5, 95, 0.11, 1_000)?
                .context("expect multiple bins after thresholding")?;

        // Expect bins 1..=8 (8 bins).
        assert_eq!(res.windows.len(), 8);

        let mut indices = Vec::new();
        let mut fracs = Vec::new();
        for w in res.windows {
            indices.push(w.idx);
            fracs.push(w.overlap_fraction as f64);
        }

        assert_eq!(indices, vec![1, 2, 3, 4, 5, 6, 7, 8]);
        for f in fracs {
            approx_eq(f, 10.0 / 90.0, 1e-6);
        }
        Ok(())
    }

    // ---------- Mixed sanity: size-mode multiple bins, sum of fractions ----------

    #[test]
    fn size_mode_overlap_fractions_sum_to_one_for_partition() -> Result<()> {
        // Fragment [3,27) len=24 across bins of 10 bp.
        // Overlaps: 7,10,7 -> fractions sum to 1.0.
        let chrom_len = 1_000;
        let mut wd_ptr = 0usize;

        let res =
            find_overlapping_windows(chrom_len, &mut wd_ptr, None, Some(10), 3, 27, 0.0, 1_000)?
                .context("partition fractions")?;

        let sum: f64 = res.windows.iter().map(|w| w.overlap_fraction as f64).sum();

        approx_eq(sum, 1.0, 1e-6);
        Ok(())
    }
}
