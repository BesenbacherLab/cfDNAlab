#[cfg(test)]
mod tests {
    use anyhow::{Context, Result};
    use cfdnalab::utils::{
        coverage::scale_genome::{
            compute_window_scaling_over_fragment, compute_window_scaling_over_overlap,
        },
        fragment::indel_counting_fragment::{IndelReadInfo, collect_fragment_with_indel_counts},
        overlaps::find_overlapping_windows,
        profiling::midpoint::midpoint_random_even_with_thread_rng,
    };

    // ------- helpers -------

    fn approx_eq(a: f64, b: f64, eps: f64) {
        assert!(
            (a - b).abs() <= eps,
            "expected ~{b}, got {a} (|Δ|={})",
            (a - b).abs()
        );
    }

    fn bed_windows(wins: &[(u64, u64)]) -> Vec<(u64, u64, u64)> {
        // Ignore original_idx; find_overlapping_windows will use scan index.
        wins.iter().map(|&(s, e)| (s, e, 0u64)).collect()
    }

    fn scaling_indices_for_fragment(
        scaling_chr: &[(u64, u64, f32)],
        frag_start: u64,
        frag_end: u64,
    ) -> Result<Vec<usize>> {
        let chrom_len = 1_000u64;
        let mut sf_ptr = 0usize;

        // Build a BED-like view of scaling bins (we only need starts/ends for overlap finding).
        let scaling_bed: Vec<(u64, u64, u64)> =
            scaling_chr.iter().map(|&(s, e, _)| (s, e, 0u64)).collect();

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
        approx_eq(w[0].1, 1.25, 1e-9); // full-window overlap -> avg scaling is the bin's weight
        approx_eq(w[0].2, 1.0, 1e-12); // overlap_fraction relative to fragment length
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
        approx_eq(w[0].1, 1.5, 1e-9); // average scaling
        approx_eq(w[0].2, 1.0, 1e-12); // full overlap of fragment with window
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

        approx_eq(*by_idx.get(&0).unwrap(), (2.0 / 6.0) * 2.0, 1e-9); // [3,5) in first window
        approx_eq(*by_idx.get(&1).unwrap(), (4.0 / 6.0) * 1.0, 1e-9); // [5,9) in second window
        approx_eq(by_idx.values().copied().sum::<f64>(), 1.3333333333, 1e-9);
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
        approx_eq(w[0].1, 8.0 / 6.0, 1e-9); // avg scaling over the overlapped span
        // Overlap_fraction is 1.0 here because the single window completely contains [2,8)
        approx_eq(w[0].2, 1.0, 1e-12);
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
        approx_eq(overlaps.windows[0].overlap_fraction as f64, 1.0, 1e-12);
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

        approx_eq(*by_idx.get(&0).unwrap(), 2.0 / 6.0, 1e-12); // [3,5) in first window
        approx_eq(*by_idx.get(&1).unwrap(), 4.0 / 6.0, 1e-12); // [5,9) in second window
        approx_eq(by_idx.values().copied().sum::<f64>(), 1.0, 1e-12);
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
            approx_eq(avg, 13.0 / 7.0, 1e-9);
            approx_eq(overlap_fraction, 1.0, 1e-12);
        }
        Ok(())
    }

    #[test]
    fn find_overlapping_windows_respects_min_overlap_threshold() -> Result<()> {
        // One window [0,100). Fragment [90,100) has 10% overlap; [1,2) has 1% overlap.
        let chrom_len = 1_000;
        let wins = bed_windows(&[(0, 100)]);
        let mut wd_ptr;

        // 1) Threshold 0.0 -> both fragments should be considered overlapping
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
                1,
                2,
                0.0,
                500,
            )?
            .context("should find tiny overlap [1,2)")?;
            assert_eq!(f2.windows.len(), 1);
        }

        // 2) Threshold 0.05 -> [1,2) should be rejected, [90,100) passes
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
                1,
                2,
                0.05,
                500,
            )?;
            assert!(f2.is_none(), "tiny overlap should be rejected at 5%");
        }

        // 3) Threshold 1.0 -> only near-full overlaps would pass; these do not
        wd_ptr = 0;
        {
            let f1 = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                0,
                100,
                1.0,
                500,
            )?;
            assert!(f1.is_some(), "exact full overlap should pass 100%");

            let f2 = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(&wins),
                None,
                90,
                100,
                1.0,
                500,
            )?;
            assert!(f2.is_none(), "partial overlap should fail 100%");
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
                m == start + len / 2 || m == start + len / 2 + 1,
                "midpoint {m} not one of the two centers"
            );
        }
    }

    #[test]
    fn indel_fragment_builder_fast_paths_and_counts() {
        // Forward and reverse on same tid, inward. Build synthetic per-read indels.
        let fwd = IndelReadInfo {
            tid: 0,
            pos: 100,
            end: 150,
            is_reverse: false,
            // One deletion in non-overlap, one within overlap
            deletions: vec![(90, 95), (120, 122)],
            // One insertion in non-overlap and one in overlap at ref pos 125
            insertions: vec![(80, 2), (125, 3)],
        };
        let rev = IndelReadInfo {
            tid: 0,
            pos: 130,
            end: 200,
            is_reverse: true,
            // One deletion in overlap overlapping [120,122] by 1 bp: [121,124)
            deletions: vec![(121, 124)],
            // Insertion at same overlap ref pos 125 but length 1 (min rules)
            insertions: vec![(125, 1)],
        };

        // 1) skip_indels = true -> None
        let r = collect_fragment_with_indel_counts(&fwd, &rev, true, true);
        assert!(r.is_none());

        // 2) count_indels = false -> Zero adjustments returned
        let r = collect_fragment_with_indel_counts(&fwd, &rev, false, false).unwrap();
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert_eq!(r.deletions_nonoverlap, 0);
        assert_eq!(r.insertions_nonoverlap, 0);
        assert_eq!(r.deletions_overlap_supported, 0);
        assert_eq!(r.insertions_overlap_supported, 0);
        assert_eq!(r.len_indel_adjusted(), r.len_ref());

        // 3) Full counting -> Check non-overlap and overlap-supported math
        let r = collect_fragment_with_indel_counts(&fwd, &rev, false, true).unwrap();
        // Fragment span is [100,200)

        // Non-overlap deletions:
        //   Forward (90,95) is entirely outside (non-overlap with mate's aligned segment), count 5 bp
        //   In overlap logic, we only split relative to the per-read aligned-overlap:
        //     Overlap on aligned segments: [max(100,130)=130, min(150,200)=150) = [130,150)
        //   Forward (120,122) has [120,130) as non-overlap => +10 bp (and [130,122) clipped to none)
        // Total deletions_nonoverlap = 5 + 10 = 15
        assert_eq!(r.deletions_nonoverlap, 15);

        // Overlap-supported deletions:
        //   Forward (120,122) contributes (clipped to overlap) [130,150) ∩ [120,122) = [130,122) -> none
        //   Reverse has (121,124) in overlap; intersection with forward's overlap part is none
        // -> 0
        assert_eq!(r.deletions_overlap_supported, 0);

        // Non-overlap insertions:
        //   Forward (80,2) outside aligned-overlap -> count 2
        //   Forward (125,3) is in overlap; do not add here (handled in overlap-supported)
        //   Reverse has only (125,1) in overlap; nothing in non-overlap
        assert_eq!(r.insertions_nonoverlap, 2);

        // Overlap-supported insertions:
        //   Both at ref pos 125 -> add min(3,1) = 1
        assert_eq!(r.insertions_overlap_supported, 1);

        // Length adjusted = ref_len + inserts_total - dels_total = 100 + (2+1) - (15+0) = 88
        assert_eq!(r.len_ref(), 100);
        assert_eq!(r.len_indel_adjusted(), 88);
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
            approx_eq(avg, 2.0, 1e-12);
            approx_eq(overlap_fraction, 1.0, 1e-12);
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
}
