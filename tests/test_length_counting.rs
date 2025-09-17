#[cfg(test)]
mod tests {
    use anyhow::{Context, Result};
    use cfdnalab::utils::{
        coverage::scale_genome::compute_scaled_window_weights, overlaps::find_overlapping_windows,
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
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr: Vec<(u64, u64, f32)> = vec![(0, 100, 1.25)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let w = compute_scaled_window_weights(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 1);
        // Full-window overlap => weight equals average scaling = 1.25
        approx_eq(w[0].1, 1.25, 1e-9);
        Ok(())
    }

    #[test]
    fn full_window_two_bins_averages_by_length() -> Result<()> {
        // Window [0,10), scaling bins: [0,5)->2.0, [5,10)->1.0, fragment [0,10)
        // Expected: (5*2 + 5*1) / 10 = 1.5
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
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr = vec![(0, 5, 2.0), (5, 10, 1.0)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let w = compute_scaled_window_weights(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 1);
        approx_eq(w[0].1, 1.5, 1e-9);
        Ok(())
    }

    #[test]
    fn two_count_windows_partition_additivity() -> Result<()> {
        // Count windows partition [0,10) into [0,5), [5,10)
        // Scaling bins: [0,5)->2.0, [5,10)->1.0
        // Fragment [3,9) crosses both windows.
        // Expected weights: w1 = (2/6)*2.0 = 0.666..., w2 = (4/6)*1.0 = 0.666..., sum=1.333...
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
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr = vec![(0, 5, 2.0), (5, 10, 1.0)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let w = compute_scaled_window_weights(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 2);

        // Map by window idx for stable assertions
        let mut by_idx = std::collections::BTreeMap::new();
        for (idx, weight) in w {
            by_idx.insert(idx, weight);
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
            1000,
        )?
        .context("count overlaps")?;

        let scaling_chr = vec![(0, 3, 1.0), (3, 6, 2.0), (6, 9, 0.5)];
        let sf_idx = scaling_indices_for_fragment(&scaling_chr, frag_start, frag_end)?;

        let w = compute_scaled_window_weights(&overlaps, &sf_idx, &scaling_chr)?;
        assert_eq!(w.len(), 1);
        approx_eq(w[0].1, 8.0 / 6.0, 1e-9);
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
            1000,
        )?
        .expect("count overlaps must be Some");

        let scaling_chr = vec![(0, 10, 1.0)];
        let sf_idx: Vec<usize> = Vec::new();

        let err = compute_scaled_window_weights(&overlaps, &sf_idx, &scaling_chr).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("scaling_bin_indices is empty"),
            "unexpected error: {msg}"
        );
        Ok(())
    }

    #[test]
    fn no_scaling_single_window_uses_overlap_fraction_full_overlap() -> Result<()> {
        // One count window fully covering the fragment -> weight == 1.0
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
            1000,
        )?
        .context("count overlaps")?;

        // No-scaling branch uses the overlap_fraction directly
        assert_eq!(overlaps.windows.len(), 1);
        approx_eq(overlaps.windows[0].overlap_fraction as f64, 1.0, 1e-12);
        Ok(())
    }

    #[test]
    fn no_scaling_two_windows_uses_overlap_fraction_partition() -> Result<()> {
        // Count windows partition [0,10) into [0,5), [5,10)
        // Fragment [3,9) -> overlaps are lengths 2 and 4 over total len 6
        // Expected weights: 2/6 and 4/6, summing to 1.0
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
            1000,
        )?
        .context("count overlaps")?;

        assert_eq!(overlaps.windows.len(), 2);

        // No-scaling branch increments by overlap_fraction
        let mut by_idx = std::collections::BTreeMap::new();
        for ow in overlaps.windows {
            by_idx.insert(ow.idx, ow.overlap_fraction as f64);
        }

        approx_eq(*by_idx.get(&0).unwrap(), 2.0 / 6.0, 1e-12); // [3,5) in first window
        approx_eq(*by_idx.get(&1).unwrap(), 4.0 / 6.0, 1e-12); // [5,9) in second window
        approx_eq(by_idx.values().copied().sum::<f64>(), 1.0, 1e-12);
        Ok(())
    }
}
