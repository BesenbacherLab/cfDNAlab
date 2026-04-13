#[cfg(test)]
mod tests_apply_scaling {
    use cfdnalab::shared::scale_genome::apply_scaling_to_coverage_in_place;

    // Assert two slices are approximately equal within eps
    fn assert_slice_eq_eps(a: &[f32], b: &[f32], eps: f32) {
        assert_eq!(
            a.len(),
            b.len(),
            "Length mismatch: {} vs {}",
            a.len(),
            b.len()
        );
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            let diff = (x - y).abs();
            assert!(
                diff <= eps,
                "Mismatch at {}: got {}, expected {}, |diff|={}",
                i,
                x,
                y,
                diff
            );
        }
    }

    #[test]
    fn scaling_mixed_bins_partial_overlap() {
        // Human verification status: unverified
        // Tile core [105, 125), len 20
        // Bins: [90,110) sf=2.0, [110,120) sf=0.0, [120,140) sf=0.5
        // Expected (multiply-by-sf):
        //   indices 0..5   *= 2.0  -> 20
        //   indices 5..15  *= 0.0  -> 0
        //   indices 15..20 *= 0.5  -> 5
        let core_start = 105u32;
        let bins = vec![(90u64, 110u64, 2.0f32), (110, 120, 0.0), (120, 140, 0.5)];
        let mut cov = vec![10.0f32; 20];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![0.0f32; 20];
        for v in &mut expected[0..5] {
            *v = 20.0;
        }
        // 5..15 stay 0.0
        for v in &mut expected[15..20] {
            *v = 5.0;
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_starts_exactly_at_bin_boundary() {
        // Human verification status: unverified
        // Tile core starts exactly at previous bin end
        // Bins: [0,100) sf=2.0, [100,200) sf=4.0
        // Tile: core_start=100, len=10 so all in second bin
        let core_start = 100u32;
        let bins = vec![(0u64, 100u64, 2.0f32), (100, 200, 4.0)];
        let mut cov = vec![8.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let expected = vec![32.0f32; 10]; // 8 * 4.0
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_zero_bin_covers_entire_tile() {
        // Human verification status: unverified
        // Entire tile lies in a zero-scaled bin
        let core_start = 60u32;
        let bins = vec![(50u64, 150u64, 0.0f32)];
        let mut cov = (0..20).map(|i| i as f32 + 1.0).collect::<Vec<_>>();
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let expected = vec![0.0f32; 20];
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_noop_on_empty_inputs() {
        // Human verification status: unverified
        // No bins or empty coverage should be a no-op
        let core_start = 100u32;

        let mut cov1 = vec![] as Vec<f32>;
        apply_scaling_to_coverage_in_place(&mut cov1, core_start, &[]);
        assert!(cov1.is_empty());

        let mut cov_with_bins = vec![] as Vec<f32>;
        apply_scaling_to_coverage_in_place(
            &mut cov_with_bins,
            core_start,
            &[(90u64, 110u64, 2.0f32)],
        );
        assert!(cov_with_bins.is_empty());

        let mut cov2 = vec![1.0f32, 2.0, 3.0];
        apply_scaling_to_coverage_in_place(&mut cov2, core_start, &[]);
        assert_slice_eq_eps(&cov2, &[1.0, 2.0, 3.0], 1e-6);
    }

    #[test]
    fn scaling_non_overlapping_bins_on_both_sides_leave_varied_tile_unchanged() {
        // Human verification status: unverified
        // Tile [500, 510), with one scaling bin entirely left and one entirely right.
        // This is the real no-op case with non-empty scaling input: the function must not touch
        // any element just because unrelated bins exist elsewhere on the chromosome.
        let core_start = 500u32;
        let bins = vec![(100u64, 200u64, 2.0f32), (520u64, 600u64, 0.25f32)];
        let mut cov = vec![1.0f32, 3.5, 2.0, 7.0, 0.5, 9.0, 4.0, 8.0, 6.0, 5.5];
        let expected = cov.clone();

        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_multiple_bins_exact_edges() {
        // Human verification status: unverified
        // Tile aligns exactly to bin edges
        // Bins: [1000,1010) sf=2, [1010,1020) sf=0, [1020,1030) sf=0.5
        let core_start = 1000u32;
        let bins = vec![
            (1000u64, 1010u64, 2.0f32),
            (1010, 1020, 0.0),
            (1020, 1030, 0.5),
        ];
        let mut cov = vec![6.0f32; 30]; // First 30 cover those bins
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![0.0f32; 30];
        for v in &mut expected[0..10] {
            *v = 12.0; // 6 * 2
        }
        for v in &mut expected[10..20] {
            *v = 0.0; // zero bin
        }
        for v in &mut expected[20..30] {
            *v = 3.0; // 6 * 0.5
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_bins_entirely_left_of_tile_noop() {
        // Human verification status: unverified
        // All bins end before the tile starts -> no effect.
        // Include one bin that touches the tile boundary exactly at `core_start` to prove the
        // implementation treats the bins as half-open and does not scale index 0 by mistake.
        let core_start = 500u32;
        let bins = vec![(100u64, 200u64, 2.0f32), (490, 500, 0.5)];
        let mut cov = vec![7.0f32, 1.0, 9.0, 3.5, 2.25, 4.0, 8.0, 6.5, 5.0, 10.0];
        let expected = cov.clone();
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_bins_entirely_right_of_tile_noop() {
        // Human verification status: unverified
        // All bins start after the tile ends -> no effect.
        // Include one bin that starts exactly at the tile end to prove the half-open convention at
        // the right boundary.
        let core_start = 100u32;
        let bins = vec![(108u64, 200u64, 2.0f32), (1000, 2000, 0.5)];
        let mut cov = vec![4.0f32, 1.5, 7.0, 2.0, 9.0, 3.0, 8.5, 6.0];
        let expected = cov.clone();
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_partial_left_overlap_only() {
        // Human verification status: unverified
        // Bin overlaps only the left edge of the tile
        // Tile: [200, 210), Bin: [195, 205) sf=3.0 -> indices 0..5 scaled by 3
        let core_start = 200u32;
        let bins = vec![(195u64, 205u64, 3.0f32)];
        let mut cov = vec![2.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![2.0f32; 10];
        for v in &mut expected[0..5] {
            *v = 6.0;
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_partial_right_overlap_only() {
        // Human verification status: unverified
        // Bin overlaps only the right edge of the tile
        // Tile: [300, 312), Bin: [308, 400) sf=0.25 -> indices 8..12 scaled by 0.25
        let core_start = 300u32;
        let bins = vec![(308u64, 400u64, 0.25f32)];
        let mut cov = vec![20.0f32; 12];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![20.0f32; 12];
        for v in &mut expected[8..12] {
            *v = 5.0;
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_many_small_bins_inside_tile() {
        // Human verification status: unverified
        // Tile: [1000, 1010). Bins: per-base alternating 2.0 and 0.5.
        let core_start = 1000u32;
        let bins = (1000u64..1010u64)
            .enumerate()
            .map(|(i, s)| {
                let sf = if i % 2 == 0 { 2.0f32 } else { 0.5f32 };
                (s, s + 1, sf)
            })
            .collect::<Vec<_>>();
        let mut cov = vec![10.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let expected: Vec<f32> = (0..10)
            .map(|i| if i % 2 == 0 { 20.0 } else { 5.0 })
            .collect();
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }
}

#[cfg(test)]
mod tests_compute_window_scaling {
    use cfdnalab::shared::{
        interval::{IndexedInterval, Interval},
        overlaps::{OverlappingWindow, OverlappingWindows, find_overlapping_windows},
        scale_genome::{compute_window_scaling_over_fragment, compute_window_scaling_over_overlap},
    };

    fn assert_f64_close(actual: f64, expected: f64, eps: f64, context: &str) {
        let diff = (actual - expected).abs();
        assert!(
            diff <= eps,
            "{context}: got {actual}, expected {expected}, |diff|={diff}"
        );
    }

    fn make_two_window_fragment_overlaps() -> anyhow::Result<OverlappingWindows> {
        // Query/fragment interval is [20,81), length 61.
        //
        // Count windows:
        // - left  window [0,50)   -> overlap [20,50) = 30 bp
        // - right window [50,100) -> overlap [50,81) = 31 bp
        //
        // So the overlap fractions relative to the fragment are:
        // - left  = 30 / 61
        // - right = 31 / 61
        let fragment = Interval::new(20, 81)?;
        let mut overlaps = OverlappingWindows::new(fragment);
        overlaps.windows.push(OverlappingWindow::new(
            0,
            Interval::new(0, 50)?,
            30.0_f32 / 61.0_f32,
        )?);
        overlaps.windows.push(OverlappingWindow::new(
            1,
            Interval::new(50, 100)?,
            31.0_f32 / 61.0_f32,
        )?);
        Ok(overlaps)
    }

    #[test]
    fn compute_window_scaling_over_fragment_uses_explicit_full_fragment_span_for_every_overlapping_window()
    -> anyhow::Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use one fragment/query interval [20,81) and two overlapping count windows:
        // - [0,50)
        // - [50,100)
        //
        // Scaling bins cover the chromosome with:
        // - [0,20):   0
        // - [20,40):  1
        // - [40,60):  1
        // - [60,80):  1
        // - [80,200): 0
        //
        // Full-fragment averaging over [20,81) is therefore:
        //   (20*1 + 20*1 + 20*1 + 1*0) / 61 = 60/61.
        //
        // The helper contract says this same full-fragment average is returned for every count
        // window that overlaps the fragment, and the reported overlap fraction is always 1.0.
        let count_overlaps = make_two_window_fragment_overlaps()?;
        let scaling_chr = vec![
            (0_u64, 20_u64, 0.0_f32),
            (20, 40, 1.0),
            (40, 60, 1.0),
            (60, 80, 1.0),
            (80, 200, 0.0),
        ];
        let scaling_bin_indices = vec![0_usize, 1, 2, 3, 4];

        // Act
        let out = compute_window_scaling_over_fragment(
            Interval::new(20, 81)?,
            &count_overlaps,
            &scaling_bin_indices,
            &scaling_chr,
        )?;

        // Assert
        assert_eq!(out.len(), 2);
        let expected_weight = 60.0_f64 / 61.0_f64;
        assert_eq!(out[0].0, 0);
        assert_f64_close(
            out[0].1,
            expected_weight,
            1e-12,
            "left full-fragment weight",
        );
        assert_f64_close(out[0].2, 1.0, 1e-12, "left full-fragment overlap fraction");
        assert_eq!(out[1].0, 1);
        assert_f64_close(
            out[1].1,
            expected_weight,
            1e-12,
            "right full-fragment weight",
        );
        assert_f64_close(out[1].2, 1.0, 1e-12, "right full-fragment overlap fraction");

        Ok(())
    }

    #[test]
    fn compute_window_scaling_over_overlap_uses_each_window_overlap_span() -> anyhow::Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same fragment/query interval [20,81), count windows, and scaling bins as the
        // previous test.
        //
        // But now the helper should average only over each window/fragment overlap:
        // - left window overlap [20,50):
        //     all 30 bp lie in scaling-factor-1 bins
        //     -> average = 1
        // - right window overlap [50,81):
        //     [50,60): 10 bp at 1
        //     [60,80): 20 bp at 1
        //     [80,81):  1 bp at 0
        //     -> average = (10 + 20 + 0) / 31 = 30/31
        //
        // The overlap fractions are still measured relative to the full fragment:
        // - left  = 30/61
        // - right = 31/61
        let count_overlaps = make_two_window_fragment_overlaps()?;
        let scaling_chr = vec![
            (0_u64, 20_u64, 0.0_f32),
            (20, 40, 1.0),
            (40, 60, 1.0),
            (60, 80, 1.0),
            (80, 200, 0.0),
        ];
        let scaling_bin_indices = vec![0_usize, 1, 2, 3, 4];

        // Act
        let out = compute_window_scaling_over_overlap(
            &count_overlaps,
            &scaling_bin_indices,
            &scaling_chr,
        )?;

        // Assert
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].0, 0);
        assert_f64_close(out[0].1, 1.0, 1e-12, "left overlap-only weight");
        assert_f64_close(
            out[0].2,
            (30.0_f32 / 61.0_f32) as f64,
            1e-7,
            "left overlap fraction",
        );
        assert_eq!(out[1].0, 1);
        assert_f64_close(
            out[1].1,
            30.0_f64 / 31.0_f64,
            1e-12,
            "right overlap-only weight",
        );
        assert_f64_close(
            out[1].2,
            (31.0_f32 / 61.0_f32) as f64,
            1e-7,
            "right overlap fraction",
        );

        Ok(())
    }

    #[test]
    fn scaling_bin_overlap_pipeline_recovers_chrom_local_indices_and_correct_fragment_weight()
    -> anyhow::Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Exercise the intended scaling pipeline directly, independent of any command:
        //
        // 1. Start from one chromosome-local scaling table:
        //      [0,20):   0
        //      [20,40):  1
        //      [40,60):  2
        //      [60,80):  4
        //      [80,200): 8
        // 2. Build an ordered BED-mode window list from that table.
        // 3. Ask the overlap finder which scaling bins touch fragment [20,81).
        // 4. Recover the overlap finder scan indices from the result.
        // 5. Feed those indices into `compute_window_scaling_over_fragment(...)`.
        //
        // In BED mode, `find_overlapping_windows(...)` returns the scan position inside the
        // supplied ordered window list, not `IndexedInterval.idx`. This still gives the correct
        // lookup key here because the window list is built in the same order as `scaling_chr`, so
        // scan positions and chromosome-local `scaling_chr` row indices are identical.
        //
        // The fragment covers:
        //   [20,40): 20 bp at weight 1
        //   [40,60): 20 bp at weight 2
        //   [60,80): 20 bp at weight 4
        //   [80,81):  1 bp at weight 8
        //
        // Therefore the full-fragment average scaling is:
        //   (20*1 + 20*2 + 20*4 + 1*8) / 61 = 148/61.
        //
        // The count windows are [0,50) and [50,100), so both overlap the fragment and must both
        // receive that same full-fragment weight.
        let fragment = Interval::new(20_u64, 81_u64)?;
        let count_overlaps = make_two_window_fragment_overlaps()?;
        let scaling_chr = vec![
            (0_u64, 20_u64, 0.0_f32),
            (20, 40, 1.0),
            (40, 60, 2.0),
            (60, 80, 4.0),
            (80, 200, 8.0),
        ];
        let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
            .iter()
            .enumerate()
            .map(|(idx, (start, end, _))| IndexedInterval::new(*start, *end, idx as u64))
            .collect::<cfdnalab::Result<_>>()?;

        // Act:
        // Recover the overlapping scaling-bin indices through the same BED-mode overlap-finder
        // path used by the commands.
        let mut scaling_ptr = 0_usize;
        let overlapping_scaling_bins = find_overlapping_windows(
            200,
            &mut scaling_ptr,
            Some(&scaling_with_bin_idx),
            None,
            fragment,
            1.0 / 100.0, // Any positive overlap; the exact denominator is not important here.
            100,
        )?
        .expect("fragment should overlap scaling bins");
        let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
            .windows
            .iter()
            .map(|window| window.idx)
            .collect();

        let per_window_scaling = compute_window_scaling_over_fragment(
            fragment,
            &count_overlaps,
            &overlapping_scaling_bin_indices,
            &scaling_chr,
        )?;

        // Assert:
        // The overlap finder must recover the touched scan positions in the ordered scaling-bin
        // list. Because that list is built in the same order as `scaling_chr`, these are also the
        // correct chromosome-local indices for indexing back into `scaling_chr`.
        assert_eq!(
            overlapping_scaling_bin_indices,
            vec![1, 2, 3, 4],
            "fragment [20,81) should touch the 2nd through 5th scaling rows in scan order"
        );

        assert_eq!(
            per_window_scaling.len(),
            2,
            "both count windows should receive the fragment-level scaling weight"
        );
        assert_eq!(
            per_window_scaling[0].0, 0,
            "left count window should retain its original window index"
        );
        assert_eq!(
            per_window_scaling[1].0, 1,
            "right count window should retain its original window index"
        );
        assert_f64_close(
            per_window_scaling[0].1,
            148.0 / 61.0,
            1e-12,
            "left window fragment-average scaling",
        );
        assert_f64_close(
            per_window_scaling[1].1,
            148.0 / 61.0,
            1e-12,
            "right window fragment-average scaling",
        );
        assert_f64_close(
            per_window_scaling[0].2,
            1.0,
            1e-12,
            "left window full-fragment overlap fraction",
        );
        assert_f64_close(
            per_window_scaling[1].2,
            1.0,
            1e-12,
            "right window full-fragment overlap fraction",
        );

        Ok(())
    }
}

#[cfg(test)]
mod tests_load_scaling_factors_tsv {
    use cfdnalab::shared::{
        bam::Contigs,
        scale_genome::{load_scaling_factors_tsv, ScalingGCMode},
    };
    use fxhash::FxHashMap;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn contigs_for_chr1(len: u32) -> Contigs {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, len));
        Contigs { contigs }
    }

    fn contigs_for_lengths(entries: &[(&str, u32)]) -> Contigs {
        let mut contigs = FxHashMap::default();
        for (idx, (chromosome, len)) in entries.iter().enumerate() {
            contigs.insert((*chromosome).to_string(), (idx as i32, *len));
        }
        Contigs { contigs }
    }

    fn write_scaling_file(contents: &str) -> anyhow::Result<NamedTempFile> {
        let mut file = NamedTempFile::new()?;
        file.write_all(contents.as_bytes())?;
        Ok(file)
    }

    #[test]
    fn load_scaling_factors_tsv_defaults_to_unknown_when_metadata_is_absent() -> anyhow::Result<()>
    {
        // The file has no metadata comments, so parsing should keep the GC mode as Unknown instead
        // of guessing "uncorrected". The two bins still fully cover chr1: [0,5) and [5,10).
        let file = write_scaling_file(
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t5\t1.25\nchr1\t5\t10\t0.75\n",
        )?;

        let loaded =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))?;

        assert_eq!(loaded.metadata.gc_mode, ScalingGCMode::Unknown);
        assert_eq!(
            loaded.bins_by_chromosome.get("chr1"),
            Some(&vec![(0, 5, 1.25_f32), (5, 10, 0.75_f32)])
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_reads_explicit_gc_mode_before_header() -> anyhow::Result<()> {
        // The file starts with one GC-mode line and one unrelated comment line, then a normal
        // header and two bins that fully cover chr1: [0,5) and [5,10). The parser should keep the
        // richer `corrected_tag` source information.
        let file = write_scaling_file(
            "# gc_mode=corrected_tag\n# generated_by=test\nchromosome\tstart\tend\tscaling_factor\nchr1\t0\t5\t1.25\nchr1\t5\t10\t0.75\n",
        )?;

        let loaded =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))?;

        assert!(
            matches!(loaded.metadata.gc_mode, ScalingGCMode::CorrectedFromTag),
            "expected corrected_tag GC mode to be preserved"
        );
        assert_eq!(
            loaded.bins_by_chromosome.get("chr1"),
            Some(&vec![(0, 5, 1.25_f32), (5, 10, 0.75_f32)])
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_invalid_gc_mode_value() -> anyhow::Result<()> {
        // `gc_mode` is enumerated metadata, so any other value should fail before row parsing.
        let file = write_scaling_file(
            "# gc_mode=maybe\nchromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t1.0\n",
        )?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("invalid metadata value should fail");

        assert!(
            err.to_string().contains("invalid value 'maybe'"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_duplicate_gc_mode_metadata() -> anyhow::Result<()> {
        // Metadata keys should be unambiguous. Two `gc_mode` lines before the header would let the
        // second one silently overwrite the first, so that input must fail during header parsing.
        let file = write_scaling_file(
            "# gc_mode=uncorrected\n# gc_mode=corrected_tag\nchromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t1.0\n",
        )?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("duplicate gc_mode metadata should fail");

        assert!(
            err.to_string().contains("duplicate scaling metadata key 'gc_mode'"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_blank_line_before_header() -> anyhow::Result<()> {
        // The format allows metadata comments before the header, but a blank line there is
        // ambiguous and should not be silently skipped when the header is required.
        let file =
            write_scaling_file("\nchromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t1.0\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("blank line before header should fail");

        assert!(
            err.to_string()
                .contains("blank lines are not allowed before the scaling TSV header"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_matches_header_case_insensitively() -> anyhow::Result<()> {
        // Column lookup should ignore case, so an uppercase header must still parse as the same
        // required four fields.
        let file =
            write_scaling_file("CHROMOSOME\tSTART\tEND\tSCALING_FACTOR\nchr1\t0\t10\t1.5\n")?;

        let loaded =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))?;

        assert_eq!(
            loaded.bins_by_chromosome.get("chr1"),
            Some(&vec![(0, 10, 1.5_f32)])
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_missing_required_header_column() -> anyhow::Result<()> {
        // Omitting `scaling_factor` from the header must fail before any row parsing starts.
        let file = write_scaling_file("chromosome\tstart\tend\nchr1\t0\t10\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("missing required header column should fail");

        assert!(
            err.to_string()
                .contains("required column 'scaling_factor' not found in header"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_short_data_row() -> anyhow::Result<()> {
        // A row missing the rightmost required field must fail with the line number and the
        // expected number of columns.
        let file =
            write_scaling_file("chromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("short row should fail");

        assert!(
            err.to_string().contains("not enough columns"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_invalid_interval() -> anyhow::Result<()> {
        // Row-level coordinate validation happens before chromosome-level contiguity checks, so a
        // zero-width interval [5,5) must fail immediately.
        let file =
            write_scaling_file("chromosome\tstart\tend\tscaling_factor\nchr1\t5\t5\t1.0\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("empty interval should fail");

        assert!(
            err.to_string().contains("invalid interval [5..5)"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_rejects_negative_scaling_factor() -> anyhow::Result<()> {
        // Scaling factors are multiplicative weights and must be finite and non-negative.
        let file =
            write_scaling_file("chromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t-1.0\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("negative scaling factor should fail");

        assert!(
            err.to_string().contains("scaling_factor must be finite and >= 0"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_sorts_requested_chromosome_bins_by_start() -> anyhow::Result<()> {
        // The loader promises sorted bins per chromosome, so an out-of-order but otherwise valid
        // file should still load as contiguous [0,5) then [5,10).
        let file = write_scaling_file(
            "chromosome\tstart\tend\tscaling_factor\nchr1\t5\t10\t0.75\nchr1\t0\t5\t1.25\n",
        )?;

        let loaded =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))?;

        assert_eq!(
            loaded.bins_by_chromosome.get("chr1"),
            Some(&vec![(0, 5, 1.25_f32), (5, 10, 0.75_f32)])
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_ignores_unrequested_chromosomes() -> anyhow::Result<()> {
        // Rows for other chromosomes should be filtered out before storage. Only the requested
        // chromosome must remain in the returned map.
        let file = write_scaling_file(
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t1.0\nchr2\t0\t5\t2.0\nchr2\t5\t10\t3.0\n",
        )?;

        let loaded =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))?;

        assert_eq!(loaded.bins_by_chromosome.len(), 1);
        assert_eq!(
            loaded.bins_by_chromosome.get("chr1"),
            Some(&vec![(0, 10, 1.0_f32)])
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_errors_when_requested_chromosome_has_no_bins() -> anyhow::Result<()>
    {
        // Filtering out other chromosomes must not silently succeed when the requested chromosome
        // ends up with no bins at all.
        let file = write_scaling_file(
            "chromosome\tstart\tend\tscaling_factor\nchr2\t0\t10\t1.0\n",
        )?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("requested chromosome without bins should fail");

        assert!(
            err.to_string()
                .contains("scaling TSV: no bins provided for chromosome 'chr1'"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_errors_when_bins_do_not_start_at_zero() -> anyhow::Result<()> {
        // Full chromosome coverage must begin at 0, so a first bin [5,10) is invalid even if it
        // would otherwise reach the contig end.
        let file =
            write_scaling_file("chromosome\tstart\tend\tscaling_factor\nchr1\t5\t10\t1.0\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("bins that start after 0 should fail");

        assert!(
            err.to_string()
                .contains("scaling TSV: bins on 'chr1' must start at 0"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_errors_when_bins_have_a_gap() -> anyhow::Result<()> {
        // Contiguous half-open bins [0,5) and [6,10) leave one uncovered base at position 5, so
        // the chromosome-level sweep must reject them as non-contiguous.
        let file = write_scaling_file(
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t5\t1.0\nchr1\t6\t10\t1.0\n",
        )?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("gapped bins should fail");

        assert!(
            err.to_string().contains("not contiguous"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_errors_when_bins_overlap() -> anyhow::Result<()> {
        // Overlapping bins [0,6) and [5,10) break the same contiguity invariant from the other
        // direction because the second bin starts before the previous one ended.
        let file = write_scaling_file(
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t6\t1.0\nchr1\t5\t10\t1.0\n",
        )?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("overlapping bins should fail");

        assert!(
            err.to_string().contains("not contiguous"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_errors_when_bins_do_not_reach_contig_end() -> anyhow::Result<()> {
        // A single bin [0,8) on a 10 bp chromosome leaves the tail uncovered, so full-coverage
        // validation must reject it.
        let file =
            write_scaling_file("chromosome\tstart\tend\tscaling_factor\nchr1\t0\t8\t1.0\n")?;

        let err =
            load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs_for_chr1(10))
                .expect_err("truncated chromosome coverage should fail");

        assert!(
            err.to_string().contains("must end at chrom_len=10 (got end=8)"),
            "unexpected error: {err}"
        );

        Ok(())
    }

    #[test]
    fn load_scaling_factors_tsv_errors_when_requested_contig_metadata_is_missing()
    -> anyhow::Result<()> {
        // The loader validates full chromosome coverage against BAM contig lengths, so requesting a
        // chromosome missing from `contigs` must fail even if the TSV rows themselves look valid.
        let file =
            write_scaling_file("chromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t1.0\n")?;
        let contigs = contigs_for_lengths(&[("chr2", 10)]);

        let err = load_scaling_factors_tsv(file.path(), &["chr1".to_string()], &contigs)
            .expect_err("missing contig metadata should fail");

        assert!(
            err.to_string().contains("missing contig info for 'chr1'"),
            "unexpected error: {err}"
        );

        Ok(())
    }
}

#[cfg(test)]
mod tests_scaling_gc_compatibility {
    use cfdnalab::shared::scale_genome::{
        ensure_scaling_gc_compatibility, scaling_gc_mode_for_run, ScalingFactorsMetadata,
        ScalingGCMode,
    };
    use std::path::Path;

    #[test]
    fn scaling_gc_compatibility_allows_unknown_scaling_gc_mode() -> anyhow::Result<()> {
        let path = Path::new("/tmp/example.scaling_factors.tsv");
        let metadata = ScalingFactorsMetadata {
            gc_mode: ScalingGCMode::Unknown,
        };

        // Missing metadata should never block use, even when the current run uses GC correction.
        ensure_scaling_gc_compatibility(path, metadata, scaling_gc_mode_for_run(true, false))?;

        Ok(())
    }

    #[test]
    fn scaling_gc_compatibility_errors_when_uncorrected_scaling_meets_corrected_run() {
        // This mismatch is now a hard error:
        // - scaling file says uncorrected coverage
        // - current run uses tag-based GC correction
        let path = Path::new("/tmp/example.scaling_factors.tsv");
        let metadata = ScalingFactorsMetadata {
            gc_mode: ScalingGCMode::Uncorrected,
        };

        let err =
            ensure_scaling_gc_compatibility(path, metadata, scaling_gc_mode_for_run(false, true))
                .expect_err("known raw-vs-corrected mismatch should error");

        assert!(
            err.to_string().contains("no GC correction"),
            "unexpected error text: {err}"
        );
        assert!(
            err.to_string().contains("via --gc-tag"),
            "unexpected error text: {err}"
        );
    }

    #[test]
    fn scaling_gc_compatibility_errors_when_corrected_scaling_meets_uncorrected_run() {
        // This is the opposite known mismatch:
        // - scaling file says corrected coverage from a tag
        // - current run applies no GC correction
        let path = Path::new("/tmp/example.scaling_factors.tsv");
        let metadata = ScalingFactorsMetadata {
            gc_mode: ScalingGCMode::CorrectedFromTag,
        };

        let err =
            ensure_scaling_gc_compatibility(path, metadata, scaling_gc_mode_for_run(false, false))
                .expect_err("known corrected-vs-raw mismatch should error");

        assert!(
            err.to_string().contains("via --gc-tag"),
            "unexpected error text: {err}"
        );
        assert!(
            err.to_string().contains("no GC correction"),
            "unexpected error text: {err}"
        );
    }

    #[test]
    fn scaling_gc_compatibility_allows_corrected_file_scaling_with_corrected_tag_run()
    -> anyhow::Result<()> {
        // File-vs-tag is allowed because both sides are still GC-corrected:
        // coverage-weights may have used `--gc-file`, while the later command uses tagged reads.
        let path = Path::new("/tmp/example.scaling_factors.tsv");
        let metadata = ScalingFactorsMetadata {
            gc_mode: ScalingGCMode::CorrectedFromFile,
        };

        ensure_scaling_gc_compatibility(path, metadata, scaling_gc_mode_for_run(false, true))?;

        Ok(())
    }

    #[test]
    fn scaling_gc_compatibility_allows_corrected_tag_scaling_with_corrected_file_run()
    -> anyhow::Result<()> {
        // The allowed file-vs-tag case must be symmetric. A scaling file built from tagged reads
        // is still compatible with a later run that uses `--gc-file`.
        let path = Path::new("/tmp/example.scaling_factors.tsv");
        let metadata = ScalingFactorsMetadata {
            gc_mode: ScalingGCMode::CorrectedFromTag,
        };

        ensure_scaling_gc_compatibility(path, metadata, scaling_gc_mode_for_run(true, false))?;

        Ok(())
    }

    #[test]
    fn scaling_gc_compatibility_allows_uncorrected_scaling_with_uncorrected_run()
    -> anyhow::Result<()> {
        // The explicit `uncorrected` mode should only reject corrected runs, not another raw run.
        let path = Path::new("/tmp/example.scaling_factors.tsv");
        let metadata = ScalingFactorsMetadata {
            gc_mode: ScalingGCMode::Uncorrected,
        };

        ensure_scaling_gc_compatibility(path, metadata, scaling_gc_mode_for_run(false, false))?;

        Ok(())
    }
}
