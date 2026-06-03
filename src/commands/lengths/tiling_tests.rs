mod tests_lengths_tiling_reducer {
    use crate::commands::lengths::counting::{LengthAxis, LengthCounts};
    use crate::commands::lengths::tiling::{
        reduce_partials_for_chr, write_cross_npy as write_cross_npy_inner,
        write_partials_npz as write_partials_npz_inner,
    };
    use crate::shared::temp_chrom_names::TempChromNameMap;
    use anyhow::Result;
    use ndarray::{Array1, Array2, ShapeBuilder};
    use ndarray_npy::NpzWriter;
    use std::sync::Arc;
    use std::{fs::File, path::PathBuf};
    use tempfile::TempDir;

    fn exact_axis(min_length: u32, max_length: u32) -> Arc<LengthAxis> {
        let edges: Vec<u32> = (min_length..=max_length + 1).collect();
        Arc::new(LengthAxis::new(edges).expect("test length axis should be valid"))
    }

    fn template_counts() -> LengthCounts {
        LengthCounts::new(exact_axis(10, 10))
    }

    fn counts_with_value(val: f64) -> LengthCounts {
        let mut lc = template_counts();
        lc.counts[0] = val;
        lc
    }

    fn expect_written_path(path: Option<PathBuf>, label: &str) -> PathBuf {
        path.unwrap_or_else(|| panic!("{label} should have been written"))
    }

    fn temp_chrom_name_map(chromosomes: &[&str]) -> TempChromNameMap {
        TempChromNameMap::from_contigs(
            &chromosomes
                .iter()
                .map(|chromosome| chromosome.to_string())
                .collect::<Vec<_>>(),
        )
        .expect("test contig temp name map should be valid")
    }

    fn write_partials_npz(
        temp_dir: &std::path::Path,
        prefix: &str,
        chr: &str,
        tile_idx: u32,
        window_idxs_chr: &[u64],
        contained_flags: &[bool],
        counts: &[LengthCounts],
    ) -> Result<Option<PathBuf>> {
        write_partials_npz_inner(
            temp_dir,
            prefix,
            chr,
            tile_idx,
            &temp_chrom_name_map(&[chr]),
            window_idxs_chr,
            contained_flags,
            counts,
        )
    }

    fn write_partials_npz_with_chromosomes(
        chromosomes: &[&str],
        temp_dir: &std::path::Path,
        prefix: &str,
        chr: &str,
        tile_idx: u32,
        window_idxs_chr: &[u64],
        contained_flags: &[bool],
        counts: &[LengthCounts],
    ) -> Result<Option<PathBuf>> {
        write_partials_npz_inner(
            temp_dir,
            prefix,
            chr,
            tile_idx,
            &temp_chrom_name_map(chromosomes),
            window_idxs_chr,
            contained_flags,
            counts,
        )
    }

    fn write_cross_npy(
        temp_dir: &std::path::Path,
        prefix: &str,
        chr: &str,
        tile_idx: u32,
        crossing_window_idxs_chr: &[u64],
    ) -> Result<Option<PathBuf>> {
        write_cross_npy_inner(
            temp_dir,
            prefix,
            chr,
            tile_idx,
            &temp_chrom_name_map(&[chr]),
            crossing_window_idxs_chr,
        )
    }

    #[test]
    fn reducer_accepts_contained_only() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(2.0)];
        let contained = vec![true];
        let partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?,
            "contained partial",
        );
        // No cross file because window is contained

        let reduced = reduce_partials_for_chr("chr1", &[partial_path], &[], 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn partial_writer_rejects_count_rows_without_matching_window_indices() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();

        let counts = vec![counts_with_value(2.0)];
        let contained = vec![true, false];

        let err = write_partials_npz(dir, "partials", "chr1", 0, &[0, 1], &contained, &counts)
            .expect_err("partial writer should reject mismatched row metadata");

        assert!(err.to_string().contains("counts length mismatch"));
    }

    #[test]
    fn reducer_counts_multiple_crossing_tiles() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![false];
        // Two tiles, both crossing the same window
        let partial_paths = vec![
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?,
                "first crossing partial",
            ),
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 1, &[0], &contained, &counts)?,
                "second crossing partial",
            ),
        ];
        let cross_paths = vec![
            expect_written_path(
                write_cross_npy(dir, "cross", "chr1", 0, &[0])?,
                "first cross",
            ),
            expect_written_path(
                write_cross_npy(dir, "cross", "chr1", 1, &[0])?,
                "second cross",
            ),
        ];

        let reduced = reduce_partials_for_chr(
            "chr1",
            partial_paths.as_slice(),
            cross_paths.as_slice(),
            1,
            &template,
        )?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_combines_contained_and_cross() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let contained_counts = vec![counts_with_value(1.0)];
        let crossing_counts = vec![counts_with_value(3.0)];
        let partial_paths = vec![
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 0, &[0], &[true], &contained_counts)?,
                "contained partial",
            ),
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 1, &[0], &[false], &crossing_counts)?,
                "crossing partial",
            ),
        ];
        let cross_paths = vec![expect_written_path(
            write_cross_npy(dir, "cross", "chr1", 1, &[0])?,
            "cross index",
        )];

        let reduced = reduce_partials_for_chr(
            "chr1",
            partial_paths.as_slice(),
            cross_paths.as_slice(),
            1,
            &template,
        )?;
        assert_eq!(reduced.len(), 1);
        // Expect 1 contained contribution and 1 crossing contribution => sum counts
        assert!((reduced[0].counts[0] - 4.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_errors_when_contribution_missing() {
        let template = template_counts();

        // No partials written -> zero contributions
        let err = reduce_partials_for_chr("chr1", &[], &[], 1, &template)
            .expect_err("should fail when contributions are missing");
        assert!(err.to_string().contains("expected 1"));
    }

    #[test]
    fn reducer_errors_on_mismatched_counts() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        // Cross file claims one contribution, but no partial exists
        let cross_path = expect_written_path(
            write_cross_npy(dir, "cross", "chr1", 0, &[0]).unwrap(),
            "cross",
        );

        let err = reduce_partials_for_chr("chr1", &[], &[cross_path], 1, &template)
            .expect_err("should fail when expected contributions not met");
        assert!(err.to_string().contains("expected 1"));
    }

    #[test]
    fn reducer_errors_on_counts_width_mismatch() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts(); // width 1

        let idxs = Array1::from(vec![0u64]);
        let contained = Array1::from(vec![1u8]);
        let counts = Array2::from_shape_vec((1, 2), vec![1.0, 0.5]).unwrap();
        let path = dir.join("partials.chr1.0.npz");
        let file = File::create(&path).unwrap();
        let mut npz = NpzWriter::new(file);
        npz.add_array("window_idx_chr", &idxs).unwrap();
        npz.add_array("contained", &contained).unwrap();
        npz.add_array("counts", &counts).unwrap();
        npz.finish().unwrap();

        let err = reduce_partials_for_chr("chr1", &[path], &[], 1, &template)
            .expect_err("should fail on counts width mismatch");
        assert!(err.to_string().contains("counts width mismatch"));
    }

    #[test]
    fn reducer_errors_on_non_contiguous_counts_rows() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = LengthCounts::new(exact_axis(10, 11)); // two-length template

        // Two rows, both targeting window 0, but stored in Fortran order so each
        // row slice is non-contiguous and should be rejected.
        let idxs = Array1::from(vec![0u64, 0u64]);
        let contained = Array1::from(vec![1u8, 1u8]);
        let counts = Array2::from_shape_vec((2, 2).f(), vec![1.0, 0.5, 2.0, 1.5]).unwrap();
        let path = dir.join("partials.chr1.0.npz");
        let file = File::create(&path).unwrap();
        let mut npz = NpzWriter::new(file);
        npz.add_array("window_idx_chr", &idxs).unwrap();
        npz.add_array("contained", &contained).unwrap();
        npz.add_array("counts", &counts).unwrap();
        npz.finish().unwrap();

        let err = reduce_partials_for_chr("chr1", &[path], &[], 1, &template)
            .expect_err("should fail on non-contiguous counts rows");
        assert!(err.to_string().contains("counts row not contiguous"));
    }

    #[test]
    fn reducer_ignores_files_from_other_chromosomes() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![true];
        let partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?,
            "chr1 partial",
        );

        // Stray files for another chromosome are not passed to the reducer because the command
        // now records exact output paths from each tile.
        write_partials_npz(dir, "partials", "chr2", 0, &[0], &contained, &counts)?;
        write_cross_npy(dir, "cross", "chr2", 0, &[0])?;

        let reduced = reduce_partials_for_chr("chr1", &[partial_path], &[], 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_uses_explicit_paths_for_overlapping_chromosome_names() -> Result<()> {
        // Human verification status: verified
        // Manual expectations:
        // - `chr1` contributes one contained count with value 1.
        // - `chr1.extra` contributes one contained count with value 5.
        // - Old dotted-substring discovery for `chr1` could also match `partials.chr1.extra.0.npz`.
        // - Passing explicit paths means the `chr1.extra` file is ignored unless the caller provides
        //   it, so the reduced `chr1` count remains 1.
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let contained = vec![true];
        let chr1_counts = vec![counts_with_value(1.0)];
        let chr1_extra_counts = vec![counts_with_value(5.0)];
        let chr1_partial_path = expect_written_path(
            write_partials_npz_with_chromosomes(
                &["chr1", "chr1.extra"],
                dir,
                "partials",
                "chr1",
                0,
                &[0],
                &contained,
                &chr1_counts,
            )?,
            "chr1 partial",
        );
        let _chr1_extra_partial_path = expect_written_path(
            write_partials_npz_with_chromosomes(
                &["chr1", "chr1.extra"],
                dir,
                "partials",
                "chr1.extra",
                0,
                &[0],
                &contained,
                &chr1_extra_counts,
            )?,
            "chr1.extra partial",
        );

        let reduced = reduce_partials_for_chr("chr1", &[chr1_partial_path], &[], 1, &template)?;

        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn write_partials_rejects_mismatched_contained() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();
        let counts = vec![template];
        let contained = vec![true, false]; // Wrong length

        let err = write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)
            .expect_err("should error on contained/idx length mismatch");
        assert!(err.to_string().contains("contained flags length mismatch"));
    }

    #[test]
    fn reducer_errors_on_out_of_bounds_partial_idx() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![false];
        // Write a partial with idx outside n_windows=1
        let partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[2], &contained, &counts).unwrap(),
            "out-of-bounds partial",
        );

        let err = reduce_partials_for_chr("chr1", &[partial_path], &[], 1, &template)
            .expect_err("should fail on out-of-bounds idx");
        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn reducer_errors_on_out_of_bounds_cross_idx() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        let cross_path = expect_written_path(
            write_cross_npy(dir, "cross", "chr1", 0, &[3]).unwrap(),
            "out-of-bounds cross",
        );
        let err = reduce_partials_for_chr("chr1", &[], &[cross_path], 1, &template)
            .expect_err("should fail on cross index out of bounds");
        assert!(err.to_string().contains("Cross index"));
    }

    #[test]
    fn reducer_separates_windows() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts0 = vec![counts_with_value(1.0)];
        let counts1 = vec![counts_with_value(2.0)];
        let contained = vec![true];

        // Window 0 contained in tile 0, window 1 contained in tile 1
        let partial_paths = vec![
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts0)?,
                "first window partial",
            ),
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 1, &[1], &contained, &counts1)?,
                "second window partial",
            ),
        ];

        let reduced = reduce_partials_for_chr("chr1", partial_paths.as_slice(), &[], 2, &template)?;
        assert_eq!(reduced.len(), 2);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        assert!((reduced[1].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn write_partials_skips_empty() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();
        let res = write_partials_npz(dir, "partials", "chr1", 0, &[], &[], &[])?;
        assert!(res.is_none());
        // Ensure reducer still errors because nothing was written
        let err = reduce_partials_for_chr("chr1", &[], &[], 1, &template)
            .expect_err("should fail when nothing written");
        assert!(err.to_string().contains("expected 1"));
        Ok(())
    }

    #[test]
    fn write_cross_skips_empty() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let res = write_cross_npy(dir, "cross", "chr1", 0, &[])?;
        assert!(res.is_none());
        Ok(())
    }
}

mod tests_lengths_tiling_helpers {
    // MOVE-MODULE-LOCAL: direct private shared tiling and fetch helper tests.

    use crate::run_like_cli::common::WindowSpec;
    use crate::shared::bam::Contigs;
    use crate::shared::interval::IndexedInterval;
    use crate::shared::tiled_run::{Tile, TileWindowSpan, build_tiles};
    use crate::shared::window_fetch::{BedFetchPolicy, fetch_span_for_tile};
    use fxhash::FxHashMap;
    use std::path::PathBuf;

    fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
        entries
            .iter()
            .map(|&(start, end, original_index)| {
                IndexedInterval::new(start, end, original_index)
                    .expect("test windows should be valid non-empty intervals")
            })
            .collect()
    }

    #[test]
    fn fetch_span_size_mode_clamps_to_halo_and_chrom() {
        // Tile: core 50-150, fetch 30-200 (halo 20 left, 50 right), chrom len 180
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 50, 150, 30, 200)
            .expect("test tile should be valid");
        let span = fetch_span_for_tile(
            &tile,
            None,
            None,
            &WindowSpec::Size(100),
            180,
            0,
            BedFetchPolicy::CandidateWindowExtent,
        )
        .expect("span expected")
        .expect("fetch span expected");
        // Window span touching core: 0..200, after halo clamp -> 30..180
        assert_eq!(span.start(), 30);
        assert_eq!(span.end(), 180);
    }

    #[test]
    fn build_tiles_aligns_to_bin_when_divisible() {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, 100u32));
        let contigs = Contigs { contigs };
        let (tiles, aligned) =
            build_tiles(&vec!["chr1".to_string()], &contigs, 30, 0, Some(10)).unwrap();
        assert!(aligned);
        // Cores should start on multiples of 10
        for t in &tiles {
            assert_eq!((t.core_start() as u64) % 10, 0);
        }
        // Expect four tiles: 0-30,30-60,60-90,90-100
        assert_eq!(tiles.len(), 4);
        assert_eq!(tiles[0].core_end(), 30);
        assert_eq!(tiles[3].core_start(), 90);
        assert_eq!(tiles[3].core_end(), 100);
    }

    #[test]
    fn build_tiles_not_aligned_when_too_few_bins() {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, 50u32));
        let contigs = Contigs { contigs };
        // With tile_bp=15 and align_bp=10, only one full 10bp bin fits,
        // and build_tiles requires at least 10 bins (k >= 10) before rounding down.
        // So alignment should be disabled and tiles keep the original 15bp size.
        let (_tiles, aligned) =
            build_tiles(&vec!["chr1".to_string()], &contigs, 15, 0, Some(10)).unwrap();
        assert!(!aligned);
    }

    #[test]
    fn fetch_span_for_tile_global_clamps_to_chrom() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 0, 50, 0, 200)
            .expect("test tile should be valid");
        let span = fetch_span_for_tile(
            &tile,
            None,
            None,
            &WindowSpec::Global,
            120,
            0,
            BedFetchPolicy::CandidateWindowExtent,
        )
        .expect("span")
        .expect("fetch span expected");
        assert_eq!(span.start(), 0);
        assert_eq!(span.end(), 120);
    }

    #[test]
    fn fetch_span_for_tile_bed_with_overlap() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 100, 160, 80, 200)
            .expect("test tile should be valid");
        let windows = indexed_windows(&[(90, 110, 0), (150, 170, 1), (250, 300, 2)]);
        let span = TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 2,
        };
        let res = fetch_span_for_tile(
            &tile,
            Some(&span),
            Some(&windows),
            &WindowSpec::Bed(PathBuf::from("dummy")),
            500,
            0,
            BedFetchPolicy::CoreOverlap,
        )
        .expect("span")
        .expect("fetch span expected");
        // min_ws=90, max_we=170, halos: left 20, right 40 -> widened to 70..210, clamped to fetch
        assert_eq!(res.start(), 80);
        assert_eq!(res.end(), 200);
    }

    #[test]
    fn fetch_span_bed_none_when_no_overlap() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 100, 150, 80, 170)
            .expect("test tile should be valid");
        // No windows overlap tile
        let windows: [IndexedInterval<u64>; 0] = [];
        let span = TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 0,
        };
        let res = fetch_span_for_tile(
            &tile,
            Some(&span),
            Some(&windows),
            &WindowSpec::Bed(PathBuf::from("dummy")),
            200,
            0,
            BedFetchPolicy::CoreOverlap,
        )
        .expect("fetch span computation should succeed");
        assert!(res.is_none());
    }

    #[test]
    fn fetch_span_size_mode_none_when_tile_right_of_chromosome() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 250, 260, 230, 270)
            .expect("test tile should be valid");
        let res = fetch_span_for_tile(
            &tile,
            None,
            None,
            &WindowSpec::Size(50),
            200,
            0,
            BedFetchPolicy::CandidateWindowExtent,
        )
        .expect("fetch span computation should succeed");
        assert!(res.is_none());
    }

    #[test]
    fn tile_constructor_rejects_empty_core() {
        let err = Tile::from_coords("chr1".to_string(), 0, 0, 100, 100, 80, 120).unwrap_err();
        assert!(format!("{err}").contains("interval end (100) must be greater than start (100)"));
    }
}
