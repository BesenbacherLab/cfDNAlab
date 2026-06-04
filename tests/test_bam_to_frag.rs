#![cfg(feature = "cmd_bam_to_frag")]

mod fixtures;

mod tests_bam_to_frag {
    // tests/bam_to_frag_integration.rs

    use anyhow::{Context, Result};
    use flate2::read::GzDecoder;
    use ndarray::array;
    use rust_htslib::bam::index;
    use rust_htslib::bam::{
        self, Format, HeaderView, Writer,
        header::Header,
        record::{Cigar, CigarString, Record},
    };
    use std::{
        collections::HashMap,
        fs::{self, File},
        io::Read as _,
        path::Path,
    };
    use tempfile::tempdir;

    // Bring your crate items into scope.
    use cfdnalab::RunOptions;
    use cfdnalab::gc_bias::GCCorrectionPackage;
    #[cfg(feature = "cmd_bam_to_bam")]
    use cfdnalab::run_like_cli::bam_to_bam::{
        BamToBamConfig, BamToBamRunResult, run_bam_to_bam as run_bam_to_bam_command,
    };
    use cfdnalab::run_like_cli::bam_to_frag::{
        BamToFragConfig, BamToFragRunResult, run_bam_to_frag,
    };
    use cfdnalab::run_like_cli::common::{ApplyGCArgFileOnly, ChromosomeArgs, IOCArgs};
    #[cfg(feature = "cmd_coverage_weights")]
    use cfdnalab::run_like_cli::coverage_weights::{
        CoverageWeightsConfig, run_coverage_weights as run_coverage_weights_command,
    };
    #[cfg(feature = "cmd_fragment_count_weights")]
    use cfdnalab::run_like_cli::fragment_count_weights::{
        FragmentCountWeightsConfig,
        run_fragment_count_weights as run_fragment_count_weights_command,
    };
    use cfdnalab::testing::{
        PairedFragmentSpec, TempBam, bam_from_fragments, single_contig_inward_pair_bam,
        twobit_from_sequences, twobit_with_single_repeating_contig,
    };
    #[cfg(feature = "cmd_bam_to_bam")]
    use cfdnalab::testing::{
        build_command_produced_gc_correction_package_for_length,
        build_command_produced_gc_correction_package_from_reference_windows,
    };
    use cfdnalab::{
        constants::GC_CORRECTION_SCHEMA_VERSION,
        reference::{ContigFootprintEntry, twobit_contig_footprint},
    };
    #[cfg(feature = "cmd_bam_to_bam")]
    use rust_htslib::bam::Read;
    #[cfg(feature = "cmd_bam_to_bam")]
    use rust_htslib::bam::record::Aux;

    fn run_bam_to_frag_for_test(cfg: &BamToFragConfig) -> Result<BamToFragRunResult> {
        run_bam_to_frag(cfg, RunOptions::new_quiet())
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    fn run_bam_to_bam_for_test(cfg: &BamToBamConfig) -> Result<BamToBamRunResult> {
        run_bam_to_bam_command(cfg, RunOptions::new_quiet())
    }

    #[cfg(feature = "cmd_coverage_weights")]
    fn run_coverage_weights(cfg: &CoverageWeightsConfig) -> Result<()> {
        run_coverage_weights_command(cfg, RunOptions::new_quiet()).map(|_| ())
    }

    #[cfg(feature = "cmd_fragment_count_weights")]
    fn run_fragment_count_weights(cfg: &FragmentCountWeightsConfig) -> Result<()> {
        run_fragment_count_weights_command(cfg, RunOptions::new_quiet()).map(|_| ())
    }

    #[test]
    fn bam_to_frag_smoke_two_chroms() -> Result<()> {
        // Temp working dir
        let work = tempdir().context("tempdir")?;
        let work_path = work.path();

        // Paths
        let bam_path = work_path.join("test.bam");
        let bai_path = work_path.join("test.bam.bai");
        let out_dir = work_path.join("out");
        fs::create_dir_all(&out_dir)?;

        // Create a tiny coordinate-sorted BAM with two chromosomes and three pairs:
        // chr1: pair A (R1 forward, R2 reverse, MAPQ=60) -> expect strand '+'
        // chr1: pair B (R1 forward, R2 reverse, MAPQ=0)  -> expect strand '+'
        // chr2: pair C (R1 reverse, R2 forward, min MAPQ=30) -> expect strand '-'
        write_test_bam(&bam_path)?;

        // Build BAI index
        index::build(
            bam_path.to_str().unwrap(),
            None,
            index::Type::Bai,
            1, // n_threads for indexing
        )
        .context("build BAI")?;
        assert!(bai_path.exists(), "BAI was not created");

        // Construct CLI config
        let ioc = IOCArgs {
            bam: bam_path.clone(),
            output_dir: out_dir.clone(),
            n_threads: 2,
        };

        // Limit to the contigs we wrote into the BAM so ChromosomeArgs resolution does not fail.
        let chromosomes = fixed_chromosome_args();

        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("sample");
        cfg.set_min_mapq(0); // keep MAPQ=0 to test min mapq behavior
        cfg.set_require_proper_pair(false); // we never require proper pair
        cfg.set_blacklist(None);

        // Run the command
        let counters = run_bam_to_frag_for_test(&cfg)
            .context("bam-to-frag run failed")?
            .counters;
        assert!(
            counters.base.counted_fragments >= 3,
            "Expected at least 3 fragments counted"
        );

        // Load merged frag file
        let frag_path = out_dir.join("sample.frag.tsv.gz");
        let recs = read_frag_gz(&frag_path)?;

        // Expect exactly 3 lines
        assert_eq!(
            recs.len(),
            3,
            "Unexpected number of fragment lines: {:?}",
            recs
        );

        // Parse into tuples for easier assertions
        // columns: chrom, start, end, mapq, strand
        let parsed: Vec<(String, u64, u64, u8, char)> = recs
            .iter()
            .map(|line| {
                let parts: Vec<&str> = line.split('\t').collect();
                assert_eq!(parts.len(), 5, "Bad line: {}", line);
                let chrom = parts[0].to_string();
                let start: u64 = parts[1].parse().unwrap();
                let end: u64 = parts[2].parse().unwrap();
                let mapq: u8 = parts[3].parse().unwrap();
                let strand: char = parts[4].chars().next().unwrap();
                (chrom, start, end, mapq, strand)
            })
            .collect();

        // Sort for stable assertions
        let mut parsed = parsed;
        parsed.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // chr1 pair A: R1 forward at 10002, R2 reverse at 10090 with 50M -> end=10140, min MQ = 60, strand '+'
        assert_eq!(parsed[0], ("chr1".into(), 10002, 10140, 60, '+'));

        // chr1 pair B: R1 forward at 10003, R2 reverse at 10087 with 50M -> end=10137, min MQ = 0, strand '+'
        assert_eq!(parsed[1], ("chr1".into(), 10003, 10137, 0, '+'));

        // chr2 pair C: R1 reverse at 20090, R2 forward at 20000 with 50M -> start=20000, end=20140, min MQ = 30, strand '-'
        assert_eq!(parsed[2], ("chr2".into(), 20000, 20140, 30, '-'));

        Ok(())
    }

    #[test]
    fn bam_to_frag_global_handles_three_chromosomes() -> Result<()> {
        let work = tempdir().context("tempdir")?;
        let bam_path = work.path().join("three_chr_global.bam");
        let out_dir = work.path().join("out_global");
        fs::create_dir_all(&out_dir)?;

        write_three_chrom_window_bam(&bam_path)?;

        index::build(bam_path.to_str().unwrap(), None, index::Type::Bai, 1).context("build BAI")?;

        let ioc = IOCArgs {
            bam: bam_path.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let mut cfg = BamToFragConfig::new(ioc, three_chromosome_args());
        cfg.set_output_prefix("three_chr_global");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);

        run_bam_to_frag_for_test(&cfg)?;

        let frag_path = out_dir.join("three_chr_global.frag.tsv.gz");
        let mut parsed = parse_frag_rows(&read_frag_gz(&frag_path)?);
        parsed.sort();

        assert_eq!(
            parsed,
            vec![
                ("chr1".to_string(), 10, 130, 60, '+'),
                ("chr1".to_string(), 400, 520, 60, '+'),
                ("chr2".to_string(), 30, 150, 60, '+'),
                ("chr2".to_string(), 430, 550, 60, '+'),
                ("chr3".to_string(), 50, 170, 60, '+'),
                ("chr3".to_string(), 460, 580, 60, '+'),
            ],
            "Global mode should keep all three chromosomes and both fragments per chromosome"
        );

        Ok(())
    }

    #[test]
    fn global_selection_matches_single_full_chromosome_bed_window() -> Result<()> {
        // Arrange:
        // `single_contig_inward_pair_bam()` contains one fragment spanning [20, 80) on chr1.
        //
        // Compare two logically equivalent selection modes:
        // - default global selection (`by_bed = None`)
        // - one BED window covering the entire chromosome [0, 200)
        //
        // Because `bam-to-frag` uses BED windows only as an inclusion filter, and the window
        // covers the entire chromosome, both runs must emit the exact same frag row:
        //   chr1 20 80 60 +
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let global_out = work.path().join("out_global_equiv");
        let bed_out = work.path().join("out_bed_equiv");
        fs::create_dir_all(&global_out)?;
        fs::create_dir_all(&bed_out)?;
        let bed_path = work.path().join("whole_chr.bed");
        fs::write(&bed_path, "chr1\t0\t200\twhole_chr\n")?;

        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };

        let mut global_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: global_out.clone(),
                n_threads: 1,
            },
            chromosomes.clone(),
        );
        global_cfg.set_output_prefix("global");
        global_cfg.set_min_mapq(0);
        global_cfg.set_require_proper_pair(false);

        let mut bed_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: bed_out.clone(),
                n_threads: 1,
            },
            chromosomes,
        );
        bed_cfg.set_output_prefix("bed");
        bed_cfg.set_min_mapq(0);
        bed_cfg.set_require_proper_pair(false);
        bed_cfg.set_by_bed(Some(bed_path));

        // Act
        run_bam_to_frag_for_test(&global_cfg)?;
        run_bam_to_frag_for_test(&bed_cfg)?;

        // Assert
        let global_rows = read_frag_gz(&global_out.join("global.frag.tsv.gz"))?;
        let bed_rows = read_frag_gz(&bed_out.join("bed.frag.tsv.gz"))?;
        let expected = vec!["chr1\t20\t80\t60\t+".to_string()];
        assert_eq!(global_rows, expected);
        assert_eq!(bed_rows, expected);

        Ok(())
    }

    #[test]
    fn by_bed_excludes_chromosomes_without_any_windows() -> Result<()> {
        // Arrange:
        // Build one fragment on each of two chromosomes:
        // - chr1: [20, 80)
        // - chr2: [120, 180)
        //
        // Then provide a BED file with a single chr1 window [0, 100).
        // In `bam-to-frag`, `--by-bed` is an inclusion filter, so chr1 should be kept and
        // chr2 should contribute nothing at all because that chromosome has no BED rows.
        let chr1_fragment = PairedFragmentSpec::new(0, 20, 60, 30).build()?;
        let mut chr2_fragment = PairedFragmentSpec::new(0, 120, 60, 30).build()?;
        chr2_fragment.forward.tid = 1;
        chr2_fragment.reverse.tid = 1;
        chr2_fragment.forward.mate_tid = Some(1);
        chr2_fragment.reverse.mate_tid = Some(1);
        let bam = bam_from_fragments(
            "bam_to_frag_missing_chr_bed_windows",
            vec![("chr1".to_string(), 300), ("chr2".to_string(), 300)],
            vec![chr1_fragment, chr2_fragment],
            Vec::new(),
        )?;

        let work = tempdir().context("tempdir")?;
        let global_out_dir = work.path().join("out_global");
        let bed_out_dir = work.path().join("out_bed");
        fs::create_dir_all(&global_out_dir)?;
        fs::create_dir_all(&bed_out_dir)?;
        let bed_path = work.path().join("chr1_only.bed");
        fs::write(&bed_path, "chr1\t0\t100\n")?;

        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        };
        let mut global_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: global_out_dir.clone(),
                n_threads: 1,
            },
            chromosomes.clone(),
        );
        global_cfg.set_output_prefix("global");
        global_cfg.set_min_mapq(0);
        global_cfg.set_require_proper_pair(false);

        let mut cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: bed_out_dir.clone(),
                n_threads: 1,
            },
            chromosomes,
        );
        cfg.set_output_prefix("bed_only_chr1");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_by_bed(Some(bed_path));

        // Act
        run_bam_to_frag_for_test(&global_cfg)?;
        run_bam_to_frag_for_test(&cfg)?;

        // Assert:
        // First prove the fixture and output path preserve both chromosomes in global mode.
        // Then show that `--by-bed` removes only the chromosome with no BED rows.
        let global_rows = read_frag_gz(&global_out_dir.join("global.frag.tsv.gz"))?;
        assert_eq!(
            global_rows,
            vec![
                "chr1\t20\t80\t60\t+".to_string(),
                "chr2\t120\t180\t60\t+".to_string(),
            ]
        );

        let bed_rows = read_frag_gz(&bed_out_dir.join("bed_only_chr1.frag.tsv.gz"))?;
        assert_eq!(bed_rows, vec!["chr1\t20\t80\t60\t+".to_string()]);

        Ok(())
    }

    #[test]
    fn chromosomes_all_follows_bam_header_order_not_lexicographic_order() -> Result<()> {
        // Arrange:
        // Build a BAM whose header order is intentionally non-lexicographic:
        //   chr2, chr10, chr1
        //
        // `bam-to-frag` resolves `--chromosomes all` through the BAM header and then
        // concatenates per-chromosome temp files in that resolved order. With one fragment
        // per chromosome, the output row order must therefore be:
        //   chr2 first, chr10 second, chr1 third
        let work = tempdir().context("tempdir")?;
        let bam_path = work.path().join("header_order.bam");
        let out_dir = work.path().join("out_header_order");
        fs::create_dir_all(&out_dir)?;

        let mut hdr = Header::new();
        hdr.push_record(
            bam::header::HeaderRecord::new(b"HD")
                .push_tag(b"VN", &"1.6")
                .push_tag(b"SO", &"coordinate"),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr2")
                .push_tag(b"LN", &1000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr10")
                .push_tag(b"LN", &1000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &1000),
        );

        let mut writer = Writer::from_path(&bam_path, &hdr, Format::Bam).context("create BAM")?;
        let header_view = HeaderView::from_header(&hdr);
        let tid_chr2 = header_view.tid(b"chr2").expect("chr2 present") as i32;
        let tid_chr10 = header_view.tid(b"chr10").expect("chr10 present") as i32;
        let tid_chr1 = header_view.tid(b"chr1").expect("chr1 present") as i32;
        let cigar = vec![Cigar::Match(40)];
        let seq = b"ACGTN".repeat(8);
        let qual = vec![30u8; 40];

        let records = vec![
            make_rec(
                b"chr2_pair",
                tid_chr2,
                10,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr2,
                50,
                true,
            ),
            make_rec(
                b"chr2_pair",
                tid_chr2,
                50,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr2,
                10,
                false,
            ),
            make_rec(
                b"chr10_pair",
                tid_chr10,
                20,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr10,
                60,
                true,
            ),
            make_rec(
                b"chr10_pair",
                tid_chr10,
                60,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr10,
                20,
                false,
            ),
            make_rec(
                b"chr1_pair",
                tid_chr1,
                30,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr1,
                70,
                true,
            ),
            make_rec(
                b"chr1_pair",
                tid_chr1,
                70,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr1,
                30,
                false,
            ),
        ];
        for record in records {
            writer.write(&record)?;
        }
        drop(writer);
        index::build(bam_path.to_str().unwrap(), None, index::Type::Bai, 1).context("build BAI")?;

        let mut cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam_path,
                output_dir: out_dir.clone(),
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["all".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_output_prefix("header_order");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);

        // Act
        run_bam_to_frag_for_test(&cfg)?;

        // Assert
        let rows = parse_frag_rows(&read_frag_gz(&out_dir.join("header_order.frag.tsv.gz"))?);
        assert_eq!(
            rows,
            vec![
                ("chr2".to_string(), 10, 90, 60, '+'),
                ("chr10".to_string(), 20, 100, 60, '+'),
                ("chr1".to_string(), 30, 110, 60, '+'),
            ]
        );

        Ok(())
    }

    #[test]
    fn explicit_chromosome_order_controls_frag_row_order() -> Result<()> {
        // Arrange:
        // Use the same intentionally non-lexicographic BAM header:
        //   chr2, chr10, chr1
        //
        // But now request an explicit chromosome order:
        //   chr1, chr2
        //
        // `bam-to-frag` documents that fragments are sorted using the chromosome order in
        // `--chromosomes`. With one fragment on each selected chromosome, the output row order
        // must therefore be:
        //   chr1 first, chr2 second
        //
        // chr10 is present in the BAM but not requested, so it must not appear at all.
        let work = tempdir().context("tempdir")?;
        let bam_path = work.path().join("explicit_order.bam");
        let out_dir = work.path().join("out_explicit_order");
        fs::create_dir_all(&out_dir)?;

        let mut hdr = Header::new();
        hdr.push_record(
            bam::header::HeaderRecord::new(b"HD")
                .push_tag(b"VN", &"1.6")
                .push_tag(b"SO", &"coordinate"),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr2")
                .push_tag(b"LN", &1000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr10")
                .push_tag(b"LN", &1000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &1000),
        );

        let mut writer = Writer::from_path(&bam_path, &hdr, Format::Bam).context("create BAM")?;
        let header_view = HeaderView::from_header(&hdr);
        let tid_chr2 = header_view.tid(b"chr2").expect("chr2 present") as i32;
        let tid_chr10 = header_view.tid(b"chr10").expect("chr10 present") as i32;
        let tid_chr1 = header_view.tid(b"chr1").expect("chr1 present") as i32;
        let cigar = vec![Cigar::Match(40)];
        let seq = b"ACGTN".repeat(8);
        let qual = vec![30u8; 40];

        let records = vec![
            make_rec(
                b"chr2_pair",
                tid_chr2,
                10,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr2,
                50,
                true,
            ),
            make_rec(
                b"chr2_pair",
                tid_chr2,
                50,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr2,
                10,
                false,
            ),
            make_rec(
                b"chr10_pair",
                tid_chr10,
                20,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr10,
                60,
                true,
            ),
            make_rec(
                b"chr10_pair",
                tid_chr10,
                60,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr10,
                20,
                false,
            ),
            make_rec(
                b"chr1_pair",
                tid_chr1,
                30,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr1,
                70,
                true,
            ),
            make_rec(
                b"chr1_pair",
                tid_chr1,
                70,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr1,
                30,
                false,
            ),
        ];
        for record in records {
            writer.write(&record)?;
        }
        drop(writer);
        index::build(bam_path.to_str().unwrap(), None, index::Type::Bai, 1).context("build BAI")?;

        let mut cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam_path,
                output_dir: out_dir.clone(),
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_output_prefix("explicit_order");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);

        // Act
        run_bam_to_frag_for_test(&cfg)?;

        // Assert
        let rows = parse_frag_rows(&read_frag_gz(&out_dir.join("explicit_order.frag.tsv.gz"))?);
        assert_eq!(
            rows,
            vec![
                ("chr1".to_string(), 30, 110, 60, '+'),
                ("chr2".to_string(), 10, 90, 60, '+'),
            ]
        );

        Ok(())
    }

    #[test]
    fn default_min_mapq_matches_explicit_zero_and_differs_from_explicit_thirty() -> Result<()> {
        // Arrange:
        // Reuse the small two-chromosome BAM fixture:
        // - pair A has min MAPQ 60
        // - pair B has min MAPQ 0
        // - pair C has min MAPQ 30
        //
        // `bam-to-frag` intentionally defaults to `min_mapq = 0`, so:
        // - default config must emit all three fragments
        // - explicit `min_mapq = 0` must match exactly
        // - explicit `min_mapq = 30` must drop only pair B
        let work = tempdir().context("tempdir")?;
        let bam_path = work.path().join("test.bam");
        let out_default = work.path().join("out_default");
        let out_zero = work.path().join("out_zero");
        let out_thirty = work.path().join("out_thirty");
        fs::create_dir_all(&out_default)?;
        fs::create_dir_all(&out_zero)?;
        fs::create_dir_all(&out_thirty)?;
        write_test_bam(&bam_path)?;
        index::build(bam_path.to_str().unwrap(), None, index::Type::Bai, 1).context("build BAI")?;

        let chromosomes = fixed_chromosome_args();
        let mut default_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam_path.clone(),
                output_dir: out_default.clone(),
                n_threads: 1,
            },
            chromosomes.clone(),
        );
        default_cfg.set_output_prefix("default");

        let mut explicit_zero_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam_path.clone(),
                output_dir: out_zero.clone(),
                n_threads: 1,
            },
            chromosomes.clone(),
        );
        explicit_zero_cfg.set_output_prefix("explicit_zero");
        explicit_zero_cfg.set_min_mapq(0);

        let mut explicit_thirty_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam_path,
                output_dir: out_thirty.clone(),
                n_threads: 1,
            },
            chromosomes,
        );
        explicit_thirty_cfg.set_output_prefix("explicit_thirty");
        explicit_thirty_cfg.set_min_mapq(30);

        // Act
        run_bam_to_frag_for_test(&default_cfg)?;
        run_bam_to_frag_for_test(&explicit_zero_cfg)?;
        run_bam_to_frag_for_test(&explicit_thirty_cfg)?;

        // Assert
        let default_rows = read_frag_gz(&out_default.join("default.frag.tsv.gz"))?;
        let explicit_zero_rows = read_frag_gz(&out_zero.join("explicit_zero.frag.tsv.gz"))?;
        let explicit_thirty_rows = read_frag_gz(&out_thirty.join("explicit_thirty.frag.tsv.gz"))?;

        let default_parsed = parse_frag_rows(&default_rows);
        let explicit_zero_parsed = parse_frag_rows(&explicit_zero_rows);
        let explicit_thirty_parsed = parse_frag_rows(&explicit_thirty_rows);

        assert_eq!(default_parsed, explicit_zero_parsed);
        assert_eq!(default_parsed.len(), 3);
        assert_eq!(explicit_thirty_parsed.len(), 2);
        assert!(
            !explicit_thirty_parsed
                .iter()
                .any(|row| row.0 == "chr1" && row.1 == 10003),
            "pair B is the only fragment below MAPQ 30 and should be removed"
        );

        Ok(())
    }

    #[test]
    fn bam_to_frag_bed_handles_three_chromosomes() -> Result<()> {
        let work = tempdir().context("tempdir")?;
        let bam_path = work.path().join("three_chr_bed.bam");
        let out_dir = work.path().join("out_bed");
        fs::create_dir_all(&out_dir)?;

        write_three_chrom_window_bam(&bam_path)?;

        index::build(bam_path.to_str().unwrap(), None, index::Type::Bai, 1).context("build BAI")?;

        let bed_path = work.path().join("three_chr_windows.bed");
        fs::write(&bed_path, "chr1\t0\t60\nchr2\t0\t80\nchr3\t40\t100\n")?;

        let ioc = IOCArgs {
            bam: bam_path.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let mut cfg = BamToFragConfig::new(ioc, three_chromosome_args());
        cfg.set_output_prefix("three_chr_bed");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_by_bed(Some(bed_path));

        run_bam_to_frag_for_test(&cfg)?;

        let frag_path = out_dir.join("three_chr_bed.frag.tsv.gz");
        let mut parsed = parse_frag_rows(&read_frag_gz(&frag_path)?);
        parsed.sort();

        assert_eq!(
            parsed,
            vec![
                ("chr1".to_string(), 10, 130, 60, '+'),
                ("chr2".to_string(), 30, 150, 60, '+'),
                ("chr3".to_string(), 50, 170, 60, '+'),
            ],
            "BED mode should keep only the fragments overlapping the per-chromosome windows"
        );

        Ok(())
    }

    #[test]
    fn bam_to_frag_gc_file_neutralize_invalid_writes_weight_one_and_keeps_row() -> Result<()> {
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out_gc_fallback");
        std::fs::create_dir_all(&out_dir)?;

        let gc_path = out_dir.join("gc_pkg.zarr");
        build_gc_package(&gc_path, 26, twobit_contig_footprint(&ref_twobit.path)?)?;

        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("gc_fallback");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: true,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            // The GC loader requires min_fragment_length > 2 * end_offset = 52.
            frag.min_fragment_length = 53;
            frag.max_fragment_length = 200;
        }

        // Manual expectations:
        // - The fixture contains one paired fragment spanning [20, 80), length 60.
        // - The GC package uses end_offset=26, leaving only 8 effective bases.
        // - The GC corrector requires at least 10 A/C/G/T bases, so the lookup fails.
        // - With `neutralize_invalid_gc=true`, `bam-to-frag` keeps the fragment and writes `gc_weight=1.0`,
        //   increments `gc_failed_fragments`, and still emits the GC column in the header.
        let counters = run_bam_to_frag_for_test(&cfg)?.counters;

        assert_eq!(counters.base.counted_fragments, 1);
        assert_eq!(counters.gc_failed_fragments, 1);

        let frag_path = out_dir.join("gc_fallback.frag.tsv.gz");
        let frag_rows = read_frag_gz(&frag_path)?;
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t1"]);

        let header_path = out_dir.join("gc_fallback.frag.header.tsv");
        let header_text = std::fs::read_to_string(&header_path)?;
        assert_eq!(
            header_text,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\n"
        );

        Ok(())
    }

    #[test]
    fn bam_to_frag_gc_file_default_behavior_skips_invalid_fragment() -> Result<()> {
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out_gc_default_skip");
        std::fs::create_dir_all(&out_dir)?;

        let gc_path = out_dir.join("gc_pkg.zarr");
        build_gc_package(&gc_path, 26, twobit_contig_footprint(&ref_twobit.path)?)?;

        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("gc_default_skip");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 53;
            frag.max_fragment_length = 200;
        }

        // Manual expectations:
        // - The one fixture fragment again fails GC lookup because only 8 effective A/C/G/T bases remain.
        // - With the default GC behavior, unusable GC information makes the fragment get skipped.
        // - The converter still writes the header with the GC column because GC correction was requested.
        let counters = run_bam_to_frag_for_test(&cfg)?.counters;

        assert_eq!(counters.base.counted_fragments, 0);
        assert_eq!(counters.gc_failed_fragments, 1);

        let frag_path = out_dir.join("gc_default_skip.frag.tsv.gz");
        let frag_rows = read_frag_gz(&frag_path)?;
        assert!(
            frag_rows.is_empty(),
            "default GC handling should skip the only fragment"
        );

        let header_path = out_dir.join("gc_default_skip.frag.header.tsv");
        let header_text = std::fs::read_to_string(&header_path)?;
        assert_eq!(
            header_text,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\n"
        );

        Ok(())
    }

    #[test]
    fn bam_to_frag_gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range()
    -> Result<()> {
        // Arrange:
        // The fixture contributes one fragment of length 60. We keep the accepted fragment length
        // range at exactly 60, then provide a GC package that only covers 10..=59.
        //
        // Because `bam-to-frag` validates the package before conversion starts, the correct
        // failure is the shared compatibility error rather than a late per-fragment lookup error.
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out_gc_range_error");
        std::fs::create_dir_all(&out_dir)?;

        let gc_path = out_dir.join("gc_pkg_short.zarr");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![10, 59],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
            correction_matrix: array![[1.0_f64]],
        };
        package.write_zarr(&gc_path)?;

        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("gc_range_error");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 60;
            frag.max_fragment_length = 60;
        }

        // Act
        let err = run_bam_to_frag_for_test(&cfg).expect_err("out-of-range GC package should fail");

        // Assert
        let msg = err.to_string();
        assert!(
            msg.contains("fragment length range [60-60] is outside the range covered by the correction package [10-59]"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn gc_file_rejects_package_with_schema_version_mismatch() -> Result<()> {
        // Arrange:
        // Build the smallest valid GC correction package shape, but make the schema version
        // incompatible. `bam-to-frag` should fail while loading the package, before writing any
        // frag rows.
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out_gc_bad_version");
        std::fs::create_dir_all(&out_dir)?;

        let gc_path = out_dir.join("gc_pkg_bad_version.zarr");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION + 1,
            end_offset: 0,
            length_edges: vec![10, 200],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
            correction_matrix: array![[1.0_f64]],
        };
        package.write_zarr(&gc_path)?;

        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("gc_bad_version");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

        // Act
        let err = run_bam_to_frag_for_test(&cfg).expect_err("schema version mismatch should fail");

        // Assert
        let msg = err.to_string();
        assert!(
            msg.contains("GC correction package schema version mismatch"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn gc_file_rejects_directory_with_clear_error() -> Result<()> {
        // Arrange:
        // Point `--gc-file` at a directory without the expected `.zarr` package suffix.
        // The command should reject that immediately during GC package loading.
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out_gc_directory_error");
        std::fs::create_dir_all(&out_dir)?;

        let gc_dir = out_dir.join("gc_package_dir");
        std::fs::create_dir_all(&gc_dir)?;

        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        };
        let chromosomes = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("gc_directory_error");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_dir.clone()),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

        // Act
        let err = run_bam_to_frag_for_test(&cfg)
            .expect_err("directory-valued GC package path should fail");

        // Assert
        let msg = err.to_string();
        assert!(
            msg.contains("GC correction package path must point to a .zarr directory"),
            "unexpected error message: {msg}"
        );
        assert!(
            msg.contains("gc_package_dir"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    #[test]
    fn bam_to_frag_and_bam_to_bam_encode_same_coverage_scaling_weight() -> Result<()> {
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let scaling_path = work.path().join("shared_scaling.tsv");
        std::fs::write(
            &scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t2\n",
        )?;

        // Manual expectations:
        // - The fixture contains one paired fragment spanning [20, 80).
        // - The scaling TSV has one chromosome-wide factor of 2.0.
        // - `bam-to-frag` averages scaling over the full fragment span, which stays 2.0.
        // - `bam-to-bam` writes the same full-fragment scaling as the `cw` tag on both mates.
        // - So the two released transformers should encode the same weight for the same fragment.

        let frag_out_dir = work.path().join("frag_out");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("scaled");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        let mut frag_scale = cfdnalab::run_like_cli::common::ScaleGenomeArgs::default();
        frag_scale.scaling_factors = Some(scaling_path.clone());
        frag_cfg.set_coverage_scaling_factors(frag_scale.scaling_factors);

        run_bam_to_frag_for_test(&frag_cfg)?;
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled.frag.tsv.gz"))?;
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t2"]);

        let bam_out = work.path().join("scaled.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        let mut bam_scale = cfdnalab::run_like_cli::common::ScaleGenomeArgs::default();
        bam_scale.scaling_factors = Some(scaling_path);
        bam_cfg.set_coverage_scaling_factors(bam_scale.scaling_factors);

        run_bam_to_bam_for_test(&bam_cfg)?;
        let mut reader = rust_htslib::bam::Reader::from_path(&bam_out)?;
        let mut cov_tags = Vec::new();
        for record in reader.records() {
            let record = record?;
            match record.aux(b"cw") {
                Ok(Aux::Float(value)) => cov_tags.push(value),
                other => panic!("expected cw float tag on every mate, got {other:?}"),
            }
        }
        assert_eq!(cov_tags, vec![2.0_f32, 2.0_f32]);

        Ok(())
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    #[test]
    fn bam_to_frag_and_bam_to_bam_encode_same_count_scaling_weight() -> Result<()> {
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let scaling_path = work.path().join("shared_count_scaling.tsv");
        std::fs::write(
            &scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t0.5\n",
        )?;

        // Manual expectations:
        // - The fixture contains one paired fragment spanning [20, 80).
        // - The count-scaling TSV has one chromosome-wide factor of 0.5.
        // - `bam-to-frag` averages count scaling over the full fragment span, which stays 0.5.
        // - `bam-to-bam` writes the same full-fragment count scaling as the `nw` tag on both mates.
        let frag_out_dir = work.path().join("frag_out");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("scaled");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_count_scaling_factors(Some(scaling_path.clone()));

        run_bam_to_frag_for_test(&frag_cfg)?;
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled.frag.tsv.gz"))?;
        let frag_header = std::fs::read_to_string(frag_out_dir.join("scaled.frag.header.tsv"))?;
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tcount_scaling_weight\n"
        );
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t0.5"]);

        let bam_out = work.path().join("scaled.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        bam_cfg.set_count_scaling_factors(Some(scaling_path));

        run_bam_to_bam_for_test(&bam_cfg)?;
        let mut reader = rust_htslib::bam::Reader::from_path(&bam_out)?;
        let mut cnt_tags = Vec::new();
        for record in reader.records() {
            let record = record?;
            match record.aux(b"nw") {
                Ok(Aux::Float(value)) => cnt_tags.push(value),
                other => panic!("expected nw float tag on every mate, got {other:?}"),
            }
        }
        assert_eq!(cnt_tags, vec![0.5_f32, 0.5_f32]);

        Ok(())
    }

    #[test]
    fn bam_to_frag_writes_coverage_and_fragment_count_scaling_columns_together() -> Result<()> {
        // Arrange:
        // Use one paired fragment spanning [20, 80) and provide two chromosome-wide scaling files:
        // - coverage scaling       -> 2.0
        // - count-based scaling -> 0.5
        //
        // `bam-to-frag` averages each scaling file over the same full fragment span, so the output
        // row should carry both constants side by side. The header must keep a stable column order:
        //   coverage_scaling_weight, count_scaling_weight
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let coverage_scaling_path = work.path().join("coverage_scaling.tsv");
        let fragment_count_scaling_path = work.path().join("fragment_count_scaling.tsv");
        std::fs::write(
            &coverage_scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t2\n",
        )?;
        std::fs::write(
            &fragment_count_scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t0.5\n",
        )?;

        let frag_out_dir = work.path().join("frag_out");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms);
        frag_cfg.set_output_prefix("scaled");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_coverage_scaling_factors(Some(coverage_scaling_path));
        frag_cfg.set_count_scaling_factors(Some(fragment_count_scaling_path));

        // Act
        run_bam_to_frag_for_test(&frag_cfg)?;

        // Assert
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled.frag.tsv.gz"))?;
        let frag_header = std::fs::read_to_string(frag_out_dir.join("scaled.frag.header.tsv"))?;
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tcoverage_scaling_weight\tcount_scaling_weight\n"
        );
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t2\t0.5"]);

        Ok(())
    }

    #[test]
    fn bam_to_frag_writes_count_scaling_column_when_only_count_scaling_is_configured() -> Result<()>
    {
        // Arrange:
        // One paired fragment spanning [20, 80) and one chromosome-wide count-scaling factor 0.5.
        //
        // Expected:
        // - header contains only `count_scaling_weight`
        // - the one fragment row carries that value
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let scaling_path = work.path().join("count_scaling.tsv");
        std::fs::write(
            &scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t0.5\n",
        )?;

        let frag_out_dir = work.path().join("frag_out");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms);
        frag_cfg.set_output_prefix("scaled");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_count_scaling_factors(Some(scaling_path));

        // Act
        run_bam_to_frag_for_test(&frag_cfg)?;

        // Assert
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled.frag.tsv.gz"))?;
        let frag_header = std::fs::read_to_string(frag_out_dir.join("scaled.frag.header.tsv"))?;
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tcount_scaling_weight\n"
        );
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t0.5"]);

        Ok(())
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    #[test]
    fn bam_to_frag_and_bam_to_bam_emit_combined_gc_scaling_and_length_metadata() -> Result<()> {
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;

        let scaling_path = work.path().join("shared_combined_scaling.tsv");
        std::fs::write(
            &scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t2\n",
        )?;

        let gc_path = work.path().join("combined_gc_pkg.zarr");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![10, 61, 100],
            gc_edges: vec![0, 51, 100],
            length_bin_frequencies: array![1.0_f64, 1.0_f64],
            reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
            correction_matrix: array![[3.0_f64, 1.0_f64], [1.0_f64, 1.0_f64]],
        };
        package.write_zarr(&gc_path)?;

        // Manual expectations:
        // - The fixture contains one paired fragment spanning [20, 80), so fragment length = 60.
        // - `simple_reference_twobit` is "ACGT" repeated; across 60 bases this gives exactly
        //   30 G/C bases, so the integer GC percentage is 50.
        // - The custom package defines:
        //   - length bin [10, 61) containing 60
        //   - GC bin [0, 51) containing 50
        //   - correction weight 3.0 in that cell
        // - The scaling TSV applies factor 2.0 over the full chromosome, so the fragment-average
        //   scaling is also 2.0.
        // - Therefore:
        //   - `bam-to-frag` should emit `chr1 20 80 60 + 3 2`
        //   - `bam-to-bam` should emit `GC=3.0`, `cw=2.0`, and `fl=60` on both mates.

        let frag_out_dir = work.path().join("frag_combined_out");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("combined");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            neutralize_invalid_gc: false,
        });
        frag_cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        let mut frag_scale = cfdnalab::run_like_cli::common::ScaleGenomeArgs::default();
        frag_scale.scaling_factors = Some(scaling_path.clone());
        frag_cfg.set_coverage_scaling_factors(frag_scale.scaling_factors);
        {
            let fragment_lengths = frag_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 100;
        }

        let bam_out = work.path().join("combined_tags.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        bam_cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        bam_cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        let mut bam_scale = cfdnalab::run_like_cli::common::ScaleGenomeArgs::default();
        bam_scale.scaling_factors = Some(scaling_path);
        bam_cfg.set_coverage_scaling_factors(bam_scale.scaling_factors);
        {
            let fragment_lengths = bam_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 100;
        }

        let frag_counters = run_bam_to_frag_for_test(&frag_cfg)?.counters;
        let bam_counters = run_bam_to_bam_for_test(&bam_cfg)?.counters;

        let frag_rows = read_frag_gz(&frag_out_dir.join("combined.frag.tsv.gz"))?;
        let frag_header = std::fs::read_to_string(frag_out_dir.join("combined.frag.header.tsv"))?;

        let mut reader = rust_htslib::bam::Reader::from_path(&bam_out)?;
        let mut observed_tags = Vec::new();
        for record in reader.records() {
            let record = record?;
            let gc = match record.aux(b"GC") {
                Ok(Aux::Float(value)) => value,
                other => panic!("expected GC float tag on every mate, got {other:?}"),
            };
            let scaling = match record.aux(b"cw") {
                Ok(Aux::Float(value)) => value,
                other => panic!("expected cw float tag on every mate, got {other:?}"),
            };
            let flen = match record.aux(b"fl") {
                Ok(Aux::U32(value)) => value,
                other => panic!("expected fl u32 tag on every mate, got {other:?}"),
            };
            observed_tags.push((gc, scaling, flen));
        }

        assert_eq!(frag_counters.base.counted_fragments, 1);
        assert_eq!(bam_counters.base.counted_fragments, 1);
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\tcoverage_scaling_weight\n"
        );
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t3\t2"]);

        assert_eq!(observed_tags.len(), 2);
        for (mate_idx, (gc, scaling, flen)) in observed_tags.into_iter().enumerate() {
            assert!(
                (gc - 3.0).abs() < 1e-6,
                "mate {mate_idx} GC tag: expected 3.0, got {gc}"
            );
            assert!(
                (scaling - 2.0).abs() < 1e-6,
                "mate {mate_idx} cw tag: expected 2.0, got {scaling}"
            );
            assert_eq!(flen, 60, "mate {mate_idx} fl tag: expected 60");
        }

        Ok(())
    }

    #[cfg(all(feature = "cmd_bam_to_bam", feature = "cmd_coverage_weights"))]
    #[test]
    fn real_coverage_weights_tsv_has_same_effect_in_bam_to_frag_and_bam_to_bam() -> Result<()> {
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;

        let weights_out_dir = work.path().join("weights_out");
        std::fs::create_dir_all(&weights_out_dir)?;
        let mut weights_cfg = CoverageWeightsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: weights_out_dir.clone(),
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
        );
        weights_cfg.set_output_prefix("coverage".to_string());
        weights_cfg.set_bin_size(40);
        weights_cfg.set_stride(20);
        weights_cfg.set_min_mapq(0);
        weights_cfg.set_require_proper_pair(false);
        {
            let frag = weights_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }
        run_coverage_weights(&weights_cfg)?;

        let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

        // Manual expectations:
        // - `coverage-weights` on the simple fixture yields stride-bin scaling factors:
        //   [0,20): 37/20
        //   [20,40): 37/45
        //   [40,60): 37/60
        //   [60,80): 37/45
        //   [80,100): 37/15
        //   remaining bins: 0
        // - The fragment written by both transformers spans [20, 80), so the full-fragment
        //   average scaling uses the three covered stride bins [20,40), [40,60), [60,80):
        //   mean = ((37/45) + (37/60) + (37/45)) / 3
        //        = (148/180 + 111/180 + 148/180) / 3
        //        = (407/180) / 3
        //        = 407/540
        // - Both transformers should therefore encode the same weight 407/540.
        let expected_weight = 407.0_f64 / 540.0_f64;

        let frag_out_dir = work.path().join("frag_real_weights");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("scaled_real");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        let mut frag_scale = cfdnalab::run_like_cli::common::ScaleGenomeArgs::default();
        frag_scale.scaling_factors = Some(scaling_path.clone());
        frag_cfg.set_coverage_scaling_factors(frag_scale.scaling_factors);

        run_bam_to_frag_for_test(&frag_cfg)?;
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled_real.frag.tsv.gz"))?;
        assert_eq!(frag_rows.len(), 1);
        let frag_parts: Vec<_> = frag_rows[0].split('\t').collect();
        assert_eq!(frag_parts.len(), 6, "expected scaling-weight frag output");
        let frag_weight: f64 = frag_parts[5].parse()?;
        assert!(
            (frag_weight - expected_weight).abs() < 1e-6,
            "bam-to-frag scaling weight: expected {expected_weight}, got {frag_weight}"
        );

        let bam_out = work.path().join("scaled_real.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        let mut bam_scale = cfdnalab::run_like_cli::common::ScaleGenomeArgs::default();
        bam_scale.scaling_factors = Some(scaling_path);
        bam_cfg.set_coverage_scaling_factors(bam_scale.scaling_factors);

        run_bam_to_bam_for_test(&bam_cfg)?;
        let mut reader = rust_htslib::bam::Reader::from_path(&bam_out)?;
        let mut cov_tags = Vec::new();
        for record in reader.records() {
            let record = record?;
            match record.aux(b"cw") {
                Ok(Aux::Float(value)) => cov_tags.push(value as f64),
                other => panic!("expected cw float tag on every mate, got {other:?}"),
            }
        }
        assert_eq!(cov_tags.len(), 2);
        for (mate_idx, value) in cov_tags.iter().enumerate() {
            assert!(
                (*value - expected_weight).abs() < 1e-6,
                "bam-to-bam cw tag for mate {mate_idx}: expected {expected_weight}, got {value}"
            );
        }

        Ok(())
    }

    #[cfg(all(feature = "cmd_bam_to_bam", feature = "cmd_fragment_count_weights"))]
    #[test]
    fn real_fragment_count_weights_tsv_has_same_effect_in_bam_to_frag_and_bam_to_bam() -> Result<()>
    {
        // Arrange:
        // Reuse the mixed-length fixture:
        // - short fragment [0,20)
        // - long  fragment [20,80)
        //
        // With bin_size == stride == 20, fragment-count-weights yields:
        // - short fragment count scaling = 0.5
        // - long  fragment count scaling = 1.5
        //
        // Both transformers average over the full fragment span, so they should encode the same
        // fragment-level count scaling for each fragment.
        //
        // We compare them at the fragment level:
        // - `bam-to-frag` writes one row per fragment
        // - `bam-to-bam` writes one tagged record per mate, so we must regroup mates by qname
        //
        // The hand-derived target values are:
        // - short -> 0.5
        // - long  -> 1.5
        let bam = mixed_length_fragment_bam("bam_to_frag_real_count_weights")?;
        let work = tempdir().context("tempdir")?;

        let weights_out_dir = work.path().join("weights_out");
        std::fs::create_dir_all(&weights_out_dir)?;
        let mut weights_cfg = FragmentCountWeightsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: weights_out_dir.clone(),
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
        );
        weights_cfg.set_output_prefix("counts".to_string());
        weights_cfg.set_bin_size(20);
        weights_cfg.set_stride(20);
        weights_cfg.set_min_mapq(0);
        weights_cfg.set_require_proper_pair(false);
        {
            let frag = weights_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 100;
        }
        run_fragment_count_weights(&weights_cfg)?;

        let scaling_path = weights_out_dir.join("counts.fragment_counts.scaling_factors.tsv");

        let frag_out_dir = work.path().join("frag_real_count_weights");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("scaled_real");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_count_scaling_factors(Some(scaling_path.clone()));
        {
            let fragment_lengths = frag_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 100;
        }

        run_bam_to_frag_for_test(&frag_cfg)?;
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled_real.frag.tsv.gz"))?;
        let mut parsed_rows: Vec<(u64, f64)> = frag_rows
            .iter()
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                assert_eq!(columns.len(), 6, "Bad line: {line}");
                Ok::<_, anyhow::Error>((columns[1].parse()?, columns[5].parse()?))
            })
            .collect::<Result<_>>()?;
        parsed_rows.sort_by(|left, right| left.0.cmp(&right.0));

        let bam_out = work.path().join("scaled_real.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        bam_cfg.set_count_scaling_factors(Some(scaling_path));
        {
            let fragment_lengths = bam_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 100;
        }

        run_bam_to_bam_for_test(&bam_cfg)?;
        let bam_fragment_cnts = read_fragment_level_cnts_from_bam(&bam_out)?;

        assert_eq!(parsed_rows.len(), 2, "expected exactly two frag rows");
        assert_eq!(
            bam_fragment_cnts.len(),
            2,
            "expected exactly two grouped BAM fragments"
        );

        // The two converters should agree on the same fragment-level value, but we still allow
        // a small tolerance because the shared scaling TSV is generated through `f32`-backed
        // smoothing/normalization before either converter consumes it.
        for ((frag_start, frag_cnt), (bam_start, bam_cnt)) in
            parsed_rows.iter().zip(bam_fragment_cnts.iter())
        {
            assert_eq!(
                frag_start, bam_start,
                "fragment starts should match between frag and BAM outputs"
            );
            assert!(
                (*frag_cnt - *bam_cnt).abs() <= 1e-6,
                "fragment start {frag_start}: bam-to-frag and bam-to-bam should encode the same fragment-level nw weight, got {frag_cnt} vs {bam_cnt}"
            );
        }

        let expected = [(0_u64, 0.5_f64), (20_u64, 1.5_f64)];
        // These are the hand-derived ideal values. The real scaling factors are produced through
        // the shared `fragment-count-weights` pipeline, whose stride-bin state and normalization
        // are currently stored as `f32`, so allow a small tolerance around the exact ideals.
        for ((start, value), (_, expected_value)) in parsed_rows.iter().zip(expected.iter()) {
            assert!(
                (*value - *expected_value).abs() <= 1e-6,
                "fragment start {start}: expected count scaling near {expected_value}, got {value}"
            );
        }

        Ok(())
    }

    #[cfg(feature = "cmd_coverage_weights")]
    #[test]
    fn real_multi_chromosome_coverage_weights_tsv_is_applied_per_chromosome_in_bam_to_frag()
    -> Result<()> {
        // Arrange:
        // Reuse the same two-chromosome fixture and hand derivation as the command-level
        // `coverage-weights` shared-global-mean test:
        //
        // chr1 fragment:
        // - span [20, 80), length 60
        // - avg-overlap profile under bin_size=40, stride=20:
        //   [1/3, 3/4, 1, 3/4, 1/4, 0, ...]
        //
        // chr2 fragment:
        // - span [20, 40), length 20
        // - avg-overlap profile:
        //   [1/3, 1/2, 1/4, 0, ...]
        //
        // Shared non-zero global mean:
        //   chr1 sum = 37/12
        //   chr2 sum = 13/12
        //   total    = 25/6 across 8 non-zero bins
        //   mean     = (25/6) / 8 = 25/48
        //
        // Inverted per-bin scaling factors are therefore:
        //   chr1 [20,40) = (25/48) / (3/4) = 25/36
        //   chr1 [40,60) = (25/48) / 1     = 25/48
        //   chr1 [60,80) = (25/48) / (3/4) = 25/36
        //   chr2 [20,40) = (25/48) / (1/2) = 25/24
        //
        // `bam-to-frag` averages scaling over the full fragment span:
        //   chr1 weight = ((25/36) + (25/48) + (25/36)) / 3 = 275/432
        //   chr2 weight = 25/24
        let mut chr2_fragment = PairedFragmentSpec::new(0, 20, 20, 10).build()?;
        chr2_fragment.forward.tid = 1;
        chr2_fragment.reverse.tid = 1;
        chr2_fragment.forward.mate_tid = Some(1);
        chr2_fragment.reverse.mate_tid = Some(1);

        let bam = bam_from_fragments(
            "bam_to_frag_real_multi_chr_scaling",
            vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
            vec![
                PairedFragmentSpec::new(0, 20, 60, 20).build()?,
                chr2_fragment,
            ],
            Vec::new(),
        )?;
        let work = tempdir().context("tempdir")?;

        let weights_out_dir = work.path().join("weights_out");
        std::fs::create_dir_all(&weights_out_dir)?;
        let mut weights_cfg = CoverageWeightsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: weights_out_dir.clone(),
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
        );
        weights_cfg.set_output_prefix("coverage".to_string());
        weights_cfg.set_bin_size(40);
        weights_cfg.set_stride(20);
        weights_cfg.set_min_mapq(0);
        weights_cfg.set_require_proper_pair(false);
        {
            let frag = weights_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }
        run_coverage_weights(&weights_cfg)?;

        let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");
        let frag_out_dir = work.path().join("frag_real_multi_chr_scaling");
        std::fs::create_dir_all(&frag_out_dir)?;
        let mut frag_cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: frag_out_dir.clone(),
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
        );
        frag_cfg.set_output_prefix("scaled_multi_chr");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_coverage_scaling_factors(Some(scaling_path));
        {
            let frag = frag_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        // Act
        let counters = run_bam_to_frag_for_test(&frag_cfg)?.counters;
        let frag_rows = read_frag_gz(&frag_out_dir.join("scaled_multi_chr.frag.tsv.gz"))?;
        let frag_header =
            std::fs::read_to_string(frag_out_dir.join("scaled_multi_chr.frag.header.tsv"))?;

        // Assert
        assert_eq!(counters.base.counted_fragments, 2);
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tcoverage_scaling_weight\n"
        );
        assert_eq!(frag_rows.len(), 2);

        let mut parsed: Vec<(String, u64, u64, u8, char, f64)> = frag_rows
            .iter()
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                assert_eq!(columns.len(), 6, "Bad line: {line}");
                Ok::<_, anyhow::Error>((
                    columns[0].to_string(),
                    columns[1].parse()?,
                    columns[2].parse()?,
                    columns[3].parse()?,
                    columns[4].chars().next().unwrap(),
                    columns[5].parse()?,
                ))
            })
            .collect::<Result<_>>()?;
        parsed.sort_by(|left, right| left.0.cmp(&right.0));

        let expected_chr1 = 275.0_f64 / 432.0_f64;
        let expected_chr2 = 25.0_f64 / 24.0_f64;

        assert_eq!(parsed[0].0, "chr1");
        assert_eq!(parsed[0].1, 20);
        assert_eq!(parsed[0].2, 80);
        assert_eq!(parsed[0].3, 60);
        assert_eq!(parsed[0].4, '+');
        assert!(
            (parsed[0].5 - expected_chr1).abs() <= 1e-6,
            "chr1 scaling: expected {expected_chr1}, got {}",
            parsed[0].5
        );

        assert_eq!(parsed[1].0, "chr2");
        assert_eq!(parsed[1].1, 20);
        assert_eq!(parsed[1].2, 40);
        assert_eq!(parsed[1].3, 60);
        assert_eq!(parsed[1].4, '+');
        assert!(
            (parsed[1].5 - expected_chr2).abs() <= 1e-6,
            "chr2 scaling: expected {expected_chr2}, got {}",
            parsed[1].5
        );

        Ok(())
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_is_neutral_in_bam_to_frag_and_bam_to_bam() -> Result<()>
    {
        let bam = single_contig_inward_pair_bam()?;
        let reference =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let gc_path = build_command_produced_gc_correction_package_for_length(
            &bam.bam,
            &reference.path,
            work.path(),
            60,
        )?;

        // Manual expectations:
        // - `simple_inward_bam` contains one fragment [20, 80), length 60.
        // - `simple_reference_twobit` is "ACGT" repeated, so over 60 bp:
        //     GC count = 30
        //     GC percentage = 50
        // - `ref-gc-bias` is run for exactly that one fragment length, and `gc-bias` is run on
        //   exactly that same fragment type over the same repeated reference.
        // - All reference mass and all sample mass therefore land in one GC-by-length cell.
        // - The resulting correction is neutral for that cell: weight 1.0.
        // - So both released converters must preserve the fragment unchanged apart from encoding
        //   the explicit neutral GC weight:
        //     `bam-to-frag`: row `chr1 20 80 60 + 1`
        //     `bam-to-bam`: both mates tagged with `GC=1.0` and `fl=60`

        let frag_out_dir = work.path().join("frag_real_gc");
        std::fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("real_gc");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            neutralize_invalid_gc: false,
        });
        frag_cfg.set_ref_2bit(Some(reference.path.clone()));
        {
            let fragment_lengths = frag_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 60;
            fragment_lengths.max_fragment_length = 60;
        }

        let bam_out = work.path().join("real_gc_tags.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        bam_cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        bam_cfg.set_ref_2bit(Some(reference.path.clone()));
        {
            let fragment_lengths = bam_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 60;
            fragment_lengths.max_fragment_length = 60;
        }

        let frag_counters = run_bam_to_frag_for_test(&frag_cfg)?.counters;
        let bam_counters = run_bam_to_bam_for_test(&bam_cfg)?.counters;

        let frag_rows = read_frag_gz(&frag_out_dir.join("real_gc.frag.tsv.gz"))?;
        let frag_header = std::fs::read_to_string(frag_out_dir.join("real_gc.frag.header.tsv"))?;

        let mut reader = rust_htslib::bam::Reader::from_path(&bam_out)?;
        let mut observed_gc_tags = Vec::new();
        let mut observed_flen_tags = Vec::new();
        for record in reader.records() {
            let record = record?;
            let gc = match record.aux(b"GC") {
                Ok(Aux::Float(value)) => value,
                other => panic!("expected GC float tag on every mate, got {other:?}"),
            };
            let flen = match record.aux(b"fl") {
                Ok(Aux::U32(value)) => value,
                other => panic!("expected fl u32 tag on every mate, got {other:?}"),
            };
            observed_gc_tags.push(gc);
            observed_flen_tags.push(flen);
        }

        assert_eq!(frag_counters.base.counted_fragments, 1);
        assert_eq!(bam_counters.base.counted_fragments, 1);
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\n"
        );
        assert_eq!(frag_rows, vec!["chr1\t20\t80\t60\t+\t1"]);
        assert_eq!(observed_gc_tags, vec![1.0_f32, 1.0_f32]);
        assert_eq!(observed_flen_tags, vec![60_u32, 60_u32]);

        Ok(())
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_changes_bam_to_frag_and_bam_to_bam_in_expected_direction()
    -> Result<()> {
        // Arrange:
        // Use the same real non-neutral producer workflow as the corresponding `gc-bias` test:
        // - Reference: chr1[0,100) all A, chr1[100,200) all C
        // - Reference windows: [0,91) and [100,191)
        // - with the `ref-gc-bias` fit rule those count starts 0..=81 and 100..=181 only, so the
        //   reference-side support is still balanced between pure-A and pure-C starts
        // - Sample BAM: one A-only 10 bp fragment and nine C-only 10 bp fragments
        //
        // The real produced GC package is hand-derived as:
        // - GC%=0   -> weight 5.0
        // - GC%=100 -> weight 5/9
        //
        // So both released converters must encode:
        // - one fragment row / mate pair with weight 5.0 for [10,20)
        // - nine fragment rows / mate pairs with weight 5/9 for [110,120) .. [190,200)
        let reference = twobit_from_sequences(
            "bam_to_frag_real_non_neutral_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100), "C".repeat(100)),
            )],
        )?;
        let starts = [10_i64, 110, 120, 130, 140, 150, 160, 170, 180, 190];
        let fragments = starts
            .into_iter()
            .map(|start| PairedFragmentSpec::new(0, start, 10, 5).build())
            .collect::<Result<Vec<_>>>()?;
        let bam = bam_from_fragments(
            "bam_to_frag_real_non_neutral_bam",
            vec![("chr1".to_string(), 200)],
            fragments,
            Vec::new(),
        )?;
        let work = tempdir().context("tempdir")?;
        let gc_path = build_command_produced_gc_correction_package_from_reference_windows(
            &bam.bam,
            &reference.path,
            work.path(),
            10,
            "chr1\t0\t91\nchr1\t100\t191\n",
            // Chromosome length 200 and fragment length 10 give:
            //   200 - 10 + 1 = 191 valid starts.
            191,
        )?;

        let frag_out_dir = work.path().join("frag_real_non_neutral");
        fs::create_dir_all(&frag_out_dir)?;
        let frag_ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: frag_out_dir.clone(),
            n_threads: 1,
        };
        let chroms = ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        };
        let mut frag_cfg = BamToFragConfig::new(frag_ioc, chroms.clone());
        frag_cfg.set_output_prefix("real_gc");
        frag_cfg.set_min_mapq(0);
        frag_cfg.set_require_proper_pair(false);
        frag_cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            neutralize_invalid_gc: false,
        });
        frag_cfg.set_ref_2bit(Some(reference.path.clone()));
        {
            let fragment_lengths = frag_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 10;
        }

        let bam_out = work.path().join("real_gc_tags.bam");
        let mut bam_cfg = BamToBamConfig::new(bam.bam.clone(), bam_out.clone(), chroms);
        bam_cfg.set_min_mapq(0);
        bam_cfg.set_require_proper_pair(false);
        bam_cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        bam_cfg.set_ref_2bit(Some(reference.path.clone()));
        {
            let fragment_lengths = bam_cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 10;
        }

        // Act
        let frag_counters = run_bam_to_frag_for_test(&frag_cfg)?.counters;
        let bam_counters = run_bam_to_bam_for_test(&bam_cfg)?.counters;

        // Assert
        let frag_rows = read_frag_gz(&frag_out_dir.join("real_gc.frag.tsv.gz"))?;
        let frag_header = fs::read_to_string(frag_out_dir.join("real_gc.frag.header.tsv"))?;
        let expected_rows = vec![
            "chr1\t10\t20\t60\t+\t5".to_string(),
            "chr1\t110\t120\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t120\t130\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t130\t140\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t140\t150\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t150\t160\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t160\t170\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t170\t180\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t180\t190\t60\t+\t0.5555555555555556".to_string(),
            "chr1\t190\t200\t60\t+\t0.5555555555555556".to_string(),
        ];

        let mut reader = rust_htslib::bam::Reader::from_path(&bam_out)?;
        let mut observed_gc_tags = Vec::new();
        let mut observed_flen_tags = Vec::new();
        for record in reader.records() {
            let record = record?;
            let gc = match record.aux(b"GC") {
                Ok(Aux::Float(value)) => value,
                other => panic!("expected GC float tag on every mate, got {other:?}"),
            };
            let flen = match record.aux(b"fl") {
                Ok(Aux::U32(value)) => value,
                other => panic!("expected fl u32 tag on every mate, got {other:?}"),
            };
            observed_gc_tags.push(gc);
            observed_flen_tags.push(flen);
        }

        assert_eq!(frag_counters.base.counted_fragments, 10);
        assert_eq!(bam_counters.base.counted_fragments, 10);
        assert_eq!(
            frag_header,
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\n"
        );
        assert_eq!(frag_rows, expected_rows);
        assert_eq!(observed_gc_tags.len(), 20);
        assert_eq!(observed_flen_tags, vec![10_u32; 20]);
        assert_eq!(observed_gc_tags[0], 5.0_f32);
        assert_eq!(observed_gc_tags[1], 5.0_f32);
        for value in observed_gc_tags.iter().skip(2) {
            assert!(
                (*value as f64 - (5.0 / 9.0)).abs() <= 1e-6,
                "expected GC tag 5/9 on C-only fragments, got {value}"
            );
        }

        Ok(())
    }

    #[test]
    fn scaling_tsv_must_cover_requested_chromosome_end_in_bam_to_frag() -> Result<()> {
        // Arrange:
        // `single_contig_inward_pair_bam()` uses chr1 length 200.
        // A scaling TSV that stops at 100 is malformed for this requested chromosome even though
        // the counted fragment itself lies inside the provided region.
        //
        // The command should therefore fail while loading scaling factors, before writing rows.
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out");
        fs::create_dir_all(&out_dir)?;
        let scaling_path = work.path().join("truncated_scaling.tsv");
        fs::write(
            &scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t100\t2.0\n",
        )?;

        let mut cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir,
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_coverage_scaling_factors(Some(scaling_path));
        cfg.set_min_mapq(0);

        // Act
        let err = run_bam_to_frag_for_test(&cfg).expect_err("truncated scaling TSV should fail");

        // Assert:
        // `bam-to-frag` also wraps the shared loader with `load scaling factors`, so inspect
        // the full error chain to reach the actual artifact-contract failure.
        let msg = format!("{err:#}");
        assert!(
            msg.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn count_scaling_tsv_must_cover_requested_chromosome_end_in_bam_to_frag() -> Result<()> {
        // Arrange:
        // Mirror the coverage-side malformed scaling regression through the separate
        // `--count-scaling-factors` path. The loader contract is identical: requested
        // chromosomes must be fully covered by contiguous bins.
        let bam = single_contig_inward_pair_bam()?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out");
        fs::create_dir_all(&out_dir)?;
        let scaling_path = work.path().join("truncated_count_scaling.tsv");
        fs::write(
            &scaling_path,
            "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t100\t2.0\n",
        )?;

        let mut cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir,
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_count_scaling_factors(Some(scaling_path));
        cfg.set_min_mapq(0);

        // Act
        let err =
            run_bam_to_frag_for_test(&cfg).expect_err("truncated count scaling TSV should fail");

        // Assert
        let msg = format!("{err:#}");
        assert!(
            msg.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn count_scaling_tsv_with_uncorrected_metadata_rejects_gc_corrected_bam_to_frag_run()
    -> Result<()> {
        // Arrange:
        // `bam-to-frag` also loads count scaling through the shared scaling-map path. When the
        // run applies `--gc-file`, an explicitly uncorrected count-scaling TSV is a known
        // incompatibility that should fail before any FRAG rows are written.
        let bam = single_contig_inward_pair_bam()?;
        let ref_twobit =
            twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
        let work = tempdir().context("tempdir")?;
        let out_dir = work.path().join("out_count_gc_mismatch");
        fs::create_dir_all(&out_dir)?;
        let scaling_path = out_dir.join("uncorrected_count_scaling.tsv");
        let gc_path = out_dir.join("gc_pkg.zarr");
        build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;
        fs::write(
            &scaling_path,
            "# gc_mode=uncorrected\nchromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t1.0\n",
        )?;

        let mut cfg = BamToFragConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir,
                n_threads: 1,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_output_prefix("count_gc_mismatch");
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_count_scaling_factors(Some(scaling_path));
        cfg.set_gc(ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

        // Act
        let err = run_bam_to_frag_for_test(&cfg)
            .expect_err("uncorrected count scaling must fail in a GC-corrected bam-to-frag run");

        // Assert
        let msg = format!("{err:#}");
        assert!(
            msg.contains("no GC correction"),
            "unexpected error message: {msg}"
        );
        assert!(
            msg.contains("via --gc-file"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    fn build_gc_package(
        path: &Path,
        end_offset: u64,
        reference_contig_footprint: Vec<ContigFootprintEntry>,
    ) -> Result<()> {
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset,
            length_edges: vec![10, 60, 200],
            gc_edges: vec![0, 50, 101],
            length_bin_frequencies: array![1.0_f64, 3.0_f64],
            reference_contig_footprint,
            correction_matrix: array![[1.0_f64, 1.0_f64], [2.0_f64, 10.0_f64]],
        };
        package.write_zarr(path)?;
        Ok(())
    }

    fn read_frag_gz(path: &Path) -> Result<Vec<String>> {
        let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mut gz = GzDecoder::new(f);
        let mut s = String::new();
        gz.read_to_string(&mut s)?;
        let lines = s.lines().map(|l| l.to_string()).collect();
        Ok(lines)
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    fn read_fragment_level_cnts_from_bam(path: &Path) -> Result<Vec<(u64, f64)>> {
        let mut reader = rust_htslib::bam::Reader::from_path(path)?;
        let mut per_qname: HashMap<String, Vec<(u64, u32, f64)>> = HashMap::new();

        for record in reader.records() {
            let record = record?;
            let qname = std::str::from_utf8(record.qname())
                .expect("synthetic BAM qnames should be valid UTF-8")
                .to_string();
            let start = record.pos() as u64;
            let flen = match record.aux(b"fl") {
                Ok(Aux::U32(value)) => value,
                other => panic!("expected fl u32 tag on every mate, got {other:?}"),
            };
            let cnt = match record.aux(b"nw") {
                Ok(Aux::Float(value)) => value as f64,
                other => panic!("expected nw float tag on every mate, got {other:?}"),
            };
            per_qname.entry(qname).or_default().push((start, flen, cnt));
        }

        let mut grouped = Vec::with_capacity(per_qname.len());
        for (qname, entries) in per_qname {
            assert_eq!(
                entries.len(),
                2,
                "expected exactly two mate records per fragment for qname {qname}"
            );

            let fragment_start = entries.iter().map(|(start, _, _)| *start).min().unwrap();
            let expected_flen = entries[0].1;
            let expected_cnt = entries[0].2;

            for (_, flen, cnt) in &entries[1..] {
                assert_eq!(
                    *flen, expected_flen,
                    "mates should agree on fl for qname {qname}"
                );
                assert!(
                    (*cnt - expected_cnt).abs() <= 1e-12,
                    "mates should carry identical nw tags for qname {qname}, got {expected_cnt} vs {cnt}"
                );
            }

            grouped.push((fragment_start, expected_cnt));
        }

        grouped.sort_by(|left, right| left.0.cmp(&right.0));
        Ok(grouped)
    }

    fn fixed_chromosome_args() -> ChromosomeArgs {
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        }
    }

    fn three_chromosome_args() -> ChromosomeArgs {
        ChromosomeArgs {
            chromosomes: Some(vec![
                "chr1".to_string(),
                "chr2".to_string(),
                "chr3".to_string(),
            ]),
            chromosomes_file: None,
        }
    }

    #[cfg(feature = "cmd_fragment_count_weights")]
    fn mixed_length_fragment_bam(name: &str) -> Result<TempBam> {
        bam_from_fragments(
            name,
            vec![("chr1".to_string(), 100)],
            vec![
                PairedFragmentSpec::new(0, 0, 20, 10).build()?,
                PairedFragmentSpec::new(0, 20, 60, 10).build()?,
            ],
            Vec::new(),
        )
    }

    /// Build a tiny BAM with two contigs and three inward-directed pairs, coordinate-sorted.
    /// All reads are primary, mapped, and paired. CIGAR is 50M.
    fn write_test_bam(path: &Path) -> Result<()> {
        // Header with two contigs
        let mut hdr = Header::new();
        hdr.push_record(
            bam::header::HeaderRecord::new(b"HD")
                .push_tag(b"VN", &"1.6")
                .push_tag(b"SO", &"coordinate"),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &100_000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr2")
                .push_tag(b"LN", &100_000),
        );

        // Create BAM writer
        let mut writer = Writer::from_path(path, &hdr, Format::Bam).context("create BAM writer")?;
        let header_view = HeaderView::from_header(&hdr);

        // Convenience closures
        let tid_chr1 = header_view.tid(b"chr1").expect("chr1 present") as i32;
        let tid_chr2 = header_view.tid(b"chr2").expect("chr2 present") as i32;

        // QNAMEs
        let q1 = b"pairA";
        let q2 = b"pairB";
        let q3 = b"pairC";

        // CIGAR 50M sequence and qual
        let cigar = vec![Cigar::Match(50)];
        let seq = b"ACGTN".repeat(10); // 50bp
        let qual = vec![30u8; 50];

        // chr1 pair A: R1 forward @10002 (MAPQ 60), R2 reverse @10090 (MAPQ 60)
        let r1_a = make_rec(
            q1, tid_chr1, 10002, false, 60, &cigar, &seq, &qual, true, tid_chr1, 10090, true,
        );
        let r2_a = make_rec(
            q1, tid_chr1, 10090, true, 60, &cigar, &seq, &qual, false, tid_chr1, 10002, false,
        );

        // chr1 pair B: R1 forward @10003 (MAPQ 0), R2 reverse @10087 (MAPQ 0)
        let r1_b = make_rec(
            q2, tid_chr1, 10003, false, 0, &cigar, &seq, &qual, true, tid_chr1, 10087, true,
        );
        let r2_b = make_rec(
            q2, tid_chr1, 10087, true, 0, &cigar, &seq, &qual, false, tid_chr1, 10003, false,
        );

        // chr2 pair C: R1 reverse @20090 (MAPQ 30), R2 forward @20000 (MAPQ 40)
        let r1_c = make_rec(
            q3, tid_chr2, 20090, true, 30, &cigar, &seq, &qual, true, tid_chr2, 20000, false,
        );
        let r2_c = make_rec(
            q3, tid_chr2, 20000, false, 40, &cigar, &seq, &qual, false, tid_chr2, 20090, true,
        );

        // Write in coordinate order
        writer.write(&r1_a)?;
        writer.write(&r1_b)?;
        writer.write(&r2_b)?;
        writer.write(&r2_a)?;
        writer.write(&r2_c)?;
        writer.write(&r1_c)?;

        Ok(())
    }

    fn write_three_chrom_window_bam(path: &Path) -> Result<()> {
        let mut hdr = Header::new();
        hdr.push_record(
            bam::header::HeaderRecord::new(b"HD")
                .push_tag(b"VN", &"1.6")
                .push_tag(b"SO", &"coordinate"),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &100_000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr2")
                .push_tag(b"LN", &100_000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr3")
                .push_tag(b"LN", &100_000),
        );

        let mut writer = Writer::from_path(path, &hdr, Format::Bam).context("create BAM writer")?;
        let header_view = HeaderView::from_header(&hdr);

        let tid_chr1 = header_view.tid(b"chr1").expect("chr1 present") as i32;
        let tid_chr2 = header_view.tid(b"chr2").expect("chr2 present") as i32;
        let tid_chr3 = header_view.tid(b"chr3").expect("chr3 present") as i32;

        let cigar = vec![Cigar::Match(40)];
        let seq = b"ACGTN".repeat(8);
        let qual = vec![30u8; 40];

        let records = vec![
            make_rec(
                b"chr1_keep",
                tid_chr1,
                10,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr1,
                90,
                true,
            ),
            make_rec(
                b"chr1_keep",
                tid_chr1,
                90,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr1,
                10,
                false,
            ),
            make_rec(
                b"chr1_drop",
                tid_chr1,
                400,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr1,
                480,
                true,
            ),
            make_rec(
                b"chr1_drop",
                tid_chr1,
                480,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr1,
                400,
                false,
            ),
            make_rec(
                b"chr2_keep",
                tid_chr2,
                30,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr2,
                110,
                true,
            ),
            make_rec(
                b"chr2_keep",
                tid_chr2,
                110,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr2,
                30,
                false,
            ),
            make_rec(
                b"chr2_drop",
                tid_chr2,
                430,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr2,
                510,
                true,
            ),
            make_rec(
                b"chr2_drop",
                tid_chr2,
                510,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr2,
                430,
                false,
            ),
            make_rec(
                b"chr3_keep",
                tid_chr3,
                50,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr3,
                130,
                true,
            ),
            make_rec(
                b"chr3_keep",
                tid_chr3,
                130,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr3,
                50,
                false,
            ),
            make_rec(
                b"chr3_drop",
                tid_chr3,
                460,
                false,
                60,
                &cigar,
                &seq,
                &qual,
                true,
                tid_chr3,
                540,
                true,
            ),
            make_rec(
                b"chr3_drop",
                tid_chr3,
                540,
                true,
                60,
                &cigar,
                &seq,
                &qual,
                false,
                tid_chr3,
                460,
                false,
            ),
        ];

        for record in records {
            writer.write(&record)?;
        }

        Ok(())
    }

    fn parse_frag_rows(rows: &[String]) -> Vec<(String, u64, u64, u8, char)> {
        rows.iter()
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                assert_eq!(columns.len(), 5, "Bad line: {line}");
                (
                    columns[0].to_string(),
                    columns[1].parse().unwrap(),
                    columns[2].parse().unwrap(),
                    columns[3].parse().unwrap(),
                    columns[4].chars().next().unwrap(),
                )
            })
            .collect()
    }

    /// Construct a paired-end record.
    ///
    /// - `qname`: read name shared by the pair
    /// - `tid`, `pos`: target and 0-based start
    /// - `is_rev`: strand
    /// - `mapq`: mapping quality
    /// - `cigar`, `seq`, `qual`: alignment, bases, quals
    /// - `is_first_in_template`: true for R1, false for R2
    /// - `mtid`, `mpos`: mate reference and 0-based pos
    /// - `mate_is_reverse`: strand flag for the mate (sets FLAG_MATE_REVERSE)
    fn make_rec(
        qname: &[u8],
        tid: i32,
        pos: i64,
        is_rev: bool,
        mapq: u8,
        cigar: &[Cigar],
        seq: &[u8],
        base_quality: &[u8],
        is_first_in_template: bool,
        mtid: i32,
        mpos: i64,
        mate_is_reverse: bool,
    ) -> Record {
        let mut rec = Record::new();
        let cigar_string = CigarString(cigar.to_vec());
        rec.set(qname, Some(&cigar_string), seq, base_quality);

        let mut flags: u16 = 0;
        flags |= 0x1; // paired
        if is_first_in_template {
            flags |= 0x40;
        } else {
            flags |= 0x80;
        }
        if is_rev {
            flags |= 0x10;
        }
        if mate_is_reverse {
            flags |= 0x20;
        }
        // No secondary/supplementary/duplicate/fail flags set.

        rec.set_tid(tid);
        rec.set_pos(pos);
        rec.set_flags(flags);
        rec.set_mapq(mapq);
        rec.set_mtid(mtid);
        rec.set_mpos(mpos);

        // Optional: template length (insert size). Sign convention is aligner-specific; not required by the code under test.
        // We skip setting isize.

        rec
    }
}
