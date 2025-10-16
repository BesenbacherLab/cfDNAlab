mod fixtures;

mod tests_fragment_kmer_command {
    use std::collections::HashMap;
    use std::path::Path;

    use crate::fixtures::{
        fragment_kmers_edge_bam, fragment_kmers_edge_reference, simple_inward_bam,
        simple_reference_twobit, write_bed, write_scaling_factors,
    };
    use anyhow::{Context, Result};
    use cfdnalab::commands::cli_common::{
        ChromosomeArgs, FragmentLengthArgs, FragmentPositionSelectionArgs, IOCArgs,
        Ref2BitRequiredArgs, ScaleGenomeArgs, WindowsArgs,
    };
    use cfdnalab::commands::fragment_kmers::config::FragmentKmersConfig;
    use cfdnalab::commands::fragment_kmers::fragment_kmers::run;
    use cfdnalab::commands::visualize_positions::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
    use cfdnalab::shared::base::make_canonical;
    use cfdnalab::shared::blacklist::BlacklistStrategy;
    use cfdnalab::shared::indel_mode::IndelMode;
    use ndarray::Array2;
    use ndarray_npy::read_npy;
    use tempfile::TempDir;

    fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
        ChromosomeArgs {
            chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
            chromosomes_file: None,
        }
    }

    #[test]
    fn counts_dinucleotides_in_global_window() -> Result<()> {
        let bam = simple_inward_bam()?;
        let reference = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;

        let mut cfg = FragmentKmersConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_output_prefix("kmers".to_string());
        cfg.set_kmer_sizes(vec![2]);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_canonical(false);

        run(&cfg)?;

        let counts_path = out_dir.path().join("kmers.k2_counts.npy");
        let motifs_path = out_dir.path().join("kmers.k2_motifs.txt");
        assert!(counts_path.exists());
        assert!(motifs_path.exists());

        let counts: Array2<f64> = read_npy(&counts_path)?;
        assert_eq!(counts.shape(), &[1, 16]);
        let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
            .lines()
            .map(|s| s.to_string())
            .collect();

        let chr1_seq = reference
            .sequence("chr1")
            .context("missing chr1 sequence in reference fixture")?;
        let start = 20usize;
        let end = 80usize; // fragment end (exclusive)
        let k = 2usize;
        let mut expected: HashMap<String, f64> = HashMap::new();
        for idx in start..=(end - k) {
            let motif = chr1_seq
                .get(idx..idx + k)
                .context("motif slice")?
                .to_string();
            *expected.entry(motif).or_insert(0.0) += 1.0;
        }

        let row = counts.row(0);
        let total: f64 = row.sum();
        assert!((total - 59.0).abs() < 1e-6);

        for (col, motif) in motif_list.iter().enumerate() {
            let expected_val = expected.get(motif).copied().unwrap_or(0.0);
            assert!(
                (row[col] - expected_val).abs() < 1e-6,
                "motif {motif} expected {expected_val} observed {}",
                row[col]
            );
        }

        Ok(())
    }

    #[test]
    fn positional_counts_restricts_starts() -> Result<()> {
        let bam = simple_inward_bam()?;
        let reference = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;

        let mut cfg = FragmentKmersConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_output_prefix("kmers_first".to_string());
        cfg.set_kmer_sizes(vec![2]);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_canonical(false);
        cfg.set_positional_counts(true);
        cfg.set_position_selection(FragmentPositionSelectionArgs {
            frame: ReferenceFrame::Left,
            positions: "1..1".to_string(),
            step: 1,
            bases_from: BasesFrom::Reference,
            mismatch_bases_from: MismatchBasesFrom::NearestRead,
        });

        run(&cfg)?;

        let counts_path = out_dir.path().join("kmers_first.k2_counts.npy");
        let motifs_path = out_dir.path().join("kmers_first.k2_motifs.txt");
        assert!(counts_path.exists());
        assert!(motifs_path.exists());

        let counts: Array2<f64> = read_npy(&counts_path)?;
        assert_eq!(counts.shape(), &[1, 16]);
        let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
            .lines()
            .map(|s| s.to_string())
            .collect();

        let chr1_seq = reference
            .sequence("chr1")
            .context("missing chr1 sequence in reference fixture")?;
        let start = 20usize;
        let k = 2usize;
        let motif = chr1_seq
            .get(start..start + k)
            .context("motif slice")?
            .to_string();

        let mut expected: HashMap<String, f64> = HashMap::new();
        expected.insert(motif.clone(), 1.0);

        let mut actual: HashMap<String, f64> = HashMap::new();
        for (motif_idx, motif) in motif_list.iter().enumerate() {
            let value = counts[[0, motif_idx]];
            if value != 0.0 {
                actual.insert(motif.clone(), value);
            }
        }

        assert_counts_close(&actual, &expected);
        Ok(())
    }

    #[test]
    fn canonical_trimers_collapse_matches_manual_counts() -> Result<()> {
        let bam = simple_inward_bam()?;
        let reference = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;

        let mut cfg = FragmentKmersConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_output_prefix("kmers".to_string());
        cfg.set_kmer_sizes(vec![3]);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_canonical(true);

        run(&cfg)?;

        let counts_path = out_dir.path().join("kmers.k3_counts.npy");
        let motifs_path = out_dir.path().join("kmers.k3_motifs.txt");
        assert!(counts_path.exists());
        assert!(motifs_path.exists());

        let counts: Array2<f64> = read_npy(&counts_path)?;
        assert_eq!(counts.shape()[0], 1);
        let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
            .lines()
            .map(|s| s.to_string())
            .collect();

        let chr1_seq = reference
            .sequence("chr1")
            .context("missing chr1 sequence in reference fixture")?;
        let start = 20usize;
        let end = 80usize;
        let k = 3usize;
        let mut expected: HashMap<String, f64> = HashMap::new();
        for idx in start..=(end - k) {
            let motif = chr1_seq
                .get(idx..idx + k)
                .context("motif slice")?
                .to_string();
            let canon = make_canonical(motif);
            *expected.entry(canon).or_insert(0.0) += 1.0;
        }

        let row = counts.row(0);
        let total: f64 = row.sum();
        assert!((total - 58.0).abs() < 1e-6);

        for motif in &motif_list {
            assert_eq!(motif, &make_canonical(motif.clone()));
        }

        for (col, motif) in motif_list.iter().enumerate() {
            let expected_val = expected.get(motif).copied().unwrap_or(0.0);
            assert!(
                (row[col] - expected_val).abs() < 1e-6,
                "motif {motif} expected {expected_val} observed {}",
                row[col]
            );
        }

        Ok(())
    }

    fn load_counts_from_output(dir: &Path, prefix: &str, k: u8) -> Result<HashMap<String, f64>> {
        let counts_path = dir.join(format!("{prefix}.k{k}_counts.npy"));
        let motifs_path = dir.join(format!("{prefix}.k{k}_motifs.txt"));
        let counts: Array2<f64> = read_npy(&counts_path)?;
        assert_eq!(
            counts.shape()[0],
            1,
            "counts matrix should have one window row"
        );
        let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
            .lines()
            .map(|s| s.to_string())
            .collect();
        let mut out = HashMap::new();
        for (idx, motif) in motif_list.iter().enumerate() {
            out.insert(motif.clone(), counts[(0, idx)]);
        }
        Ok(out)
    }

    fn assert_counts_close(actual: &HashMap<String, f64>, expected: &HashMap<String, f64>) {
        for (motif, exp) in expected {
            let obs = actual.get(motif).copied().unwrap_or(0.0);
            assert!(
                (obs - exp).abs() < 1e-6,
                "motif {motif} expected {exp} observed {obs}"
            );
        }
        for (motif, obs) in actual {
            let exp = expected.get(motif).copied().unwrap_or(0.0);
            assert!(
                (obs - exp).abs() < 1e-6,
                "motif {motif} expected {exp} observed {obs}"
            );
        }
    }

    #[test]
    fn complex_edge_cases_respect_scaling_and_blacklists() -> Result<()> {
        let bam = fragment_kmers_edge_bam()?;
        let reference = fragment_kmers_edge_reference()?;
        let out_dir = TempDir::new()?;
        let chromosomes = ["chr1"];

        let fragment_lengths = FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 1000,
        };

        let mut cfg_base = FragmentKmersConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            base_chromosomes(&chromosomes),
        );
        cfg_base.set_output_prefix("edge_base".to_string());
        cfg_base.set_kmer_sizes(vec![2]);
        cfg_base.set_indel_mode(IndelMode::Adjust);
        cfg_base.set_min_mapq(0);
        cfg_base.set_require_proper_pair(false);
        cfg_base.set_canonical(false);
        cfg_base.set_ignore_gap(true);
        cfg_base.set_positional_counts(true);
        cfg_base.set_position_selection(FragmentPositionSelectionArgs {
            frame: ReferenceFrame::Left,
            positions: "2..-1".to_string(),
            step: 1,
            bases_from: BasesFrom::Reference,
            mismatch_bases_from: MismatchBasesFrom::NearestRead,
        });
        {
            let fl = cfg_base.fragment_lengths_mut();
            fl.min_fragment_length = fragment_lengths.min_fragment_length;
            fl.max_fragment_length = fragment_lengths.max_fragment_length;
        }

        run(&cfg_base)?;

        let observed_base = load_counts_from_output(out_dir.path(), "edge_base", 2)?;
        let expected_base: HashMap<String, f64> = vec![
            ("AA", 2.0),
            ("AC", 5.0),
            ("AG", 1.0),
            ("CC", 5.0),
            ("CG", 4.0),
            ("CT", 3.0),
            ("GA", 1.0),
            ("GG", 1.0),
            ("GT", 3.0),
            ("TA", 5.0),
            ("TG", 1.0),
            ("TT", 2.0),
        ]
        .into_iter()
        .map(|(m, c)| (m.to_string(), c))
        .collect();
        // Sanity check: manually computed expectations for edge-case fragments.
        assert_counts_close(&observed_base, &expected_base);

        let blacklist_path = out_dir.path().join("mask.bed");
        write_bed(
            &blacklist_path,
            &[("chr1", 9, 11, "mask"), ("chr1", 22, 23, "mask")],
        )?;
        let scaling_path = out_dir.path().join("scaling.tsv");
        write_scaling_factors(
            &scaling_path,
            &[
                ("chr1", 0, 6, 1.0),
                ("chr1", 6, 8, 0.0),
                ("chr1", 8, 20, 1.5),
                ("chr1", 20, 40, 0.5),
            ],
        )?;

        let mut cfg_scaled = cfg_base.clone();
        cfg_scaled.set_output_prefix("edge_scaled".to_string());
        cfg_scaled.blacklist = Some(vec![blacklist_path.clone()]);
        cfg_scaled.blacklist_strategy = BlacklistStrategy::Proportion(1.0);
        let mut scale_args = ScaleGenomeArgs::default();
        scale_args.scaling_factors = Some(scaling_path.clone());
        cfg_scaled.set_scale_genome(scale_args);

        run(&cfg_scaled)?;

        let observed_scaled = load_counts_from_output(out_dir.path(), "edge_scaled", 2)?;

        let expected_scaled: HashMap<String, f64> = vec![
            ("AA", 3.0),
            ("AC", 3.0),
            ("CC", 4.5),
            ("CG", 4.0),
            ("CT", 1.5),
            ("GA", 1.0),
            ("GG", 1.5),
            ("GT", 3.0),
            ("TA", 4.0),
            ("TG", 1.0),
            ("TT", 0.5),
        ]
        .into_iter()
        .map(|(m, c)| (m.to_string(), c))
        .collect();

        assert_counts_close(&observed_scaled, &expected_scaled);

        Ok(())
    }
}

#[cfg(test)]
mod tests_fragment_kmers_tiling {
    use anyhow::Result;
    use cfdnalab::{
        commands::fragment_kmers::tiling::{
            TileKmerCountEntry, TileWindowCounts, merge_tile_counts,
        },
        shared::kmers::kmer_codec::{KmerOrientation, KmerSpec, build_kmer_specs},
    };

    fn code_for_motif(spec: &KmerSpec, motif: &str) -> u64 {
        let limit = 5u64.pow(spec.k as u32);
        for code in 0..limit {
            if spec.decode_kmer(code) == motif {
                return code;
            }
        }
        panic!("motif {} not encodable", motif);
    }

    #[test]
    fn merge_tile_counts_merges_two_tiles() -> Result<()> {
        let kmer_specs = build_kmer_specs(&[3])?;
        let spec3 = &kmer_specs[&3];
        let code_aaa = code_for_motif(spec3, "AAA");

        let payload_a = vec![TileWindowCounts {
            original_idx: 0,
            entries: vec![TileKmerCountEntry {
                k: 3,
                code: code_aaa,
                orientation: KmerOrientation::Forward,
                value: 1.5,
            }],
        }];

        let payload_b = vec![TileWindowCounts {
            original_idx: 0,
            entries: vec![TileKmerCountEntry {
                k: 3,
                code: code_aaa,
                orientation: KmerOrientation::Forward,
                value: 2.0,
            }],
        }];

        let merged = merge_tile_counts(vec![payload_a, payload_b], 1, &kmer_specs)?;
        assert_eq!(merged.len(), 1);
        let window_counts = merged[0].counts.get(&3).unwrap();
        let value = window_counts.get("AAA").copied().unwrap_or_default();
        assert!((value - 3.5).abs() < 1e-9);
        Ok(())
    }

    #[test]
    fn merge_tile_counts_merges_three_tiles() -> Result<()> {
        let kmer_specs = build_kmer_specs(&[3])?;
        let spec3 = &kmer_specs[&3];
        let code_aaa = code_for_motif(spec3, "AAA");
        let code_aac = code_for_motif(spec3, "AAC");

        let payload_1 = vec![
            TileWindowCounts {
                original_idx: 0,
                entries: vec![TileKmerCountEntry {
                    k: 3,
                    code: code_aaa,
                    orientation: KmerOrientation::Forward,
                    value: 1.0,
                }],
            },
            TileWindowCounts {
                original_idx: 1,
                entries: vec![TileKmerCountEntry {
                    k: 3,
                    code: code_aac,
                    orientation: KmerOrientation::Forward,
                    value: 2.0,
                }],
            },
        ];

        let payload_2 = vec![
            TileWindowCounts {
                original_idx: 0,
                entries: vec![TileKmerCountEntry {
                    k: 3,
                    code: code_aaa,
                    orientation: KmerOrientation::Forward,
                    value: 3.0,
                }],
            },
            TileWindowCounts {
                original_idx: 2,
                entries: vec![TileKmerCountEntry {
                    k: 3,
                    code: code_aaa,
                    orientation: KmerOrientation::Forward,
                    value: 5.0,
                }],
            },
        ];

        let payload_3 = vec![
            TileWindowCounts {
                original_idx: 0,
                entries: vec![TileKmerCountEntry {
                    k: 3,
                    code: code_aaa,
                    orientation: KmerOrientation::Forward,
                    value: 0.5,
                }],
            },
            TileWindowCounts {
                original_idx: 1,
                entries: vec![TileKmerCountEntry {
                    k: 3,
                    code: code_aac,
                    orientation: KmerOrientation::Forward,
                    value: 1.5,
                }],
            },
        ];

        let merged = merge_tile_counts(vec![payload_1, payload_2, payload_3], 3, &kmer_specs)?;
        assert_eq!(merged.len(), 3);

        let win0 = merged[0].counts.get(&3).unwrap();
        assert!((win0.get("AAA").copied().unwrap_or_default() - 4.5).abs() < 1e-9);

        let win1 = merged[1].counts.get(&3).unwrap();
        assert!((win1.get("AAC").copied().unwrap_or_default() - 3.5).abs() < 1e-9);

        let win2 = merged[2].counts.get(&3).unwrap();
        assert_eq!(win2.len(), 1);
        assert!((win2.get("AAA").copied().unwrap_or_default() - 5.0).abs() < 1e-9);
        Ok(())
    }

    #[test]
    fn merge_tile_counts_rejects_out_of_range_indices() {
        let kmer_specs = build_kmer_specs(&[3]).expect("build specs");
        let spec3 = &kmer_specs[&3];
        let code_aaa = code_for_motif(spec3, "AAA");

        let payload = vec![TileWindowCounts {
            original_idx: 5,
            entries: vec![TileKmerCountEntry {
                k: 3,
                code: code_aaa,
                orientation: KmerOrientation::Forward,
                value: 1.0,
            }],
        }];

        let result = merge_tile_counts(vec![payload], 2, &kmer_specs);
        assert!(result.is_err());
    }
}

mod tests_fragment_kmer_positions {
    use std::num::NonZeroUsize;

    use cfdnalab::{
        commands::{
            fragment_kmers::{
                fragment_kmers::count_kmers_at_positions,
                positions::{PositionSelection, PositionSelectionCache},
            },
            visualize_positions::{LinearRange, PositionsSpec, ReferenceFrame},
        },
        shared::{
            fragment::segment_kmer_fragment::FragmentWithKmerSegments,
            kmers::kmer_codec::{
                Kmer, KmerCodes, KmerOrientation, KmerSpec, build_kmer_specs,
                build_left_aligned_codes_per_k,
            },
        },
    };
    use fxhash::FxHashMap;
    use smallvec::smallvec;

    #[test]
    fn given_left_frame_when_counting_then_collects_forward_kmers() {
        let seq = b"ACGTAC";
        let context = TestContext::new(seq);
        let cache = build_cache(
            ReferenceFrame::Left,
            PositionsSpec::Linear(LinearRange::All),
            context.fragment.len(),
        );
        let selections = cache
            .offsets(context.fragment.len())
            .expect("left frame offsets");

        let mut counts = FxHashMap::default();
        count_kmers_at_positions(
            &context.fragment,
            selections,
            &context.positional_codes_by_k,
            &context.kmer_specs,
            &mut counts,
            None,
            0,
            context.fragment.len(),
        );

        let expected = expected_counts(
            selections,
            context.fragment.len() as usize,
            &context.positional_codes_by_k[&context.k],
            context.k as usize,
        );
        assert_eq!(counts, expected);
        assert!(
            counts
                .keys()
                .all(|kmer| matches!(kmer.orientation, KmerOrientation::Forward))
        );
    }

    #[test]
    fn given_right_frame_when_counting_then_collects_reverse_kmers() {
        let seq = b"ACGTAC";
        let context = TestContext::new(seq);
        let cache = build_cache(
            ReferenceFrame::Right,
            PositionsSpec::Linear(LinearRange::All),
            context.fragment.len(),
        );
        let selections = cache
            .offsets(context.fragment.len())
            .expect("right frame offsets");

        let mut counts = FxHashMap::default();
        count_kmers_at_positions(
            &context.fragment,
            selections,
            &context.positional_codes_by_k,
            &context.kmer_specs,
            &mut counts,
            None,
            0,
            context.fragment.len(),
        );

        let expected = expected_counts(
            selections,
            context.fragment.len() as usize,
            &context.positional_codes_by_k[&context.k],
            context.k as usize,
        );
        assert_eq!(counts, expected);
        assert!(
            counts
                .keys()
                .all(|kmer| matches!(kmer.orientation, KmerOrientation::Reverse))
        );
    }

    #[test]
    fn given_per_end_frame_when_counting_then_collects_both_orientations() {
        let seq = b"ACGTAC";
        let context = TestContext::new(seq);
        let cache = build_cache(
            ReferenceFrame::PerEnd,
            PositionsSpec::Linear(LinearRange::All),
            context.fragment.len(),
        );
        let selections = cache
            .offsets(context.fragment.len())
            .expect("per-end offsets");

        let mut counts = FxHashMap::default();
        count_kmers_at_positions(
            &context.fragment,
            selections,
            &context.positional_codes_by_k,
            &context.kmer_specs,
            &mut counts,
            None,
            0,
            context.fragment.len(),
        );

        let expected = expected_counts(
            selections,
            context.fragment.len() as usize,
            &context.positional_codes_by_k[&context.k],
            context.k as usize,
        );
        assert_eq!(counts, expected);

        assert!(
            counts
                .keys()
                .any(|kmer| kmer.orientation == KmerOrientation::Forward)
        );
        assert!(
            counts
                .keys()
                .any(|kmer| kmer.orientation == KmerOrientation::Reverse)
        );
    }

    #[test]
    fn given_nearest_frame_when_counting_then_splits_orientations_by_half() {
        let seq = b"ACGTAC";
        let context = TestContext::new(seq);
        let cache = build_cache(
            ReferenceFrame::Nearest,
            PositionsSpec::Nearest(cfdnalab::commands::visualize_positions::NearestRange::All),
            context.fragment.len(),
        );
        let selections = cache
            .offsets(context.fragment.len())
            .expect("nearest offsets");

        let mut counts = FxHashMap::default();
        count_kmers_at_positions(
            &context.fragment,
            selections,
            &context.positional_codes_by_k,
            &context.kmer_specs,
            &mut counts,
            None,
            0,
            context.fragment.len(),
        );

        let expected = expected_counts(
            selections,
            context.fragment.len() as usize,
            &context.positional_codes_by_k[&context.k],
            context.k as usize,
        );
        assert_eq!(counts, expected);

        assert!(
            counts
                .keys()
                .any(|kmer| kmer.orientation == KmerOrientation::Forward)
        );
        assert!(
            counts
                .keys()
                .any(|kmer| kmer.orientation == KmerOrientation::Reverse)
        );
    }

    struct TestContext {
        fragment: FragmentWithKmerSegments,
        kmer_specs: FxHashMap<u8, KmerSpec>,
        positional_codes_by_k: FxHashMap<u8, KmerCodes>,
        k: u8,
    }

    impl TestContext {
        fn new(seq: &[u8]) -> Self {
            let k_values = [3u8];
            let kmer_specs = build_kmer_specs(&k_values).expect("kmer specs");
            let positional_codes_by_k = build_left_aligned_codes_per_k(seq, &kmer_specs);
            let fragment_len = seq.len() as u32;
            let fragment = FragmentWithKmerSegments {
                tid: 0,
                start: 0,
                end: fragment_len,
                segments: smallvec![(0, fragment_len)],
            };

            Self {
                fragment,
                kmer_specs,
                positional_codes_by_k,
                k: k_values[0],
            }
        }
    }

    fn build_cache(
        frame: ReferenceFrame,
        positions: PositionsSpec,
        length: u32,
    ) -> PositionSelectionCache {
        PositionSelectionCache::new(
            frame,
            &positions,
            NonZeroUsize::new(1).expect("non-zero step"),
            length,
            length,
        )
        .expect("build selection cache")
    }

    fn expected_counts(
        selections: &[PositionSelection],
        fragment_len: usize,
        codes: &KmerCodes,
        k: usize,
    ) -> FxHashMap<Kmer, f64> {
        let mut expected = FxHashMap::default();
        for selection in selections {
            let offset = selection.offset() as usize;
            match selection.orientation() {
                cfdnalab::commands::fragment_kmers::positions::PositionOrientation::Forward => {
                    if offset + k > fragment_len {
                        continue;
                    }
                    let code = codes.get(offset);
                    *expected
                        .entry(Kmer {
                            k: k as u8,
                            code,
                            orientation: KmerOrientation::Forward,
                        })
                        .or_insert(0.0) += 1.0;
                }
                cfdnalab::commands::fragment_kmers::positions::PositionOrientation::Reverse => {
                    if offset + 1 < k || offset >= fragment_len {
                        continue;
                    }
                    let start = offset + 1 - k;
                    if start + k > fragment_len {
                        continue;
                    }
                    let code = codes.get(start);
                    *expected
                        .entry(Kmer {
                            k: k as u8,
                            code,
                            orientation: KmerOrientation::Reverse,
                        })
                        .or_insert(0.0) += 1.0;
                }
            }
        }
        expected
    }
}
