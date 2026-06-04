#![cfg(feature = "cmd_transitions")]

mod fixtures;

mod transitions_command_tests {
    use crate::fixtures::single_position_selection;
    use anyhow::{Context, Result};
    use cfdnalab::RunOptions;
    use cfdnalab::run_like_cli::common::{
        ChromosomeArgs, IOCArgs, Ref2BitRequiredArgs, ReferenceFrame,
    };
    use cfdnalab::run_like_cli::transitions::{TransitionsConfig, run_transitions};
    use cfdnalab::testing::{
        Cigar, FragmentSpec, ReadSpec, bam_from_fragments, twobit_from_sequences,
    };
    use ndarray::{Array3, s};
    use ndarray_npy::read_npy;
    use std::collections::HashMap;
    use std::fs;
    use tempfile::TempDir;

    fn run_transitions_quiet(cfg: &TransitionsConfig) -> Result<()> {
        run_transitions(cfg, RunOptions::new_quiet()).map(|_| ())
    }

    fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
        ChromosomeArgs {
            chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
            chromosomes_file: None,
        }
    }

    /// Build two inward-facing fragments whose left starts differ by four bases.
    ///
    /// Layout (reference positions in brackets):
    /// - Fragment A forward read: [0..4), reverse read: [20..24)
    /// - Fragment B forward read: [4..8), reverse read: [24..28)
    ///
    /// By spacing the fragments four bases apart, the first three offsets (0, 1, 2) sample
    /// different dinucleotide transitions from each fragment, which leads to distinct
    /// conditional probabilities that we can reason about analytically.
    fn synthetic_fragments() -> Vec<FragmentSpec> {
        const FLAG_FIRST_MATE: u16 = 0x40;
        const FLAG_SECOND_MATE: u16 = 0x80;
        const FLAG_PROPER_PAIR: u16 = 0x2;
        const FLAG_MATE_REVERSE: u16 = 0x20;

        fn make_fragment(start: i64, mate_start: i64) -> FragmentSpec {
            let read_len = 4u32;
            let insert = mate_start - start + read_len as i64;
            FragmentSpec {
                forward: ReadSpec {
                    tid: 0,
                    pos: start,
                    cigar: vec![Cigar::Match(read_len)],
                    seq: vec![b'A'; read_len as usize],
                    base_quality: 30,
                    is_reverse: false,
                    mapq: 60,
                    flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                    mate_tid: Some(0),
                    mate_pos: Some(mate_start),
                    insert_size: insert,
                },
                reverse: ReadSpec {
                    tid: 0,
                    pos: mate_start,
                    cigar: vec![Cigar::Match(read_len)],
                    seq: vec![b'T'; read_len as usize],
                    base_quality: 30,
                    is_reverse: true,
                    mapq: 60,
                    flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                    mate_tid: Some(0),
                    mate_pos: Some(start),
                    insert_size: -insert,
                },
            }
        }

        vec![make_fragment(0, 20), make_fragment(4, 24)]
    }

    #[test]
    fn run_transitions_produces_expected_frequencies() -> Result<()> {
        // Reference repeats "AACGAGTTACGA" so the motifs at successive offsets are predictable
        // Reference is a repetition of "AACGAGTTACGA". Motifs encountered by the fragments:
        // Fragment A offsets 0,1,2 -> "AA", "AC", "CG"
        // Fragment B offsets 0,1,2 -> (shifted by four bases) -> "AG", "GT", "TT"
        // These expectations underpin the assertions later in the test.
        let reference_seq = "AACGAGTTACGAACGAGTTACGAACGAGTTACGA";
        let reference = twobit_from_sequences(
            "transitions_manual",
            vec![("chr1".to_string(), reference_seq.to_string())],
        )?;
        let chroms = vec![("chr1".to_string(), reference_seq.len() as u32)];
        let bam = bam_from_fragments(
            "transitions_manual_bam",
            chroms,
            synthetic_fragments(),
            Vec::new(),
        )?;
        let out_dir = TempDir::new()?;

        let mut cfg = TransitionsConfig::new(
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
        cfg.set_output_prefix("transitions_simple".to_string());
        cfg.set_orders(vec![1]);
        cfg.set_canonical(false);
        cfg.set_save_sparse(false);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_position_selection(single_position_selection(ReferenceFrame::Left, "1..4", 1));
        {
            let lengths = cfg.fragment_lengths_mut();
            // The synthetic fragments span [0, 24) and [4, 28), so both have length 24.
            lengths.min_fragment_length = 24;
            lengths.max_fragment_length = 24;
        }

        run_transitions_quiet(&cfg)?;

        let prefix = cfg.shared_args.output_prefix.trim();
        let freqs_path = out_dir.path().join(format!("{prefix}.k2_left_freqs.npy"));
        let motifs_path = out_dir.path().join(format!("{prefix}.k2_left_motifs.txt"));
        let positions_path = out_dir.path().join(format!("{prefix}.left_positions.txt"));

        assert!(
            freqs_path.exists(),
            "expected freqs at {}",
            freqs_path.display()
        );
        assert!(
            motifs_path.exists(),
            "expected motifs at {}",
            motifs_path.display()
        );
        assert!(
            positions_path.exists(),
            "expected positions at {}",
            positions_path.display()
        );

        // Ensure temporary working directories are cleaned up after the run
        for entry in fs::read_dir(out_dir.path())? {
            let entry = entry?;
            if entry.file_type()?.is_dir() {
                assert!(
                    !entry.file_name().to_string_lossy().starts_with("tmp."),
                    "temporary directory {} should be removed",
                    entry.path().display()
                );
            }
        }

        let freqs: Array3<f64> = read_npy(&freqs_path)?;
        assert_eq!(freqs.shape()[0], 1);

        let motifs: Vec<String> = fs::read_to_string(&motifs_path)?
            .lines()
            .map(|line| line.to_string())
            .collect();

        // `fragment-kmers` stores offsets as zero-based integers. With "1..4" we expect
        // offsets {0,1,2} in the dense output. Record the index of each offset so we can
        // slice the frequency array with explicit intent
        let positions: Vec<i32> = fs::read_to_string(&positions_path)?
            .lines()
            .map(|line| line.parse::<i32>().expect("position offset"))
            .collect();
        let offset0_idx = positions
            .iter()
            .position(|&p| p == 0)
            .expect("offset 0 should be present for position range 1..4");
        let offset1_idx = positions
            .iter()
            .position(|&p| p == 1)
            .expect("offset 1 should be present for position range 1..4");
        let offset2_idx = positions
            .iter()
            .position(|&p| p == 2)
            .expect("offset 2 should be present for position range 1..4");

        // Forward fragments begin at positions 0 and 4, so offsets 0..2 observe two different dinucleotides each:
        // offset 0 -> {"AA", "AG"}, offset 1 -> {"AC", "GT"}, offset 2 -> {"CG", "TT"}.
        let reference_seq = reference
            .sequence("chr1")
            .context("missing chr1 sequence in reference fixture")?;

        let motif_offset0_fragment1 = &reference_seq[0..2];
        let motif_offset0_fragment2 = &reference_seq[4..6];
        let motif_offset1_fragment1 = &reference_seq[1..3];
        let motif_offset1_fragment2 = &reference_seq[5..7];
        let motif_offset2_fragment1 = &reference_seq[2..4];
        let motif_offset2_fragment2 = &reference_seq[6..8];

        // Build a map from motif -> column index for fast lookups when we project the slices.
        // Precompute motif -> column index so assertions can address probabilities by name.
        let motif_idx: HashMap<&str, usize> = motifs
            .iter()
            .enumerate()
            .map(|(idx, motif)| (motif.as_str(), idx))
            .collect();

        let slice_offset0 = freqs.slice(s![0, offset0_idx, ..]);
        // Expected conditional probabilities (prefix-normalised):
        // offset 0 -> { "AA": 0.5, "AG": 0.5 }
        // offset 1 -> { "AC": 1.0 (prefix 'A'), "GT": 1.0 (prefix 'G') }
        // offset 2 -> { "CG": 1.0 (prefix 'C'), "TT": 1.0 (prefix 'T') }

        // Offset 0 sees two starts (positions 0 and 4); both motifs share prefix 'A', so the
        // conditional probability splits evenly between them.
        for (motif, expected) in [
            (motif_offset0_fragment1, 0.5),
            (motif_offset0_fragment2, 0.5),
        ] {
            let idx = *motif_idx
                .get(motif)
                .with_context(|| format!("missing motif {motif} at offset 0"))?;
            let observed = slice_offset0[idx];
            assert!(
                (observed - expected).abs() < 1e-9,
                "offset 0 motif {motif} expected {expected} observed {observed}"
            );
        }
        for (&motif, &idx) in &motif_idx {
            if motif == motif_offset0_fragment1 || motif == motif_offset0_fragment2 {
                continue;
            }
            // Every other motif should retain zero probability mass at this offset.
            let observed = slice_offset0[idx];
            assert!(
                observed.abs() < 1e-9,
                "offset 0 motif {motif} should be 0, observed {observed}"
            );
        }

        let slice_offset1 = freqs.slice(s![0, offset1_idx, ..]);
        // Offset 1 corresponds to positions 1 and 5, yielding "AC" (prefix 'A') and "GT" (prefix 'G').
        // Because compute_transition_frequencies normalises within each prefix bucket, both motifs
        // retain probability 1.0 for their respective prefixes.
        for (motif, expected) in [
            (motif_offset1_fragment1, 1.0),
            (motif_offset1_fragment2, 1.0),
        ] {
            let idx = *motif_idx
                .get(motif)
                .with_context(|| format!("missing motif {motif} at offset 1"))?;
            let observed = slice_offset1[idx];
            assert!(
                (observed - expected).abs() < 1e-9,
                "offset 1 motif {motif} expected {expected} observed {observed}"
            );
        }
        for (&motif, &idx) in &motif_idx {
            if motif == motif_offset1_fragment1 || motif == motif_offset1_fragment2 {
                continue;
            }
            // Only "AC" and "GT" should be present at offset 1.
            let observed = slice_offset1[idx];
            assert!(
                observed.abs() < 1e-9,
                "offset 1 motif {motif} should be 0, observed {observed}"
            );
        }

        let slice_offset2 = freqs.slice(s![0, offset2_idx, ..]);
        // Offset 2 produces "CG" (prefix 'C') and "TT" (prefix 'T'), again leading to per-prefix
        // probabilities of 1.0 because each prefix only observes a single continuation.
        for (motif, expected) in [
            (motif_offset2_fragment1, 1.0),
            (motif_offset2_fragment2, 1.0),
        ] {
            let idx = *motif_idx
                .get(motif)
                .with_context(|| format!("missing motif {motif} at offset 2"))?;
            let observed = slice_offset2[idx];
            assert!(
                (observed - expected).abs() < 1e-9,
                "offset 2 motif {motif} expected {expected} observed {observed}"
            );
        }
        for (&motif, &idx) in &motif_idx {
            if motif == motif_offset2_fragment1 || motif == motif_offset2_fragment2 {
                continue;
            }
            // Offset 2 is restricted to motifs "CG" and "TT".
            let observed = slice_offset2[idx];
            assert!(
                observed.abs() < 1e-9,
                "offset 2 motif {motif} should be 0, observed {observed}"
            );
        }

        Ok(())
    }
}
