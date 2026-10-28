#![cfg(feature = "cmd_ends")]

//! Molecular side selection through the public runner and output loader.

use anyhow::Result;
use cfdnalab::{
    RunOptions,
    output_loaders::{EndMotifStorageMode, EndsOutput, load_ends_output},
    run_like_cli::{
        common::{BlacklistStrategy, ChromosomeArgs, IOCArgs, IndelMotifFilterPolicy},
        ends::{EndSideSelection, EndsConfig, KmerSource, run_ends},
    },
    testing::{
        Bed4Row, Cigar, FragmentSpec, PairedFragmentSpec, TempBam, TempTwoBit,
        bam::paired_fragment, bam_from_fragments, twobit_from_sequences, write_bed4,
    },
};
use std::{collections::BTreeMap, fs, path::Path};
use tempfile::TempDir;

/// Make one 40 bp fragment with independently recognizable motifs at both ends.
///
/// BAM stores bases in reference orientation. Left read bases begin ACGA, and right read
/// bases end TGCC. Swapping mate-role flags changes molecular orientation without moving
/// either alignment, which distinguishes molecular side selection from genomic left/right.
fn side_fragment(original_is_reverse: bool) -> Result<FragmentSpec> {
    let mut fragment = paired_fragment(&PairedFragmentSpec::new(0, 40, 40, 20))?;
    fragment.forward.seq[..4].copy_from_slice(b"ACGA");
    fragment.reverse.seq[16..].copy_from_slice(b"TGCC");
    if original_is_reverse {
        fragment.forward.flags = (fragment.forward.flags & !0xc0) | 0x80;
        fragment.reverse.flags = (fragment.reverse.flags & !0xc0) | 0x40;
    }
    Ok(fragment)
}

/// Reference alleles differ from the reads, so reference-backed extraction has its own oracle.
///
/// At the left boundary: AAGC | CTGA. At the right boundary: GACT | TCAG.
fn side_reference() -> Result<TempTwoBit> {
    let mut sequence = vec![b'A'; 200];
    sequence[36..44].copy_from_slice(b"AAGCCTGA");
    sequence[76..84].copy_from_slice(b"GACTTCAG");
    twobit_from_sequences(
        "side_reference",
        vec![("chr1".into(), String::from_utf8(sequence)?)],
    )
}

fn side_bam(fragment: FragmentSpec) -> Result<TempBam> {
    bam_from_fragments(
        "side_pair",
        vec![("chr1".into(), 200)],
        vec![fragment],
        vec![],
    )
}

fn side_config(bam: &Path, output_dir: &Path, inside: usize, outside: usize) -> EndsConfig {
    let mut config = EndsConfig::new(
        IOCArgs {
            bam: bam.into(),
            output_dir: output_dir.into(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".into()]),
            chromosomes_file: None,
        },
        inside,
        outside,
    );
    config.output_prefix = "sample".into();
    // The fragment starts in [0,50), while its right endpoint lies in the neighboring tile
    config.tile_size = 50;
    config
}

fn assert_global_counts(output: &EndsOutput, expected: &[(&str, f64)]) {
    assert_eq!(output.row_count(), 1);
    let observed: BTreeMap<&str, f64> = output
        .motif_labels()
        .iter()
        .enumerate()
        .filter_map(|(column, motif)| {
            let count = output.count(0, column).expect("valid output cell");
            (count != 0.0).then_some((motif.as_str(), count))
        })
        .collect();
    assert_eq!(observed, expected.iter().copied().collect());
}

#[test]
fn side_default_and_explicit_both_preserve_pooled_counts() -> Result<()> {
    // Both retains left ACGA and reverse-complements right TGCC to GGCA
    // The right outside TCAG becomes CTGA. Mate-role swaps must not change pooled results
    let reference = side_reference()?;
    for original_is_reverse in [false, true] {
        let bam = side_bam(side_fragment(original_is_reverse)?)?;
        for explicit_both in [false, true] {
            let output_dir = TempDir::new()?;
            let mut config = side_config(&bam.bam, output_dir.path(), 4, 4);
            assert_eq!(config.side, EndSideSelection::Both);
            if explicit_both {
                config.set_side(EndSideSelection::Both);
            }
            config.ref_2bit = Some(reference.path.clone());
            let result = run_ends(&config, RunOptions::new_quiet())?;
            assert_global_counts(
                &load_ends_output(result.motif_output_path)?,
                &[("AAGC_ACGA", 1.0), ("CTGA_GGCA", 1.0)],
            );
            assert_eq!(result.counters.base.counted_fragments, 1);
            assert_eq!(result.counters.counted_motifs, 2);
        }
    }
    Ok(())
}

#[test]
fn side_paired_motifs_follow_original_strand_in_dense_and_sparse_outputs() -> Result<()> {
    // Reverse-strand EM3 is RC(ACGA)=TCGT, not the genomic right-end sequence
    // Unequal halves use left outside GC / inside ACG and right outside TC / inside GCC
    let cases = [
        (false, EndSideSelection::FivePrime, 4, 0, "_ACGA"),
        (false, EndSideSelection::ThreePrime, 4, 0, "_TGCC"),
        (true, EndSideSelection::FivePrime, 4, 0, "_GGCA"),
        (true, EndSideSelection::ThreePrime, 4, 0, "_TCGT"),
        (false, EndSideSelection::FivePrime, 0, 4, "AAGC_"),
        (false, EndSideSelection::ThreePrime, 0, 4, "TCAG_"),
        (true, EndSideSelection::FivePrime, 0, 4, "CTGA_"),
        (true, EndSideSelection::ThreePrime, 0, 4, "GCTT_"),
        (false, EndSideSelection::FivePrime, 3, 2, "GC_ACG"),
        (false, EndSideSelection::ThreePrime, 3, 2, "TC_GCC"),
        (true, EndSideSelection::FivePrime, 3, 2, "GA_GGC"),
        (true, EndSideSelection::ThreePrime, 3, 2, "GC_CGT"),
    ];
    let reference = side_reference()?;
    for (original_is_reverse, side, inside, outside, expected) in cases {
        let bam = side_bam(side_fragment(original_is_reverse)?)?;
        for dense in [false, true] {
            let output_dir = TempDir::new()?;
            let mut config = side_config(&bam.bam, output_dir.path(), inside, outside);
            config.side = side;
            config.all_motifs = dense;
            if outside > 0 {
                config.ref_2bit = Some(reference.path.clone());
            }
            let result = run_ends(&config, RunOptions::new_quiet())?;
            let output = load_ends_output(result.motif_output_path)?;
            assert_global_counts(&output, &[(expected, 1.0)]);
            assert_eq!(
                output.storage_mode(),
                if dense {
                    EndMotifStorageMode::Dense
                } else {
                    EndMotifStorageMode::SparseCoo
                }
            );
            assert_eq!(
                output.motif_count(),
                if dense {
                    4_usize.pow((inside + outside) as u32)
                } else {
                    1
                }
            );
            assert_eq!(result.counters.base.counted_fragments, 1);
            assert_eq!(result.counters.counted_motifs, 1);
            let settings: serde_json::Value = serde_json::from_str(&fs::read_to_string(
                output_dir.path().join("sample.end_settings.json"),
            )?)?;
            assert_eq!(settings["side"], side.as_str());
        }
    }
    Ok(())
}

#[test]
fn side_three_prime_reference_backed_halves_keep_outside_inside_order() -> Result<()> {
    // Forward original: outside TCAG, reference inside GACT
    // Reverse original: RC(AAGC)=GCTT outside and RC(CTGA)=TCAG inside
    let reference = side_reference()?;
    for (original_is_reverse, expected) in [(false, "TCAG_GACT"), (true, "GCTT_TCAG")] {
        let bam = side_bam(side_fragment(original_is_reverse)?)?;
        let output_dir = TempDir::new()?;
        let mut config = side_config(&bam.bam, output_dir.path(), 4, 4);
        config.side = EndSideSelection::ThreePrime;
        config.source_inside = KmerSource::Reference;
        config.ref_2bit = Some(reference.path.clone());
        let result = run_ends(&config, RunOptions::new_quiet())?;
        assert_global_counts(
            &load_ends_output(result.motif_output_path)?,
            &[(expected, 1.0)],
        );
    }
    Ok(())
}

#[test]
fn side_unpaired_full_fragment_reads_use_their_own_strand() -> Result<()> {
    // A full-fragment read has the same stored ends ACGA...TGCC on either alignment strand
    // Forward molecular ends are ACGA/TGCC, reverse molecular ends are GGCA/TCGT
    for original_is_reverse in [false, true] {
        let mut read = paired_fragment(&PairedFragmentSpec::new(0, 40, 40, 40))?.forward;
        read.seq[..4].copy_from_slice(b"ACGA");
        read.seq[36..].copy_from_slice(b"TGCC");
        read.is_reverse = original_is_reverse;
        read.flags = 0;
        let bam = bam_from_fragments(
            "side_single",
            vec![("chr1".into(), 200)],
            vec![],
            vec![read],
        )?;
        for (side, expected) in [
            (
                EndSideSelection::FivePrime,
                if original_is_reverse {
                    "_GGCA"
                } else {
                    "_ACGA"
                },
            ),
            (
                EndSideSelection::ThreePrime,
                if original_is_reverse {
                    "_TCGT"
                } else {
                    "_TGCC"
                },
            ),
        ] {
            let output_dir = TempDir::new()?;
            let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
            config.unpaired.reads_are_fragments = true;
            config.side = side;
            let result = run_ends(&config, RunOptions::new_quiet())?;
            assert_global_counts(
                &load_ends_output(result.motif_output_path)?,
                &[(expected, 1.0)],
            );
            assert_eq!(result.counters.counted_motifs, 1);
        }
    }
    Ok(())
}

#[test]
fn side_specific_counting_errors_on_ambiguous_mate_roles_but_both_keeps_counting() -> Result<()> {
    for (left_role, right_role) in [
        (0, 0x80),
        (0x40, 0),
        (0x40, 0x40),
        (0x80, 0x80),
        (0xc0, 0x80),
    ] {
        let mut fragment = side_fragment(false)?;
        fragment.forward.flags = (fragment.forward.flags & !0xc0) | left_role;
        fragment.reverse.flags = (fragment.reverse.flags & !0xc0) | right_role;
        let bam = side_bam(fragment)?;
        for side in [
            EndSideSelection::Both,
            EndSideSelection::FivePrime,
            EndSideSelection::ThreePrime,
        ] {
            let output_dir = TempDir::new()?;
            let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
            config.side = side;
            let result = run_ends(&config, RunOptions::new_quiet());
            if side == EndSideSelection::Both {
                assert_global_counts(
                    &load_ends_output(result?.motif_output_path)?,
                    &[("_ACGA", 1.0), ("_GGCA", 1.0)],
                );
            } else {
                let error = format!("{:#}", result.expect_err("ambiguous strand must fail"));
                assert!(
                    error.contains("requires exactly one R1 and one R2"),
                    "{error}"
                );
                assert!(error.contains(side.as_str()), "{error}");
            }
        }
    }
    Ok(())
}

#[test]
fn side_selected_and_grouped_motifs_match_the_three_prime_orientation() -> Result<()> {
    // The reverse original's 3' end is GCTT_TCGT. Its uncomplemented AAGC_ACGA
    // must not count, even when listed as a separate motif or target group
    let bam = side_bam(side_fragment(true)?)?;
    let reference = side_reference()?;
    for grouped in [false, true] {
        for dense in [false, true] {
            let output_dir = TempDir::new()?;
            let motifs_file = output_dir.path().join("motifs.tsv");
            fs::write(
                &motifs_file,
                if grouped {
                    "GCTT_TCGT\tkept\nAAGC_ACGA\tother\n"
                } else {
                    "GCTT_TCGT\nAAGC_ACGA\n"
                },
            )?;
            let mut config = side_config(&bam.bam, output_dir.path(), 4, 4);
            config.side = EndSideSelection::ThreePrime;
            config.ref_2bit = Some(reference.path.clone());
            config.motifs_file = Some(motifs_file);
            config.all_motifs = dense;
            let result = run_ends(&config, RunOptions::new_quiet())?;
            let output = load_ends_output(result.motif_output_path)?;
            assert_global_counts(
                &output,
                &[(if grouped { "kept" } else { "GCTT_TCGT" }, 1.0)],
            );
            assert_eq!(output.motif_count(), if dense { 2 } else { 1 });
            assert_eq!(result.counters.counted_motifs, 1);
        }
    }
    Ok(())
}

#[test]
fn side_endpoint_assignment_uses_selected_coordinates_across_tiles_and_groups() -> Result<()> {
    // The left endpoint is base 40 and the right endpoint is base 79
    // For bins of size 50 these belong to rows 0 and 1. The same row order is
    // used by the alphabetically named BED groups left/right
    for original_is_reverse in [false, true] {
        let bam = side_bam(side_fragment(original_is_reverse)?)?;
        for grouped in [false, true] {
            let output_dir = TempDir::new()?;
            let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
            config.side = EndSideSelection::ThreePrime;
            if grouped {
                let windows = output_dir.path().join("sites.bed");
                write_bed4(
                    &windows,
                    &[
                        Bed4Row::new("chr1", 40, 41, "left"),
                        Bed4Row::new("chr1", 79, 80, "right"),
                    ],
                )?;
                config.windows.by_grouped_bed = Some(windows);
            } else {
                config.windows.by_size = Some(50);
            }
            let result = run_ends(&config, RunOptions::new_quiet())?;
            let output = load_ends_output(result.motif_output_path)?;
            assert_eq!(output.row_count(), if grouped { 2 } else { 4 });
            let motif = if original_is_reverse {
                "_TCGT"
            } else {
                "_TGCC"
            };
            let selected_row = if original_is_reverse { 0 } else { 1 };
            assert_eq!(output.motif_labels(), &[motif]);
            for row in 0..output.row_count() {
                assert_eq!(
                    output.count_for_motif(row, motif)?,
                    Some(if row == selected_row { 1.0 } else { 0.0 })
                );
            }
            assert_eq!(result.counters.base.counted_fragments, 1);
            assert_eq!(result.counters.counted_motifs, 1);
        }
    }
    Ok(())
}

#[test]
fn side_selection_preserves_fragment_wide_base_quality_filters() -> Result<()> {
    // Unselected R1 has quality 10 and selected R2 has quality 40
    // The full-fragment mean is 25, so a fragment threshold of 30 still rejects it
    // An end-only threshold of 30 removes R1 but retains the requested R2 end
    let mut fragment = side_fragment(false)?;
    fragment.forward.base_quality = 10;
    let bam = side_bam(fragment)?;
    for (filter, expected_count) in [("mean in fragment >= 30", 0), ("mean in end >= 30", 1)] {
        let output_dir = TempDir::new()?;
        let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
        config.side = EndSideSelection::ThreePrime;
        config.bq_filter = vec![filter.parse().expect("valid base quality filter")];
        let result = run_ends(&config, RunOptions::new_quiet())?;
        let output = load_ends_output(result.motif_output_path)?;
        assert_global_counts(
            &output,
            if expected_count == 0 {
                &[]
            } else {
                &[("_TGCC", 1.0)]
            },
        );
        assert_eq!(result.counters.base.counted_fragments, expected_count);
        assert_eq!(result.counters.counted_motifs, expected_count);
    }
    Ok(())
}

#[test]
fn side_selection_respects_selected_end_clipping() -> Result<()> {
    for affected_end_is_right in [false, true] {
        let mut fragment = side_fragment(false)?;
        if affected_end_is_right {
            fragment.reverse.cigar = vec![Cigar::Match(16), Cigar::SoftClip(4)];
        } else {
            fragment.forward.cigar = vec![Cigar::SoftClip(4), Cigar::Match(16)];
        }
        let bam = side_bam(fragment)?;
        let output_dir = TempDir::new()?;
        let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
        config.side = EndSideSelection::ThreePrime;
        let result = run_ends(&config, RunOptions::new_quiet())?;
        assert_global_counts(
            &load_ends_output(result.motif_output_path)?,
            if affected_end_is_right {
                &[]
            } else {
                &[("_TGCC", 1.0)]
            },
        );
        assert_eq!(
            result.counters.counted_motifs,
            if affected_end_is_right { 0 } else { 1 }
        );
    }

    Ok(())
}

#[test]
fn side_selection_preserves_fragment_wide_indel_filters() -> Result<()> {
    // The unselected left end has an insertion within its first four bases
    // skip-affected-fragment must still reject the whole fragment
    let mut fragment = side_fragment(false)?;
    fragment.forward.cigar = vec![Cigar::Match(3), Cigar::Ins(1), Cigar::Match(16)];
    let bam = side_bam(fragment)?;
    let output_dir = TempDir::new()?;
    let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
    config.side = EndSideSelection::ThreePrime;
    config.indel_filter = IndelMotifFilterPolicy::SkipAffectedFragment;
    let result = run_ends(&config, RunOptions::new_quiet())?;
    assert_global_counts(&load_ends_output(result.motif_output_path)?, &[]);
    assert_eq!(result.counters.counted_motifs, 0);
    Ok(())
}

#[test]
fn side_selected_end_blacklist_validation_does_not_count_the_other_end() -> Result<()> {
    let bam = side_bam(side_fragment(false)?)?;
    let reference = side_reference()?;
    for (masked_start, expected_count) in [(40, 1), (79, 0)] {
        let output_dir = TempDir::new()?;
        let blacklist = output_dir.path().join("blacklist.bed");
        write_bed4(
            &blacklist,
            &[Bed4Row::new("chr1", masked_start, masked_start + 1, "mask")],
        )?;
        let mut config = side_config(&bam.bam, output_dir.path(), 4, 0);
        config.side = EndSideSelection::ThreePrime;
        config.ref_2bit = Some(reference.path.clone());
        config.blacklist = Some(vec![blacklist]);
        config.blacklist_strategy = BlacklistStrategy::All;
        let result = run_ends(&config, RunOptions::new_quiet())?;
        assert_global_counts(
            &load_ends_output(result.motif_output_path)?,
            if expected_count == 0 {
                &[]
            } else {
                &[("_TGCC", 1.0)]
            },
        );
        assert_eq!(result.counters.counted_motifs, expected_count);
    }
    Ok(())
}

#[cfg(feature = "cli")]
#[test]
fn side_cli_accepts_only_documented_values_and_defaults_to_both() -> Result<()> {
    let bam = side_bam(side_fragment(false)?)?;
    for side in [
        None,
        Some("both"),
        Some("5p"),
        Some("3p"),
        Some("left"),
        Some("right"),
        Some("3P"),
    ] {
        let output_dir = TempDir::new()?;
        let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_cfdna"));
        command
            .arg("ends")
            .arg("--bam")
            .arg(&bam.bam)
            .arg("--output-dir")
            .arg(output_dir.path())
            .args([
                "--output-prefix",
                "sample",
                "--chromosomes",
                "chr1",
                "--k-inside",
                "4",
                "--k-outside",
                "0",
                "--n-threads",
                "1",
            ]);
        if let Some(side) = side {
            command.args(["--side", side]);
        }
        let output = command.output()?;
        if matches!(side, Some("left" | "right" | "3P")) {
            assert!(!output.status.success());
            assert!(String::from_utf8_lossy(&output.stderr).contains("invalid value"));
        } else {
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            let counts = load_ends_output(output_dir.path().join("sample.end_motifs.zarr"))?;
            match side {
                Some("5p") => assert_global_counts(&counts, &[("_ACGA", 1.0)]),
                Some("3p") => assert_global_counts(&counts, &[("_TGCC", 1.0)]),
                _ => assert_global_counts(&counts, &[("_ACGA", 1.0), ("_GGCA", 1.0)]),
            }
        }
    }
    Ok(())
}

#[test]
fn side_selected_long_motifs_use_original_orientation_in_the_subspace_codec() -> Result<()> {
    // Thirty inside bases exceed the full radix-5 codec limit
    // The reverse original's 3' inside reference is RC(CTGA + 26 A) = 26 T + TCAG
    // Its outside half is RC(AAGC) = GCTT. The uncomplemented label must not count
    let bam = side_bam(side_fragment(true)?)?;
    let reference = side_reference()?;
    let output_dir = TempDir::new()?;
    let expected = format!("GCTT_{}TCAG", "T".repeat(26));
    let uncomplemented = format!("AAGC_CTGA{}", "A".repeat(26));
    let motifs_file = output_dir.path().join("long_motifs.tsv");
    fs::write(&motifs_file, format!("{expected}\n{uncomplemented}\n"))?;
    let mut config = side_config(&bam.bam, output_dir.path(), 30, 4);
    config.side = EndSideSelection::ThreePrime;
    config.source_inside = KmerSource::Reference;
    config.ref_2bit = Some(reference.path.clone());
    config.motifs_file = Some(motifs_file);
    let result = run_ends(&config, RunOptions::new_quiet())?;
    assert_global_counts(
        &load_ends_output(result.motif_output_path)?,
        &[(expected.as_str(), 1.0)],
    );
    assert_eq!(result.counters.counted_motifs, 1);
    Ok(())
}
