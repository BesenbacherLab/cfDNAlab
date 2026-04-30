#![cfg(feature = "cmd_fragment_kmers")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ApplyGCArgs, ChromosomeArgs, IOCArgs, WindowsArgs};
use cfdnalab::commands::fragment_kmers::{config::FragmentKmersConfig, fragment_kmers::run};
use cfdnalab::shared::io::dot_join;
use fixtures::{
    ReadSpec, bam_from_specs, late_origin_gc_reference_sequence, simple_inward_bam,
    simple_reference_twobit, twobit_from_sequences, write_bed, write_two_bin_gc_package,
};
use ndarray::Array2;
use ndarray_npy::read_npy;
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|chr| chr.to_string()).collect()),
        chromosomes_file: None,
    }
}

#[test]
fn bed_windowed_runs_write_prefixed_bins_tsv_with_exact_blacklisted_fractions() -> Result<()> {
    // Arrange:
    // - `simple_inward_bam()` gives one 60 bp fragment on chr1 spanning [20,80).
    // - The two BED windows are [10,20) and [20,30).
    // - The blacklist interval [15,20) overlaps only the first window for 5 of its 10 bases.
    // - With a non-empty output prefix, the bins metadata should follow the same prefixed filename
    //   contract as the primary count outputs.
    let bam = simple_inward_bam()?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 20, "left"), ("chr1", 20, 30, "right")],
    )?;
    write_bed(&blacklist_bed, &[("chr1", 15, 20, "masked")])?;

    let mut cfg = FragmentKmersConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("sampleA".to_string());
    cfg.set_kmer_sizes(vec![1]);
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Act
    run(&cfg)?;
    let bins_tsv =
        std::fs::read_to_string(out_dir.path().join(dot_join(&["sampleA", "bins.tsv"])))?;

    // Assert
    assert_eq!(
        bins_tsv,
        concat!(
            "chrom\tstart\tend\tblacklisted_fraction\n",
            "chr1\t10\t20\t0.5\n",
            "chr1\t20\t30\t0\n"
        )
    );
    assert!(!out_dir.path().join("bins.tsv").exists());
    Ok(())
}

#[test]
fn gc_file_late_tile_window_uses_reference_coordinates_after_fetch_narrowing() -> Result<()> {
    // Arrange:
    // - The one unpaired fragment spans [900,961), while the BED window starts far from tile
    //   origin 0.
    // - The reference is shorter than the BAM chromosome, but long enough for the narrowed
    //   window-derived fetch span. Reading the full tile reference would overrun the reference.
    // - The correct fragment interval [900,961) is all C, so it lands in the high-GC correction
    //   bin with weight 7.0. Using prefix-local origin 0 would see A-only sequence instead.
    // - With k=1 over the selected reference span, the 61 selected bases are all C.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1_500)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 900,
            cigar: vec![('M', 61)],
            seq: vec![b'A'; 61],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        "fragment_kmers_late_tile_gc_origin",
    )?;
    let reference = twobit_from_sequences(
        "fragment_kmers_late_tile_gc_origin_ref",
        vec![("chr1".to_string(), late_origin_gc_reference_sequence())],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("late_window.bed");
    let gc_path = out_dir.path().join("two_bin_gc_package.npz");
    write_bed(&bed_path, &[("chr1", 900, 961, "late")])?;
    write_two_bin_gc_package(&gc_path, 61, 2.0, 7.0)?;

    let mut cfg = FragmentKmersConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("kmers".to_string());
    cfg.set_kmer_sizes(vec![1]);
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
    });
    cfg.set_min_mapq(0);
    cfg.shared_args.unpaired.reads_are_fragments = true;
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 61;
        frag.max_fragment_length = 61;
    }
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });

    // Act
    run(&cfg)?;
    let counts: Array2<f64> = read_npy(out_dir.path().join("kmers.k1_counts.npy"))?;
    let motifs = std::fs::read_to_string(out_dir.path().join("kmers.k1_motifs.txt"))?;

    // Assert
    assert_eq!(counts.shape(), &[1, 4]);
    for (motif, expected) in [("A", 0.0), ("C", 427.0), ("G", 0.0), ("T", 0.0)] {
        let column = motifs
            .lines()
            .position(|observed| observed == motif)
            .expect("motif should be present");
        assert_eq!(
            counts[(0, column)],
            expected,
            "unexpected count for {motif}"
        );
    }
    assert_eq!(counts.sum(), 427.0);
    Ok(())
}

// #![cfg(feature = "cmd_fragment_kmers_tests")]

// mod fixtures;

// mod tests_fragment_kmer_command {
//     use std::collections::HashMap;
//     use std::path::Path;

//     use crate::fixtures::{
//         FragmentSpec, ReadSpec, bam_from_specs, fragment_kmers_edge_bam,
//         fragment_kmers_edge_reference, simple_inward_bam, simple_reference_twobit,
//         single_position_selection, twobit_from_sequences, write_bed, write_scaling_factors,
//     };
//     use anyhow::{Context, Result, bail};
//     use cfdnalab::commands::cli_common::{
//         ChromosomeArgs, FragmentLengthArgs, IOCArgs, Ref2BitRequiredArgs, ScaleGenomeArgs,
//         WindowsArgs,
//     };
//     use cfdnalab::commands::fragment_kmers::config::FragmentKmersConfig;
//     use cfdnalab::commands::fragment_kmers::fragment_kmers::run;
//     use cfdnalab::commands::visualize_positions::ReferenceFrame;
//     use cfdnalab::shared::base::make_canonical;
//     use cfdnalab::shared::blacklist::BlacklistStrategy;
//     use cfdnalab::shared::indel_mode::IndelMode;
//     use ndarray::{Array2, Array3};
//     use ndarray_npy::read_npy;
//     use tempfile::TempDir;

//     pub(crate) fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
//         ChromosomeArgs {
//             chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
//             chromosomes_file: None,
//         }
//     }

//     pub(crate) fn build_revcomp_assets(
//         left: &str,
//     ) -> Result<(
//         crate::fixtures::TwoBitFixture,
//         crate::fixtures::BamFixture,
//         crate::fixtures::TwoBitFixture,
//         crate::fixtures::BamFixture,
//         u32,
//         u32,
//     )> {
//         let right = revcomp(left);
//         let combined = format!("{left}{right}");

//         let left_ref =
//             twobit_from_sequences("revcomp_left", vec![("chr1".to_string(), left.to_string())])?;
//         let combined_ref =
//             twobit_from_sequences("revcomp_combined", vec![("chr1".to_string(), combined)])?;

//         let left_len = left.len() as u32;
//         if left_len % 2 != 0 {
//             bail!("left reference length must be even to fold symmetrically");
//         }
//         let read_len_left = left_len / 2;
//         assert!(
//             read_len_left > 0 && read_len_left * 2 == left_len,
//             "left reference length must be even and >= 2"
//         );

//         let left_fragment = make_fragment_pair(0, 0, read_len_left);
//         let left_bam = bam_from_specs(
//             vec![("chr1".to_string(), left_len)],
//             vec![left_fragment],
//             Vec::new(),
//             "revcomp_left",
//         )?;

//         let read_len_combined = left_len;
//         let combined_fragments = vec![make_fragment_pair(0, 0, read_len_combined)];
//         let combined_bam = bam_from_specs(
//             vec![("chr1".to_string(), left_len * 2)],
//             combined_fragments,
//             Vec::new(),
//             "revcomp_combined",
//         )?;

//         Ok((
//             left_ref,
//             left_bam,
//             combined_ref,
//             combined_bam,
//             left_len,
//             left_len * 2,
//         ))
//     }

//     pub(crate) fn revcomp(seq: &str) -> String {
//         seq.chars()
//             .rev()
//             .map(|c| match c {
//                 'A' | 'a' => 'T',
//                 'C' | 'c' => 'G',
//                 'G' | 'g' => 'C',
//                 'T' | 't' => 'A',
//                 other => panic!("unexpected base {other}"),
//             })
//             .collect()
//     }

//     pub(crate) fn make_fragment_pair(tid: usize, start: i64, read_len: u32) -> FragmentSpec {
//         const FLAG_FIRST_MATE: u16 = 0x40;
//         const FLAG_SECOND_MATE: u16 = 0x80;
//         const FLAG_PROPER_PAIR: u16 = 0x2;
//         const FLAG_MATE_REVERSE: u16 = 0x20;

//         let insert_size = (read_len as i64) * 2;
//         FragmentSpec {
//             forward: ReadSpec {
//                 tid,
//                 pos: start,
//                 cigar: vec![('M', read_len)],
//                 seq: vec![b'A'; read_len as usize],
//                 qual: 40,
//                 is_reverse: false,
//                 mapq: 60,
//                 flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
//                 mate_tid: Some(tid),
//                 mate_pos: Some(start + read_len as i64),
//                 insert_size,
//             },
//             reverse: ReadSpec {
//                 tid,
//                 pos: start + read_len as i64,
//                 cigar: vec![('M', read_len)],
//                 seq: vec![b'T'; read_len as usize],
//                 qual: 40,
//                 is_reverse: true,
//                 mapq: 60,
//                 flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
//                 mate_tid: Some(tid),
//                 mate_pos: Some(start),
//                 insert_size: -insert_size,
//             },
//         }
//     }

//     pub(crate) fn manual_kmer_counts(seq: &str, k: usize) -> HashMap<String, f64> {
//         let mut counts = HashMap::new();
//         if seq.len() < k {
//             return counts;
//         }
//         for idx in 0..=(seq.len() - k) {
//             let motif = &seq[idx..idx + k];
//             *counts.entry(motif.to_string()).or_insert(0.0) += 1.0;
//         }
//         counts
//     }

//     pub(crate) fn manual_offset_counts(
//         seq: &str,
//         k: usize,
//         offsets: &[usize],
//     ) -> HashMap<String, f64> {
//         let mut counts = HashMap::new();
//         for &offset in offsets {
//             let motif = &seq[offset..offset + k];
//             *counts.entry(motif.to_string()).or_insert(0.0) += 1.0;
//         }
//         counts
//     }

//     pub(crate) fn load_positional_group_counts(
//         dir: &Path,
//         prefix: &str,
//         k: u8,
//         group: &str,
//     ) -> Result<HashMap<String, f64>> {
//         let counts_path = dir.join(format!("{prefix}.k{k}_{group}_counts.npy"));
//         let motifs_path = dir.join(format!("{prefix}.k{k}_{group}_motifs.txt"));
//         let counts: Array3<f64> = read_npy(&counts_path)?;
//         let motifs: Vec<String> = std::fs::read_to_string(&motifs_path)?
//             .lines()
//             .map(|line| line.to_string())
//             .collect();
//         let mut totals = vec![0.0f64; motifs.len()];
//         for window_idx in 0..counts.shape()[0] {
//             for pos_idx in 0..counts.shape()[1] {
//                 for motif_idx in 0..counts.shape()[2] {
//                     totals[motif_idx] += counts[[window_idx, pos_idx, motif_idx]];
//                 }
//             }
//         }
//         let mut out = HashMap::new();
//         for (motif, total) in motifs.into_iter().zip(totals.into_iter()) {
//             out.insert(motif, total);
//         }
//         Ok(out)
//     }

//     pub(crate) fn assert_count_map_matches(
//         observed: &HashMap<String, f64>,
//         expected: &HashMap<String, f64>,
//         context: &str,
//     ) {
//         for (motif, exp) in expected {
//             let obs = observed.get(motif).copied().unwrap_or_default();
//             assert!(
//                 (obs - exp).abs() < 1e-6,
//                 "{context}: motif {motif} expected {exp} observed {obs}"
//             );
//         }
//         for (motif, obs) in observed {
//             if !expected.contains_key(motif) {
//                 assert!(
//                     obs.abs() < 1e-6,
//                     "{context}: unexpected non-zero motif {motif}: {obs}"
//                 );
//             }
//         }
//     }

//     #[test]
//     fn counts_dinucleotides_in_global_window() -> Result<()> {
//         let bam = crate::fixtures::simple_inward_bam()?;
//         let reference = crate::fixtures::simple_reference_twobit()?;
//         let out_dir = TempDir::new()?;

//         let mut cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: out_dir.path().to_path_buf(),
//                 n_threads: 2,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         cfg.set_output_prefix("kmers".to_string());
//         cfg.set_kmer_sizes(vec![2]);
//         cfg.set_windows(WindowsArgs::default());
//         cfg.set_indel_mode(IndelMode::Ignore);
//         cfg.set_min_mapq(0);
//         cfg.set_require_proper_pair(false);
//         cfg.set_canonical(false);

//         run(&cfg)?;

//         let counts_path = out_dir.path().join("kmers.k2_counts.npy");
//         let motifs_path = out_dir.path().join("kmers.k2_motifs.txt");
//         assert!(counts_path.exists());
//         assert!(motifs_path.exists());

//         let counts: Array2<f64> = read_npy(&counts_path)?;
//         assert_eq!(counts.shape(), &[1, 16]);
//         let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
//             .lines()
//             .map(|s| s.to_string())
//             .collect();

//         let chr1_seq = reference
//             .sequence("chr1")
//             .context("missing chr1 sequence in reference fixture")?;
//         let start = 20usize;
//         let end = 80usize; // fragment end (exclusive)
//         let k = 2usize;
//         let mut expected: HashMap<String, f64> = HashMap::new();
//         for idx in start..=(end - k) {
//             let motif = chr1_seq
//                 .get(idx..idx + k)
//                 .context("motif slice")?
//                 .to_string();
//             *expected.entry(motif).or_insert(0.0) += 1.0;
//         }

//         let row = counts.row(0);
//         let total: f64 = row.sum();
//         assert!((total - 59.0).abs() < 1e-6);

//         for (col, motif) in motif_list.iter().enumerate() {
//             let expected_val = expected.get(motif).copied().unwrap_or(0.0);
//             assert!(
//                 (row[col] - expected_val).abs() < 1e-6,
//                 "motif {motif} expected {expected_val} observed {}",
//                 row[col]
//             );
//         }

//         Ok(())
//     }

//     #[test]
//     fn positional_counts_restricts_starts() -> Result<()> {
//         let bam = simple_inward_bam()?;
//         let reference = simple_reference_twobit()?;
//         let out_dir = TempDir::new()?;

//         let mut cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: out_dir.path().to_path_buf(),
//                 n_threads: 2,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         cfg.set_output_prefix("kmers_first".to_string());
//         cfg.set_kmer_sizes(vec![2]);
//         cfg.set_windows(WindowsArgs::default());
//         cfg.set_indel_mode(IndelMode::Ignore);
//         cfg.set_min_mapq(0);
//         cfg.set_require_proper_pair(false);
//         cfg.set_canonical(false);
//         cfg.set_positional_counts(true);
//         cfg.set_position_selection(single_position_selection(ReferenceFrame::Left, "1..1", 1));

//         run(&cfg)?;

//         let counts_path = out_dir.path().join("kmers_first.k2_left_counts.npy");
//         let motifs_path = out_dir.path().join("kmers_first.k2_left_motifs.txt");
//         let positions_path = out_dir.path().join("kmers_first.left_positions.txt");
//         assert!(counts_path.exists());
//         assert!(motifs_path.exists());
//         assert!(positions_path.exists());

//         let counts: Array3<f64> = read_npy(&counts_path)?;
//         assert_eq!(counts.shape(), &[1, 1, 16]);
//         let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
//             .lines()
//             .map(|s| s.to_string())
//             .collect();
//         let positions: Vec<i32> = std::fs::read_to_string(&positions_path)?
//             .lines()
//             .map(|line| line.parse::<i32>().expect("position"))
//             .collect();
//         assert_eq!(positions, vec![0]);

//         let chr1_seq = reference
//             .sequence("chr1")
//             .context("missing chr1 sequence in reference fixture")?;
//         let start = 20usize;
//         let k = 2usize;
//         let motif = chr1_seq
//             .get(start..start + k)
//             .context("motif slice")?
//             .to_string();

//         let mut expected: HashMap<String, f64> = HashMap::new();
//         expected.insert(motif.clone(), 1.0);

//         let mut actual: HashMap<String, f64> = HashMap::new();
//         for (motif_idx, motif) in motif_list.iter().enumerate() {
//             let value = counts[[0, 0, motif_idx]];
//             if value != 0.0 {
//                 actual.insert(motif.clone(), value);
//             }
//         }

//         assert_counts_close(&actual, &expected);
//         Ok(())
//     }

//     #[test]
//     fn canonical_trimers_collapse_matches_manual_counts() -> Result<()> {
//         let bam = simple_inward_bam()?;
//         let reference = simple_reference_twobit()?;
//         let out_dir = TempDir::new()?;

//         let mut cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: out_dir.path().to_path_buf(),
//                 n_threads: 2,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         cfg.set_output_prefix("kmers".to_string());
//         cfg.set_kmer_sizes(vec![3]);
//         cfg.set_windows(WindowsArgs::default());
//         cfg.set_indel_mode(IndelMode::Ignore);
//         cfg.set_min_mapq(0);
//         cfg.set_require_proper_pair(false);
//         cfg.set_canonical(true);

//         run(&cfg)?;

//         let counts_path = out_dir.path().join("kmers.k3_counts.npy");
//         let motifs_path = out_dir.path().join("kmers.k3_motifs.txt");
//         assert!(counts_path.exists());
//         assert!(motifs_path.exists());

//         let counts: Array2<f64> = read_npy(&counts_path)?;
//         assert_eq!(counts.shape()[0], 1);
//         let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
//             .lines()
//             .map(|s| s.to_string())
//             .collect();

//         let chr1_seq = reference
//             .sequence("chr1")
//             .context("missing chr1 sequence in reference fixture")?;
//         let start = 20usize;
//         let end = 80usize;
//         let k = 3usize;
//         let mut expected: HashMap<String, f64> = HashMap::new();
//         for idx in start..=(end - k) {
//             let motif = chr1_seq
//                 .get(idx..idx + k)
//                 .context("motif slice")?
//                 .to_string();
//             let canon = make_canonical(motif);
//             *expected.entry(canon).or_insert(0.0) += 1.0;
//         }

//         let row = counts.row(0);
//         let total: f64 = row.sum();
//         assert!((total - 58.0).abs() < 1e-6);

//         for motif in &motif_list {
//             assert_eq!(motif, &make_canonical(motif.clone()));
//         }

//         for (col, motif) in motif_list.iter().enumerate() {
//             let expected_val = expected.get(motif).copied().unwrap_or(0.0);
//             assert!(
//                 (row[col] - expected_val).abs() < 1e-6,
//                 "motif {motif} expected {expected_val} observed {}",
//                 row[col]
//             );
//         }

//         Ok(())
//     }

//     pub(crate) fn load_counts_from_output(
//         dir: &Path,
//         prefix: &str,
//         k: u8,
//     ) -> Result<HashMap<String, f64>> {
//         let dense_path = dir.join(format!("{prefix}.k{k}_counts.npy"));
//         if dense_path.exists() {
//             let motifs_path = dir.join(format!("{prefix}.k{k}_motifs.txt"));
//             let counts: Array2<f64> = read_npy(&dense_path)?;
//             assert_eq!(
//                 counts.shape()[0],
//                 1,
//                 "counts matrix should have one window row"
//             );
//             let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
//                 .lines()
//                 .map(|s| s.to_string())
//                 .collect();
//             let mut out = HashMap::new();
//             for (idx, motif) in motif_list.iter().enumerate() {
//                 out.insert(motif.clone(), counts[(0, idx)]);
//             }
//             return Ok(out);
//         }

//         // Positional output is split per-group (left/right/mid). Aggregate counts over windows and positions.
//         let mut aggregates: HashMap<String, f64> = HashMap::new();
//         let groups = ["left", "right", "mid"];
//         for group in groups {
//             let counts_path = dir.join(format!("{prefix}.k{k}_{group}_counts.npy"));
//             if !counts_path.exists() {
//                 continue;
//             }
//             let motifs_path = dir.join(format!("{prefix}.k{k}_{group}_motifs.txt"));
//             let counts: Array3<f64> = read_npy(&counts_path)?;
//             let motif_list: Vec<String> = std::fs::read_to_string(&motifs_path)?
//                 .lines()
//                 .map(|s| s.to_string())
//                 .collect();

//             let mut totals = vec![0.0f64; motif_list.len()];
//             for window_idx in 0..counts.shape()[0] {
//                 for pos_idx in 0..counts.shape()[1] {
//                     for motif_idx in 0..counts.shape()[2] {
//                         totals[motif_idx] += counts[(window_idx, pos_idx, motif_idx)];
//                     }
//                 }
//             }

//             for (motif, total) in motif_list.iter().zip(totals.into_iter()) {
//                 *aggregates.entry(motif.clone()).or_insert(0.0) += total;
//             }
//         }

//         if aggregates.is_empty() {
//             bail!(
//                 "no counts files found for prefix '{}' and k {} in {}",
//                 prefix,
//                 k,
//                 dir.display()
//             );
//         }

//         Ok(aggregates)
//     }

//     fn assert_counts_close(actual: &HashMap<String, f64>, expected: &HashMap<String, f64>) {
//         for (motif, exp) in expected {
//             let obs = actual.get(motif).copied().unwrap_or(0.0);
//             assert!(
//                 (obs - exp).abs() < 1e-6,
//                 "motif {motif} expected {exp} observed {obs}"
//             );
//         }
//         for (motif, obs) in actual {
//             let exp = expected.get(motif).copied().unwrap_or(0.0);
//             assert!(
//                 (obs - exp).abs() < 1e-6,
//                 "motif {motif} expected {exp} observed {obs}"
//             );
//         }
//     }

//     #[test]
//     fn complex_edge_cases_left_frame_respect_scaling_and_blacklists() -> Result<()> {
//         let bam = fragment_kmers_edge_bam()?;
//         let reference = fragment_kmers_edge_reference()?;
//         let out_dir = TempDir::new()?;
//         let chromosomes = ["chr1"];

//         let fragment_lengths = FragmentLengthArgs {
//             min_fragment_length: 10,
//             max_fragment_length: 1000,
//         };

//         let mut cfg_base = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: out_dir.path().to_path_buf(),
//                 n_threads: 2,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&chromosomes),
//         );
//         cfg_base.set_output_prefix("edge_base".to_string());
//         cfg_base.set_kmer_sizes(vec![2]);
//         cfg_base.set_indel_mode(IndelMode::Adjust);
//         cfg_base.set_min_mapq(0);
//         cfg_base.set_require_proper_pair(false);
//         cfg_base.set_canonical(false);
//         cfg_base.set_ignore_gap(true);
//         cfg_base.set_positional_counts(true);
//         cfg_base.set_position_selection(single_position_selection(ReferenceFrame::Left, "2..", 1));
//         {
//             let fl = cfg_base.fragment_lengths_mut();
//             fl.min_fragment_length = fragment_lengths.min_fragment_length;
//             fl.max_fragment_length = fragment_lengths.max_fragment_length;
//         }

//         run(&cfg_base)?;

//         // Explaining expectations:
//         //
//         // Reference:
//         // 0 A 1 C 2 G 3 T 4 G 5 A 6 C 7 C 8 T 9 T
//         // 10 A 11 G 12 G 13 C 14 T 15 A 16 A 17 C 18 C 19 G
//         // 20 T 21 A 22 C 23 G 24 T 25 T 26 A 27 G 28 C 29 C
//         // 30 G 31 A 32 T 33 T 34 A 35 C 36 A 37 A 38 G 39 T
//         //
//         // frame = Left, positions = “2..”, so skip the first base of each fragment (offset 0) and count forward 2-mers only.
//         // insertions/deletions split segments. A 2-mer must be fully inside a single contiguous segment (no crossing the I/D boundary)
//         //
//         // # Fragment 1
//         //
//         // forward: 0..10 (10M)
//         //
//         // reverse: 14..24 (10M)
//         //
//         // segments: [0..10) and [14..24)
//         //
//         // Skip first left base (abs 0); k=2 allowed starts:
//         //
//         // [0..10): starts 1..=8 -> CG, GT, TG, GA, AC, CC, CT, TT
//         //
//         // [14..24): starts 14..=22 -> TA, AA, AC, CC, CG, GT, TA, AC, CG
//         //
//         // Counts from F1
//         // AA 1, AC 3, CC 2, CG 3, CT 1, GA 1, GT 2, TA 2, TG 1, TT 1
//         //
//         // # Fragment 2 (has 4M 1I 4M on forward)
//         //
//         // start..end: 5..21
//         //
//         // forward read splits the reference into [5..9) and [9..13) because of the insertion (I consumes read, not reference).
//         //
//         // reverse read adds [13..21).
//         //
//         // segments: [5..9), [9..13), [13..21)
//         //
//         // Skip first left base (abs 5). k=2 allowed starts:
//         //
//         // [5..9): 6..=7 -> CC, CT (note: TT at 8 is NOT allowed, last start = 7)
//         //
//         // [9..13): 9..=11 -> TA, AG, GG (GC at 12 is NOT allowed, last start = 11)
//         //
//         // [13..21): 13..=19 -> CT, TA, AA, AC, CC, CG, GT
//         //
//         // Counts from F2
//         // AA 1, AC 1, CC 2, CG 1, CT 2, AG 1, GG 1, GT 1, TA 2
//         //
//         // # Fragment 3 (has 3M 1D 5M on forward)
//         //
//         // start..end: 16..27
//         //
//         // Deletion consumes reference -> gap at [19..20).
//         //
//         // reverse read 20..27.
//         //
//         // segments: [16..19) and [20..27)
//         //
//         // Skip first left base (abs 16). k=2 allowed starts:
//         //
//         // [16..19): 17 -> CC
//         //
//         // [20..27): 20..=25 -> TA, AC, CG, GT, TT, TA
//         //
//         // Counts from F3
//         // AC 1, CC 1, CG 1, GT 1, TA 2, TT 1

//         let observed_base = load_counts_from_output(out_dir.path(), "edge_base", 2)?;
//         println!("{:?}", observed_base);
//         let expected_base: HashMap<String, f64> = vec![
//             ("AA", 2.0),
//             ("AC", 5.0),
//             ("AG", 1.0),
//             ("CC", 5.0),
//             ("CG", 5.0),
//             ("CT", 3.0),
//             ("GA", 1.0),
//             ("GG", 1.0),
//             ("GT", 4.0),
//             ("TA", 6.0),
//             ("TG", 1.0),
//             ("TT", 2.0),
//         ]
//         .into_iter()
//         .map(|(m, c)| (m.to_string(), c))
//         .collect();
//         // Sanity check: manually computed expectations for edge-case fragments.
//         assert_counts_close(&observed_base, &expected_base);

//         println!("Next setup");

//         let blacklist_path = out_dir.path().join("mask.bed");
//         write_bed(
//             &blacklist_path,
//             &[("chr1", 9, 11, "mask"), ("chr1", 22, 23, "mask")],
//         )?;
//         let scaling_path = out_dir.path().join("scaling.tsv");
//         write_scaling_factors(
//             &scaling_path,
//             &[
//                 ("chr1", 0, 6, 1.0),
//                 ("chr1", 6, 8, 0.0),
//                 ("chr1", 8, 20, 1.5),
//                 ("chr1", 20, 40, 0.5),
//             ],
//         )?;

//         let mut cfg_scaled = cfg_base.clone();
//         cfg_scaled.set_output_prefix("edge_scaled".to_string());
//         cfg_scaled.set_blacklist(Some(vec![blacklist_path.clone()]));
//         cfg_scaled.set_blacklist_strategy(BlacklistStrategy::Proportion(1.0));
//         let mut scale_args = ScaleGenomeArgs::default();
//         scale_args.scaling_factors = Some(scaling_path.clone());
//         cfg_scaled.set_scale_genome(scale_args);

//         run(&cfg_scaled)?;

//         // Explaining expectations:
//         //
//         // Blacklist (N-mask on reference): positions [9,11) => {(8),9,10} and [22,23) => {(21),22}
//         // Scaling: [0,6) -> 1.0, [6,8) -> 0.0 (also N-masked), [8,20) -> 1.5, [20,40) -> 0.5
//         //
//         // A start is valid iff both bases (start and start+1) are not N-masked
//         // Weighting = scaling weight at the start base
//         //
//         // per-fragment contributions (motif -> sum of weights)
//         //
//         // # Fragment 1
//         //
//         // starts kept (0-indexed): 1,2,3,4,5,14,15,16,17,18,19,20
//         //
//         // yields:
//         // CG 1.0(@1)+1.5(@18)=2.5;
//         // GT 1.0(@2)+1.5(@19)=2.5;
//         // TG 1.0(@3);
//         // GA 1.0(@4);
//         // TA 1.5(@14)+0.5(@20)=2.0;
//         // AA 1.5(@15);
//         // AC 1.5(@16);
//         // CC 1.5(@17)
//         //
//         // # Fragment 2
//         //
//         // starts kept: 11,13,14,15,16,17,18,19
//         // (note: 6,7 masked by scaling=0; 9,10 masked by blacklist)
//         //
//         // yields:
//         // GG 1.5(@11); CT 1.5(@13); TA 1.5(@14); AA 1.5(@15);
//         // AC 1.5(@16); CC 1.5(@17); CG 1.5(@18); GT 1.5(@19)
//         //
//         // # Fragment 3
//         //
//         // starts kept: 17,20,23,24,25
//         // (22 excluded due to blacklist)
//         //
//         // yields:
//         // CC 1.5(@17); TA 0.5(@20)+0.5(@25)=1.0; GT 0.5(@23); TT 0.5(@24)

//         let observed_scaled = load_counts_from_output(out_dir.path(), "edge_scaled", 2)?;

//         let expected_scaled: HashMap<String, f64> = vec![
//             ("AA", 3.0),
//             ("AC", 3.0),
//             ("CC", 4.5),
//             ("CG", 4.0),
//             ("CT", 1.5),
//             ("GA", 1.0),
//             ("GG", 1.5),
//             ("GT", 4.5),
//             ("TA", 4.5),
//             ("TG", 1.0),
//             ("TT", 0.5),
//         ]
//         .into_iter()
//         .map(|(m, c)| (m.to_string(), c))
//         .collect();

//         assert_counts_close(&observed_scaled, &expected_scaled);

//         Ok(())
//     }

//     #[test]
//     fn complex_edge_cases_right_frame_respect_scaling_and_blacklists() -> Result<()> {
//         let bam = fragment_kmers_edge_bam()?;
//         let reference = fragment_kmers_edge_reference()?;
//         let out_dir = TempDir::new()?;
//         let chromosomes = ["chr1"];

//         let fragment_lengths = FragmentLengthArgs {
//             min_fragment_length: 10,
//             max_fragment_length: 1000,
//         };

//         let mut cfg_base = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: out_dir.path().to_path_buf(),
//                 n_threads: 2,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&chromosomes),
//         );
//         cfg_base.set_output_prefix("edge_base_right".to_string());
//         cfg_base.set_kmer_sizes(vec![2]);
//         cfg_base.set_indel_mode(IndelMode::Adjust);
//         cfg_base.set_min_mapq(0);
//         cfg_base.set_require_proper_pair(false);
//         cfg_base.set_canonical(false);
//         cfg_base.set_ignore_gap(true);
//         cfg_base.set_positional_counts(true);
//         cfg_base.set_position_selection(single_position_selection(ReferenceFrame::Right, "2..", 1));
//         {
//             let fl = cfg_base.fragment_lengths_mut();
//             fl.min_fragment_length = fragment_lengths.min_fragment_length;
//             fl.max_fragment_length = fragment_lengths.max_fragment_length;
//         }

//         // Base (no blacklist/scaling)
//         run(&cfg_base)?;
//         let observed_base = load_counts_from_output(out_dir.path(), "edge_base_right", 2)?;
//         // Hand-derived (reverse anchoring, terminal-base weight=1.0 everywhere, no masking)
//         let expected_base: HashMap<String, f64> = vec![
//             ("AA", 2.0),
//             ("AC", 3.0),
//             ("AG", 3.0),
//             ("AT", 0.0),
//             ("CA", 1.0),
//             ("CC", 1.0),
//             ("CG", 4.0),
//             ("CT", 1.0),
//             ("GA", 0.0),
//             ("GC", 0.0),
//             ("GG", 5.0),
//             ("GT", 8.0),
//             ("TA", 5.0),
//             ("TC", 1.0),
//             ("TG", 0.0),
//             ("TT", 2.0),
//         ]
//         .into_iter()
//         .map(|(m, c)| (m.to_string(), c))
//         .collect();
//         assert_counts_close(&observed_base, &expected_base);

//         // Blacklist + scaling scenario

//         let blacklist_path = out_dir.path().join("mask.bed");
//         write_bed(
//             &blacklist_path,
//             &[("chr1", 9, 11, "mask"), ("chr1", 22, 23, "mask")],
//         )?;
//         let scaling_path = out_dir.path().join("scaling.tsv");
//         write_scaling_factors(
//             &scaling_path,
//             &[
//                 ("chr1", 0, 6, 1.0),
//                 ("chr1", 6, 8, 0.0), // also N-masked
//                 ("chr1", 8, 20, 1.5),
//                 ("chr1", 20, 40, 0.5),
//             ],
//         )?;

//         let mut cfg_scaled = cfg_base.clone();
//         cfg_scaled.set_output_prefix("edge_scaled_right".to_string());
//         cfg_scaled.set_blacklist(Some(vec![blacklist_path.clone()]));
//         cfg_scaled.set_blacklist_strategy(BlacklistStrategy::Proportion(1.0));
//         let mut scale_args = ScaleGenomeArgs::default();
//         scale_args.scaling_factors = Some(scaling_path.clone());
//         cfg_scaled.set_scale_genome(scale_args);

//         run(&cfg_scaled)?;
//         let observed_scaled = load_counts_from_output(out_dir.path(), "edge_scaled_right", 2)?;

//         // Hand-derived with masking (N at 6,7,9,10,22) and terminal-base weights:
//         // weight(p) = 1.0 for p∈[0,6), 0.0 for p∈[6,8), 1.5 for p∈[8,20), 0.5 for p∈[20,40).
//         // reverse k-mers use terminal base index 'p' for weighting and span [p-1, p].
//         // Consequently, anchors at p=6,7,9,10 are invalid (touch masked bases), p=23 is invalid (touches 22),
//         // but p=21 remains valid (spans 20–21 and does not touch 22).
//         let expected_scaled: HashMap<String, f64> = vec![
//             ("GT", 5.5),
//             ("GG", 4.5),
//             ("CG", 4.0),
//             ("TA", 4.0),
//             ("TT", 3.0),
//             ("AC", 2.0),
//             ("AG", 1.5),
//             ("CC", 1.5),
//             ("CA", 1.0),
//             ("TC", 1.0),
//             ("AA", 0.5),
//         ]
//         .into_iter()
//         .map(|(m, c)| (m.to_string(), c))
//         .collect();
//         assert_counts_close(&observed_scaled, &expected_scaled);

//         Ok(())
//     }
// }

// #[cfg(test)]
// mod tests_fragment_kmers_tiling {
//     use anyhow::Result;
//     use cfdnalab::{
//         commands::fragment_kmers::{
//             positions::PositionGroup,
//             tiling::{TileKmerCountEntry, TileWindowCounts, merge_tile_counts},
//         },
//         shared::kmers::kmer_codec::{KmerSpec, build_kmer_specs},
//     };

//     fn code_for_motif(spec: &KmerSpec, motif: &str) -> u64 {
//         let limit = 5u64.pow(spec.k as u32);
//         for code in 0..limit {
//             if spec.decode_kmer(code) == motif {
//                 return code;
//             }
//         }
//         panic!("motif {} not encodable", motif);
//     }

//     #[test]
//     fn merge_tile_counts_merges_two_tiles() -> Result<()> {
//         let kmer_specs = build_kmer_specs(&[3])?;
//         let spec3 = &kmer_specs[&3];
//         let code_aaa = code_for_motif(spec3, "AAA");

//         let payload_a = vec![TileWindowCounts {
//             original_idx: 0,
//             entries: vec![TileKmerCountEntry {
//                 k: 3,
//                 code: code_aaa,
//                 position: None,
//                 group: PositionGroup::Left,
//                 value: 1.5,
//             }],
//         }];

//         let payload_b = vec![TileWindowCounts {
//             original_idx: 0,
//             entries: vec![TileKmerCountEntry {
//                 k: 3,
//                 code: code_aaa,
//                 position: None,
//                 group: PositionGroup::Left,
//                 value: 2.0,
//             }],
//         }];

//         let merged = merge_tile_counts(vec![payload_a, payload_b], 1, &kmer_specs)?;
//         assert_eq!(merged.len(), 1);
//         let window_counts = merged[0].counts.get(&3).unwrap();
//         let value = window_counts.get("AAA").copied().unwrap_or_default();
//         assert!((value - 3.5).abs() < 1e-9);
//         Ok(())
//     }

//     #[test]
//     fn merge_tile_counts_merges_three_tiles() -> Result<()> {
//         let kmer_specs = build_kmer_specs(&[3])?;
//         let spec3 = &kmer_specs[&3];
//         let code_aaa = code_for_motif(spec3, "AAA");
//         let code_aac = code_for_motif(spec3, "AAC");

//         let payload_1 = vec![
//             TileWindowCounts {
//                 original_idx: 0,
//                 entries: vec![TileKmerCountEntry {
//                     k: 3,
//                     code: code_aaa,
//                     position: None,
//                     group: PositionGroup::Left,
//                     value: 1.0,
//                 }],
//             },
//             TileWindowCounts {
//                 original_idx: 1,
//                 entries: vec![TileKmerCountEntry {
//                     k: 3,
//                     code: code_aac,
//                     position: None,
//                     group: PositionGroup::Left,
//                     value: 2.0,
//                 }],
//             },
//         ];

//         let payload_2 = vec![
//             TileWindowCounts {
//                 original_idx: 0,
//                 entries: vec![TileKmerCountEntry {
//                     k: 3,
//                     code: code_aaa,
//                     position: None,
//                     group: PositionGroup::Left,
//                     value: 3.0,
//                 }],
//             },
//             TileWindowCounts {
//                 original_idx: 2,
//                 entries: vec![TileKmerCountEntry {
//                     k: 3,
//                     code: code_aaa,
//                     position: None,
//                     group: PositionGroup::Left,
//                     value: 5.0,
//                 }],
//             },
//         ];

//         let payload_3 = vec![
//             TileWindowCounts {
//                 original_idx: 0,
//                 entries: vec![TileKmerCountEntry {
//                     k: 3,
//                     code: code_aaa,
//                     position: None,
//                     group: PositionGroup::Left,
//                     value: 0.5,
//                 }],
//             },
//             TileWindowCounts {
//                 original_idx: 1,
//                 entries: vec![TileKmerCountEntry {
//                     k: 3,
//                     code: code_aac,
//                     position: None,
//                     group: PositionGroup::Left,
//                     value: 1.5,
//                 }],
//             },
//         ];

//         let merged = merge_tile_counts(vec![payload_1, payload_2, payload_3], 3, &kmer_specs)?;
//         assert_eq!(merged.len(), 3);

//         let win0 = merged[0].counts.get(&3).unwrap();
//         assert!((win0.get("AAA").copied().unwrap_or_default() - 4.5).abs() < 1e-9);

//         let win1 = merged[1].counts.get(&3).unwrap();
//         assert!((win1.get("AAC").copied().unwrap_or_default() - 3.5).abs() < 1e-9);

//         let win2 = merged[2].counts.get(&3).unwrap();
//         assert_eq!(win2.len(), 1);
//         assert!((win2.get("AAA").copied().unwrap_or_default() - 5.0).abs() < 1e-9);
//         Ok(())
//     }

//     #[test]
//     fn merge_tile_counts_rejects_out_of_range_indices() {
//         let kmer_specs = build_kmer_specs(&[3]).expect("build specs");
//         let spec3 = &kmer_specs[&3];
//         let code_aaa = code_for_motif(spec3, "AAA");

//         let payload = vec![TileWindowCounts {
//             original_idx: 5,
//             entries: vec![TileKmerCountEntry {
//                 k: 3,
//                 code: code_aaa,
//                 position: None,
//                 group: PositionGroup::Left,
//                 value: 1.0,
//             }],
//         }];

//         let result = merge_tile_counts(vec![payload], 2, &kmer_specs);
//         assert!(result.is_err());
//     }
// }

// mod tests_fragment_kmer_positions {
//     use std::num::NonZeroUsize;

//     use cfdnalab::{
//         commands::{
//             fragment_kmers::{
//                 fragment_kmers::count_kmers_at_positions,
//                 positions::{PositionGroup, PositionSelection, PositionSelectionCache},
//                 tiling::CountKey,
//             },
//             visualize_positions::{LinearRange, PositionsSpec, ReferenceFrame},
//         },
//         shared::{
//             fragment::segment_kmer_fragment::FragmentWithKmerSegments,
//             kmers::kmer_codec::{
//                 KmerCodes, KmerOrientation, KmerSpec, build_kmer_specs,
//                 build_left_aligned_codes_per_k,
//             },
//         },
//     };
//     use fxhash::FxHashMap;
//     use smallvec::smallvec;

//     #[test]
//     fn given_left_frame_when_counting_then_collects_forward_kmers() {
//         let seq = b"ACGTAC";
//         let context = TestContext::new(seq);
//         let cache = build_cache(
//             ReferenceFrame::Left,
//             PositionsSpec::Linear(LinearRange::All),
//             context.fragment.len(),
//         );
//         let selections = cache
//             .offsets(context.fragment.len())
//             .expect("left frame offsets");

//         let mut counts: FxHashMap<CountKey, f64> = FxHashMap::default();
//         count_kmers_at_positions(
//             &context.fragment,
//             selections,
//             true,
//             &context.positional_codes_by_k,
//             &context.kmer_specs,
//             &mut counts,
//             None,
//             0,
//             context.fragment.len(),
//             ReferenceFrame::Left,
//         );

//         let expected = expected_counts(
//             selections,
//             context.fragment.len() as usize,
//             &context.positional_codes_by_k[&context.k],
//             context.k as usize,
//         );
//         assert_eq!(counts, expected);
//         assert!(counts.keys().all(|key| key.group == PositionGroup::Left));
//         assert!(
//             counts
//                 .keys()
//                 .all(|key| matches!(key.orientation(), KmerOrientation::Forward))
//         );
//     }

//     #[test]
//     fn given_right_frame_when_counting_then_collects_reverse_kmers() {
//         let seq = b"ACGTAC";
//         let context = TestContext::new(seq);
//         let cache = build_cache(
//             ReferenceFrame::Right,
//             PositionsSpec::Linear(LinearRange::All),
//             context.fragment.len(),
//         );
//         let selections = cache
//             .offsets(context.fragment.len())
//             .expect("right frame offsets");

//         let mut counts: FxHashMap<CountKey, f64> = FxHashMap::default();
//         count_kmers_at_positions(
//             &context.fragment,
//             selections,
//             true,
//             &context.positional_codes_by_k,
//             &context.kmer_specs,
//             &mut counts,
//             None,
//             0,
//             context.fragment.len(),
//             ReferenceFrame::Right,
//         );

//         let expected = expected_counts(
//             selections,
//             context.fragment.len() as usize,
//             &context.positional_codes_by_k[&context.k],
//             context.k as usize,
//         );
//         assert_eq!(counts, expected);
//         assert!(
//             counts
//                 .keys()
//                 .all(|key| matches!(key.orientation(), KmerOrientation::Reverse))
//         );
//         assert!(counts.keys().all(|key| key.group == PositionGroup::Right));
//     }

//     #[test]
//     fn given_per_end_frame_when_counting_then_collects_both_orientations() {
//         let seq = b"ACGTAC";
//         let context = TestContext::new(seq);
//         let cache = build_cache(
//             ReferenceFrame::PerEnd,
//             PositionsSpec::Linear(LinearRange::All),
//             context.fragment.len(),
//         );
//         let selections = cache
//             .offsets(context.fragment.len())
//             .expect("per-end offsets");

//         let mut counts: FxHashMap<CountKey, f64> = FxHashMap::default();
//         count_kmers_at_positions(
//             &context.fragment,
//             selections,
//             true,
//             &context.positional_codes_by_k,
//             &context.kmer_specs,
//             &mut counts,
//             None,
//             0,
//             context.fragment.len(),
//             ReferenceFrame::PerEnd,
//         );

//         let expected = expected_counts(
//             selections,
//             context.fragment.len() as usize,
//             &context.positional_codes_by_k[&context.k],
//             context.k as usize,
//         );
//         assert_eq!(counts, expected);

//         assert!(
//             counts
//                 .keys()
//                 .any(|key| key.orientation() == KmerOrientation::Forward)
//         );
//         assert!(
//             counts
//                 .keys()
//                 .any(|key| key.orientation() == KmerOrientation::Reverse)
//         );
//         assert!(counts.keys().any(|key| key.group == PositionGroup::Left));
//         assert!(counts.keys().any(|key| key.group == PositionGroup::Right));
//         assert!(counts.keys().any(|key| key.group == PositionGroup::Left));
//         assert!(counts.keys().any(|key| key.group == PositionGroup::Right));
//     }

//     #[test]
//     fn given_nearest_frame_when_counting_then_splits_orientations_by_half() {
//         let seq = b"ACGTAC";
//         let context = TestContext::new(seq);
//         let cache = build_cache(
//             ReferenceFrame::Nearest,
//             PositionsSpec::Nearest(cfdnalab::commands::visualize_positions::NearestRange::All),
//             context.fragment.len(),
//         );
//         let selections = cache
//             .offsets(context.fragment.len())
//             .expect("nearest offsets");

//         let mut counts: FxHashMap<CountKey, f64> = FxHashMap::default();
//         count_kmers_at_positions(
//             &context.fragment,
//             selections,
//             true,
//             &context.positional_codes_by_k,
//             &context.kmer_specs,
//             &mut counts,
//             None,
//             0,
//             context.fragment.len(),
//             ReferenceFrame::Nearest,
//         );

//         let expected = expected_counts_nearest(
//             selections,
//             context.fragment.len() as usize,
//             &context.positional_codes_by_k[&context.k],
//             context.k as usize,
//         );
//         assert_eq!(counts, expected);

//         assert!(
//             counts
//                 .keys()
//                 .any(|key| key.orientation() == KmerOrientation::Forward)
//         );
//         assert!(
//             counts
//                 .keys()
//                 .any(|key| key.orientation() == KmerOrientation::Reverse)
//         );
//     }

//     struct TestContext {
//         fragment: FragmentWithKmerSegments,
//         kmer_specs: FxHashMap<u8, KmerSpec>,
//         positional_codes_by_k: FxHashMap<u8, KmerCodes>,
//         k: u8,
//     }

//     impl TestContext {
//         fn new(seq: &[u8]) -> Self {
//             let k_values = [3u8];
//             let kmer_specs = build_kmer_specs(&k_values).expect("kmer specs");
//             let positional_codes_by_k = build_left_aligned_codes_per_k(seq, &kmer_specs);
//             let fragment_len = seq.len() as u32;
//             let fragment = FragmentWithKmerSegments {
//                 tid: 0,
//                 start: 0,
//                 end: fragment_len,
//                 segments: smallvec![(0, fragment_len)],
//                 gc_tag: Default::default(),
//             };

//             Self {
//                 fragment,
//                 kmer_specs,
//                 positional_codes_by_k,
//                 k: k_values[0],
//             }
//         }
//     }

//     fn build_cache(
//         frame: ReferenceFrame,
//         positions: PositionsSpec,
//         length: u32,
//     ) -> PositionSelectionCache {
//         PositionSelectionCache::new(
//             frame,
//             &positions,
//             NonZeroUsize::new(1).expect("non-zero step"),
//             length,
//             length,
//         )
//         .expect("build selection cache")
//     }

//     fn expected_counts(
//         selections: &[PositionSelection],
//         fragment_len: usize,
//         codes: &KmerCodes,
//         k: usize,
//     ) -> FxHashMap<CountKey, f64> {
//         let mut expected = FxHashMap::default();
//         for selection in selections {
//             let offset = selection.offset() as usize;
//             match selection.orientation() {
//                 cfdnalab::commands::fragment_kmers::positions::PositionOrientation::Forward => {
//                     if offset + k > fragment_len {
//                         continue;
//                     }
//                     let code = codes.get(offset);
//                     let key = CountKey {
//                         k: k as u8,
//                         code,
//                         position: Some(selection.offset() as i32),
//                         group: selection.group(),
//                     };
//                     *expected.entry(key).or_insert(0.0) += 1.0;
//                 }
//                 cfdnalab::commands::fragment_kmers::positions::PositionOrientation::Reverse => {
//                     if offset + 1 < k || offset >= fragment_len {
//                         continue;
//                     }
//                     let start = offset + 1 - k;
//                     if start + k > fragment_len {
//                         continue;
//                     }
//                     let code = codes.get(start);
//                     let key = CountKey {
//                         k: k as u8,
//                         code,
//                         position: Some(selection.offset() as i32),
//                         group: selection.group(),
//                     };
//                     *expected.entry(key).or_insert(0.0) += 1.0;
//                 }
//             }
//         }
//         expected
//     }

//     fn expected_counts_nearest(
//         selections: &[PositionSelection],
//         fragment_len: usize,
//         codes: &KmerCodes,
//         k: usize,
//     ) -> FxHashMap<CountKey, f64> {
//         let mut expected = FxHashMap::default();
//         let len = fragment_len as u64;
//         let k_span = k as u64;
//         let half = len / 2; // floor

//         // Midpoint guards that mirror count_kmers_at_positions(… ReferenceFrame::Nearest …)
//         let (left_max_start, right_min_anchor) = if (len % 2) == 1 {
//             // Odd: exclude the physical midpoint at `half`
//             // Forward start <= half - k
//             // Reverse anchor (offset) >= half + k
//             (half.saturating_sub(k_span), half.saturating_add(k_span))
//         } else {
//             // Even: choose base nearest each side's start
//             // Forward start <= (L/2) - k
//             // Reverse anchor (offset) >= (L/2) + (k-1)
//             (
//                 half.saturating_sub(k_span),
//                 half.saturating_add(k_span.saturating_sub(1)),
//             )
//         };

//         for selection in selections {
//             let offset = selection.offset() as u64;
//             match selection.orientation() {
//                 cfdnalab::commands::fragment_kmers::positions::PositionOrientation::Forward => {
//                     // Reject if it crosses midpoint
//                     if offset > left_max_start {
//                         continue;
//                     }
//                     // Usual bounds
//                     if offset.saturating_add(k_span) > len {
//                         continue;
//                     }
//                     let start = offset as usize;
//                     let code = codes.get(start);
//                     let key = CountKey {
//                         k: k as u8,
//                         code,
//                         position: Some(selection.offset() as i32),
//                         group: selection.group(),
//                     };
//                     *expected.entry(key).or_insert(0.0) += 1.0;
//                 }
//                 cfdnalab::commands::fragment_kmers::positions::PositionOrientation::Reverse => {
//                     // Anchor must be far enough right so k-mer starts on/right of the right half
//                     if offset < right_min_anchor {
//                         continue;
//                     }
//                     if offset + 1 < k_span || offset >= len {
//                         continue;
//                     }
//                     let start = (offset + 1 - k_span) as usize;
//                     if (start as u64).saturating_add(k_span) > len {
//                         continue;
//                     }
//                     let code = codes.get(start);
//                     let key = CountKey {
//                         k: k as u8,
//                         code,
//                         position: Some(selection.offset() as i32),
//                         group: selection.group(),
//                     };
//                     *expected.entry(key).or_insert(0.0) += 1.0;
//                 }
//             }
//         }
//         expected
//     }
// }
// mod revcomp_tests {
//     use crate::fixtures::{simple_inward_bam, simple_reference_twobit, single_position_selection};
//     use crate::tests_fragment_kmer_command::{
//         assert_count_map_matches, base_chromosomes, build_revcomp_assets, load_counts_from_output,
//         load_positional_group_counts, manual_kmer_counts, manual_offset_counts,
//     };
//     use anyhow::Result;
//     use cfdnalab::commands::cli_common::{IOCArgs, Ref2BitRequiredArgs};
//     use cfdnalab::commands::fragment_kmers::config::FragmentKmersConfig;
//     use cfdnalab::commands::fragment_kmers::fragment_kmers::{run, run_inner};
//     use cfdnalab::commands::visualize_positions::ReferenceFrame;
//     use cfdnalab::shared::fragment::segment_kmer_fragment::FragmentWithKmerSegments;
//     use cfdnalab::shared::fragment_iterators::fragments_with_kmer_segments_from_bam;
//     use cfdnalab::shared::read::default_include_read_paired_end;
//     use tempfile::TempDir;

//     #[test]
//     fn nearest_counts_double_on_revcomp_reference() -> Result<()> {
//         let left_seq = "AGTACGCT";
//         let k = 2u8;
//         let expected = manual_kmer_counts(left_seq, k as usize);
//         let (
//             left_ref,
//             left_bam,
//             combined_ref,
//             combined_bam,
//             left_fragment_len,
//             combined_fragment_len,
//         ) = build_revcomp_assets(left_seq)?;

//         // Baseline: left reference only, left frame
//         let baseline_dir = TempDir::new()?;
//         let mut baseline_cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: left_bam.bam.clone(),
//                 output_dir: baseline_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: left_ref.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         baseline_cfg.set_output_prefix("left_only".to_string());
//         baseline_cfg.set_kmer_sizes(vec![k]);
//         baseline_cfg.set_min_mapq(0);
//         baseline_cfg.set_require_proper_pair(false);
//         baseline_cfg.set_ignore_gap(true);
//         baseline_cfg.set_canonical(false);
//         baseline_cfg.set_position_selection(single_position_selection(
//             ReferenceFrame::Left,
//             "..",
//             1,
//         ));
//         {
//             let lengths = baseline_cfg.fragment_lengths_mut();
//             lengths.min_fragment_length = left_fragment_len;
//             lengths.max_fragment_length = left_fragment_len;
//         }
//         run(&baseline_cfg)?;
//         let baseline_counts = load_counts_from_output(baseline_dir.path(), "left_only", k)?;
//         assert_count_map_matches(&baseline_counts, &expected, "baseline left-only");

//         // Combined reference: nearest frame should double every count
//         let combined_dir = TempDir::new()?;
//         let mut combined_cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: combined_bam.bam.clone(),
//                 output_dir: combined_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: combined_ref.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         combined_cfg.set_output_prefix("nearest".to_string());
//         combined_cfg.set_kmer_sizes(vec![k]);
//         combined_cfg.set_min_mapq(0);
//         combined_cfg.set_require_proper_pair(false);
//         combined_cfg.set_ignore_gap(true);
//         combined_cfg.set_canonical(false);
//         combined_cfg.set_position_selection(single_position_selection(
//             ReferenceFrame::Nearest,
//             "..",
//             1,
//         ));
//         {
//             let lengths = combined_cfg.fragment_lengths_mut();
//             lengths.min_fragment_length = combined_fragment_len;
//             lengths.max_fragment_length = combined_fragment_len;
//         }
//         run(&combined_cfg)?;
//         let combined_counts = load_counts_from_output(combined_dir.path(), "nearest", k)?;
//         for (motif, value) in &expected {
//             let observed = combined_counts.get(motif).copied().unwrap_or_default();
//             assert!(
//                 (observed - value * 2.0).abs() < 1e-6,
//                 "combined nearest global: motif {motif} expected {} observed {observed}",
//                 value * 2.0
//             );
//         }

//         Ok(())
//     }

//     #[test]
//     fn nearest_positional_counts_reflect_revcomp_symmetry() -> Result<()> {
//         let left_seq = "AGTACGCT";
//         let k = 2u8;
//         let expected_offsets = manual_offset_counts(left_seq, k as usize, &[0, 1, 2]);
//         let (
//             left_ref,
//             left_bam,
//             combined_ref,
//             combined_bam,
//             left_fragment_len,
//             combined_fragment_len,
//         ) = build_revcomp_assets(left_seq)?;

//         // Baseline positional counts using left frame
//         let baseline_dir = TempDir::new()?;
//         let mut baseline_cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: left_bam.bam.clone(),
//                 output_dir: baseline_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: left_ref.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         baseline_cfg.set_output_prefix("left_pos".to_string());
//         baseline_cfg.set_kmer_sizes(vec![k]);
//         baseline_cfg.set_min_mapq(0);
//         baseline_cfg.set_require_proper_pair(false);
//         baseline_cfg.set_ignore_gap(true);
//         baseline_cfg.set_canonical(false);
//         baseline_cfg.set_positional_counts(true);
//         baseline_cfg.set_position_selection(single_position_selection(
//             ReferenceFrame::Left,
//             "1..3",
//             1,
//         ));
//         {
//             let lengths = baseline_cfg.fragment_lengths_mut();
//             lengths.min_fragment_length = left_fragment_len;
//             lengths.max_fragment_length = left_fragment_len;
//         }
//         run(&baseline_cfg)?;
//         let baseline_left =
//             load_positional_group_counts(baseline_dir.path(), "left_pos", k, "left")?;
//         assert_count_map_matches(
//             &baseline_left,
//             &expected_offsets,
//             "baseline positional left",
//         );

//         // Combined positional counts with nearest frame
//         let combined_dir = TempDir::new()?;
//         let mut combined_cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: combined_bam.bam.clone(),
//                 output_dir: combined_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: combined_ref.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         combined_cfg.set_output_prefix("nearest_pos".to_string());
//         combined_cfg.set_kmer_sizes(vec![k]);
//         combined_cfg.set_min_mapq(0);
//         combined_cfg.set_require_proper_pair(false);
//         combined_cfg.set_ignore_gap(true);
//         combined_cfg.set_canonical(false);
//         combined_cfg.set_positional_counts(true);
//         combined_cfg.set_position_selection(single_position_selection(
//             ReferenceFrame::Nearest,
//             "1..3",
//             1,
//         ));
//         {
//             let lengths = combined_cfg.fragment_lengths_mut();
//             lengths.min_fragment_length = combined_fragment_len;
//             lengths.max_fragment_length = combined_fragment_len;
//         }
//         run(&combined_cfg)?;
//         let nearest_left =
//             load_positional_group_counts(combined_dir.path(), "nearest_pos", k, "left")?;
//         let nearest_right =
//             load_positional_group_counts(combined_dir.path(), "nearest_pos", k, "right")?;

//         assert_count_map_matches(&nearest_left, &expected_offsets, "nearest positional left");
//         assert_count_map_matches(
//             &nearest_right,
//             &expected_offsets,
//             "nearest positional right",
//         );
//         for (motif, expected) in &expected_offsets {
//             let total = nearest_left.get(motif).copied().unwrap_or_default()
//                 + nearest_right.get(motif).copied().unwrap_or_default();
//             assert!(
//                 (total - expected * 2.0).abs() < 1e-6,
//                 "motif {motif} expected {} combined observed {total}",
//                 expected * 2.0
//             );
//         }

//         Ok(())
//     }

//     #[test]
//     fn per_end_half_positions_are_balanced() -> Result<()> {
//         let left_seq = "AGTACGCT";
//         let k = 2u8;
//         let expected_half = manual_offset_counts(left_seq, k as usize, &[0, 1, 2, 3]);
//         let (
//             _left_ref,
//             _left_bam,
//             combined_ref,
//             combined_bam,
//             _left_fragment_len,
//             combined_fragment_len,
//         ) = build_revcomp_assets(left_seq)?;

//         let out_dir = TempDir::new()?;
//         let mut cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: combined_bam.bam.clone(),
//                 output_dir: out_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: combined_ref.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         cfg.set_output_prefix("per_end".to_string());
//         cfg.set_kmer_sizes(vec![k]);
//         cfg.set_min_mapq(0);
//         cfg.set_require_proper_pair(false);
//         cfg.set_ignore_gap(true);
//         cfg.set_canonical(false);
//         cfg.set_positional_counts(true);
//         cfg.set_position_selection(single_position_selection(ReferenceFrame::PerEnd, "..4", 1));
//         {
//             let lengths = cfg.fragment_lengths_mut();
//             lengths.min_fragment_length = combined_fragment_len;
//             lengths.max_fragment_length = combined_fragment_len;
//         }
//         run(&cfg)?;

//         let per_end_left = load_positional_group_counts(out_dir.path(), "per_end", k, "left")?;
//         let per_end_right = load_positional_group_counts(out_dir.path(), "per_end", k, "right")?;
//         assert_count_map_matches(&per_end_left, &expected_half, "per-end left");
//         assert_count_map_matches(&per_end_right, &expected_half, "per-end right");

//         Ok(())
//     }

//     #[test]
//     fn run_inner_positions_match_counts() -> Result<()> {
//         use ndarray::{Array2, Array3};
//         use rust_htslib::bam::{Read, Reader};
//         use std::collections::{HashMap, HashSet};

//         let bam = simple_inward_bam()?;
//         let reference = simple_reference_twobit()?;

//         let positional_dir = TempDir::new()?;
//         let counts_dir = TempDir::new()?;

//         let k_sizes = vec![2u8];
//         let k = k_sizes[0];
//         let positional_prefix = "pos_reconstruct";
//         let counts_prefix = "agg_reconstruct";

//         let configure = |cfg: &mut FragmentKmersConfig, output_prefix: &str, positional: bool| {
//             cfg.set_output_prefix(output_prefix.to_string());
//             cfg.set_kmer_sizes(k_sizes.clone());
//             cfg.set_min_mapq(0);
//             cfg.set_require_proper_pair(false);
//             cfg.set_ignore_gap(false);
//             cfg.set_canonical(false);
//             cfg.set_positional_counts(positional);
//             cfg.set_position_selection(single_position_selection(ReferenceFrame::Nearest, "..", 1));
//             let lengths = cfg.fragment_lengths_mut();
//             lengths.min_fragment_length = 20;
//             lengths.max_fragment_length = 120;
//         };

//         let mut positional_cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: positional_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         configure(&mut positional_cfg, positional_prefix, true);

//         let mut aggregate_cfg = FragmentKmersConfig::new(
//             IOCArgs {
//                 bam: bam.bam.clone(),
//                 output_dir: counts_dir.path().to_path_buf(),
//                 n_threads: 1,
//             },
//             Ref2BitRequiredArgs {
//                 ref_2bit: reference.path.clone(),
//             },
//             base_chromosomes(&["chr1"]),
//         );
//         configure(&mut aggregate_cfg, counts_prefix, false);

//         run_inner(&aggregate_cfg)?;
//         run_inner(&positional_cfg)?;

//         let aggregate_counts_path = counts_dir
//             .path()
//             .join(format!("{counts_prefix}.k{k}_counts.npy"));
//         let aggregate_motifs_path = counts_dir
//             .path()
//             .join(format!("{counts_prefix}.k{k}_motifs.txt"));

//         let aggregate_counts: Array2<f64> = ndarray_npy::read_npy(&aggregate_counts_path)?;
//         let aggregate_motifs: Vec<String> = std::fs::read_to_string(&aggregate_motifs_path)?
//             .lines()
//             .map(|line| line.to_string())
//             .collect();
//         assert_eq!(
//             aggregate_counts.shape()[1],
//             aggregate_motifs.len(),
//             "motif list should align with aggregate matrix columns"
//         );

//         let n_windows = aggregate_counts.shape()[0];

//         let total_aggregate: f64 = aggregate_counts.iter().copied().sum();
//         let mut aggregate_totals_by_motif: HashMap<String, f64> = HashMap::new();
//         for (motif_idx, motif) in aggregate_motifs.iter().enumerate() {
//             let total = (0..n_windows)
//                 .map(|w| aggregate_counts[[w, motif_idx]])
//                 .sum::<f64>();
//             aggregate_totals_by_motif.insert(motif.clone(), total);
//         }

//         let mut groups: Vec<String> = std::fs::read_dir(positional_dir.path())?
//             .filter_map(|entry| {
//                 let entry = entry.ok()?;
//                 let name = entry.file_name();
//                 let name = name.to_string_lossy();
//                 let prefix_with_dot = format!("{positional_prefix}.");
//                 if let Some(rest) = name.strip_prefix(&prefix_with_dot) {
//                     if let Some(group) = rest.strip_suffix("_positions.txt") {
//                         return Some(group.to_string());
//                     }
//                 }
//                 None
//             })
//             .collect();
//         groups.sort();
//         groups.dedup();

//         let mut positional_totals_by_motif: HashMap<String, f64> = HashMap::new();
//         let mut positional_total = 0.0f64;
//         let mut processed_groups = 0usize;
//         for group in groups {
//             let counts_path = positional_dir
//                 .path()
//                 .join(format!("{positional_prefix}.k{k}_{group}_counts.npy"));
//             if !counts_path.exists() {
//                 continue;
//             }
//             processed_groups += 1;

//             let offsets_path = positional_dir
//                 .path()
//                 .join(format!("{positional_prefix}.{group}_positions.txt"));
//             let offsets: Vec<i32> = std::fs::read_to_string(&offsets_path)?
//                 .lines()
//                 .map(|line| line.trim().parse::<i32>().expect("offset"))
//                 .collect();
//             assert!(
//                 !offsets.is_empty(),
//                 "expected positional metadata for group {group}"
//             );

//             let counts: Array3<f64> = ndarray_npy::read_npy(&counts_path)?;
//             assert_eq!(
//                 counts.shape()[1],
//                 offsets.len(),
//                 "axis 1 should mirror stored positions for group {group}"
//             );
//             assert_eq!(
//                 counts.shape()[0],
//                 n_windows,
//                 "window axis must match aggregate counts"
//             );

//             let motifs_path = positional_dir
//                 .path()
//                 .join(format!("{positional_prefix}.k{k}_{group}_motifs.txt"));
//             let group_motifs: Vec<String> = std::fs::read_to_string(&motifs_path)?
//                 .lines()
//                 .map(|line| line.to_string())
//                 .collect();
//             assert_eq!(
//                 counts.shape()[2],
//                 group_motifs.len(),
//                 "motif axis should mirror motif list for group {group}"
//             );

//             let mut observed_offsets = HashSet::new();
//             for (pos_idx, offset) in offsets.iter().enumerate() {
//                 let mut has_signal = false;
//                 for window_idx in 0..n_windows {
//                     for (motif_idx, motif) in group_motifs.iter().enumerate() {
//                         let value = counts[[window_idx, pos_idx, motif_idx]];
//                         if value > 0.0 {
//                             has_signal = true;
//                             *positional_totals_by_motif
//                                 .entry(motif.clone())
//                                 .or_insert(0.0) += value;
//                             positional_total += value;
//                         }
//                     }
//                 }
//                 if has_signal {
//                     observed_offsets.insert(*offset);
//                 }
//             }

//             let expected_offsets: HashSet<i32> = offsets.into_iter().collect();
//             assert_eq!(
//                 observed_offsets, expected_offsets,
//                 "positions reconstructed from counts differed for group {group}"
//             );
//         }
//         assert!(
//             processed_groups > 0,
//             "expected at least one positional group with counts"
//         );

//         use anyhow::Error;
//         let require_proper_pair = positional_cfg.shared_args.require_proper_pair;
//         let min_mapq = positional_cfg.shared_args.min_mapq;
//         let indel_mode = positional_cfg.shared_args.indel_mode;
//         let include_gap = !positional_cfg.shared_args.ignore_gap;
//         let length_filter = positional_cfg.shared_args.fragment_lengths.clone();

//         let include_read = move |rec: &rust_htslib::bam::Record| {
//             default_include_read_paired_end(rec, require_proper_pair, min_mapq)
//         };
//         let fragment_filter =
//             move |fragment: &FragmentWithKmerSegments| length_filter.contains(fragment.len());

//         let mut reader = Reader::from_path(&bam.bam)?;
//         let fragments: Vec<FragmentWithKmerSegments> = fragments_with_kmer_segments_from_bam(
//             reader.records().map(|r| r.map_err(Error::from)),
//             include_read,
//             indel_mode,
//             include_gap,
//             0,
//             None,
//             fragment_filter,
//             false,
//         )
//         .collect::<Result<Vec<_>>>()?;
//         assert!(
//             !fragments.is_empty(),
//             "expected at least one fragment after filtering"
//         );

//         let expected_total: f64 = fragments
//             .iter()
//             .map(|fragment| expected_nearest_positions_for_length(fragment.len(), k as u32) as f64)
//             .sum();
//         assert!(
//             expected_total > 0.0,
//             "expected at least one counted position"
//         );

//         let tolerance = 1e-6f64;
//         assert!(
//             (total_aggregate - positional_total).abs() < tolerance,
//             "aggregate total {} mismatched positional total {}",
//             total_aggregate,
//             positional_total
//         );
//         assert!(
//             (total_aggregate - expected_total).abs() < tolerance,
//             "aggregate total {} mismatched expected positional count {}",
//             total_aggregate,
//             expected_total
//         );
//         assert!(
//             (positional_total - expected_total).abs() < tolerance,
//             "positional total {} mismatched expected positional count {}",
//             positional_total,
//             expected_total
//         );

//         for (motif, aggregate) in &aggregate_totals_by_motif {
//             let positional = positional_totals_by_motif
//                 .get(motif)
//                 .copied()
//                 .unwrap_or_default();
//             assert!(
//                 (aggregate - positional).abs() < tolerance,
//                 "motif {motif} aggregate {aggregate} positional {positional}"
//             );
//         }
//         for (motif, positional) in &positional_totals_by_motif {
//             let aggregate = aggregate_totals_by_motif
//                 .get(motif)
//                 .copied()
//                 .unwrap_or_default();
//             assert!(
//                 (aggregate - positional).abs() < tolerance,
//                 "motif {motif} aggregate {aggregate} positional {positional}"
//             );
//         }

//         Ok(())
//     }

//     fn expected_nearest_positions_for_length(length: u32, k: u32) -> u64 {
//         if k == 0 || length < k {
//             return 0;
//         }
//         let even_length = if length % 2 == 1 { length - 1 } else { length };
//         even_length.saturating_sub(2 * (k - 1)) as u64
//     }

//     #[test]
//     fn global_windowing_handles_three_chromosomes() -> Result<()> {
//         // Planned regression:
//         // - Build one simple fragment per chromosome.
//         // - Run the global window path.
//         // - Assert that the output keeps chromosome order and counts all three chromosomes.
//         // This should mirror the three-chromosome smoke tests added to the active suites.
//         todo!("fragment_kmers suite is intentionally commented out until the command works again");
//     }

//     #[test]
//     fn by_size_windowing_handles_three_chromosomes() -> Result<()> {
//         // Planned regression:
//         // - Use one fixed-size window per chromosome.
//         // - Assert that all three chromosomes produce output rows and that window indexing
//         //   does not silently assume chromosome-local dense indices.
//         todo!("fragment_kmers suite is intentionally commented out until the command works again");
//     }

//     #[test]
//     fn by_bed_windowing_handles_three_chromosomes() -> Result<()> {
//         // Planned regression:
//         // - Load one BED window per chromosome.
//         // - Assert that later chromosomes survive the reducer path instead of being dropped
//         //   by a hidden original-index assumption.
//         todo!("fragment_kmers suite is intentionally commented out until the command works again");
//     }
// }
