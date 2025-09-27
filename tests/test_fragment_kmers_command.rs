mod fixtures;

use std::collections::HashMap;
use std::str;

use anyhow::{Context, Result};
use cfdnalab::cli_common::{ChromosomeArgs, IOCArgs, Ref2BitRequiredArgs, WindowsArgs};
use cfdnalab::fragment_kmers::{FragmentKmersConfig, run};
use cfdnalab::utils::base::make_canonical;
use cfdnalab::utils::indel_mode::IndelMode;
use fixtures::{simple_inward_bam, simple_reference_twobit};
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

    run(cfg)?;

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

    run(cfg)?;

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
