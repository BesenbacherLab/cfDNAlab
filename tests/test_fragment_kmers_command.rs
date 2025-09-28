mod fixtures;

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use cfdnalab::cli_common::{
    ChromosomeArgs, FragmentLengthArgs, IOCArgs, Ref2BitRequiredArgs, ScaleGenomeArgs, WindowsArgs,
};
use cfdnalab::fragment_kmers::{FragmentKmersConfig, run};
use cfdnalab::utils::base::make_canonical;
use cfdnalab::utils::blacklist::BlacklistStrategy;
use cfdnalab::utils::indel_mode::IndelMode;
use fixtures::{
    fragment_kmers_edge_bam, fragment_kmers_edge_reference, simple_inward_bam,
    simple_reference_twobit, write_bed, write_scaling_factors,
};
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
    cfg_base.set_end_offset(1);
    {
        let fl = cfg_base.fragment_lengths_mut();
        fl.min_fragment_length = fragment_lengths.min_fragment_length;
        fl.max_fragment_length = fragment_lengths.max_fragment_length;
    }

    run(cfg_base.clone())?;

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

    run(cfg_scaled.clone())?;

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
