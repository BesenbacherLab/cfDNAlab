mod fixtures;

use anyhow::Result;
use cfdnalab::cli_common::{ChromosomeArgs, IOCArgs};
use cfdnalab::normalize_genome::{NormalizeGenomeConfig, run};
use fixtures::simple_inward_bam;
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
}

#[test]
fn coverage_scaling_written_with_expected_ranges() -> Result<()> {
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = NormalizeGenomeConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_bin_size(40);
    cfg.set_stride(20);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    run(cfg)?;

    let tsv_path = out_dir.path().join("coverage_scaling_factors.tsv");
    assert!(tsv_path.exists());
    let content = std::fs::read_to_string(&tsv_path)?;
    let mut lines = content.lines();
    let header = lines.next().unwrap_or("");
    assert!(header.contains("chromosome"));
    let mut saw_zero = false;
    let mut saw_non_zero = false;
    for line in lines {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 6);
        let start: u64 = parts[1].parse().unwrap();
        let scaling: f64 = parts[5].parse().unwrap();
        if scaling == 0.0 {
            saw_zero = true;
        }
        if start >= 20 && start < 80 && scaling > 0.0 {
            saw_non_zero = true;
        }
    }
    assert!(saw_zero, "expected uncovered stride bin with scaling 0");
    assert!(
        saw_non_zero,
        "expected covered stride bin with positive scaling"
    );

    Ok(())
}

#[test]
fn check_bin_sizes_rejects_invalid_stride() {
    let mut cfg = NormalizeGenomeConfig::new(
        IOCArgs {
            bam: std::path::PathBuf::new(),
            output_dir: std::path::PathBuf::new(),
            n_threads: 1,
        },
        ChromosomeArgs::default(),
    );
    cfg.set_bin_size(30);
    cfg.set_stride(40);
    let err = cfg.check_bin_sizes().unwrap_err();
    assert!(format!("{err}").contains("stride"));
}
