mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{AssignToWindowArgs, ChromosomeArgs, IOCArgs, WindowsArgs};
use cfdnalab::commands::lengths::config::LengthsConfig;
use cfdnalab::commands::lengths::lengths::run;
use cfdnalab::shared::indel_mode::IndelMode;
use fixtures::simple_inward_bam;
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
fn counts_reference_lengths_global_window() -> Result<()> {
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = LengthsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_indel_mode(IndelMode::Ignore);
    cfg.set_windows(WindowsArgs::default());
    cfg.set_window_assignment(AssignToWindowArgs::default());
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    run(&cfg)?;

    let npy_path = out_dir.path().join("all_length_counts.npy");
    assert!(npy_path.exists());
    let arr: Array2<f64> = read_npy(&npy_path)?;
    assert_eq!(arr.shape(), &[1, 191]);
    let len60_idx = 60 - 10; // min_fragment_length
    assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
    assert_eq!(arr[(0, len60_idx - 1)], 0.0);

    Ok(())
}
