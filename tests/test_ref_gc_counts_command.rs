#![cfg(feature = "cmd_ref_gc_counts")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, FragmentLengthArgs, WindowsArgs};
use cfdnalab::commands::ref_gc_counts::{config::RefGCCountsConfig, ref_gc_counts::run};
use fixtures::{simple_reference_twobit, write_bed};
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|chr| chr.to_string()).collect()),
        chromosomes_file: None,
    }
}

#[test]
fn bed_windowed_runs_write_ref_gc_bins_tsv_with_exact_blacklisted_fractions() -> Result<()> {
    // Arrange:
    // - The two BED windows are [10,20) and [20,30).
    // - The blacklist interval [15,20) overlaps only the first window for 5 of its 10 bases.
    // - `ref-gc-counts` should therefore persist:
    //     chr1  10  20  0.5
    //     chr1  20  30  0
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 20, "left"), ("chr1", 20, 30, "right")],
    )?;
    write_bed(&blacklist_bed, &[("chr1", 15, 20, "masked")])?;

    let cfg = RefGCCountsConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        n_threads: 1,
        n_positions: 8,
        seed: Some(7),
        windows: WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        },
        chromosomes: base_chromosomes(&["chr1"]),
        blacklist: Some(vec![blacklist_bed]),
        fragment_lengths: FragmentLengthArgs {
            min_fragment_length: 20,
            max_fragment_length: 20,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
    };

    // Act
    run(&cfg)?;
    let bins_tsv = std::fs::read_to_string(out_dir.path().join("ref_gc_bins.tsv"))?;

    // Assert
    assert_eq!(
        bins_tsv,
        concat!(
            "chrom\tstart\tend\tblacklisted_fraction\n",
            "chr1\t10\t20\t0.5\n",
            "chr1\t20\t30\t0\n"
        )
    );
    Ok(())
}
