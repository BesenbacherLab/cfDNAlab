//! Helpers for GC correction package test inputs.
//!
//! These helpers create GC correction packages for tests that need to exercise
//! downstream command behavior with a real package file. The hand-authored
//! package writers are small and deterministic. The command-produced helpers
//! run `ref-gc-bias` and `gc-bias` with quiet options and require both
//! `cmd_gc_bias` and `cmd_ref_gc_bias`.
//!
//! This module deliberately exposes GC correction package helpers, not private
//! reference-GC loaders. If a test needs to inspect a reference-GC package, add
//! a test-facing summary type instead of publishing command internals.

#[cfg(feature = "cmd_gc_bias")]
use anyhow::{Result, ensure};
#[cfg(feature = "cmd_gc_bias")]
use std::path::Path;
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use std::path::PathBuf;
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use tempfile::TempDir;

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use crate::{
    RunOptions,
    commands::{
        cli_common::{
            ChromosomeArgs, FragmentLengthArgs, GCWindowsArgs, IOCArgs, LoggingArgs,
            Ref2BitRequiredArgs,
        },
        gc_bias::{
            config::{GCConfig, OutlierMethodArg},
            gc_bias::run_gc_bias,
        },
        ref_gc_bias::{
            config::{RefGCBiasConfig, RefGCWindowsArgs},
            ref_gc_bias::run_ref_gc_bias,
        },
    },
    reference::twobit_contig_lengths,
};

/// Write a hand-authored one-length-bin, one-GC-bin correction package.
///
/// Requires the `testing` and `cmd_gc_bias` cargo features.
///
/// The package covers exactly one half-open fragment length bin:
/// `[fragment_length, fragment_length + 1)`. It has one GC bin covering
/// `[0, 101)`, so every GC percentage maps to the same multiplicative weight.
/// The correction matrix has shape `(1, 1)` and contains `weight`.
///
/// Use this when a test needs a valid GC correction package where the expected
/// correction factor is the same for every fragment. The package is
/// hand-authored, not produced by running `gc-bias`, so it is best for tests of
/// downstream package consumption rather than tests of GC-bias estimation.
///
/// The package is written to `path` as a Zarr store. `end_offset` is 0,
/// `length_bin_frequencies` is `[1.0]`, and `reference_contig_footprint` is
/// empty.
#[cfg(feature = "cmd_gc_bias")]
pub fn write_constant_gc_correction_package(
    path: &Path,
    fragment_length: u32,
    weight: f64,
) -> Result<()> {
    let package = crate::gc_bias::GCCorrectionPackage {
        version: crate::constants::GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![fragment_length, fragment_length + 1],
        gc_edges: vec![0, 101],
        length_bin_frequencies: ndarray::array![1.0_f64],
        reference_contig_footprint: Vec::new(),
        correction_matrix: ndarray::array![[weight]],
    };
    package.write_zarr(path)?;
    Ok(())
}

/// Write a hand-authored unit-weight correction package for one fragment length.
///
/// Requires the `testing` and `cmd_gc_bias` cargo features.
///
/// This is the explicit all-weights-are-1 wrapper. The package covers exactly
/// one half-open fragment length bin: `[fragment_length, fragment_length + 1)`.
/// It has one GC bin covering `[0, 101)`, a `(1, 1)` correction matrix
/// containing `1.0`, `length_bin_frequencies = [1.0]`, `end_offset = 0`, and
/// an empty `reference_contig_footprint`.
///
/// Use this when a test needs a valid package that should leave every fragment
/// weight unchanged. It is hand-authored, not produced by running `gc-bias`.
#[cfg(feature = "cmd_gc_bias")]
pub fn write_unit_gc_correction_package(path: &Path, fragment_length: u32) -> Result<()> {
    write_constant_gc_correction_package(path, fragment_length, 1.0)
}

/// Write a hand-authored unit-weight correction package for a fragment length range.
///
/// Requires the `testing` and `cmd_gc_bias` cargo features.
///
/// This helper creates one length bin per integer fragment length in the
/// inclusive range `min_fragment_length..=max_fragment_length`. GC still has a
/// single bin `[0, 101)`. Every correction weight is `1.0`, so the package is
/// neutral by construction.
///
/// Length-bin frequencies are uniform across the generated length bins.
/// `end_offset` is 0 and `reference_contig_footprint` is empty. The package is
/// hand-authored, not produced by running `gc-bias`.
#[cfg(feature = "cmd_gc_bias")]
pub fn write_unit_gc_correction_package_for_range(
    path: &Path,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<()> {
    ensure!(
        min_fragment_length <= max_fragment_length,
        "minimum fragment length ({min_fragment_length}) must be <= maximum fragment length ({max_fragment_length})"
    );
    let final_edge = max_fragment_length
        .checked_add(1)
        .ok_or_else(|| anyhow::anyhow!("maximum fragment length must be less than u32::MAX"))?;
    let length_edges: Vec<u32> = (min_fragment_length..=final_edge).collect();
    let length_bin_count = length_edges.len() - 1;
    let length_bin_frequency = 1.0 / length_bin_count as f64;
    let package = crate::gc_bias::GCCorrectionPackage {
        version: crate::constants::GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges,
        gc_edges: vec![0, 101],
        length_bin_frequencies: ndarray::Array1::from_elem(length_bin_count, length_bin_frequency),
        reference_contig_footprint: Vec::new(),
        correction_matrix: ndarray::Array2::from_elem((length_bin_count, 1), 1.0),
    };
    package.write_zarr(path)?;
    Ok(())
}

/// Write a hand-authored one-length-bin, two-GC-bin correction package.
///
/// Requires the `testing` and `cmd_gc_bias` cargo features.
///
/// The package covers exactly one half-open fragment length bin:
/// `[fragment_length, fragment_length + 1)`. GC values below 51 use
/// `low_gc_weight`, and GC values from 51 through 100 use `high_gc_weight`.
/// The correction matrix has shape `(1, 2)` and contains those two weights in
/// GC-bin order.
///
/// `reference_contig_footprint` is written into the package unchanged. Pass an
/// empty vector when the test only needs correction weights and not footprint
/// propagation.
///
/// The package is hand-authored and written to `path` as a Zarr store.
/// `end_offset` is 0 and `length_bin_frequencies` is `[1.0]`.
#[cfg(feature = "cmd_gc_bias")]
pub fn write_two_bin_gc_correction_package(
    path: &Path,
    fragment_length: u32,
    low_gc_weight: f64,
    high_gc_weight: f64,
    reference_contig_footprint: Vec<crate::reference::ContigFootprintEntry>,
) -> Result<()> {
    let package = crate::gc_bias::GCCorrectionPackage {
        version: crate::constants::GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![fragment_length, fragment_length + 1],
        gc_edges: vec![0, 51, 101],
        length_bin_frequencies: ndarray::array![1.0_f64],
        reference_contig_footprint,
        correction_matrix: ndarray::array![[low_gc_weight, high_gc_weight]],
    };
    package.write_zarr(path)?;
    Ok(())
}

/// Build a command-produced GC correction package for one fragment length.
///
/// Requires the `testing`, `cmd_gc_bias`, and `cmd_ref_gc_bias` cargo features.
///
/// This runs the real `ref-gc-bias` then `gc-bias` producer chain. The helper
/// is intended for tests that need a valid package artifact produced through
/// the same path as command output.
///
/// The reference side uses global windows on `chr1`, a deterministic seed, and
/// at most 100 sampled start positions. The generated package is written under
/// `out_dir`, and the returned path points to `gc_bias_correction.zarr`.
///
/// `fragment_length` is used as both minimum and maximum fragment length. The
/// BAM and reference are expected to contain a `chr1` contig because the helper
/// configures both producer commands for that chromosome.
///
/// This helper does not promise unit weights. The resulting correction matrix
/// depends on the supplied BAM, reference, and producer settings.
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
pub fn build_command_produced_gc_correction_package_for_length(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    fragment_length: u32,
) -> Result<PathBuf> {
    build_command_produced_gc_correction_package_for_range(
        bam_path,
        reference_path,
        out_dir,
        fragment_length,
        fragment_length,
    )
}

/// Build a command-produced GC correction package for a fragment length range.
///
/// Requires the `testing`, `cmd_gc_bias`, and `cmd_ref_gc_bias` cargo features.
///
/// This is the range version of
/// `build_command_produced_gc_correction_package_for_length`. It is useful when
/// a test BAM contains sparse observed lengths but the correction package must
/// honestly cover a broader configured range.
///
/// The helper computes a valid sampled-position count from the `chr1` reference
/// length and `max_fragment_length`. It fails if the reference has no valid
/// start positions for the requested range.
///
/// The helper runs `ref-gc-bias` with global windows, `seed = 7`, at most 100
/// sampled positions, no interpolation, no smoothing, `end_offset = 0`, and
/// one thread. It then runs `gc-bias` with global windows, one thread, minimum
/// MAPQ 0, no extreme-GC-bin or short-length-bin outlier handling, minimum
/// length-bin mass 0, and minimum length-bin width 1.
///
/// The package is written below `out_dir` in a deterministic subdirectory named
/// from the requested length range. The returned path points to the generated
/// `gc_bias_correction.zarr` store.
///
/// This helper does not promise unit weights. The resulting correction matrix
/// depends on the supplied BAM, reference, and producer settings.
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
pub fn build_command_produced_gc_correction_package_for_range(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<PathBuf> {
    let chromosomes = vec!["chr1".to_string()];
    let chrom_lengths = twobit_contig_lengths(reference_path, &chromosomes)?;
    let total_possible_starts: usize = chrom_lengths
        .values()
        .map(|&chrom_len| {
            chrom_len
                .checked_sub(max_fragment_length as usize)
                .map(|remaining| remaining + 1)
                .unwrap_or(0)
        })
        .sum();
    let n_positions = total_possible_starts.min(100);
    ensure!(
        n_positions > 0,
        "GC package helper has no valid reference start positions for fragment length range {}-{}",
        min_fragment_length,
        max_fragment_length
    );

    let ref_gc_dir = TempDir::new()?;
    let ref_cfg = RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs {
            ref_2bit: reference_path.to_path_buf(),
        },
        output_dir: ref_gc_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: base_chromosomes(&["chr1"]),
        blacklist: None,
        fragment_lengths: FragmentLengthArgs {
            min_fragment_length,
            max_fragment_length,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };
    run_ref_gc_bias(&ref_cfg, RunOptions::new_quiet())?;

    let gc_out_dir = out_dir.join(format!(
        "command_produced_gc_correction_len_{}-{}",
        min_fragment_length, max_fragment_length
    ));
    std::fs::create_dir_all(&gc_out_dir)?;
    let mut gc_cfg = GCConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: gc_out_dir.clone(),
            n_threads: 1,
        },
        reference_path.to_path_buf(),
        ref_gc_dir.path().join("ref_gc_package.zarr"),
        base_chromosomes(&["chr1"]),
    );
    configure_gc_bias_common(&mut gc_cfg);
    gc_cfg.set_min_length_bin_mass(0.0);
    gc_cfg.set_min_length_bin_width(1);
    run_gc_bias(&gc_cfg, RunOptions::new_quiet())?;

    Ok(gc_out_dir.join("gc_bias_correction.zarr"))
}

/// Build a command-produced GC correction package from caller-supplied reference windows.
///
/// Requires the `testing`, `cmd_gc_bias`, and `cmd_ref_gc_bias` cargo features.
///
/// Use this when the reference-side counts and resulting weights must be
/// derivable from explicit windows at the test site.
///
/// `reference_windows_bed` is written to a temporary BED file and passed to
/// `ref-gc-bias`. `n_positions` is caller-controlled so tests can keep the
/// sampled-start arithmetic visible in the expected-value derivation.
///
/// The helper configures both producer commands for `chr1`, uses `fragment_length`
/// as both minimum and maximum fragment length, runs with one thread, disables
/// interpolation and smoothing in the reference package, and uses no outlier
/// handling in the observed package. The reference-side seed is 23.
///
/// The package is written below `out_dir` in a deterministic subdirectory named
/// from `fragment_length`. The returned path points to the generated
/// `gc_bias_correction.zarr` store.
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
pub fn build_command_produced_gc_correction_package_from_reference_windows(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    fragment_length: u32,
    reference_windows_bed: &str,
    n_positions: usize,
) -> Result<PathBuf> {
    build_command_produced_gc_correction_package_from_reference_windows_for_range(
        bam_path,
        reference_path,
        out_dir,
        fragment_length,
        fragment_length,
        reference_windows_bed,
        n_positions,
    )
}

/// Build a command-produced GC correction package from reference windows and a length range.
///
/// Requires the `testing`, `cmd_gc_bias`, and `cmd_ref_gc_bias` cargo features.
///
/// This is the range version of
/// `build_command_produced_gc_correction_package_from_reference_windows`. It
/// writes `reference_windows_bed` to a temporary BED file, runs `ref-gc-bias`
/// for `min_fragment_length..=max_fragment_length`, then runs `gc-bias` against
/// the generated reference package.
///
/// The helper configures both producer commands for `chr1`, runs with one
/// thread, disables interpolation and smoothing in the reference package, and
/// uses no outlier handling in the observed package. The reference-side seed is
/// 23.
///
/// The package is written below `out_dir` in a deterministic subdirectory named
/// from the requested length range. The returned path points to the generated
/// `gc_bias_correction.zarr` store.
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
pub fn build_command_produced_gc_correction_package_from_reference_windows_for_range(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    min_fragment_length: u32,
    max_fragment_length: u32,
    reference_windows_bed: &str,
    n_positions: usize,
) -> Result<PathBuf> {
    ensure!(
        min_fragment_length <= max_fragment_length,
        "minimum fragment length ({min_fragment_length}) must be <= maximum fragment length ({max_fragment_length})"
    );
    let ref_gc_dir = TempDir::new()?;
    let bed_path = ref_gc_dir.path().join("reference_windows.bed");
    std::fs::write(&bed_path, reference_windows_bed)?;
    let ref_cfg = RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs {
            ref_2bit: reference_path.to_path_buf(),
        },
        output_dir: ref_gc_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions,
        seed: Some(23),
        windows: RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: base_chromosomes(&["chr1"]),
        blacklist: None,
        fragment_lengths: FragmentLengthArgs {
            min_fragment_length,
            max_fragment_length,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };
    run_ref_gc_bias(&ref_cfg, RunOptions::new_quiet())?;

    let gc_out_dir = out_dir.join(format!(
        "command_produced_gc_correction_from_reference_windows_len_{}-{}",
        min_fragment_length, max_fragment_length
    ));
    std::fs::create_dir_all(&gc_out_dir)?;
    let mut gc_cfg = GCConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: gc_out_dir.clone(),
            n_threads: 1,
        },
        reference_path.to_path_buf(),
        ref_gc_dir.path().join("ref_gc_package.zarr"),
        base_chromosomes(&["chr1"]),
    );
    configure_gc_bias_common(&mut gc_cfg);
    gc_cfg.set_min_gc_bin_mass(1.0);
    gc_cfg.set_min_length_bin_mass(0.0);
    gc_cfg.set_min_length_bin_width(1);
    run_gc_bias(&gc_cfg, RunOptions::new_quiet())?;

    Ok(gc_out_dir.join("gc_bias_correction.zarr"))
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
fn base_chromosomes(chromosome_names: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(
            chromosome_names
                .iter()
                .map(|name| name.to_string())
                .collect(),
        ),
        chromosomes_file: None,
    }
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
fn configure_gc_bias_common(gc_cfg: &mut GCConfig) {
    gc_cfg.set_min_mapq(0);
    gc_cfg.set_tile_size(1_000_000);
    gc_cfg.set_min_window_acgt_pct(0);
    gc_cfg.set_num_extreme_gc_bins(0);
    gc_cfg.set_num_short_length_bins(0);
    gc_cfg.outlier_method = OutlierMethodArg::None;
    gc_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
}
