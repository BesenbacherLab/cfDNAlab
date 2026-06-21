use super::*;
use crate::commands::ref_gc_bias::zarr::{ReferenceGCZarrPackage, write_reference_gc_package_zarr};
use crate::run_like_cli::{
    common::{
        ChromosomeArgs, FragmentLengthArgs, GCWindowsArgs, IOCArgs, LoggingArgs,
        Ref2BitRequiredArgs,
    },
    gc_bias::{GCConfig, OutlierMethodArg, OutlierScopeArg, run_gc_bias as run_gc_bias_command},
    ref_gc_bias::{RefGCBiasConfig, run_ref_gc_bias as run_ref_gc_bias_command},
};
use crate::testing::bam::paired_fragment as build_paired_fragment;
use crate::testing::{
    FragmentSpec, PairedFragmentSpec, TempBam, TempTwoBit, bam_from_fragments,
    bam_from_fragments_with_record_indexed_names, single_contig_inward_pair_bam,
    twobit_from_sequences, twobit_with_single_repeating_contig,
};
use crate::{
    commands::gc_bias::{
        binning::{BinnedAxis, bins_from_edges, compute_bin_edges},
        counting::{GCCounts, build_gc_prefixes, gc_percent_widths},
        load_reference_bias::{ReferenceGCMetadata, load_reference_gc_data},
        outliers::{
            OutlierAction, OutlierRule, OutlierScope, OutlierStats, apply_outliers_to_matrix,
            interpolated_quantile, outlier_bounds,
        },
        support_masking::build_extreme_bins_support_mask,
    },
    shared::interval::Interval,
};
use anyhow::Context;
use anyhow::Result;
use fxhash::FxHashMap;
use ndarray::array;
use ndarray_npy::read_npy;
use tempfile::{TempDir, tempdir};

/* Helpers */

const GC_COMMAND_F64_TOL: f64 = 1e-6;

fn run_gc_bias(config: &GCConfig) -> Result<()> {
    run_gc_bias_command(config, RunOptions::new_quiet()).map(|_| ())
}

fn run_ref_gc_bias(config: &RefGCBiasConfig) -> Result<()> {
    run_ref_gc_bias_command(config, RunOptions::new_quiet()).map(|_| ())
}

fn paired_fragment(start: i64, fragment_length: i64, read_length: i64) -> FragmentSpec {
    build_paired_fragment(&PairedFragmentSpec::new(
        0,
        start,
        fragment_length,
        read_length,
    ))
    .expect("paired fragment spec in GC-bias test fixture should be valid")
}

fn simple_reference_twobit() -> Result<TempTwoBit> {
    twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)
}

fn assert_gc_command_close(actual: f64, expected: f64, context: &str) {
    // The outlier helpers estimate quantiles/bounds in `f32` and only then write the matrix
    // back as `f64`, so the stable contract here is "matches the hand-derived value within the
    // command's float precision", not bit-exact `f64` arithmetic on the ideal decimal.
    assert!(
        (actual - expected).abs() <= GC_COMMAND_F64_TOL,
        "{context}: expected {expected}, got {actual}"
    );
}

fn make_gc_bias_cfg(
    bam_path: &std::path::Path,
    reference_path: &std::path::Path,
    ref_gc_dir: &std::path::Path,
    output_dir: &std::path::Path,
) -> GCConfig {
    let ioc = IOCArgs {
        bam: bam_path.to_path_buf(),
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
    };
    let mut cfg = GCConfig::new(
        ioc,
        reference_path.to_path_buf(),
        ref_gc_dir.join("ref_gc_package.zarr"),
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    cfg.set_min_mapq(0);
    cfg.set_tile_size(1_000_000);
    cfg.set_min_window_acgt_pct(0);
    cfg.set_save_intermediates(true);
    cfg
}

fn write_reference_package_for_single_length(
    reference_path: &std::path::Path,
    out_dir: &TempDir,
    fragment_length: u32,
    end_offset: u8,
) -> Result<()> {
    let cfg = RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs {
            ref_2bit: reference_path.to_path_buf(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        // These tests use small synthetic references, for example:
        // - `simple_reference_twobit()` is 256 bp
        // - the smallest custom GC-bias fixture here is 200 bp
        //
        // `ref-gc-bias` samples fragment starts from the set of valid start positions,
        // so `n_positions` must stay below that count:
        //   valid_starts = chrom_len - fragment_length + 1
        //
        // Using 100 keeps the helper valid for all current fixtures while still exercising
        // the full producer -> consumer path. The exact number is not part of the behavior
        // under test in this file.
        n_positions: 100,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: FragmentLengthArgs {
            min_fragment_length: fragment_length,
            max_fragment_length: fragment_length,
        },
        end_offset,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };
    run_ref_gc_bias(&cfg)
}

struct ReferencePackageFixture {
    version: Vec<u32>,
    skip_interpolation: Vec<bool>,
    smoothing_radius: Vec<u32>,
    smoothing_sigma: Vec<f64>,
    skip_smoothing: Vec<bool>,
    chromosomes: Vec<String>,
    length_range: [u32; 2],
    end_offset: u32,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
}

impl Default for ReferencePackageFixture {
    fn default() -> Self {
        Self {
            version: vec![GC_CORRECTION_SCHEMA_VERSION],
            skip_interpolation: vec![false],
            smoothing_radius: vec![2],
            smoothing_sigma: vec![0.55],
            skip_smoothing: vec![true],
            chromosomes: vec!["chr1".to_string()],
            length_range: [30, 31],
            end_offset: 10,
            reference_contig_footprint: Vec::new(),
        }
    }
}

fn write_reference_gc_package_fixture(
    out_dir: &std::path::Path,
    fixture: ReferencePackageFixture,
) -> Result<()> {
    let counts = array![[1.0_f64, 2.0_f64], [3.0_f64, 4.0_f64]];
    let support_unobservables = array![[true, false], [true, true]];
    let support_outliers = array![[true, true], [false, true]];
    let gc_percent_widths = array![[10_u16, 20_u16], [30_u16, 40_u16]];

    write_reference_zarr_arrays(
        out_dir,
        &counts,
        &support_unobservables,
        &support_outliers,
        &gc_percent_widths,
        fixture.length_range[0],
        fixture.length_range[1],
        fixture.end_offset,
        fixture.skip_interpolation.first().copied().unwrap_or(false),
        fixture.smoothing_radius.first().copied().unwrap_or(2) as u8,
        fixture.smoothing_sigma.first().copied().unwrap_or(0.55),
        fixture.skip_smoothing.first().copied().unwrap_or(true),
        &fixture.chromosomes,
        &fixture.reference_contig_footprint,
    )?;
    let package_path = out_dir.join("ref_gc_package.zarr");
    if fixture.version != vec![GC_CORRECTION_SCHEMA_VERSION] {
        tamper_reference_gc_root_attribute(
            &package_path,
            "cfdnalab_schema_version",
            serde_json::json!(fixture.version[0]),
        )?;
    }
    if fixture.skip_smoothing.len() != 1 {
        tamper_reference_gc_root_attribute(
            &package_path,
            "skip_smoothing",
            serde_json::json!(fixture.skip_smoothing),
        )?;
    }
    Ok(())
}

fn write_reference_gc_package_with_count_row_mismatch(out_dir: &std::path::Path) -> Result<()> {
    let counts = array![[1.0_f64, 2.0_f64], [3.0_f64, 4.0_f64], [5.0_f64, 6.0_f64]];
    let support_unobservables = array![[true, false], [true, true], [false, true]];
    let support_outliers = array![[true, true], [false, true], [true, false]];
    let gc_percent_widths = array![[10_u16, 20_u16], [30_u16, 40_u16], [50_u16, 60_u16]];

    write_reference_zarr_arrays(
        out_dir,
        &counts,
        &support_unobservables,
        &support_outliers,
        &gc_percent_widths,
        30,
        32,
        10,
        false,
        2,
        0.55,
        true,
        &["chr1".to_string()],
        &[],
    )?;
    tamper_reference_gc_array_shape(&out_dir.join("ref_gc_package.zarr"), "counts", &[2, 2])
}

fn write_reference_gc_package_with_shape_mismatch(out_dir: &std::path::Path) -> Result<()> {
    write_reference_gc_package_fixture(out_dir, ReferencePackageFixture::default())?;
    let package_path = out_dir.join("ref_gc_package.zarr");
    tamper_reference_gc_array_shape(&package_path, "support_mask_unobservables", &[1, 2])?;
    tamper_reference_gc_array_shape(&package_path, "support_mask_outliers", &[1, 2])
}

fn write_two_bin_reference_gc_package(
    out_dir: &std::path::Path,
    length_range: (u32, u32),
    chromosomes: &[&str],
    reference_contig_footprint: Vec<ContigFootprintEntry>,
) -> Result<()> {
    let n_lengths = (length_range.1 - length_range.0 + 1) as usize;
    let mut counts = ndarray::Array2::<f64>::zeros((n_lengths, 101));
    let mut support_outliers = ndarray::Array2::<bool>::from_elem((n_lengths, 101), false);
    let gc_percent_widths = gc_percent_widths(length_range.0 as usize, length_range.1 as usize, 0);
    let support_unobservables = gc_percent_widths.mapv(|width| width > 0);

    for row_idx in 0..n_lengths {
        counts[(row_idx, 0)] = 1.0;
        counts[(row_idx, 100)] = 1.0;
        support_outliers[(row_idx, 0)] = true;
        support_outliers[(row_idx, 100)] = true;
    }

    let chromosomes: Vec<String> = chromosomes.iter().map(|name| (*name).to_string()).collect();
    write_reference_zarr_arrays(
        out_dir,
        &counts,
        &support_unobservables,
        &support_outliers,
        &gc_percent_widths,
        length_range.0,
        length_range.1,
        0,
        true,
        2,
        0.55,
        true,
        &chromosomes,
        &reference_contig_footprint,
    )
}

#[allow(clippy::too_many_arguments)]
fn write_reference_zarr_arrays(
    out_dir: &std::path::Path,
    counts: &ndarray::Array2<f64>,
    support_unobservables: &ndarray::Array2<bool>,
    support_outliers: &ndarray::Array2<bool>,
    gc_percent_widths: &ndarray::Array2<u16>,
    min_fragment_length: u32,
    max_fragment_length: u32,
    end_offset: u32,
    skip_interpolation: bool,
    smoothing_radius: u8,
    smoothing_sigma: f64,
    skip_smoothing: bool,
    chromosomes: &[String],
    reference_contig_footprint: &[ContigFootprintEntry],
) -> Result<()> {
    let writer_end_offset = u8::try_from(end_offset).unwrap_or(0);
    write_reference_gc_package_zarr(
        &out_dir.join("ref_gc_package.zarr"),
        ReferenceGCZarrPackage {
            counts,
            support_unobservables,
            support_outliers,
            gc_percent_widths,
            length_min: min_fragment_length as usize,
            length_max: max_fragment_length as usize,
            end_offset: writer_end_offset,
            skip_interpolation,
            smoothing_radius,
            smoothing_sigma,
            skip_smoothing,
            chromosomes,
            reference_contig_footprint,
        },
    )?;
    if u8::try_from(end_offset).is_err() {
        tamper_reference_gc_root_attribute(
            &out_dir.join("ref_gc_package.zarr"),
            "end_offset",
            serde_json::json!(end_offset),
        )?;
    }
    Ok(())
}

fn tamper_reference_gc_root_attribute(
    package_path: &std::path::Path,
    field_name: &str,
    value: serde_json::Value,
) -> Result<()> {
    let metadata_path = package_path.join("zarr.json");
    let mut metadata: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&metadata_path)?)?;
    metadata
        .get_mut("attributes")
        .and_then(serde_json::Value::as_object_mut)
        .context("reference GC Zarr root metadata should have attributes")?
        .insert(field_name.to_string(), value);
    std::fs::write(&metadata_path, serde_json::to_vec_pretty(&metadata)?)?;
    Ok(())
}

fn tamper_reference_gc_array_shape(
    package_path: &std::path::Path,
    array_name: &str,
    shape: &[usize],
) -> Result<()> {
    let metadata_path = package_path.join(array_name).join("zarr.json");
    let mut metadata: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&metadata_path)?)?;
    *metadata
        .get_mut("shape")
        .with_context(|| format!("{array_name} Zarr metadata should have a shape field"))? =
        serde_json::json!(shape);
    std::fs::write(&metadata_path, serde_json::to_vec_pretty(&metadata)?)?;
    Ok(())
}

fn overwrite_reference_gc_length_axis(
    package_path: &std::path::Path,
    values: &[i32],
) -> Result<()> {
    let store = std::sync::Arc::new(zarrs::filesystem::FilesystemStore::new(package_path)?);
    let array = zarrs::array::Array::open(store, "/length")?;
    array
        .store_chunk(&[0], values)
        .context("overwrite malformed reference GC length axis")?;
    Ok(())
}

fn write_balanced_two_length_reference_gc_package(
    out_dir: &std::path::Path,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
) -> Result<()> {
    // Hand-built but still realistic reference package for run-level outlier tests.
    //
    // The package covers two fragment lengths, 10 and 11. We intentionally place reference
    // mass only at GC% 0 and GC% 100 for both rows, because the paired sample fixture also
    // lives only in those two extreme classes.
    //
    // But the metadata must still look like a real package:
    // - `support_mask_unobservables` follows the true theoretical GC%-reachability geometry
    //   for lengths 10 and 11
    // - `gc_percent_widths` uses the real rounding widths for those lengths
    // - `support_mask_outliers` marks only the empirically populated bins (0 and 100)
    //
    // That keeps the reference-side normalization denominator restricted to the two populated
    // bins while still preserving realistic theoretical metadata.
    //
    // With that support mask, once `gc-bias` bins the GC axis into `[0]` and `[1..100]`, the
    // reference-side binned rows are perfectly balanced:
    //   length 10 -> [1, 1]
    //   length 11 -> [1, 1]
    //
    // That keeps the run-level derivation simple: after per-length normalization the raw
    // correction matrix is driven entirely by the sample BAM's within-row imbalance.
    let mut counts = ndarray::Array2::<f64>::zeros((2, 101));
    counts[(0, 0)] = 1.0;
    counts[(0, 100)] = 1.0;
    counts[(1, 0)] = 1.0;
    counts[(1, 100)] = 1.0;

    let gc_percent_widths = gc_percent_widths(10, 11, 0);
    let support_unobservables = gc_percent_widths.mapv(|width| width > 0);
    let mut support_outliers = ndarray::Array2::<bool>::from_elem((2, 101), false);
    support_outliers[(0, 0)] = true;
    support_outliers[(0, 100)] = true;
    support_outliers[(1, 0)] = true;
    support_outliers[(1, 100)] = true;

    write_reference_zarr_arrays(
        out_dir,
        &counts,
        &support_unobservables,
        &support_outliers,
        &gc_percent_widths,
        10,
        11,
        0,
        true,
        2,
        0.55,
        true,
        &["chr1".to_string()],
        &reference_contig_footprint,
    )
}

fn write_three_bin_reference_gc_package(
    out_dir: &std::path::Path,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
) -> Result<()> {
    // One fragment length, with reference mass only at GC% 0, 50, and 100.
    //
    // As above, keep the metadata realistic:
    // - theoretical support and width correction follow the true length-10 rounding geometry
    // - empirical outlier support is restricted to the three populated GC bins
    //
    // That keeps the reference-side normalization easy to reason about for GC binning tests:
    // every empirically supported GC point starts with the same mass 1.0.
    let mut counts = ndarray::Array2::<f64>::zeros((1, 101));
    counts[(0, 0)] = 1.0;
    counts[(0, 50)] = 1.0;
    counts[(0, 100)] = 1.0;

    let gc_percent_widths = gc_percent_widths(10, 10, 0);
    let support_unobservables = gc_percent_widths.mapv(|width| width > 0);
    let mut support_outliers = ndarray::Array2::<bool>::from_elem((1, 101), false);
    support_outliers[(0, 0)] = true;
    support_outliers[(0, 50)] = true;
    support_outliers[(0, 100)] = true;

    write_reference_zarr_arrays(
        out_dir,
        &counts,
        &support_unobservables,
        &support_outliers,
        &gc_percent_widths,
        10,
        10,
        0,
        true,
        2,
        0.55,
        true,
        &["chr1".to_string()],
        &reference_contig_footprint,
    )
}

fn make_two_length_outlier_fixture() -> Result<(TempTwoBit, TempBam)> {
    let reference = twobit_from_sequences(
        "gc_bias_two_length_outlier_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;

    // Length-10 row:
    // - one pure-A fragment  -> GC%=0
    // - nine pure-C fragments -> GC%=100
    //
    // Length-11 row:
    // - five pure-A fragments  -> GC%=0
    // - five pure-C fragments -> GC%=100
    //
    // Later, after global mean-scaling and GC binning into `[0]` and `[1..100]`, the sample
    // binned rows are proportional to:
    //   length 10 -> [1, 9] -> normalized to [0.2, 1.8]
    //   length 11 -> [5, 5] -> normalized to [1.0, 1.0]
    let mut fragments = Vec::new();
    fragments.push(paired_fragment(10, 10, 5));
    for start in [110_i64, 120, 130, 140, 150, 160, 170, 180, 190] {
        fragments.push(paired_fragment(start, 10, 5));
    }
    for start in [20_i64, 30, 40, 50, 60] {
        fragments.push(paired_fragment(start, 11, 5));
    }
    for start in [120_i64, 130, 140, 150, 160] {
        fragments.push(paired_fragment(start, 11, 5));
    }

    // Several length-10 and length-11 fragments deliberately share the same left start. Use
    // strict BAM identity here so those stacked molecules do not collapse onto one qname.
    let bam = bam_from_fragments_with_record_indexed_names(
        "gc_bias_two_length_outlier_bam",
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
    )?;
    Ok((reference, bam))
}

fn make_two_length_low_mass_tail_fixture() -> Result<(TempTwoBit, TempBam)> {
    let reference = twobit_from_sequences(
        "gc_bias_two_length_low_mass_tail_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;

    // Length-10 row:
    // - one pure-A fragment  -> GC%=0
    // - nine pure-C fragments -> GC%=100
    //
    // Length-11 row:
    // - one pure-A fragment -> GC%=0
    // - one pure-C fragment -> GC%=100
    //
    // So before any length binning the per-row normalized correction rows are:
    //   length 10 -> [0.2, 1.8]
    //   length 11 -> [1.0, 1.0]
    //
    // But the total row masses are deliberately unequal:
    //   length 10 -> 10
    //   length 11 -> 2
    //
    // That makes length 11 a clean "low-mass tail" row for testing greedy length binning by
    // percentage mass.
    let mut fragments = Vec::new();
    fragments.push(paired_fragment(10, 10, 5));
    for start in [110_i64, 120, 130, 140, 150, 160, 170, 180, 190] {
        fragments.push(paired_fragment(start, 10, 5));
    }
    fragments.push(paired_fragment(20, 11, 5));
    fragments.push(paired_fragment(120, 11, 5));

    // The low-mass tail fixture reuses start 120 across two distinct fragment lengths, so
    // each synthetic fragment needs its own qname.
    let bam = bam_from_fragments_with_record_indexed_names(
        "gc_bias_two_length_low_mass_tail_bam",
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
    )?;
    Ok((reference, bam))
}

fn make_three_gc_bin_fixture() -> Result<(TempTwoBit, TempBam)> {
    let reference = twobit_from_sequences(
        "gc_bias_three_gc_bin_reference",
        vec![(
            "chr1".to_string(),
            format!(
                "{}{}{}",
                "A".repeat(40),
                "CCCCCAAAAA".repeat(4),
                "C".repeat(80)
            ),
        )],
    )?;

    // One fragment length only: 10 bp.
    //
    // Chosen starts land in three exact GC classes:
    // - start 10  -> AAAAAAAAAA   -> GC%=0
    // - start 40  -> CCCCCAAAAA   -> GC%=50
    // - start 100 -> CCCCCCCCCC   -> GC%=100
    //
    // Counts are deliberately imbalanced:
    //   GC%=0   -> 1 fragment
    //   GC%=50  -> 5 fragments
    //   GC%=100 -> 9 fragments
    let mut fragments = Vec::new();
    fragments.push(paired_fragment(10, 10, 5));
    for _ in 0..5 {
        fragments.push(paired_fragment(40, 10, 5));
    }
    for _ in 0..9 {
        fragments.push(paired_fragment(100, 10, 5));
    }

    // This fixture intentionally stacks five fragments at GC%=50 and nine at GC%=100. Use
    // strict identity so repeated starts still represent repeated molecules.
    let bam = bam_from_fragments_with_record_indexed_names(
        "gc_bias_three_gc_bin_bam",
        vec![("chr1".to_string(), 160)],
        fragments,
        Vec::new(),
    )?;
    Ok((reference, bam))
}

/* Tests */

#[test]
fn get_fragment_gc_uses_sequence_interval_as_prefix_origin() -> Result<()> {
    // Manual derivation:
    // - Prefixes are built from the loaded reference slice [900,961), not from chromosome
    //   origin 0.
    // - The sequence slice is 61 C bases, so fragment [900,961) has 61 GC bases.
    // - A local-origin bug would either ask the 61 bp prefix for [900,961) or otherwise fail
    //   to count the loaded slice as the fragment interval.
    let prefixes = build_gc_prefixes(&[b'C'; 61]);
    let fragment_interval = Interval::new(900_u64, 961_u64)?;
    let sequence_interval = Interval::new(900_u64, 961_u64)?;

    let gc_count = get_fragment_gc(fragment_interval, sequence_interval, 0, &prefixes, 0.0)?;

    assert_eq!(gc_count, Some(61));
    Ok(())
}

#[test]
fn get_fragment_gc_returns_none_when_fragment_is_outside_loaded_sequence() -> Result<()> {
    // Manual derivation:
    // - Prefixes cover only [900,961).
    // - Fragment [961,1022) is a valid reference interval, but its contracted GC window is
    //   completely outside the loaded sequence, so this is a legitimate missing correction
    //   rather than an indexing error.
    let prefixes = build_gc_prefixes(&[b'C'; 61]);
    let fragment_interval = Interval::new(961_u64, 1022_u64)?;
    let sequence_interval = Interval::new(900_u64, 961_u64)?;

    let gc_count = get_fragment_gc(fragment_interval, sequence_interval, 0, &prefixes, 0.0)?;

    assert_eq!(gc_count, None);
    Ok(())
}

#[test]
fn masks_extreme_gc_bins_per_side_in_square_matrix() {
    // Arrange: 6x6 matrix with two extreme GC bins on each side.
    let expected = array![
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
    ];

    // Act: build the support mask after binning.
    let mask = build_extreme_bins_support_mask((6, 6), 2, 0);

    // Assert: the central two GC bins remain supported across all lengths.
    assert_eq!(mask, expected);
}

#[test]
fn masks_shortest_length_bins_in_matrix() {
    // Arrange: 5x4 matrix with one shortest length bin masked.
    let expected = array![
        [false, false, false, false],
        [true, true, true, true],
        [true, true, true, true],
        [true, true, true, true],
        [true, true, true, true],
    ];

    // Act: build the support mask after binning.
    let mask = build_extreme_bins_support_mask((5, 4), 0, 1);

    // Assert: the central three length bins remain supported across all GC bins.
    assert_eq!(mask, expected);
}

#[test]
fn interpolates_masked_short_length_row() -> Result<()> {
    // Arrange: first length row is masked; other rows are supported.
    let mut matrix = array![
        [0.0_f64, 0.0_f64],
        [2.0_f64, 2.0_f64],
        [4.0_f64, 4.0_f64],
        [6.0_f64, 6.0_f64],
    ];
    let mask = build_extreme_bins_support_mask((4, 2), 0, 1);

    // Act: interpolate masked bins.
    interpolate_masked_corrections(&mut matrix, &mask)?;

    // Assert:
    // - the masked first row is filled from the nearest supported row
    // - the supported rows remain unchanged
    let expected = array![
        [2.0_f64, 2.0_f64],
        [2.0_f64, 2.0_f64],
        [4.0_f64, 4.0_f64],
        [6.0_f64, 6.0_f64],
    ];
    assert_eq!(matrix, expected);
    Ok(())
}

#[test]
fn round_trips_bins_to_edges_and_back() {
    // Arrange: build a simple BinnedAxis where bins group indices as [0-1], [2-4], and [5-7].
    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();
    let bins: [Vec<usize>; 3] = [vec![0, 1], vec![2, 3, 4], vec![5, 6, 7]];
    for (bin_idx, indices) in bins.iter().enumerate() {
        bin_to_indices.insert(bin_idx, indices.clone());
        for &idx in indices {
            index_to_bin.insert(idx, bin_idx);
        }
    }
    let axis = BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins: 3,
    };

    // Act: compute edges then reconstruct the bins.
    let edges = compute_bin_edges(&axis, 0, 7).expect("edges should be computed");
    let reconstructed_axis = bins_from_edges(edges.as_slice()).expect("rebuild should work");

    // Assert: the derived edges match the expected bin boundaries, and the reconstructed
    // axis matches the original bin layout.
    assert_eq!(edges, vec![0, 2, 5, 7]);
    assert_eq!(reconstructed_axis.num_bins, axis.num_bins);
    assert_eq!(reconstructed_axis.bin_to_indices, axis.bin_to_indices);
    assert_eq!(reconstructed_axis.index_to_bin, axis.index_to_bin);
}

#[test]
fn apply_outliers_per_length_winsorizes_rows() {
    let mut matrix = array![[1.0_f64, 2.0_f64, 100.0_f64], [1.0_f64, 5.0_f64, 6.0_f64]];
    let mask = array![[true, true, true], [true, true, true]];

    let stats = apply_outliers_to_matrix(
        &mut matrix,
        Some(&mask),
        OutlierScope::PerLength,
        OutlierRule::Quantile {
            lower: 0.0,
            upper: 0.5,
        },
        OutlierAction::Winsorize,
    );

    assert_eq!(matrix[[0, 0]], 1.0);
    assert_eq!(matrix[[0, 1]], 2.0);
    assert_eq!(matrix[[0, 2]], 2.0); // Clamped
    assert_eq!(matrix[[1, 0]], 1.0);
    assert_eq!(matrix[[1, 1]], 5.0);
    assert_eq!(matrix[[1, 2]], 5.0); // Clamped
    assert_eq!(
        stats,
        OutlierStats {
            total_examined: 6,
            total_outliers_handled: 2,
            unsupported_examined: 0,
            unsupported_outliers_handled: 0,
            hard_clamped: 0
        }
    );
}

#[test]
fn quantile_outliers_symmetry_clamps_extremes() {
    let mut matrix = array![[1.0_f64, 1.0_f64, 100.0_f64]];

    apply_outliers_to_matrix(
        &mut matrix,
        None,
        OutlierScope::Global,
        OutlierRule::Quantile {
            lower: 0.25,
            upper: 0.75,
        },
        OutlierAction::Winsorize,
    );

    assert!((matrix[[0, 0]] - 1.0).abs() < 1e-6);
    assert!((matrix[[0, 1]] - 1.0).abs() < 1e-6);
    assert!((matrix[[0, 2]] - 50.5).abs() < 1e-6);
}

#[test]
fn masked_cells_are_clamped_but_not_counted() {
    let mut matrix = array![[1.0_f64, 2.0_f64, 100.0_f64]];
    let mask = array![[true, true, false]];

    let stats = apply_outliers_to_matrix(
        &mut matrix,
        Some(&mask),
        OutlierScope::Global,
        OutlierRule::TukeyIqr { k: 1.0 },
        OutlierAction::Winsorize,
    );

    assert!((matrix[[0, 0]] - 1.0).abs() < 1e-6);
    assert!((matrix[[0, 1]] - 2.0).abs() < 1e-6);
    assert!((matrix[[0, 2]] - 2.25).abs() < 1e-6); // Unsupported cell still clamped
    assert_eq!(
        stats,
        OutlierStats {
            total_examined: 2,
            total_outliers_handled: 0,
            unsupported_examined: 1,
            unsupported_outliers_handled: 1,
            hard_clamped: 0
        }
    );
}

#[test]
fn interpolated_quantile_weights_neighbors_by_offset() {
    // Arrange
    let values = vec![0.0_f32, 10.0_f32, 20.0_f32, 30.0_f32, 40.0_f32];

    // Act
    let p_0 = interpolated_quantile(&values, 0.0);
    let p_05 = interpolated_quantile(&values, 0.5);
    let p_06 = interpolated_quantile(&values, 0.6);
    let p_08 = interpolated_quantile(&values, 0.8);
    let p_1 = interpolated_quantile(&values, 1.0);

    // Assert
    assert!((p_0 - 0.0).abs() < 1e-6);
    assert!((p_05 - 20.0).abs() < 1e-6);
    assert!((p_06 - 24.0).abs() < 1e-6); // 40% from 20 to 30
    assert!((p_08 - 32.0).abs() < 1e-6); // 20% from 30 to 40
    assert!((p_1 - 40.0).abs() < 1e-6);
}

#[test]
fn quantile_bounds_interpolate_between_indices() {
    // Arrange: Percentiles fall between indices, so bounds should blend neighbors
    let values = vec![0.0_f32, 10.0_f32, 20.0_f32, 30.0_f32, 40.0_f32];

    // Act: compute bounds for percentiles that require interpolation.
    let bounds = outlier_bounds(
        &values,
        OutlierRule::Quantile {
            lower: 0.6,
            upper: 0.8,
        },
    )
    .expect("quantile bounds should exist");

    // Assert: 0.6 is 40% from element 2 (20) to 3 (30); 0.8 is 20% from 3 (30) to 4 (40)
    assert!((bounds.0 - 24.0).abs() < 1e-6);
    assert!((bounds.1 - 32.0).abs() < 1e-6);
}

#[test]
fn iqr_outliers_per_length_clamps_high_values() {
    let mut matrix = array![[1.0_f64, 2.0_f64, 8.0_f64]];

    apply_outliers_to_matrix(
        &mut matrix,
        None,
        OutlierScope::PerLength,
        OutlierRule::TukeyIqr { k: 0.5 },
        OutlierAction::Winsorize,
    );

    assert!((matrix[[0, 0]] - 1.0).abs() < 1e-6);
    assert!((matrix[[0, 1]] - 2.0).abs() < 1e-6);
    assert!((matrix[[0, 2]] - 6.75).abs() < 1e-6);
}

#[test]
fn stddev_outliers_global_clamps_tail() {
    let mut matrix = array![[1.0_f64, 1.0_f64, 10.0_f64]];

    apply_outliers_to_matrix(
        &mut matrix,
        None,
        OutlierScope::Global,
        OutlierRule::StdDev { k: 1.0 },
        OutlierAction::Winsorize,
    );

    assert!((matrix[[0, 2]] - 8.2426405).abs() < 1e-5);
}

#[test]
fn mad_outliers_symmetrically_clamp() {
    let mut matrix = array![[1.0_f64, 2.0_f64, 3.0_f64, 9.0_f64]];

    apply_outliers_to_matrix(
        &mut matrix,
        None,
        OutlierScope::Global,
        OutlierRule::Mad { k: 1.0 },
        OutlierAction::Winsorize,
    );

    assert!((matrix[[0, 0]] - 1.0174).abs() < 1e-4);
    assert!((matrix[[0, 1]] - 2.0).abs() < 1e-6);
    assert!((matrix[[0, 2]] - 3.0).abs() < 1e-6);
    assert!((matrix[[0, 3]] - 3.9826).abs() < 1e-4);
}

#[test]
fn per_length_scope_differs_from_global() {
    let mut matrix = array![[1.0_f64, 100.0_f64], [1.0_f64, 1.0_f64]];

    apply_outliers_to_matrix(
        &mut matrix,
        None,
        OutlierScope::PerLength,
        OutlierRule::Quantile {
            lower: 0.25,
            upper: 0.75,
        },
        OutlierAction::Winsorize,
    );

    assert!((matrix[[0, 0]] - 25.75).abs() < 1e-6);
    assert!((matrix[[0, 1]] - 75.25).abs() < 1e-6);
    assert!((matrix[[1, 0]] - 1.0).abs() < 1e-6);
    assert!((matrix[[1, 1]] - 1.0).abs() < 1e-6);
}

#[test]
fn should_use_effective_length_when_binning_to_gc_percent_with_end_offset() {
    // Arrange: one 30bp fragment with 20 GC bases after trimming 5bp from each end
    let mut counts = GCCounts::new(30, 30, 5, (0, 0)).expect("counts init");
    counts.incr(30, 20);

    // Act
    let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");

    // Assert: value lands in the 100% bin, not in the 67% bin (which used full length)
    assert_eq!(grid[(0, 100)], 1.0);
    assert_eq!(grid[(0, 67)], 0.0);
}

#[test]
fn should_not_smooth_into_gc_counts_beyond_effective_length() {
    // Arrange: length=6, end_offset=2 -> effective length is 2bp, so gc>2 is unreachable.
    let mut counts = GCCounts::new(6, 6, 2, (0, 0)).expect("counts init");
    counts.set(6, 2, 10.0);

    // Act: smooth only the reachable portion of the row.
    counts
        .smooth_length_rows_in_place(1.0, 1)
        .expect("smoothing should succeed for valid sigma and radius");

    // Assert: unreachable GC counts are absent and storage matches the effective length.
    assert!(counts.get(6, 3).is_none());
    assert_eq!(counts.borrow_raw_counts().len(), 3);
}

#[test]
fn should_place_gc_counts_in_matching_percent_bins() {
    // Arrange: one length row with distinct weights per GC count.
    let mut counts = GCCounts::new(10, 10, 0, (0, 0)).expect("counts init");
    for gc in 0..=10 {
        counts.set(10, gc, (gc + 1) as f64); // unique weight per bin
    }

    // Act
    let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");
    let row = grid.row(0);

    // Assert: each GC count lands in its integer percent bin.
    for gc in 0..=10 {
        let pct_bin = (gc * 10) as usize;
        assert!(
            (row[pct_bin] - (gc + 1) as f64).abs() < 1e-12,
            "gc {} expected at pct {}, got {}",
            gc,
            pct_bin,
            row[pct_bin]
        );
    }
}

#[test]
fn should_round_half_up_for_fractional_percentages() {
    // Arrange: length=3 has fractional percentages for gc=1 and gc=2.
    let mut counts = GCCounts::new(3, 3, 0, (0, 0)).expect("counts init");
    counts.set(3, 1, 2.0); // 33.3...% -> 33 via half-up
    counts.set(3, 2, 3.0); // 66.6...% -> 67 via half-up

    // Act
    let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");
    let row = grid.row(0);

    // Assert: derive the half-up bins explicitly
    // calculate_gc_bin does round_half_up(100 * gc / effective_length) via (100 * gc + len/2) / len
    // Effective length is 3 (no end trimming)
    // gc=1 -> (100 * 1 + 3/2) / 3 = (100 + 1) / 3 = 33
    // gc=2 -> (100 * 2 + 3/2) / 3 = (200 + 1) / 3 = 67
    // Mass must land only in those bins
    for (idx, &val) in row.iter().enumerate() {
        match idx {
            33 => assert!(
                (val - 2.0).abs() < 1e-12,
                "bin {} expected 2.0, got {}",
                idx,
                val
            ),
            67 => assert!(
                (val - 3.0).abs() < 1e-12,
                "bin {} expected 3.0, got {}",
                idx,
                val
            ),
            _ => assert!(val.abs() < 1e-12, "bin {} expected 0, got {}", idx, val),
        }
    }
}

#[test]
fn should_propagate_acgt_totals_and_length_metadata() {
    // Arrange
    let mut counts = GCCounts::new(5, 6, 1, (8, 12)).expect("counts init");
    counts.set(5, 2, 1.0);
    counts.set(6, 3, 2.0);

    // Act
    let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");

    // Assert: shapes match the two length bins and 101 GC bins
    assert_eq!(grid.nrows(), 2);
    assert_eq!(grid.ncols(), 101);

    let row_len5 = grid.row(0);
    let row_len6 = grid.row(1);
    // Derivation with end offsets
    // End offset is 1 so effective length = length - 2
    // calculate_gc_bin uses (100 * gc + eff_len/2) / eff_len
    // len5 -> eff3: gc=2 gives (100 * 2 + 3/2) / 3 = (200 + 1) / 3 = 67
    // len6 -> eff4: gc=3 gives (100 * 3 + 4/2) / 4 = (300 + 2) / 4 = 75
    assert!((row_len5[67] - 1.0).abs() < 1e-12);
    assert!((row_len6[75] - 2.0).abs() < 1e-12);
}

#[test]
fn reports_offsets_based_on_effective_length() {
    // length_min=3, length_max=5, end_offset=1 -> effective lengths: 1,2,3
    let counts = GCCounts::new(3, 5, 1, (0, 0)).expect("init counts");

    let bounds_len3 = counts.length_bounds(3).expect("len3 bounds");
    let bounds_len4 = counts.length_bounds(4).expect("len4 bounds");
    let bounds_len5 = counts.length_bounds(5).expect("len5 bounds");

    assert_eq!(bounds_len3, (0, 2)); // size 2 for effective len 1 (gc 0..1)
    assert_eq!(bounds_len4, (2, 5)); // size 3 for effective len 2 (gc 0..2)
    assert_eq!(bounds_len5, (5, 9)); // size 4 for effective len 3 (gc 0..3)

    // Verify the slice lengths match the effective length + 1
    assert_eq!(bounds_len3.1 - bounds_len3.0, 2);
    assert_eq!(bounds_len4.1 - bounds_len4.0, 3);
    assert_eq!(bounds_len5.1 - bounds_len5.0, 4);
}

#[test]
fn row_bounds_errors_outside_length_range() {
    let counts = GCCounts::new(10, 12, 0, (0, 0)).expect("init counts");
    assert!(counts.length_bounds(9).is_err());
    assert!(counts.length_bounds(13).is_err());
}

#[test]
fn leaves_zero_rows_untouched_in_mean_scaling() {
    // Arrange: first length row has no mass; second has values that should be mean-scaled.
    let counts = array![[0.0, 0.0], [2.0, 4.0]];
    let mask = array![[true, true], [true, true]];

    // Act
    let scaled = mean_scale_per_length_array(&counts, 0.0, Some(&mask));

    // Assert: empty row stays zero; non-empty row divides by its mean (3.0).
    assert!(
        scaled.row(0).iter().all(|&value| value == 0.0),
        "zero row should remain zero after scaling"
    );
    assert!((scaled[(1, 0)] - 2.0 / 3.0).abs() < 1e-12);
    assert!((scaled[(1, 1)] - 4.0 / 3.0).abs() < 1e-12);
}

#[test]
fn save_intermediates_writes_expected_sequence_and_mean_scaled_average_counts() -> Result<()> {
    // Arrange:
    // Use a single global window and a reference package that already disables smoothing and
    // interpolation. In that configuration `gc-bias` should save exactly six intermediate
    // arrays:
    //   0 avg_cfdna_counts
    //   1 normalized_avg_cfdna_counts
    //   2 binned_ref_counts
    //   3 binned_cfdna_counts
    //   4 normalized_binned_cfdna_counts
    //   5 normalized_binned_ref_counts
    //
    // The strongest low-level coherence check in this branch is the first normalization step:
    // `normalized_avg_cfdna_counts` must equal `avg_cfdna_counts / supported_mean`, where the
    // mean is taken only over the reference outlier-support mask.
    let bam = single_contig_inward_pair_bam()?;
    let reference = simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;

    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_save_intermediates(true);

    // Act
    run_gc_bias(&cfg)?;

    // Assert:
    // No interpolation/smoothing intermediates should exist for this reference package, so the
    // numbering must stay dense across exactly six saved arrays.
    let mut intermediate_files: Vec<String> = std::fs::read_dir(out_dir.path())?
        .filter_map(|entry| entry.ok())
        .filter_map(|entry| {
            let name = entry.file_name().into_string().ok()?;
            if name.starts_with("gc_bias.") && name.ends_with(".npy") {
                Some(name)
            } else {
                None
            }
        })
        .collect();
    intermediate_files.sort();
    assert_eq!(
        intermediate_files,
        vec![
            "gc_bias.avg_cfdna_counts.0.npy".to_string(),
            "gc_bias.binned_cfdna_counts.3.npy".to_string(),
            "gc_bias.binned_ref_counts.2.npy".to_string(),
            "gc_bias.normalized_avg_cfdna_counts.1.npy".to_string(),
            "gc_bias.normalized_binned_cfdna_counts.4.npy".to_string(),
            "gc_bias.normalized_binned_ref_counts.5.npy".to_string(),
        ]
    );

    let avg_counts: ndarray::Array2<f64> =
        read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
    let normalized_avg: ndarray::Array2<f64> = read_npy(
        out_dir
            .path()
            .join("gc_bias.normalized_avg_cfdna_counts.1.npy"),
    )?;
    let reference_data = load_reference_gc_data(&ref_gc_dir.path().join("ref_gc_package.zarr"))?;

    // The support mask defines exactly which cells contribute to the mean-scaling denominator.
    let mut supported_sum = 0.0_f64;
    let mut supported_count = 0usize;
    for (value, supported) in avg_counts
        .iter()
        .zip(reference_data.outliers_support_mask.iter())
    {
        if *supported {
            supported_sum += *value;
            supported_count += 1;
        }
    }
    assert!(
        supported_count > 0,
        "fixture must have supported reference bins"
    );
    let supported_mean = supported_sum / supported_count as f64;
    assert!(
        supported_mean > 0.0,
        "supported mean must be positive for mean scaling"
    );

    for ((row_idx, col_idx), avg_value) in avg_counts.indexed_iter() {
        let expected = *avg_value / supported_mean;
        let actual = normalized_avg[(row_idx, col_idx)];
        assert!(
            (actual - expected).abs() < 1e-12,
            "normalized avg mismatch at ({row_idx}, {col_idx}): expected {expected}, got {actual}"
        );
    }

    Ok(())
}

#[test]
fn rejects_correction_package_components_with_invalid_weights() -> Result<()> {
    // Arrange: one length bin spanning 30..=31 and one GC bin spanning 0..=100.
    // The package writer should reject invalid final correction weights before an NPZ can be
    // written, and the error should identify the offending bin.
    let length_bins = bins_from_edges(&[30, 31])?;
    let gc_bins = bins_from_edges(&[0, 100])?;
    let reference_metadata = ReferenceGCMetadata {
        min_fragment_length: 30,
        max_fragment_length: 31,
        end_offset: 10,
        chromosomes: vec!["chr1".to_string()],
        reference_contig_footprint: Vec::new(),
        skip_interpolation: false,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
    };

    for invalid_weight in [f64::NAN, f64::INFINITY, -0.25] {
        // Act
        let error = GCCorrectionPackage::from_components(
            GC_CORRECTION_SCHEMA_VERSION,
            &length_bins,
            &gc_bins,
            array![[invalid_weight]],
            array![1.0_f64],
            &reference_metadata,
        )
        .expect_err("invalid correction weight should fail package construction");

        // Assert
        let message = error.to_string();
        assert!(
            message.contains("GC correction matrix contains invalid weight"),
            "unexpected error message: {message}"
        );
        assert!(
            message.contains("length bin 0 [30-31], GC bin 0 [0-100]"),
            "unexpected error message: {message}"
        );
    }

    Ok(())
}

#[test]
fn multi_chromosome_cross_tile_windows_match_hand_derived_counts_in_each_window_mode() -> Result<()>
{
    // Arrange:
    // Build asymmetric chromosomes with 100 bp logical windows:
    // - chr1 is all A, so every length-10 fragment contributes GC% 0
    // - chr2 is all C, so every length-10 fragment contributes GC% 100
    // - chr1 has one counted window
    // - chr2 has two counted windows
    //
    // The intended behavior is independent of tile layout. We use `tile_size = 75`, so all
    // 100 bp fixed-size and BED windows cross a tile boundary:
    // - [0,100) crosses tiles [0,75) and [75,150)
    // - [100,200) crosses tiles [75,150) and [150,200)
    //
    // This same fixture is run through all public windowing modes:
    // - global: one chromosome-wide window per chromosome, no per-window mean scaling
    // - by-size: one 100 bp window on chr1 and two 100 bp windows on chr2
    // - by-BED: the same three 100 bp windows, explicitly listed
    //
    // Hand-derived expected saved `avg_cfdna_counts`:
    // - global mode keeps raw counts, so the chr1 fragment gives GC% 0 = 1 and the two chr2
    //   fragments give GC% 100 = 2
    // - by-size and by-BED scale each pure counted window to 11 at the observed GC bin,
    //   because a 10 bp fragment has 11 reachable GC-count states, then average across three
    //   counted windows:
    //     GC% 0   -> 11 / 3
    //     GC% 100 -> 22 / 3
    let reference = twobit_from_sequences(
        "gc_bias_multi_chr_cross_tile_reference",
        vec![
            ("chr1".to_string(), "A".repeat(100)),
            ("chr2".to_string(), "C".repeat(200)),
        ],
    )?;

    let fragment_on_chromosome = |tid: usize, start: i64| {
        let mut fragment = paired_fragment(start, 10, 5);
        fragment.forward.tid = tid;
        fragment.reverse.tid = tid;
        fragment.forward.mate_tid = Some(tid);
        fragment.reverse.mate_tid = Some(tid);
        fragment
    };
    let bam = bam_from_fragments(
        "gc_bias_multi_chr_cross_tile_bam",
        vec![("chr1".to_string(), 100), ("chr2".to_string(), 200)],
        vec![
            fragment_on_chromosome(0, 10),
            fragment_on_chromosome(1, 10),
            fragment_on_chromosome(1, 110),
        ],
        Vec::new(),
    )?;

    let ref_gc_dir = TempDir::new()?;
    write_two_bin_reference_gc_package(
        ref_gc_dir.path(),
        (10, 10),
        &["chr1", "chr2"],
        twobit_contig_footprint(&reference.path)?,
    )?;

    let bed_dir = TempDir::new()?;
    let bed_path = bed_dir.path().join("windows.bed");
    std::fs::write(&bed_path, "chr1\t0\t100\nchr2\t0\t100\nchr2\t100\t200\n")?;

    let cases = vec![
        (
            "global",
            GCWindowsArgs {
                by_size: None,
                by_bed: None,
                global: true,
            },
            1.0,
            2.0,
        ),
        (
            "by-size",
            GCWindowsArgs {
                by_size: Some(100),
                by_bed: None,
                global: false,
            },
            11.0 / 3.0,
            22.0 / 3.0,
        ),
        (
            "by-BED",
            GCWindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
                global: false,
            },
            11.0 / 3.0,
            22.0 / 3.0,
        ),
    ];

    for (case_name, windows, expected_gc0, expected_gc100) in cases {
        let out_dir = TempDir::new()?;
        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        };
        let mut cfg = GCConfig::new(
            ioc,
            reference.path.clone(),
            ref_gc_dir.path().join("ref_gc_package.zarr"),
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_windows(windows);
        cfg.set_min_mapq(0);
        cfg.set_tile_size(75);
        cfg.set_min_window_acgt_pct(0);
        cfg.set_save_intermediates(true);

        // Act
        run_gc_bias(&cfg).with_context(|| {
            format!("{case_name} multi-chromosome cross-tile gc-bias run failed")
        })?;

        // Assert
        let avg_counts: ndarray::Array2<f64> =
            read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(avg_counts.dim(), (1, 101), "{case_name}");
        for (gc_pct, &value) in avg_counts.row(0).iter().enumerate() {
            match gc_pct {
                0 => assert!(
                    (value - expected_gc0).abs() < 1e-12,
                    "{case_name}: expected {expected_gc0} at GC% 0, got {value}"
                ),
                100 => assert!(
                    (value - expected_gc100).abs() < 1e-12,
                    "{case_name}: expected {expected_gc100} at GC% 100, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "{case_name}: expected no mass outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }
    }

    Ok(())
}

#[test]
fn gc_bias_run_rejects_reference_package_with_non_scalar_metadata_array() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let reference = simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_gc_package_fixture(
        ref_gc_dir.path(),
        ReferencePackageFixture {
            skip_smoothing: vec![true, false],
            ..ReferencePackageFixture::default()
        },
    )?;
    let out_dir = TempDir::new()?;
    let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

    // Manual expectations:
    // - The reference package is malformed before the command ever starts sample counting:
    //   `skip_smoothing` is written as a length-2 array instead of a scalar metadata field.
    // - `gc-bias` loads the reference package at the start of `run()`, so the correct
    //   behavior is an immediate loader failure with the scalar-shape guardrail message.
    let err = run_gc_bias(&cfg).expect_err("non-scalar reference metadata should fail");
    let msg = err.to_string();
    assert!(
        msg.contains("skip_smoothing") && msg.contains("must be a bool"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn gc_bias_run_rejects_reference_package_with_schema_version_mismatch() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let reference = simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_gc_package_fixture(
        ref_gc_dir.path(),
        ReferencePackageFixture {
            version: vec![GC_CORRECTION_SCHEMA_VERSION + 1],
            ..ReferencePackageFixture::default()
        },
    )?;
    let out_dir = TempDir::new()?;
    let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

    // Manual expectations:
    // - The package schema version is intentionally incompatible.
    // - `gc-bias` must fail while loading the reference-GC artifact, before producing any
    //   sample-side intermediates or correction output.
    let err = run_gc_bias(&cfg).expect_err("schema version mismatch should fail");
    let msg = err.to_string();
    assert!(
        msg.contains("Reference GC package schema version mismatch"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn gc_bias_run_rejects_reference_package_with_different_chromosomes() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let reference = simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_gc_package_fixture(
        ref_gc_dir.path(),
        ReferencePackageFixture {
            chromosomes: vec!["chr2".to_string()],
            ..ReferencePackageFixture::default()
        },
    )?;
    let out_dir = TempDir::new()?;
    let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

    // Manual expectations:
    // - The run selects `chr1` from the BAM through `make_gc_bias_cfg()`.
    // - The hand-written reference package claims it was built for `chr2`.
    // - `gc-bias` must reject this before using the reference counts for correction.
    let err = run_gc_bias(&cfg).expect_err("chromosome-mismatched reference package should fail");
    let msg = err.to_string();
    assert!(
        msg.contains("built for chromosomes [chr2]") && msg.contains("selected [chr1]"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn gc_bias_run_rejects_by_size_smaller_than_reference_max_fragment_length() -> Result<()> {
    // Human verification status: verified
    // Manual expectations:
    // - The reference package fixture covers fragment lengths 30..=31, so gc-bias inherits
    //   max_fragment_length = 31.
    // - A fixed window size of 30 cannot preserve the fixed-size two-buffer counting invariant.
    // - The command should fail directly after loading the reference package, before creating
    //   the output directory.
    let bam = single_contig_inward_pair_bam()?;
    let reference = simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_gc_package_fixture(ref_gc_dir.path(), ReferencePackageFixture::default())?;
    let output_parent = TempDir::new()?;
    let output_dir = output_parent.path().join("not_created");
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), &output_dir);
    cfg.set_windows(GCWindowsArgs {
        by_size: Some(30),
        by_bed: None,
        global: false,
    });

    let err = run_gc_bias(&cfg).expect_err("too-small fixed window should fail");
    let msg = err.to_string();

    assert!(
        msg.contains("--by-size (30) must be >= max fragment length from --ref-gc-file (31)"),
        "unexpected error message: {msg}"
    );
    assert!(
        !output_dir.exists(),
        "fixed-window validation should fail before creating output_dir"
    );

    Ok(())
}

#[test]
fn loads_versioned_reference_gc_package() -> Result<()> {
    // Arrange: write a minimal reference package with the current schema version and scalar
    // metadata fields.
    let tmp = tempdir()?;
    write_reference_gc_package_fixture(tmp.path(), ReferencePackageFixture::default())?;

    // Act
    let loaded = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))?;

    // Assert: arrays and scalar metadata survive round-trip exactly.
    assert_eq!(
        loaded.counts,
        array![[1.0_f64, 2.0_f64], [3.0_f64, 4.0_f64]]
    );
    assert_eq!(
        loaded.unobservables_support_mask,
        array![[true, false], [true, true]]
    );
    assert_eq!(
        loaded.outliers_support_mask,
        array![[true, true], [false, true]]
    );
    assert_eq!(
        loaded.gc_percent_widths,
        array![[10_u16, 20_u16], [30_u16, 40_u16]]
    );
    assert_eq!(loaded.metadata.min_fragment_length, 30);
    assert_eq!(loaded.metadata.max_fragment_length, 31);
    assert_eq!(loaded.metadata.end_offset, 10);
    assert_eq!(loaded.metadata.chromosomes, vec!["chr1".to_string()]);
    assert_eq!(loaded.metadata.reference_contig_footprint, Vec::new());
    assert!(!loaded.metadata.skip_interpolation);
    assert_eq!(loaded.metadata.smoothing_radius, 2);
    assert!((loaded.metadata.smoothing_sigma - 0.55).abs() < 1e-12);
    assert!(loaded.metadata.skip_smoothing);
    Ok(())
}

#[test]
fn rejects_reference_gc_package_with_non_scalar_metadata_array() -> Result<()> {
    // Arrange: `skip_smoothing` is written with two values. This should fail cleanly instead of
    // indexing `[0]` and panicking.
    let tmp = tempdir()?;
    write_reference_gc_package_fixture(
        tmp.path(),
        ReferencePackageFixture {
            skip_smoothing: vec![true, false],
            ..ReferencePackageFixture::default()
        },
    )?;

    // Act
    let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))
        .expect_err("expected scalar-length error");

    // Assert
    assert!(
        error.to_string().contains("skip_smoothing")
            && error.to_string().contains("must be a bool")
    );
    Ok(())
}

#[test]
fn rejects_reference_gc_package_with_schema_version_mismatch() -> Result<()> {
    // Arrange: same package shape, but an incompatible version number.
    let tmp = tempdir()?;
    write_reference_gc_package_fixture(
        tmp.path(),
        ReferencePackageFixture {
            version: vec![GC_CORRECTION_SCHEMA_VERSION + 1],
            ..ReferencePackageFixture::default()
        },
    )?;

    // Act
    let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))
        .expect_err("expected schema version mismatch");

    // Assert
    assert!(
        error
            .to_string()
            .contains("Reference GC package schema version mismatch")
    );
    Ok(())
}

#[test]
fn rejects_reference_gc_package_with_row_count_mismatched_to_length_axis() -> Result<()> {
    // Arrange: the package has two count rows, but the Zarr length axis names three concrete
    // fragment length rows: 30, 31, and 32.
    let tmp = tempdir()?;
    write_reference_gc_package_with_count_row_mismatch(tmp.path())?;

    // Act
    let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))
        .expect_err("expected row-count mismatch");

    // Assert
    assert!(
        error
            .to_string()
            .contains("row count 2 does not match length axis [30, 32] (expected 3)")
    );
    Ok(())
}

#[test]
fn rejects_reference_gc_package_with_inverted_length_axis() -> Result<()> {
    // Arrange: `[31, 30]` cannot name an ordered inclusive set of length rows.
    let tmp = tempdir()?;
    write_reference_gc_package_fixture(tmp.path(), ReferencePackageFixture::default())?;
    overwrite_reference_gc_length_axis(&tmp.path().join("ref_gc_package.zarr"), &[31, 30])?;

    // Act
    let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))
        .expect_err("expected inverted length axis error");

    // Assert
    assert!(
        error
            .to_string()
            .contains("length axis must contain contiguous integer values")
    );
    Ok(())
}

#[test]
fn rejects_reference_gc_package_with_out_of_range_end_offset() -> Result<()> {
    // Arrange: the writer stores `end_offset` as u32, but the command metadata model uses u8.
    // Values above 255 must not be truncated.
    let tmp = tempdir()?;
    write_reference_gc_package_fixture(
        tmp.path(),
        ReferencePackageFixture {
            end_offset: 300,
            ..ReferencePackageFixture::default()
        },
    )?;

    // Act
    let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))
        .expect_err("expected out-of-range end_offset error");

    // Assert
    assert!(
        error
            .to_string()
            .contains("end_offset in reference GC package must fit in u8")
    );
    Ok(())
}

#[test]
fn rejects_reference_gc_package_with_too_short_effective_minimum_length() -> Result<()> {
    // Arrange: `min_fragment_length = 30` and `end_offset = 11` leaves only 8 bp after
    // trimming both ends. `ref-gc-bias` refuses to write such a package, so the loader should
    // reject one if it appears on disk.
    let tmp = tempdir()?;
    write_reference_gc_package_fixture(
        tmp.path(),
        ReferencePackageFixture {
            end_offset: 11,
            ..ReferencePackageFixture::default()
        },
    )?;

    // Act
    let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.zarr"))
        .expect_err("expected invalid effective minimum length");

    // Assert
    let expected_message = format!(
        "min_fragment_length (30) - 2 * end_offset (11) must be >= {}",
        MIN_ACGT_BASES_FOR_GC_FRACTION
    );
    assert!(
        error.to_string().contains(&expected_message),
        "unexpected error message: {error}"
    );
    Ok(())
}

#[test]
fn gc_bias_run_rejects_reference_package_with_incompatible_support_mask_shape() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let reference = simple_reference_twobit()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_gc_package_with_shape_mismatch(ref_gc_dir.path())?;
    let out_dir = TempDir::new()?;
    let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

    // Manual expectations:
    // - The reference package is syntactically present, but the support masks are the wrong
    //   shape for the count matrix:
    //     counts                     = (2, 2)
    //     support_mask_unobservables = (1, 2)
    //     support_mask_outliers      = (1, 2)
    // - `gc-bias` must reject this immediately while loading the artifact rather than trying
    //   to continue with inconsistent masking semantics.
    let err = run_gc_bias(&cfg).expect_err("shape-mismatched reference package should fail");
    let msg = err.to_string();
    assert!(
        msg.contains("Reference counts") && msg.contains("incompatible shapes"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn quantile_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
    // Arrange:
    // Use the synthetic two-length reference package and BAM fixture defined above.
    //
    // Reference package after GC binning:
    // - GC bins are `[0]` and `[1..100]`
    // - both length rows are balanced, so normalized reference rows are:
    //     length 10 -> [1.0, 1.0]
    //     length 11 -> [1.0, 1.0]
    //
    // Sample BAM after the same GC binning:
    // - length 10 raw row is [1, 9], so normalized row is [0.2, 1.8]
    // - length 11 raw row is [5, 5], so normalized row is [1.0, 1.0]
    //
    // Therefore, before outlier handling, the raw correction matrix is:
    //   [[0.2, 1.8],
    //    [1.0, 1.0]]
    //
    // `--outlier-method none`:
    // - no winsorization
    // - no hard clamp, since all values already lie inside [0.1, 10]
    // - inversion gives:
    //     [[5.0, 5/9],
    //      [1.0, 1.0]]
    //
    // `--outlier-method quantile --outlier-scope per-length --outlier-quantiles 0.25,0.75`:
    // - length-10 row sorted values are [0.2, 1.8]
    // - with the command's linear interpolation:
    //     Q25 = 0.2 + 0.25 * (1.8 - 0.2) = 0.6
    //     Q75 = 0.2 + 0.75 * (1.8 - 0.2) = 1.4
    // - winsorized length-10 row becomes [0.6, 1.4]
    // - length-11 row stays [1.0, 1.0]
    // - inversion gives:
    //     [[5/3, 5/7],
    //      [1.0, 1.0]]
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_none = TempDir::new()?;
    let out_quantile = TempDir::new()?;

    let mut none_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        out_none.path(),
    );
    none_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    none_cfg.set_min_length_bin_mass(0.0);
    none_cfg.set_min_length_bin_width(1);
    none_cfg.set_min_gc_bin_mass(1.0);
    none_cfg.set_num_extreme_gc_bins(0);
    none_cfg.set_num_short_length_bins(0);
    none_cfg.outlier_method = OutlierMethodArg::None;

    let mut quantile_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        out_quantile.path(),
    );
    quantile_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    quantile_cfg.set_min_length_bin_mass(0.0);
    quantile_cfg.set_min_length_bin_width(1);
    quantile_cfg.set_min_gc_bin_mass(1.0);
    quantile_cfg.set_num_extreme_gc_bins(0);
    quantile_cfg.set_num_short_length_bins(0);
    quantile_cfg.outlier_method = OutlierMethodArg::Quantile;
    quantile_cfg.outlier_scope = OutlierScopeArg::PerLength;
    quantile_cfg.outlier_quantiles = vec![0.25, 0.75];

    // Act
    run_gc_bias(&none_cfg)?;
    run_gc_bias(&quantile_cfg)?;

    // Assert
    let package_none =
        GCCorrectionPackage::from_file(out_none.path().join("gc_bias_correction.zarr"))?;
    let package_quantile =
        GCCorrectionPackage::from_file(out_quantile.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(package_none.correction_matrix.dim(), (2, 2));
    assert_eq!(package_quantile.correction_matrix.dim(), (2, 2));
    assert_eq!(package_none.gc_edges, vec![0, 1, 100]);
    assert_eq!(package_quantile.gc_edges, vec![0, 1, 100]);
    assert_eq!(package_none.length_bin_frequencies.len(), 2);
    assert_eq!(package_quantile.length_bin_frequencies.len(), 2);
    assert!((package_none.length_bin_frequencies[0] - 0.5).abs() < 1e-12);
    assert!((package_none.length_bin_frequencies[1] - 0.5).abs() < 1e-12);
    assert!((package_quantile.length_bin_frequencies[0] - 0.5).abs() < 1e-12);
    assert!((package_quantile.length_bin_frequencies[1] - 0.5).abs() < 1e-12);

    assert!((package_none.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
    assert!((package_none.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
    assert!((package_none.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
    assert!((package_none.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

    assert_gc_command_close(
        package_quantile.correction_matrix[(0, 0)],
        5.0 / 3.0,
        "quantile row0 col0",
    );
    assert_gc_command_close(
        package_quantile.correction_matrix[(0, 1)],
        5.0 / 7.0,
        "quantile row0 col1",
    );
    assert_gc_command_close(
        package_quantile.correction_matrix[(1, 0)],
        1.0,
        "quantile row1 col0",
    );
    assert_gc_command_close(
        package_quantile.correction_matrix[(1, 1)],
        1.0,
        "quantile row1 col1",
    );

    Ok(())
}

#[test]
fn quantile_outlier_scope_global_differs_from_per_length_in_real_command() -> Result<()> {
    // Arrange:
    // Reuse the same raw correction matrix derivation as the previous test:
    //   [[0.2, 1.8],
    //    [1.0, 1.0]]
    //
    // With `quantile` and explicit `Q25/Q75 = 0.25/0.75`:
    //
    // Per-length scope:
    // - length 10 row [0.2, 1.8] -> [0.6, 1.4]
    // - length 11 row [1.0, 1.0] -> [1.0, 1.0]
    // - final weights:
    //     [[5/3, 5/7],
    //      [1.0, 1.0]]
    //
    // Global scope:
    // - full sorted matrix values are [0.2, 1.0, 1.0, 1.8]
    // - linear interpolation gives:
    //     Q25 = 0.2 + 0.75 * (1.0 - 0.2) = 0.8
    //     Q75 = 1.0 + 0.25 * (1.8 - 1.0) = 1.2
    // - winsorized matrix becomes:
    //     [[0.8, 1.2],
    //      [1.0, 1.0]]
    // - final weights are therefore:
    //     [[1.25, 5/6],
    //      [1.0, 1.0]]
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_per_length = TempDir::new()?;
    let out_global = TempDir::new()?;

    let make_cfg = |output_dir: &std::path::Path| {
        let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), output_dir);
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = OutlierMethodArg::Quantile;
        cfg.outlier_quantiles = vec![0.25, 0.75];
        cfg
    };

    let mut per_length_cfg = make_cfg(out_per_length.path());
    per_length_cfg.outlier_scope = OutlierScopeArg::PerLength;

    let mut global_cfg = make_cfg(out_global.path());
    global_cfg.outlier_scope = OutlierScopeArg::Global;

    // Act
    run_gc_bias(&per_length_cfg)?;
    run_gc_bias(&global_cfg)?;

    // Assert
    let package_per_length =
        GCCorrectionPackage::from_file(out_per_length.path().join("gc_bias_correction.zarr"))?;
    let package_global =
        GCCorrectionPackage::from_file(out_global.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(package_per_length.correction_matrix.dim(), (2, 2));
    assert_eq!(package_global.correction_matrix.dim(), (2, 2));

    assert_gc_command_close(
        package_per_length.correction_matrix[(0, 0)],
        5.0 / 3.0,
        "per-length quantile row0 col0",
    );
    assert_gc_command_close(
        package_per_length.correction_matrix[(0, 1)],
        5.0 / 7.0,
        "per-length quantile row0 col1",
    );
    assert_gc_command_close(
        package_per_length.correction_matrix[(1, 0)],
        1.0,
        "per-length quantile row1 col0",
    );
    assert_gc_command_close(
        package_per_length.correction_matrix[(1, 1)],
        1.0,
        "per-length quantile row1 col1",
    );

    assert_gc_command_close(
        package_global.correction_matrix[(0, 0)],
        1.25,
        "global quantile row0 col0",
    );
    assert_gc_command_close(
        package_global.correction_matrix[(0, 1)],
        5.0 / 6.0,
        "global quantile row0 col1",
    );
    assert_gc_command_close(
        package_global.correction_matrix[(1, 0)],
        1.0,
        "global quantile row1 col0",
    );
    assert_gc_command_close(
        package_global.correction_matrix[(1, 1)],
        1.0,
        "global quantile row1 col1",
    );

    Ok(())
}

#[test]
fn iqr_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
    // Arrange:
    // Reuse the same raw correction matrix as the quantile tests:
    //   [[0.2, 1.8],
    //    [1.0, 1.0]]
    //
    // We now use `--outlier-method iqr --outlier-scope per-length --outlier-k 0.25`.
    //
    // For the skewed length-10 row [0.2, 1.8]:
    // - Q1 = 0.6
    // - Q3 = 1.4
    // - IQR = 0.8
    // - Tukey bounds with k=0.25 are:
    //     lower = 0.6 - 0.25 * 0.8 = 0.4
    //     upper = 1.4 + 0.25 * 0.8 = 1.6
    // - the row is winsorized to [0.4, 1.6]
    //
    // The balanced length-11 row [1.0, 1.0] has IQR=0 and therefore stays [1.0, 1.0].
    //
    // After the command's final inversion step, the package must store:
    //   [[1 / 0.4, 1 / 1.6],
    //    [1 / 1.0, 1 / 1.0]]
    // =
    //   [[2.5, 0.625],
    //    [1.0, 1.0]]
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_min_length_bin_mass(0.0);
    cfg.set_min_length_bin_width(1);
    cfg.set_min_gc_bin_mass(1.0);
    cfg.set_num_extreme_gc_bins(0);
    cfg.set_num_short_length_bins(0);
    cfg.outlier_method = OutlierMethodArg::Iqr;
    cfg.outlier_scope = OutlierScopeArg::PerLength;
    cfg.outlier_k = 0.25;

    // Act
    run_gc_bias(&cfg)?;

    // Assert
    let package = GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
    assert_eq!(package.correction_matrix.dim(), (2, 2));
    assert_eq!(package.gc_edges, vec![0, 1, 100]);

    assert_gc_command_close(package.correction_matrix[(0, 0)], 2.5, "iqr row0 col0");
    assert_gc_command_close(package.correction_matrix[(0, 1)], 0.625, "iqr row0 col1");
    assert_gc_command_close(package.correction_matrix[(1, 0)], 1.0, "iqr row1 col0");
    assert_gc_command_close(package.correction_matrix[(1, 1)], 1.0, "iqr row1 col1");

    Ok(())
}

#[test]
fn stddev_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
    // Arrange:
    // Reuse the same raw correction matrix as the other run-level outlier tests:
    //   [[0.2, 1.8],
    //    [1.0, 1.0]]
    //
    // We now use `--outlier-method stddev --outlier-scope per-length --outlier-k 0.6`.
    //
    // For the skewed length-10 row [0.2, 1.8]:
    // - mean = (0.2 + 1.8) / 2 = 1.0
    // - sd   = sqrt((0.8^2 + 0.8^2) / 2) = 0.8
    // - bounds are:
    //     lower = 1.0 - 0.6 * 0.8 = 0.52 = 13/25
    //     upper = 1.0 + 0.6 * 0.8 = 1.48 = 37/25
    // - winsorized row becomes [0.52, 1.48]
    //
    // The balanced length-11 row [1.0, 1.0] has sd=0 and therefore stays [1.0, 1.0].
    //
    // After the command's final inversion step, the package must store:
    //   [[1 / (13/25), 1 / (37/25)],
    //    [1 / 1.0,     1 / 1.0]]
    // =
    //   [[25/13, 25/37],
    //    [1.0,   1.0]]
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_min_length_bin_mass(0.0);
    cfg.set_min_length_bin_width(1);
    cfg.set_min_gc_bin_mass(1.0);
    cfg.set_num_extreme_gc_bins(0);
    cfg.set_num_short_length_bins(0);
    cfg.outlier_method = OutlierMethodArg::Stddev;
    cfg.outlier_scope = OutlierScopeArg::PerLength;
    cfg.outlier_k = 0.6;

    // Act
    run_gc_bias(&cfg)?;

    // Assert
    let package = GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
    assert_eq!(package.correction_matrix.dim(), (2, 2));
    assert_eq!(package.gc_edges, vec![0, 1, 100]);

    assert_gc_command_close(
        package.correction_matrix[(0, 0)],
        25.0 / 13.0,
        "stddev row0 col0",
    );
    assert_gc_command_close(
        package.correction_matrix[(0, 1)],
        25.0 / 37.0,
        "stddev row0 col1",
    );
    assert_gc_command_close(package.correction_matrix[(1, 0)], 1.0, "stddev row1 col0");
    assert_gc_command_close(package.correction_matrix[(1, 1)], 1.0, "stddev row1 col1");

    Ok(())
}

#[test]
fn mad_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
    // Arrange:
    // Reuse the same raw correction matrix as the other run-level outlier tests:
    //   [[0.2, 1.8],
    //    [1.0, 1.0]]
    //
    // We now use `--outlier-method mad --outlier-scope per-length --outlier-k 0.5`.
    //
    // For the skewed length-10 row [0.2, 1.8]:
    // - median = 1.0
    // - absolute deviations are [0.8, 0.8]
    // - the implementation scales MAD by 1.4826, so:
    //     scaled_mad = 1.4826 * 0.8 = 1.18608
    // - bounds are:
    //     lower = 1.0 - 0.5 * 1.18608 = 0.40696
    //     upper = 1.0 + 0.5 * 1.18608 = 1.59304
    // - winsorized row becomes [0.40696, 1.59304]
    //
    // The balanced length-11 row [1.0, 1.0] has zero MAD and therefore stays [1.0, 1.0].
    //
    // After the command's final inversion step, the package must store:
    //   [[1 / 0.40696, 1 / 1.59304],
    //    [1 / 1.0,     1 / 1.0]]
    // =
    //   [[2.457244200903038, 0.627730626600084],
    //    [1.0,               1.0]]
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_min_length_bin_mass(0.0);
    cfg.set_min_length_bin_width(1);
    cfg.set_min_gc_bin_mass(1.0);
    cfg.set_num_extreme_gc_bins(0);
    cfg.set_num_short_length_bins(0);
    cfg.outlier_method = OutlierMethodArg::Mad;
    cfg.outlier_scope = OutlierScopeArg::PerLength;
    cfg.outlier_k = 0.5;

    // Act
    run_gc_bias(&cfg)?;

    // Assert
    let package = GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
    assert_eq!(package.correction_matrix.dim(), (2, 2));
    assert_eq!(package.gc_edges, vec![0, 1, 100]);

    assert_gc_command_close(
        package.correction_matrix[(0, 0)],
        2.457244200903038,
        "mad row0 col0",
    );
    assert_gc_command_close(
        package.correction_matrix[(0, 1)],
        0.627730626600084,
        "mad row0 col1",
    );
    assert_gc_command_close(package.correction_matrix[(1, 0)], 1.0, "mad row1 col0");
    assert_gc_command_close(package.correction_matrix[(1, 1)], 1.0, "mad row1 col1");

    Ok(())
}

#[test]
fn hard_clamp_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
    // Arrange:
    // Use a hand-built reference package that supports only GC% 0 and GC% 100 for length 10.
    // That keeps the mean-scaling denominator restricted to the two relevant cells.
    //
    // Sample BAM:
    // - one fragment in the A-only region  -> GC%=0
    // - 999 fragments in the C-only region -> GC%=100
    //
    // Before GC binning, global raw counts on the supported cells are [1, 999], so:
    // - supported mean = (1 + 999) / 2 = 500
    // - normalized supported counts = [0.002, 1.998]
    //
    // With `min_gc_bin_mass = 0.09%`, the total normalized mass is exactly 2.0, so:
    // - min bin mass = 0.0018
    // - the first GC bin `[0]` already exceeds threshold on its own
    // - the second bin is `[1..100]`
    //
    // The reference package is balanced across those same two bins, so after reference-side
    // normalization the raw correction row is:
    //   [0.002, 1.998]
    //
    // Outlier handling is disabled, so only the hard safety clamp applies:
    // - low cell 0.002 is clamped up to 0.1
    // - high cell 1.998 is unchanged
    // giving:
    //   [0.1, 1.998]
    //
    // The command then re-normalizes the row to mean 1.0 before inversion:
    // - mean = (0.1 + 1.998) / 2 = 1.049
    // - normalized row = [0.1 / 1.049, 1.998 / 1.049]
    // - final multiplicative weights after inversion are:
    //     [1.049 / 0.1, 1.049 / 1.998]
    //   = [10.49, 0.525025025025025]
    //
    // This is the important contract: the hard clamp happens *before* the final re-centering,
    // so the written package can end up slightly outside the nominal clamp range afterwards.
    let reference = twobit_from_sequences(
        "gc_bias_hard_clamp_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;

    let mut fragments = Vec::new();
    fragments.push(paired_fragment(10, 10, 5));
    for _ in 0..999 {
        fragments.push(paired_fragment(120, 10, 5));
    }
    // The hard-clamp setup stacks 999 fragments at the same C-only start, so it must opt
    // into unique qnames to keep those molecules distinct in the paired-end parser.
    let bam = bam_from_fragments_with_record_indexed_names(
        "gc_bias_hard_clamp_bam",
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
    )?;

    let ref_gc_dir = TempDir::new()?;
    write_two_bin_reference_gc_package(
        ref_gc_dir.path(),
        (10, 10),
        &["chr1"],
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_min_length_bin_mass(0.0);
    cfg.set_min_length_bin_width(1);
    // Stay slightly below the exact 0.002 boundary so the GC bin split is not sensitive to
    // floating-point equality at the threshold.
    cfg.set_min_gc_bin_mass(0.09);
    cfg.set_num_extreme_gc_bins(0);
    cfg.set_num_short_length_bins(0);
    cfg.outlier_method = OutlierMethodArg::None;

    // Act
    run_gc_bias(&cfg)?;

    // Assert
    let package = GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
    assert_eq!(package.correction_matrix.dim(), (1, 2));
    assert_eq!(package.gc_edges, vec![0, 1, 100]);
    assert!((package.correction_matrix[(0, 0)] - 10.49).abs() < 1e-12);
    assert!((package.correction_matrix[(0, 1)] - 0.525025025025025).abs() < 1e-12);

    Ok(())
}

#[test]
fn min_length_bin_width_merges_two_lengths_into_one_binned_correction_row() -> Result<()> {
    // Arrange:
    // Reuse the same two-length fixture as the run-level outlier tests. Before any length
    // binning, the normalized cfDNA rows are exactly:
    //   length 10 -> [0.2, 1.8]
    //   length 11 -> [1.0, 1.0]
    //
    // and the balanced handcrafted reference package gives:
    //   reference rows -> [1.0, 1.0] for both lengths
    //
    // With:
    // - `min_length_bin_mass = 0.0`
    // - `min_length_bin_width = 2`
    // the greedy length binning must merge the two adjacent lengths into one bin, because a
    // bin cannot close until it has width >= 2.
    //
    // The merged binned rows are then the simple mean of the two source rows:
    //   merged cfDNA row = ([0.2, 1.8] + [1.0, 1.0]) / 2 = [0.6, 1.4]
    //   merged ref row   = ([1.0, 1.0] + [1.0, 1.0]) / 2 = [1.0, 1.0]
    //
    // So the raw merged correction row is [0.6, 1.4]. Its mean is already 1.0, and with
    // outlier handling disabled no other transform changes it before inversion.
    //
    // The final written correction row must therefore be:
    //   [1 / 0.6, 1 / 1.4] = [5/3, 5/7]
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_min_length_bin_mass(0.0);
    cfg.set_min_length_bin_width(2);
    cfg.set_min_gc_bin_mass(1.0);
    cfg.set_num_extreme_gc_bins(0);
    cfg.set_num_short_length_bins(0);
    cfg.outlier_method = OutlierMethodArg::None;

    // Act
    run_gc_bias(&cfg)?;

    // Assert
    let package = GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
    assert_eq!(package.correction_matrix.dim(), (1, 2));
    assert_eq!(package.length_edges, vec![10, 11]);
    assert_eq!(package.gc_edges, vec![0, 1, 100]);
    assert_eq!(package.length_bin_frequencies.len(), 1);
    assert!((package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
    assert!((package.correction_matrix[(0, 0)] - (5.0 / 3.0)).abs() < 1e-12);
    assert!((package.correction_matrix[(0, 1)] - (5.0 / 7.0)).abs() < 1e-12);

    Ok(())
}

#[test]
fn num_short_length_bins_neutralizes_the_shortest_length_row_in_real_command() -> Result<()> {
    // Arrange:
    // Reuse the same two-length fixture as the other run-level `gc-bias` tests.
    //
    // Baseline, with no short-length masking:
    //   length 10 -> raw correction [0.2, 1.8] -> final weights [5, 5/9]
    //   length 11 -> raw correction [1.0, 1.0] -> final weights [1, 1]
    //
    // Now set `num_short_length_bins = 1`. After length binning there are exactly two
    // length rows, so the shortest row is the entire length-10 row.
    //
    // The support mask then becomes:
    //   row 0 (length 10): unsupported everywhere
    //   row 1 (length 11): supported everywhere
    //
    // The pipeline contract is:
    // 1. Unsupported entries in the normalized cfDNA and reference matrices are set to 1.0.
    // 2. The raw correction row for the masked shortest length therefore becomes [1, 1].
    // 3. Re-centering and inversion keep [1, 1] unchanged.
    //
    // So the shortest row must change from the informative baseline `[5, 5/9]` to the
    // neutral row `[1, 1]`, while the longer row stays `[1, 1]`.
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let baseline_out = TempDir::new()?;
    let masked_out = TempDir::new()?;

    let mut baseline_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        baseline_out.path(),
    );
    baseline_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    baseline_cfg.set_min_length_bin_mass(0.0);
    baseline_cfg.set_min_length_bin_width(1);
    baseline_cfg.set_min_gc_bin_mass(1.0);
    baseline_cfg.set_num_extreme_gc_bins(0);
    baseline_cfg.set_num_short_length_bins(0);
    baseline_cfg.outlier_method = OutlierMethodArg::None;

    let mut masked_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        masked_out.path(),
    );
    masked_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    masked_cfg.set_min_length_bin_mass(0.0);
    masked_cfg.set_min_length_bin_width(1);
    masked_cfg.set_min_gc_bin_mass(1.0);
    masked_cfg.set_num_extreme_gc_bins(0);
    masked_cfg.set_num_short_length_bins(1);
    masked_cfg.outlier_method = OutlierMethodArg::None;

    // Act
    run_gc_bias(&baseline_cfg)?;
    run_gc_bias(&masked_cfg)?;

    // Assert
    let baseline_package =
        GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.zarr"))?;
    let masked_package =
        GCCorrectionPackage::from_file(masked_out.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(baseline_package.correction_matrix.dim(), (2, 2));
    assert_eq!(masked_package.correction_matrix.dim(), (2, 2));
    assert_eq!(baseline_package.gc_edges, vec![0, 1, 100]);
    assert_eq!(masked_package.gc_edges, vec![0, 1, 100]);

    assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

    assert!((masked_package.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);
    assert!((masked_package.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);
    assert!((masked_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
    assert!((masked_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

    Ok(())
}

#[test]
fn num_extreme_gc_bins_neutralizes_a_two_bin_gc_axis_in_real_command() -> Result<()> {
    // Arrange:
    // Reuse the same two-length fixture again. With `min_gc_bin_mass = 1.0`, the sparse
    // counts collapse to exactly two GC bins:
    //   left  bin = GC% 0
    //   right bin = GC% 1..100
    //
    // Baseline, with `num_extreme_gc_bins = 0`:
    //   length 10 -> [5, 5/9]
    //   length 11 -> [1, 1]
    //
    // Now set `num_extreme_gc_bins = 1`. On a 2-bin GC axis, "one extreme bin from each side"
    // masks both columns:
    //   leftmost column  -> masked
    //   rightmost column -> masked
    //
    // So every correction cell is unsupported. The pipeline then sets all masked normalized
    // cfDNA/reference counts to 1.0 before division, yielding a raw correction matrix of all
    // ones. Re-centering and inversion keep it at all ones.
    //
    // This is an important boundary contract: on a 2-bin GC axis, one extreme bin per side
    // completely neutralizes the matrix rather than leaving any informative GC correction.
    let (reference, bam) = make_two_length_outlier_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let baseline_out = TempDir::new()?;
    let masked_out = TempDir::new()?;

    let mut baseline_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        baseline_out.path(),
    );
    baseline_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    baseline_cfg.set_min_length_bin_mass(0.0);
    baseline_cfg.set_min_length_bin_width(1);
    baseline_cfg.set_min_gc_bin_mass(1.0);
    baseline_cfg.set_num_extreme_gc_bins(0);
    baseline_cfg.set_num_short_length_bins(0);
    baseline_cfg.outlier_method = OutlierMethodArg::None;

    let mut masked_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        masked_out.path(),
    );
    masked_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    masked_cfg.set_min_length_bin_mass(0.0);
    masked_cfg.set_min_length_bin_width(1);
    masked_cfg.set_min_gc_bin_mass(1.0);
    masked_cfg.set_num_extreme_gc_bins(1);
    masked_cfg.set_num_short_length_bins(0);
    masked_cfg.outlier_method = OutlierMethodArg::None;

    // Act
    run_gc_bias(&baseline_cfg)?;
    run_gc_bias(&masked_cfg)?;

    // Assert
    let baseline_package =
        GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.zarr"))?;
    let masked_package =
        GCCorrectionPackage::from_file(masked_out.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(baseline_package.correction_matrix.dim(), (2, 2));
    assert_eq!(masked_package.correction_matrix.dim(), (2, 2));
    assert_eq!(baseline_package.gc_edges, vec![0, 1, 100]);
    assert_eq!(masked_package.gc_edges, vec![0, 1, 100]);

    assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

    for value in masked_package.correction_matrix.iter() {
        assert!((*value - 1.0).abs() < 1e-12);
    }

    Ok(())
}

#[test]
fn min_length_bin_mass_merges_a_sparse_tail_length_into_the_previous_bin() -> Result<()> {
    // Arrange:
    // Use a two-length fixture where the shorter length has mass 10 and the longer tail length
    // has mass 2.
    //
    // Baseline, with `min_length_bin_mass = 0` and width 1:
    //   length 10 -> [0.2, 1.8] -> final weights [5, 5/9]
    //   length 11 -> [1.0, 1.0] -> final weights [1, 1]
    //
    // Now set `min_length_bin_mass = 20%`.
    // Total binned mass is 12, so the minimum bin mass is:
    //   12 * 0.20 = 2.4
    //
    // Greedy binning over the length axis then behaves as follows:
    // - row for length 10 has mass 10, so it can close a bin by itself
    // - row for length 11 has mass 2, so it cannot form its own bin
    // - the final underweight tail row is therefore appended to the previous bin
    //
    // The important implementation detail is that greedy *length* merging happens before the
    // later per-row mean-scaling step in `gc-bias`.
    //
    // In this fixture the globally normalized GC rows are:
    //   length 10 -> [1/3, 3]
    //   length 11 -> [1/3, 1/3]
    //
    // So after greedy length merging, the single merged row is the arithmetic mean of those
    // pre-row-normalized rows:
    //   ([1/3, 3] + [1/3, 1/3]) / 2 = [1/3, 5/3]
    //
    // The reference side is still balanced after the same merge:
    //   [1, 1]
    //
    // Per-row mean scaling then keeps the raw correction row at:
    //   [1/3, 5/3]
    //
    // With no outlier handling, the final multiplicative row is therefore:
    //   [1 / (1/3), 1 / (5/3)] = [3, 3/5]
    let (reference, bam) = make_two_length_low_mass_tail_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_balanced_two_length_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let baseline_out = TempDir::new()?;
    let merged_out = TempDir::new()?;

    let mut baseline_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        baseline_out.path(),
    );
    baseline_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    baseline_cfg.set_min_length_bin_mass(0.0);
    baseline_cfg.set_min_length_bin_width(1);
    baseline_cfg.set_min_gc_bin_mass(1.0);
    baseline_cfg.set_num_extreme_gc_bins(0);
    baseline_cfg.set_num_short_length_bins(0);
    baseline_cfg.outlier_method = OutlierMethodArg::None;

    let mut merged_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        merged_out.path(),
    );
    merged_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    merged_cfg.set_min_length_bin_mass(20.0);
    merged_cfg.set_min_length_bin_width(1);
    merged_cfg.set_min_gc_bin_mass(1.0);
    merged_cfg.set_num_extreme_gc_bins(0);
    merged_cfg.set_num_short_length_bins(0);
    merged_cfg.outlier_method = OutlierMethodArg::None;

    // Act
    run_gc_bias(&baseline_cfg)?;
    run_gc_bias(&merged_cfg)?;

    // Assert
    let baseline_package =
        GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.zarr"))?;
    let merged_package =
        GCCorrectionPackage::from_file(merged_out.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(baseline_package.correction_matrix.dim(), (2, 2));
    assert_eq!(merged_package.correction_matrix.dim(), (1, 2));
    assert_eq!(merged_package.length_edges, vec![10, 11]);
    assert_eq!(merged_package.gc_edges, vec![0, 1, 100]);
    assert_eq!(merged_package.length_bin_frequencies.len(), 1);
    assert!((merged_package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);

    assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

    assert!((merged_package.correction_matrix[(0, 0)] - 3.0).abs() < 1e-12);
    assert!((merged_package.correction_matrix[(0, 1)] - (3.0 / 5.0)).abs() < 1e-12);

    Ok(())
}

#[test]
fn min_gc_bin_mass_greedily_merges_sparse_gc_tail_bins_in_real_command() -> Result<()> {
    // Arrange:
    // Use a one-length fixture with exact GC-class masses:
    //   GC%=0   -> 1
    //   GC%=50  -> 5
    //   GC%=100 -> 9
    //
    // The handcrafted reference package marks those same three GC points as equally likely.
    //
    // Baseline, with `min_gc_bin_mass = 1%`:
    // - total sample mass = 15
    // - min bin mass = 0.15
    // - greedy GC binning therefore closes bins at:
    //     [0], [1..50], [51..100]
    // - sample binned row is [1, 5, 9]
    // - reference binned row is [1, 1, 1]
    // - normalized correction row is already [0.2, 1.0, 1.8]
    // - final weights are therefore [5, 1, 5/9]
    //
    // Now set `min_gc_bin_mass = 25%`:
    // - min bin mass = 15 * 0.25 = 3.75
    // - the first sparse GC point (mass 1 at GC%=0) cannot close a bin by itself
    // - it gets merged with the next non-zero point at GC%=50
    // - the resulting GC bins are:
    //     [0..50], [51..100]
    //
    // So the binned rows become:
    //   sample    -> [1 + 5, 9] = [6, 9]
    //   reference -> [1 + 1, 1] = [2, 1]
    //
    // Per-row normalization gives:
    //   sample    -> [0.8, 1.2]
    //   reference -> [4/3, 2/3]
    //
    // Raw correction is therefore:
    //   [0.8 / (4/3), 1.2 / (2/3)] = [0.6, 1.8]
    //
    // Re-centering that row to mean 1.0 gives [0.5, 1.5], and inversion yields:
    //   [2, 2/3]
    let (reference, bam) = make_three_gc_bin_fixture()?;
    let ref_gc_dir = TempDir::new()?;
    write_three_bin_reference_gc_package(
        ref_gc_dir.path(),
        twobit_contig_footprint(&reference.path)?,
    )?;

    let baseline_out = TempDir::new()?;
    let merged_out = TempDir::new()?;

    let mut baseline_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        baseline_out.path(),
    );
    baseline_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    baseline_cfg.set_min_length_bin_mass(0.0);
    baseline_cfg.set_min_length_bin_width(1);
    baseline_cfg.set_min_gc_bin_mass(1.0);
    baseline_cfg.set_num_extreme_gc_bins(0);
    baseline_cfg.set_num_short_length_bins(0);
    baseline_cfg.outlier_method = OutlierMethodArg::None;

    let mut merged_cfg = make_gc_bias_cfg(
        &bam.bam,
        &reference.path,
        ref_gc_dir.path(),
        merged_out.path(),
    );
    merged_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    merged_cfg.set_min_length_bin_mass(0.0);
    merged_cfg.set_min_length_bin_width(1);
    merged_cfg.set_min_gc_bin_mass(25.0);
    merged_cfg.set_num_extreme_gc_bins(0);
    merged_cfg.set_num_short_length_bins(0);
    merged_cfg.outlier_method = OutlierMethodArg::None;

    // Act
    run_gc_bias(&baseline_cfg)?;
    run_gc_bias(&merged_cfg)?;

    // Assert
    let baseline_package =
        GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.zarr"))?;
    let merged_package =
        GCCorrectionPackage::from_file(merged_out.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(baseline_package.correction_matrix.dim(), (1, 3));
    assert_eq!(baseline_package.gc_edges, vec![0, 1, 51, 100]);
    assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);
    assert!((baseline_package.correction_matrix[(0, 2)] - (5.0 / 9.0)).abs() < 1e-12);

    assert_eq!(merged_package.correction_matrix.dim(), (1, 2));
    assert_eq!(merged_package.gc_edges, vec![0, 51, 100]);
    assert!((merged_package.correction_matrix[(0, 0)] - 2.0).abs() < 1e-12);
    assert!((merged_package.correction_matrix[(0, 1)] - (2.0 / 3.0)).abs() < 1e-12);

    Ok(())
}

#[test]
fn gc_bias_transfers_reference_gc_package_footprint_to_correction_package() -> Result<()> {
    // Human verification status: verified
    // Manual expectations:
    // - `ref-gc-bias` writes the 2bit contig footprint into the reference package.
    // - `gc-bias` validates that footprint against its current `--ref-2bit`.
    // - The final GC correction package should carry forward that validated footprint exactly.
    let reference = simple_reference_twobit()?;
    let bam = single_contig_inward_pair_bam()?;
    let ref_gc_dir = TempDir::new()?;
    write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;
    let out_dir = TempDir::new()?;
    let mut cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
    cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
    cfg.set_num_extreme_gc_bins(0);
    cfg.set_num_short_length_bins(0);
    cfg.outlier_method = OutlierMethodArg::None;

    run_gc_bias(&cfg)?;

    let expected_footprint = twobit_contig_footprint(&reference.path)?;
    let reference_data = load_reference_gc_data(&ref_gc_dir.path().join("ref_gc_package.zarr"))?;
    let correction_package =
        GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;

    assert_eq!(
        reference_data.metadata.reference_contig_footprint,
        expected_footprint
    );
    assert_eq!(
        correction_package.reference_contig_footprint,
        expected_footprint
    );

    Ok(())
}

mod test_gc_tag_values {
    #[cfg(feature = "cmd_gc_bias")]
    use crate::commands::gc_bias::{
        correct::GCCorrector, counting::build_gc_prefixes, package::GCCorrectionPackage,
    };
    use rust_htslib::bam::record::{Aux, Record};

    use crate::shared::base::ZEROISH_F32_TOLERANCE;
    #[cfg(feature = "cmd_gc_bias")]
    use crate::shared::constants::GC_CORRECTION_SCHEMA_VERSION;
    use crate::shared::gc_tag::{
        ClassifiedGCTagWeight, GCTagValue, MIN_REASONABLE_GC_WEIGHT, combine_gc_tag_values,
        read_gc_tag_from_record,
    };
    #[cfg(feature = "cmd_gc_bias")]
    use crate::shared::interval::Interval;
    #[cfg(feature = "cmd_gc_bias")]
    use ndarray::array;

    #[test]
    fn gc_tag_values_follow_supported_range_and_zero_snap_rules() {
        // Arrange: start with a sane weight
        let mut rec_ok = Record::new();
        rec_ok.push_aux(b"GC", Aux::Float(2.5)).expect("set GC tag");
        let ok = read_gc_tag_from_record(&rec_ok, b"GC");

        // Assert: valid weight passes through
        assert_eq!(ok.weight, Some(2.5));
        assert!(!ok.was_missing);
        assert!(!ok.had_invalid);
        assert!(!ok.was_out_of_range);

        // Arrange: record carrying a wildly high weight that should be treated as invalid
        let mut rec_high = Record::new();
        rec_high
            .push_aux(b"GC", Aux::Float(1.1e3))
            .expect("set GC tag");
        let high = read_gc_tag_from_record(&rec_high, b"GC");

        // Assert: extreme values are rejected to avoid runaway coverage
        assert!(high.weight.is_none());
        assert!(high.had_invalid);
        assert!(high.was_out_of_range);

        // Arrange: meaningfully negative values are invalid, not zero-snapped.
        let mut rec_neg = Record::new();
        rec_neg
            .push_aux(b"GC", Aux::Float(-3.0))
            .expect("set GC tag");
        let neg = read_gc_tag_from_record(&rec_neg, b"GC");
        assert!(neg.weight.is_none());
        assert!(neg.had_invalid);
        assert!(neg.was_out_of_range);

        // Arrange: NaN should be invalid but not counted as out-of-range
        let mut rec_nan = Record::new();
        rec_nan
            .push_aux(b"GC", Aux::Float(f32::NAN))
            .expect("set GC tag");
        let nan = read_gc_tag_from_record(&rec_nan, b"GC");

        assert!(nan.weight.is_none());
        assert!(nan.had_invalid);
        assert!(!nan.was_out_of_range);

        // Arrange: tiny positive values near zero are snapped to zero.
        let mut rec_tiny = Record::new();
        rec_tiny
            .push_aux(b"GC", Aux::Float(ZEROISH_F32_TOLERANCE))
            .expect("set GC tag");
        let tiny = read_gc_tag_from_record(&rec_tiny, b"GC");
        assert_eq!(tiny.weight, Some(0.0));
        assert!(!tiny.had_invalid);
        assert!(!tiny.was_out_of_range);

        // Arrange: the zero-snap window is symmetric around zero.
        let mut rec_tiny_negative = Record::new();
        rec_tiny_negative
            .push_aux(b"GC", Aux::Float(-ZEROISH_F32_TOLERANCE))
            .expect("set GC tag");
        let tiny_negative = read_gc_tag_from_record(&rec_tiny_negative, b"GC");
        assert_eq!(tiny_negative.weight, Some(0.0));
        assert!(!tiny_negative.had_invalid);
        assert!(!tiny_negative.was_out_of_range);

        // Arrange: positive values below the minimum supported GC weight are invalid.
        let mut rec_low = Record::new();
        rec_low
            .push_aux(b"GC", Aux::Float(MIN_REASONABLE_GC_WEIGHT / 10.0))
            .expect("set GC tag");
        let low = read_gc_tag_from_record(&rec_low, b"GC");
        assert!(low.weight.is_none());
        assert!(low.had_invalid);
        assert!(low.was_out_of_range);
    }

    #[test]
    fn gc_tag_values_just_below_minimum_supported_weight_are_invalid() {
        // Arrange: choose the nearest representable f32 below the supported lower bound.
        // This is a stronger boundary than "/10" because it proves the exact cutoff behavior.
        let just_below_min = f32::from_bits(MIN_REASONABLE_GC_WEIGHT.to_bits() - 1);
        let mut rec = Record::new();
        rec.push_aux(b"GC", Aux::Float(just_below_min))
            .expect("set GC tag");

        // Act
        let observed = read_gc_tag_from_record(&rec, b"GC");

        // Assert: values below 1e-3 remain invalid even when they are only one f32 step lower.
        assert!(observed.weight.is_none());
        assert!(observed.had_invalid);
        assert!(observed.was_out_of_range);
    }

    #[test]
    fn missing_gc_tag_is_reported_separately() {
        let rec = Record::new();
        let missing = read_gc_tag_from_record(&rec, b"GC");
        assert!(missing.weight.is_none());
        assert!(missing.was_missing);
        assert!(!missing.had_invalid);
        assert!(!missing.was_out_of_range);
    }

    #[test]
    fn combining_valid_weights_averages_before_final_range_check() {
        let mut rec_a = Record::new();
        rec_a.push_aux(b"GC", Aux::Float(2.0)).expect("set GC tag");
        let mut rec_b = Record::new();
        rec_b.push_aux(b"GC", Aux::Float(4.0)).expect("set GC tag");

        let a = read_gc_tag_from_record(&rec_a, b"GC");
        let b = read_gc_tag_from_record(&rec_b, b"GC");
        let combined = combine_gc_tag_values(&a, &b);

        assert_eq!(combined.weight, Some(3.0));
        assert!(!combined.had_invalid);
        assert!(!combined.was_out_of_range);
    }

    #[test]
    fn combining_paired_tags_reuses_single_usable_mate_and_keeps_zero_precedence() {
        let mut rec_zero = Record::new();
        rec_zero
            .push_aux(b"GC", Aux::Float(0.0))
            .expect("set GC tag");
        let mut rec_valid = Record::new();
        rec_valid
            .push_aux(b"GC", Aux::Float(4.0))
            .expect("set GC tag");

        let zero = read_gc_tag_from_record(&rec_zero, b"GC");
        let valid = read_gc_tag_from_record(&rec_valid, b"GC");
        let zero_combined = combine_gc_tag_values(&zero, &valid);
        assert_eq!(zero_combined.weight, Some(0.0));
        assert!(!zero_combined.had_invalid);

        let missing = read_gc_tag_from_record(&Record::new(), b"GC");
        let missing_combined = combine_gc_tag_values(&valid, &missing);
        assert_eq!(missing_combined.weight, Some(4.0));
        assert!(!missing_combined.was_missing);
        assert!(!missing_combined.had_invalid);
        assert!(!missing_combined.was_out_of_range);

        let mut rec_low = Record::new();
        rec_low
            .push_aux(b"GC", Aux::Float(MIN_REASONABLE_GC_WEIGHT / 10.0))
            .expect("set GC tag");
        let low = read_gc_tag_from_record(&rec_low, b"GC");
        let invalid_combined = combine_gc_tag_values(&valid, &low);
        assert!(invalid_combined.weight.is_none());
        assert!(invalid_combined.had_invalid);
        assert!(invalid_combined.was_out_of_range);
    }

    #[test]
    fn gc_tag_classify_exposes_one_explicit_state() {
        assert_eq!(
            GCTagValue {
                weight: Some(2.5),
                was_missing: false,
                had_invalid: false,
                was_out_of_range: false,
            }
            .classify()
            .expect("valid classification"),
            ClassifiedGCTagWeight::Usable(2.5)
        );
        assert_eq!(
            GCTagValue::missing()
                .classify()
                .expect("missing classification"),
            ClassifiedGCTagWeight::Missing
        );
        assert_eq!(
            GCTagValue {
                weight: None,
                was_missing: false,
                had_invalid: true,
                was_out_of_range: true,
            }
            .classify()
            .expect("invalid classification"),
            ClassifiedGCTagWeight::Invalid { out_of_range: true }
        );
    }

    #[test]
    fn gc_tag_classify_rejects_inconsistent_internal_state() {
        let err = GCTagValue {
            weight: None,
            was_missing: false,
            had_invalid: false,
            was_out_of_range: false,
        }
        .classify()
        .expect_err("inconsistent state should error");

        assert!(
            err.to_string().contains("inconsistent GC tag state"),
            "unexpected error: {err}"
        );
    }

    #[cfg(feature = "cmd_gc_bias")]
    #[test]
    fn gc_file_weights_follow_the_same_sanity_rules() {
        let prefixes = build_gc_prefixes(b"AAAAAAAAAA");
        let interval = Interval::new(0_u64, 10_u64).expect("valid interval");
        let scenarios = [
            ("negative_below_snap_window_is_unusable", -3.0_f64, None),
            (
                "tiny_negative_becomes_zero",
                -(ZEROISH_F32_TOLERANCE as f64),
                Some(0.0_f64),
            ),
            (
                "tiny_positive_becomes_zero",
                ZEROISH_F32_TOLERANCE as f64,
                Some(0.0_f64),
            ),
            (
                "too_small_positive_is_unusable",
                (MIN_REASONABLE_GC_WEIGHT / 10.0) as f64,
                None,
            ),
            ("too_large_positive_is_unusable", 1.1e3_f64, None),
        ];

        for (name, weight, expected) in scenarios {
            let package = GCCorrectionPackage {
                version: GC_CORRECTION_SCHEMA_VERSION,
                end_offset: 0,
                length_edges: vec![10, 11],
                gc_edges: vec![0, 101],
                length_bin_frequencies: array![1.0_f64],
                reference_contig_footprint: Vec::new(),
                correction_matrix: array![[weight]],
            };
            let corrector = GCCorrector::from_package(&package).expect("build corrector");

            let observed = corrector
                .correct_fragment(interval, &prefixes)
                .expect("correct fragment");
            assert_eq!(observed, expected, "unexpected sanitized weight for {name}");
        }
    }

    #[cfg(feature = "cmd_gc_bias")]
    #[test]
    fn gc_file_weights_just_below_minimum_supported_weight_are_invalid() {
        // Arrange: use the nearest representable f64 below the exact f64 threshold used by the
        // sanitizer so this checks the boundary itself, not just a clearly too-small value.
        let just_below_min = f64::from_bits((MIN_REASONABLE_GC_WEIGHT as f64).to_bits() - 1);
        let prefixes = build_gc_prefixes(b"AAAAAAAAAA");
        let interval = Interval::new(0_u64, 10_u64).expect("valid interval");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![10, 11],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            reference_contig_footprint: Vec::new(),
            correction_matrix: array![[just_below_min]],
        };
        let corrector = GCCorrector::from_package(&package).expect("build corrector");

        // Act
        let observed = corrector
            .correct_fragment(interval, &prefixes)
            .expect("correct fragment");

        // Assert: the GC-file path rejects values that fall just below the accepted range.
        assert_eq!(observed, None);
    }
}

mod test_fragment_iterator_gc_tags {
    use crate::{
        shared::{
            fragment::{
                minimal_fragment::Fragment, segment_fragment::FragmentWithSegments,
            },
            fragment_iterators::{
                fragments_from_bam, fragments_with_segments_from_bam,
            },
            gc_tag::GCTagValue,
        },
    };
    #[cfg(feature = "cmd_ends")]
    use crate::{
        commands::ends::config_structs::{ClipStrategy, KmerSource},
        shared::{
            fragment::ends_fragment::FragmentWithEnds,
            fragment_iterators::fragments_with_ends_from_bam,
            indel_mode::IndelMotifFilterPolicy,
        },
    };
    #[cfg(feature = "cmd_fragment_kmers")]
    use crate::shared::{
        fragment::segment_kmer_fragment::FragmentWithKmerSegments,
        fragment_iterators::fragments_with_kmer_segments_from_bam,
        indel_mode::IndelMode,
    };
    use anyhow::Result;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};

    fn assert_valid_gc_tag(observed: GCTagValue, expected_weight: f32) {
        assert_eq!(observed.weight, Some(expected_weight));
        assert!(!observed.was_missing);
        assert!(!observed.had_invalid);
        assert!(!observed.was_out_of_range);
    }

    fn make_record(
        qname: &[u8],
        tid: i32,
        pos: i64,
        is_reverse: bool,
        seq_len: usize,
        gc_weight: f32,
    ) -> Record {
        let mut record = Record::new();
        record.set_tid(tid);
        record.set_pos(pos);
        record.set_flags(if is_reverse { 0x11 } else { 0x1 });
        record.set_mapq(60);

        let cigar = CigarString(vec![Cigar::Match(seq_len as u32)]);
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        record.set(qname, Some(&cigar), &seq, &qual);
        record
            .push_aux(b"GC", Aux::Float(gc_weight))
            .expect("set GC tag");

        record
    }

    fn first_fragment(iter: impl Iterator<Item = Result<Fragment>>) -> Fragment {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    fn first_segment_fragment(
        iter: impl Iterator<Item = Result<FragmentWithSegments>>,
    ) -> FragmentWithSegments {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    #[cfg(feature = "cmd_fragment_kmers")]
    fn first_kmer_segment_fragment(
        iter: impl Iterator<Item = Result<FragmentWithKmerSegments>>,
    ) -> FragmentWithKmerSegments {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    #[cfg(feature = "cmd_ends")]
    fn first_end_fragment(
        iter: impl Iterator<Item = Result<FragmentWithEnds>>,
    ) -> FragmentWithEnds {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    #[test]
    fn basic_fragment_iterator_paired_uses_configured_gc_tag() {
        // Arrange: two mates with GC weights 2 and 4 should average to 3 on the fragment.
        let qname = b"pair_basic";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act: build fragments through the same iterator used by basic-fragment commands.
        let fragment = first_fragment(fragments_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            Some(b"GC"),
            |_fragment: &Fragment| true,
            false,
        ));

        // Assert: the configured GC tag is preserved and combined at fragment level.
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[test]
    fn basic_fragment_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange: a single read-as-fragment should keep its own GC-tag value.
        let record = make_record(b"single_basic", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_fragment(fragments_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            Some(b"GC"),
            |_fragment: &Fragment| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }

    #[test]
    fn segment_fragment_iterator_paired_uses_configured_gc_tag() {
        // Arrange
        let qname = b"pair_segments";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act
        let fragment = first_segment_fragment(fragments_with_segments_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            1,
            true,
            Some(b"GC"),
            |_fragment: &FragmentWithSegments| true,
            false,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[test]
    fn segment_fragment_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange
        let record = make_record(b"single_segments", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_segment_fragment(fragments_with_segments_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            1,
            true,
            Some(b"GC"),
            |_fragment: &FragmentWithSegments| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }

    #[cfg(feature = "cmd_fragment_kmers")]
    #[test]
    fn kmer_segment_iterator_paired_uses_configured_gc_tag() {
        // Arrange
        let qname = b"pair_kmers";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act
        let fragment = first_kmer_segment_fragment(fragments_with_kmer_segments_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            IndelMode::Ignore,
            true,
            0,
            Some(b"GC"),
            |_fragment: &FragmentWithKmerSegments| true,
            false,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[cfg(feature = "cmd_fragment_kmers")]
    #[test]
    fn kmer_segment_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange
        let record = make_record(b"single_kmers", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_kmer_segment_fragment(fragments_with_kmer_segments_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            IndelMode::Ignore,
            true,
            0,
            Some(b"GC"),
            |_fragment: &FragmentWithKmerSegments| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }

    #[cfg(feature = "cmd_ends")]
    #[test]
    fn ends_iterator_paired_uses_configured_gc_tag() {
        // Arrange
        let qname = b"pair_ends";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act
        let fragment = first_end_fragment(fragments_with_ends_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            u32::MAX,
            &[],
            Some(b"GC"),
            |_fragment: &FragmentWithEnds| true,
            false,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[cfg(feature = "cmd_ends")]
    #[test]
    fn ends_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange
        let record = make_record(b"single_ends", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_end_fragment(fragments_with_ends_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            u32::MAX,
            &[],
            Some(b"GC"),
            |_fragment: &FragmentWithEnds| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }
}
