#![cfg(feature = "cmd_midpoints")]

// KEEP-IN-TESTS: all active tests in this file cover midpoint/profile command output, errors, or artifacts.

mod fixtures;

use anyhow::Result;
use cfdnalab::RunOptions;
#[cfg(feature = "cli")]
use anyhow::{Context, bail};
#[cfg(feature = "cmd_bam_to_bam")]
use cfdnalab::run_like_cli::bam_to_bam::{
    BamToBamConfig, run_bam_to_bam as run_bam_to_bam_command,
};
use cfdnalab::run_like_cli::common::{ApplyGCArgs, ChromosomeArgs, IOCArgs, ScaleGenomeArgs};
#[cfg(feature = "cmd_coverage_weights")]
use cfdnalab::run_like_cli::coverage_weights::{
    CoverageWeightsConfig, run_coverage_weights as run_coverage_weights_command,
};
use cfdnalab::gc_bias::GCCorrectionPackage;
use cfdnalab::run_like_cli::midpoints::{
    MidpointSmoothing, MidpointsConfig, run_midpoints as run_midpoints_command,
};
use cfdnalab::{
    blacklist::BlacklistStrategy,
    constants::GC_CORRECTION_SCHEMA_VERSION,
    reference::{ContigFootprintEntry, twobit_contig_footprint},
};
use fixtures::{
    FragmentSpec, ReadSpec, bam_from_specs, bam_from_specs_strict_identity,
    build_real_neutral_gc_package, build_real_neutral_gc_package_for_range,
    build_real_non_neutral_gc_package, late_origin_gc_reference_sequence,
    read_midpoint_zarr_counts, read_midpoint_zarr_i32_1d, simple_reference_twobit,
    twobit_from_sequences, write_bed, write_two_bin_gc_package,
};
use ndarray::Array3;
use ndarray::array;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read, Reader};
use serde_json::Value;
use std::collections::HashMap;
#[cfg(feature = "cli")]
use std::ffi::OsStr;
#[cfg(feature = "cli")]
use std::path::Path;
use std::path::PathBuf;
#[cfg(feature = "cli")]
use std::process::Command;
use tempfile::TempDir;

fn run(cfg: &MidpointsConfig) -> Result<()> {
    run_midpoints_command(cfg, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(feature = "cmd_bam_to_bam")]
fn run_bam_to_bam(cfg: &BamToBamConfig) -> Result<()> {
    run_bam_to_bam_command(cfg, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(feature = "cmd_coverage_weights")]
fn run_coverage_weights(cfg: &CoverageWeightsConfig) -> Result<()> {
    run_coverage_weights_command(cfg, RunOptions::new_quiet()).map(|_| ())
}

const MIDPOINT_F32_TOL: f32 = 1e-5;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
}

fn base_midpoints_config_for_length_bins() -> MidpointsConfig {
    MidpointsConfig::new(
        IOCArgs {
            bam: PathBuf::from("dummy.bam"),
            output_dir: PathBuf::from("out"),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        PathBuf::from("intervals.bed"),
    )
}

fn paired_fragment_on_tid(
    tid: usize,
    start: i64,
    fragment_len: i64,
    read_len: i64,
) -> FragmentSpec {
    paired_fragment_on_tid_with_mapq(tid, start, fragment_len, read_len, 60)
}

fn paired_fragment_on_tid_with_mapq(
    tid: usize,
    start: i64,
    fragment_len: i64,
    read_len: i64,
    mapq: u8,
) -> FragmentSpec {
    const FLAG_FIRST_MATE: u16 = 0x40;
    const FLAG_SECOND_MATE: u16 = 0x80;
    const FLAG_PROPER_PAIR: u16 = 0x2;
    const FLAG_MATE_REVERSE: u16 = 0x20;

    let reverse_start = start + fragment_len - read_len;
    let insert_size = fragment_len;
    FragmentSpec {
        forward: ReadSpec {
            tid,
            pos: start,
            cigar: vec![('M', read_len as u32)],
            seq: vec![b'A'; read_len as usize],
            qual: 40,
            is_reverse: false,
            mapq,
            flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
            mate_tid: Some(tid),
            mate_pos: Some(reverse_start),
            insert_size,
        },
        reverse: ReadSpec {
            tid,
            pos: reverse_start,
            cigar: vec![('M', read_len as u32)],
            seq: vec![b'T'; read_len as usize],
            qual: 40,
            is_reverse: true,
            mapq,
            flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
            mate_tid: Some(tid),
            mate_pos: Some(start),
            insert_size: -insert_size,
        },
    }
}

fn single_read_fragment_bam(
    name: &str,
    fragment_start: i64,
    fragment_len: u32,
) -> Result<fixtures::BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: fragment_start,
            cigar: vec![('M', fragment_len)],
            seq: vec![b'A'; fragment_len as usize],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        name,
    )
}

fn assert_midpoint_profile_row_matches(row: &[f32], expected_weight_at_pos5: f32, context: &str) {
    // `midpoints` accumulates profile mass in `Vec<f32>` and writes the final NPY as `f32`.
    // So the scientifically correct contract here is the stored `f32` value, not ideal `f64`
    // arithmetic carried all the way through the test expectation.
    for (position, value) in row.iter().enumerate() {
        let expected = if position == 5 {
            expected_weight_at_pos5
        } else {
            0.0
        };
        assert!(
            (value - expected).abs() <= MIDPOINT_F32_TOL,
            "{context}: unexpected midpoint weight at position {position}: expected {expected}, got {value}"
        );
    }
}

#[derive(Debug)]
struct TaggedBamFixture {
    _tempdir: TempDir,
    bam: PathBuf,
}

fn build_bai_for_test_bam(bam_path: &std::path::Path) -> Result<()> {
    let bai_path = bam_path.with_extension("bam.bai");
    bam::index::build(bam_path, None, bam::index::Type::Bai, 1)?;
    let target = bam_path.with_extension("bai");
    if bai_path.exists() {
        std::fs::rename(&bai_path, &target)?;
    }
    Ok(())
}

fn bam_with_gc_tags(
    base_bam: &std::path::Path,
    name: &str,
    tags: &[Option<f32>],
) -> Result<TaggedBamFixture> {
    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));

    let mut reader = Reader::from_path(base_bam)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;

    for (record_index, record_result) in reader.records().enumerate() {
        let mut record = record_result?;
        if let Some(Some(tag_value)) = tags.get(record_index) {
            record.push_aux(b"GC", Aux::Float(*tag_value))?;
        }
        writer.write(&record)?;
    }

    drop(writer);
    build_bai_for_test_bam(&bam_path)?;

    Ok(TaggedBamFixture {
        _tempdir: tempdir,
        bam: bam_path,
    })
}

fn read_group_index_map(path: &std::path::Path) -> Result<HashMap<String, usize>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = HashMap::new();
    for line in text.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let idx = fields.next().unwrap().parse::<usize>()?;
        let name = fields.next().unwrap().to_string();
        out.insert(name, idx);
    }
    Ok(out)
}

fn read_group_index_eligible_intervals(path: &std::path::Path) -> Result<HashMap<String, usize>> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header = lines.next().unwrap_or("");
    assert_eq!(
        header, "group_idx\tgroup_name\teligible_intervals",
        "midpoint group index must expose eligible profile intervals per group"
    );

    let mut out = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            fields.len(),
            3,
            "midpoint group index row should have group_idx, group_name, and eligible_intervals"
        );
        let group_name = fields[1].to_string();
        let eligible_intervals = fields[2].parse::<usize>()?;
        out.insert(group_name, eligible_intervals);
    }
    Ok(out)
}

#[derive(Debug)]
struct MidpointsAxisContractFixture {
    _interval_dir: TempDir,
    bam: fixtures::BamFixture,
    intervals_path: PathBuf,
}

fn midpoint_axis_contract_fixture() -> Result<MidpointsAxisContractFixture> {
    let high_tile_offset = 1_000_000_i64;
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 1_100_000),
            ("chr2".to_string(), 1_100_000),
        ],
        vec![
            // Positive: chr1 groupA, length bin [20,60), midpoint 50, position 10.
            paired_fragment_on_tid(0, 40, 21, 10),
            // Positive: chr1 groupA, length bin [60,120), midpoint 1,000,205, position 25.
            paired_fragment_on_tid(0, high_tile_offset + 168, 75, 20),
            // Positive: chr2 groupB, length bin [20,60), midpoint 45, position 25.
            paired_fragment_on_tid(1, 35, 21, 10),
            // Positive: chr2 groupB, length bin [60,120), midpoint 1,000,080, position 20.
            paired_fragment_on_tid(1, high_tile_offset + 43, 75, 20),
            // Positive: chr1 groupC, length bin [20,60), midpoint 310, position 10.
            paired_fragment_on_tid(0, 300, 21, 10),
            // Positive: chr2 groupC, length bin [20,60), midpoint 310, position 10.
            paired_fragment_on_tid(1, 300, 21, 10),
            // Negative: midpoint sits inside a groupA window, but length 19 is below [20,60).
            paired_fragment_on_tid(0, 56, 19, 10),
            // Negative: midpoint sits inside a groupB window, but length 120 is the final
            // exclusive edge of [60,120).
            paired_fragment_on_tid(1, high_tile_offset + 20, 120, 20),
            // Negative: midpoint sits inside a groupA window, but MAPQ is below the configured
            // threshold.
            paired_fragment_on_tid_with_mapq(0, high_tile_offset + 180, 21, 10, 20),
            // Negative: accepted length and MAPQ, but the midpoint does not overlap any BED row.
            paired_fragment_on_tid(1, 500_000, 21, 10),
        ],
        Vec::new(),
        "midpoints_axis_contract",
    )?;

    let interval_dir = TempDir::new()?;
    let intervals_path = interval_dir.path().join("windows.bed");
    write_bed(
        &intervals_path,
        &[
            ("chr1", 40, 80, "groupA"),
            ("chr1", 300, 340, "groupC"),
            ("chr1", 1_000_180, 1_000_220, "groupA"),
            ("chr2", 20, 60, "groupB"),
            ("chr2", 300, 340, "groupC"),
            ("chr2", 1_000_060, 1_000_100, "groupB"),
        ],
    )?;

    Ok(MidpointsAxisContractFixture {
        _interval_dir: interval_dir,
        bam,
        intervals_path,
    })
}

fn expected_midpoint_axis_contract_counts() -> Array3<f32> {
    let mut expected = Array3::<f32>::zeros((3, 2, 40));
    expected[[0, 0, 10]] = 1.0;
    expected[[0, 1, 25]] = 1.0;
    expected[[1, 0, 10]] = 2.0;
    expected[[2, 0, 25]] = 1.0;
    expected[[2, 1, 20]] = 1.0;
    expected
}

fn midpoint_axis_contract_config(
    fixture: &MidpointsAxisContractFixture,
    output_dir: &std::path::Path,
    n_threads: usize,
    output_prefix: &str,
) -> MidpointsConfig {
    let mut config = MidpointsConfig::new(
        IOCArgs {
            bam: fixture.bam.bam.clone(),
            output_dir: output_dir.to_path_buf(),
            n_threads,
        },
        base_chromosomes(&["chr1", "chr2"]),
        fixture.intervals_path.clone(),
    );
    config.set_output_prefix(output_prefix);
    config.set_length_bins(vec![20, 60, 120]);
    config.set_smoothing(MidpointSmoothing::None);
    config.set_tile_size(1_000_000);
    config.set_min_mapq(30);
    config.set_require_proper_pair(false);
    config.set_scale_genome(ScaleGenomeArgs::default());
    config
}

#[cfg(feature = "cli")]
fn binary_name(base_name: &str) -> String {
    if cfg!(windows) {
        format!("{base_name}.exe")
    } else {
        base_name.to_string()
    }
}

#[cfg(feature = "cli")]
fn cfdna_bin_path() -> Result<PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_cfdna") {
        return Ok(PathBuf::from(path));
    }

    let current_exe = std::env::current_exe().context("failed to read current test binary path")?;
    let deps_dir = current_exe
        .parent()
        .context("failed to derive deps directory from current test binary path")?;
    let target_dir = deps_dir
        .parent()
        .context("failed to derive target directory from deps path")?;

    let direct_path = target_dir.join(binary_name("cfdna"));
    if direct_path.is_file() {
        return Ok(direct_path);
    }

    let mut hashed_candidates = Vec::new();
    for entry in std::fs::read_dir(deps_dir).with_context(|| {
        format!(
            "failed to list candidate binaries in {}",
            deps_dir.display()
        )
    })? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let file_name = match path.file_name().and_then(OsStr::to_str) {
            Some(name) => name,
            None => continue,
        };
        let extension = path.extension().and_then(OsStr::to_str);
        let looks_like_hashed_binary = file_name.starts_with("cfdna-");
        let is_makefile_dep = extension == Some("d");
        if looks_like_hashed_binary && !is_makefile_dep {
            hashed_candidates.push(path);
        }
    }
    hashed_candidates.sort_by_key(|path| {
        std::fs::metadata(path)
            .and_then(|meta| meta.modified())
            .ok()
            .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
    });
    if let Some(path) = hashed_candidates.into_iter().last() {
        return Ok(path);
    }

    bail!(
        "Could not locate cfdna binary. Tried CARGO_BIN_EXE_cfdna, {}, and hashed binaries under {}",
        direct_path.display(),
        deps_dir.display()
    );
}

#[cfg(feature = "cli")]
fn command_output(command_name: &str, args: &[&str]) -> Result<std::process::Output> {
    Command::new(cfdna_bin_path()?)
        .arg(command_name)
        .args(args)
        .output()
        .with_context(|| format!("failed running cfdna {command_name} {}", args.join(" ")))
}

#[cfg(feature = "cli")]
fn assert_success_with_logs(output: &std::process::Output, command_desc: &str) {
    let stdout_text = String::from_utf8_lossy(&output.stdout);
    let stderr_text = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "Expected {command_desc} to succeed.\nstdout:\n{stdout_text}\nstderr:\n{stderr_text}"
    );
}

#[cfg(feature = "cli")]
fn path_text(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

fn write_minimal_gc_package_excluding_length_61(
    path: &std::path::Path,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
) -> Result<()> {
    // Smallest possible valid GC package that only covers fragment lengths 10..=60 and a single
    // GC bin spanning 0..=100.
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 60],
        gc_edges: vec![0, 101],
        correction_matrix: array![[1.0_f64]],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint,
    };
    package.write_zarr(path)?;
    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
fn make_simple_coverage_weights_config(
    out_dir: &std::path::Path,
    bam: &std::path::Path,
) -> CoverageWeightsConfig {
    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.to_path_buf(),
            output_dir: out_dir.to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_bin_size(20);
    cfg.set_stride(20);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_output_prefix("coverage".to_string());
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }
    cfg
}

#[test]
fn length_bin_range_spec_matches_brace_expansion_edges() -> Result<()> {
    // Arrange: Hand-derived expected edges for 100..220 with step 10.
    // The end is an edge (not a counted length), so we expect:
    // 100, 110, 120, ..., 220.
    let expected_edges = vec![
        100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220,
    ];

    let mut edge_list_config = base_midpoints_config_for_length_bins();
    edge_list_config.set_length_bins(expected_edges.clone());

    let mut range_spec_config = base_midpoints_config_for_length_bins();
    range_spec_config.set_length_bins_spec("100:220:10");

    // Act
    let edges_from_edge_list = edge_list_config.resolve_length_bins()?;
    let edges_from_range_spec = range_spec_config.resolve_length_bins()?;

    // Assert
    assert_eq!(edges_from_edge_list, expected_edges);
    assert_eq!(edges_from_range_spec, expected_edges);
    assert_eq!(edges_from_edge_list, edges_from_range_spec);

    Ok(())
}

#[test]
fn midpoints_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero() -> Result<()>
{
    // Arrange:
    // Count one group over one 11 bp window [45, 56). Use three identical 61 bp fragments
    // with midpoint exactly 50, so each accepted fragment contributes one count at profile
    // position 50 - 45 = 5.
    //
    // MAPQ setup:
    // - fragment A: MAPQ 60
    // - fragment B: MAPQ 0
    // - fragment C: MAPQ 30
    //
    // Therefore:
    // - default `min_mapq = 30`: counts A and C -> total mass 2 at position 5
    // - explicit `min_mapq = 30`: identical to default
    // - explicit `min_mapq = 0`: counts A, B, and C -> total mass 3 at position 5
    let fragment_with_mapq = |mapq: u8| -> FragmentSpec {
        let mut fragment = paired_fragment_on_tid(0, 20, 61, 20);
        fragment.forward.mapq = mapq;
        fragment.reverse.mapq = mapq;
        fragment
    };
    // All three fragments share the same start, so the strict fixture is required to keep them as
    // three distinct molecules instead of one reused qname.
    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 200)],
        vec![
            fragment_with_mapq(60),
            fragment_with_mapq(0),
            fragment_with_mapq(30),
        ],
        Vec::new(),
        "midpoints_default_min_mapq",
    )?;
    let temp = TempDir::new()?;
    let intervals = temp.path().join("sites.bed");
    write_bed(&intervals, &[("chr1", 45, 56, "groupA")])?;
    let out_default = TempDir::new()?;
    let out_thirty = TempDir::new()?;
    let out_zero = TempDir::new()?;

    let make_cfg = |out_dir: &std::path::Path, prefix: &str| {
        let mut cfg = MidpointsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
            intervals.clone(),
        );
        cfg.set_output_prefix(prefix);
        cfg.set_length_bins(vec![61, 62]);
        cfg.set_smoothing(MidpointSmoothing::None);
        cfg.set_require_proper_pair(false);
        cfg
    };

    let default_cfg = make_cfg(out_default.path(), "default");
    let mut explicit_thirty_cfg = make_cfg(out_thirty.path(), "explicit_thirty");
    explicit_thirty_cfg.set_min_mapq(30);
    let mut explicit_zero_cfg = make_cfg(out_zero.path(), "explicit_zero");
    explicit_zero_cfg.set_min_mapq(0);

    // Act
    run(&default_cfg)?;
    run(&explicit_thirty_cfg)?;
    run(&explicit_zero_cfg)?;

    // Assert
    let read_profiles = |dir: &TempDir, prefix: &str| -> Result<Array3<f32>> {
        let counts_path = dir.path().join(format!("{prefix}.midpoint_profiles.zarr"));
        read_midpoint_zarr_counts(&counts_path).map_err(Into::into)
    };

    let default_arr = read_profiles(&out_default, "default")?;
    let explicit_thirty_arr = read_profiles(&out_thirty, "explicit_thirty")?;
    let explicit_zero_arr = read_profiles(&out_zero, "explicit_zero")?;

    assert_eq!(default_arr.shape(), &[1, 1, 11]);
    assert_eq!(default_arr, explicit_thirty_arr);
    assert_eq!(default_arr[[0, 0, 5]], 2.0);
    assert_eq!(default_arr.sum(), 2.0);

    assert_eq!(explicit_zero_arr.shape(), &[1, 1, 11]);
    assert_eq!(explicit_zero_arr[[0, 0, 5]], 3.0);
    assert_eq!(explicit_zero_arr.sum(), 3.0);

    Ok(())
}

#[test]
fn unpaired_single_read_matches_paired_midpoint_profile_for_same_span() -> Result<()> {
    // Arrange:
    // Compare two representations of the same physical fragment span [20, 81):
    // - paired fragment of length 61
    // - one unpaired read with aligned span [20, 81)
    //
    // We use an odd fragment length so midpoint placement is deterministic:
    //   midpoint = 20 + floor(61 / 2) = 50
    // For one window [45, 56), that lands at profile position:
    //   50 - 45 = 5
    //
    // Both representations must therefore produce the same 3D midpoint profile array with a
    // single count at [group=0, length_bin=0, position=5].
    let paired_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_paired_parity",
    )?;
    let unpaired_bam = single_read_fragment_bam("midpoints_unpaired_parity", 20, 61)?;
    let paired_out = TempDir::new()?;
    let unpaired_out = TempDir::new()?;
    let intervals = paired_out.path().join("sites.bed");
    write_bed(&intervals, &[("chr1", 45, 56, "groupA")])?;

    let make_cfg = |bam_path: &std::path::Path, out_dir: &std::path::Path, unpaired: bool| {
        let mut cfg = MidpointsConfig::new(
            IOCArgs {
                bam: bam_path.to_path_buf(),
                output_dir: out_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
            intervals.clone(),
        );
        cfg.set_output_prefix("sites");
        cfg.set_length_bins(vec![61, 62]);
        cfg.set_smoothing(MidpointSmoothing::None);
        cfg.set_tile_size(1_000);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.unpaired.reads_are_fragments = unpaired;
        cfg
    };

    let paired_cfg = make_cfg(&paired_bam.bam, paired_out.path(), false);
    let unpaired_cfg = make_cfg(&unpaired_bam.bam, unpaired_out.path(), true);

    // Act
    run(&paired_cfg)?;
    run(&unpaired_cfg)?;

    // Assert
    let paired_arr: Array3<f32> =
        read_midpoint_zarr_counts(paired_out.path().join("sites.midpoint_profiles.zarr"))?;
    let unpaired_arr: Array3<f32> =
        read_midpoint_zarr_counts(unpaired_out.path().join("sites.midpoint_profiles.zarr"))?;

    assert_eq!(paired_arr, unpaired_arr);
    assert_eq!(paired_arr.shape(), &[1, 1, 11]);
    assert_eq!(paired_arr[[0, 0, 5]], 1.0);
    assert_eq!(paired_arr.sum(), 1.0);

    Ok(())
}

#[test]
fn bed_sites_mixed_core_and_halo_rows_keep_only_the_core_midpoint_count_across_tile_sizes()
-> Result<()> {
    // Arrange:
    // - one unpaired fragment span [5,16), so the deterministic midpoint is 10
    // - BED site [10,11) is the true core-overlap site and should receive one count at position 0
    // - BED site [22,23) is downstream and should stay zero
    // - with tile_size=10 the two sites fall in different tiles; with tile_size=1000 they do not
    // - the final grouped midpoint profiles must therefore be identical across tile sizes
    let bam = single_read_fragment_bam("midpoints_mixed_core_and_halo_rows", 5, 11)?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let temp = TempDir::new()?;
        let bed_path = temp.path().join(format!("sites_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 10, 11, "group_core"),
                ("chr1", 22, 23, "group_halo"),
            ],
        )?;

        let mut cfg = MidpointsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: temp.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
            bed_path,
        );
        cfg.set_output_prefix("sites");
        cfg.set_length_bins(vec![11, 12]);
        cfg.set_smoothing(MidpointSmoothing::None);
        cfg.set_tile_size(tile_size);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.unpaired.reads_are_fragments = true;

        run(&cfg)?;

        let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
        let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
        let map_path = temp.path().join("sites.group_index.tsv");
        let group_to_idx = read_group_index_map(&map_path)?;

        outputs.push((arr, group_to_idx));
    }

    let (small_tile_arr, small_tile_groups) = &outputs[0];
    let (large_tile_arr, large_tile_groups) = &outputs[1];

    assert_eq!(small_tile_groups, large_tile_groups);
    assert_eq!(small_tile_arr, large_tile_arr);
    assert_eq!(small_tile_arr.shape(), &[2, 1, 1]);
    assert_eq!(small_tile_arr[[small_tile_groups["group_core"], 0, 0]], 1.0);
    assert_eq!(small_tile_arr[[small_tile_groups["group_halo"], 0, 0]], 0.0);
    assert_eq!(small_tile_arr.sum(), 1.0);

    Ok(())
}

#[test]
fn later_tile_site_keeps_midpoint_count_when_window_span_starts_after_zero() -> Result<()> {
    // Arrange:
    // - Two one-base sites are present on the chromosome, so the later tile's cached span starts
    //   at chromosome-wide window index 1.
    // - The unpaired fragment spans [37,48), has length 11, and has deterministic midpoint 42.
    // - With tile_size=20, the midpoint belongs to tile [40,60), where the tile-local window list
    //   contains only [42,43). The window scan must therefore start at local index 0.
    // - With tile_size=1000, both sites are in one tile. The final profiles must match the
    //   multi-tile run exactly.
    let bam = single_read_fragment_bam("midpoints_later_tile_window_span_pointer", 37, 11)?;
    let tile_sizes = [20_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let temp = TempDir::new()?;
        let bed_path = temp.path().join(format!("sites_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 10, 11, "group_early"),
                ("chr1", 42, 43, "group_late"),
            ],
        )?;

        let mut cfg = MidpointsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: temp.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
            bed_path,
        );
        cfg.set_output_prefix("sites");
        cfg.set_length_bins(vec![11, 12]);
        cfg.set_smoothing(MidpointSmoothing::None);
        cfg.set_tile_size(tile_size);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.unpaired.reads_are_fragments = true;

        run(&cfg)?;

        let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
        let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
        let map_path = temp.path().join("sites.group_index.tsv");
        let group_to_idx = read_group_index_map(&map_path)?;

        outputs.push((arr, group_to_idx));
    }

    let (small_tile_arr, small_tile_groups) = &outputs[0];
    let (large_tile_arr, large_tile_groups) = &outputs[1];

    assert_eq!(small_tile_groups, large_tile_groups);
    assert_eq!(small_tile_arr, large_tile_arr);
    assert_eq!(small_tile_arr.shape(), &[2, 1, 1]);
    assert_eq!(
        small_tile_arr[[small_tile_groups["group_early"], 0, 0]],
        0.0
    );
    assert_eq!(small_tile_arr[[small_tile_groups["group_late"], 0, 0]], 1.0);
    assert_eq!(small_tile_arr.sum(), 1.0);

    Ok(())
}

#[test]
fn core_overlap_bed_site_is_kept_for_midpoints() -> Result<()> {
    // Arrange:
    // - one unpaired fragment span [5,16), so the deterministic midpoint is 10
    // - BED site [10,11) overlaps that midpoint and must receive one count
    let bam = single_read_fragment_bam("midpoints_core_site", 5, 11)?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("sites.bed");
    write_bed(&bed_path, &[("chr1", 10, 11, "group_core")])?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![11, 12]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(10);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.unpaired.reads_are_fragments = true;

    run(&cfg)?;

    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;
    let group_to_idx = read_group_index_map(&temp.path().join("sites.group_index.tsv"))?;
    assert_eq!(arr.shape(), &[1, 1, 1]);
    assert_eq!(arr[[group_to_idx["group_core"], 0, 0]], 1.0);
    assert_eq!(arr.sum(), 1.0);
    Ok(())
}

#[test]
fn stranded_bed_intervals_mirror_midpoint_positions_in_command_output() -> Result<()> {
    // Arrange:
    // - One unpaired read-as-fragment spans [6,17), so its midpoint is 11.
    // - The same genomic interval [10,15) is supplied twice, once as + and once as -.
    // - Forward position is 11 - 10 = 1.
    // - Reverse position mirrors within the half-open interval: (15 - 1) - 11 = 3.
    // - Because both intervals overlap the same midpoint, each group receives one count at its
    //   strand-oriented position.
    let bam = single_read_fragment_bam("midpoints_stranded_profile_mirror", 6, 11)?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("stranded_sites.bed");
    std::fs::write(
        &bed_path,
        "chr1\t10\t15\tgroup_plus\t0\t+\nchr1\t10\t15\tgroup_minus\t0\t-\n",
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![11, 12]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.unpaired.reads_are_fragments = true;

    // Act
    run(&cfg)?;

    // Assert
    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;
    let group_to_idx = read_group_index_map(&temp.path().join("sites.group_index.tsv"))?;
    let group_plus = group_to_idx["group_plus"];
    let group_minus = group_to_idx["group_minus"];

    assert_eq!(arr.shape(), &[2, 1, 5]);
    assert_eq!(
        arr.slice(ndarray::s![group_plus, 0, ..]).to_vec(),
        vec![0.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        arr.slice(ndarray::s![group_minus, 0, ..]).to_vec(),
        vec![0.0, 0.0, 0.0, 1.0, 0.0]
    );
    assert_eq!(arr.sum(), 2.0);

    Ok(())
}

#[test]
fn even_length_midpoint_tie_counts_exactly_one_of_two_adjacent_edge_windows() -> Result<()> {
    // Arrange:
    // One even-length fragment spans [40,50), so `midpoints` randomizes the tie and places the
    // midpoint at either 44 or 45.
    //
    // Put two 1 bp windows exactly on those two candidate midpoint bases:
    //   groupA -> [44,45)
    //   groupB -> [45,46)
    //
    // Exactly one of them must receive the count, never both and never neither.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 100)],
        vec![paired_fragment_on_tid(0, 40, 10, 5)],
        Vec::new(),
        "midpoints_even_tie_two_edge_windows",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[("chr1", 44, 45, "groupA"), ("chr1", 45, 46, "groupB")],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![10, 11]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Act
    run(&cfg)?;

    // Assert
    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;
    let group_to_idx = read_group_index_map(&temp.path().join("sites.group_index.tsv"))?;
    assert_eq!(arr.shape(), &[2, 1, 1]);
    assert_eq!(arr.sum(), 1.0);

    let group_a = arr[[group_to_idx["groupA"], 0, 0]];
    let group_b = arr[[group_to_idx["groupB"], 0, 0]];
    let is_valid_one_hot = (group_a == 1.0 && group_b == 0.0) || (group_a == 0.0 && group_b == 1.0);
    assert!(
        is_valid_one_hot,
        "even-length midpoint tie must count exactly one adjacent edge window, got groupA={group_a}, groupB={group_b}"
    );

    Ok(())
}

#[test]
fn blacklist_midpoint_filtering_checks_both_centers_for_even_fragments() -> Result<()> {
    // Arrange:
    // Use the same even-length fragment [40,50), whose placement midpoint is randomized between
    // 44 and 45. Count it against one 2 bp window [44,46), so without blacklist the profile row
    // must be one-hot:
    //   [1,0] if midpoint=44
    //   [0,1] if midpoint=45
    //
    // Now choose blacklist strategy `Midpoint`.
    // The shared blacklist helper is intentionally conservative for even-length fragments:
    // either central base can blacklist the fragment.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 100)],
        vec![paired_fragment_on_tid(0, 40, 10, 5)],
        Vec::new(),
        "midpoints_blacklist_midpoint_central_base_contract",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    let left_center_blacklist_path = temp.path().join("blacklist_left_center.bed");
    let right_center_blacklist_path = temp.path().join("blacklist_right_center.bed");
    write_bed(&bed_path, &[("chr1", 44, 46, "groupA")])?;
    std::fs::write(&left_center_blacklist_path, "chr1\t44\t45\n")?;
    std::fs::write(&right_center_blacklist_path, "chr1\t45\t46\n")?;

    let make_cfg = |output_dir: &std::path::Path, blacklist: Option<Vec<PathBuf>>| {
        let mut cfg = MidpointsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: output_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
            bed_path.clone(),
        );
        cfg.set_output_prefix("sites");
        cfg.set_length_bins(vec![10, 11]);
        cfg.set_smoothing(MidpointSmoothing::None);
        cfg.set_tile_size(1_000);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scale_genome(ScaleGenomeArgs::default());
        cfg.blacklist = blacklist;
        cfg.blacklist_strategy = cfdnalab::blacklist::BlacklistStrategy::Midpoint;
        // This test targets fragment-level midpoint blacklist behavior. The interval prefilter
        // would correctly remove this tiny test window before the fragment filter runs
        cfg.set_keep_blacklisted_intervals(true);
        cfg
    };

    let baseline_out = temp.path().join("baseline");
    let left_blacklisted_out = temp.path().join("blacklisted_left_center");
    let right_blacklisted_out = temp.path().join("blacklisted_right_center");
    std::fs::create_dir_all(&baseline_out)?;
    std::fs::create_dir_all(&left_blacklisted_out)?;
    std::fs::create_dir_all(&right_blacklisted_out)?;

    // Act
    run(&make_cfg(&baseline_out, None))?;
    run(&make_cfg(
        &left_blacklisted_out,
        Some(vec![left_center_blacklist_path]),
    ))?;
    run(&make_cfg(
        &right_blacklisted_out,
        Some(vec![right_center_blacklist_path]),
    ))?;

    // Assert
    let baseline_arr: Array3<f32> =
        read_midpoint_zarr_counts(baseline_out.join("sites.midpoint_profiles.zarr"))?;
    assert_eq!(baseline_arr.shape(), &[1, 1, 2]);
    let baseline_row = baseline_arr.slice(ndarray::s![0, 0, ..]).to_vec();
    let baseline_is_valid_one_hot =
        baseline_row == vec![1.0, 0.0] || baseline_row == vec![0.0, 1.0];
    assert!(
        baseline_is_valid_one_hot,
        "without blacklist the even-length midpoint must land at exactly one central base, got {:?}",
        baseline_row
    );

    for (case_name, output_dir) in [
        ("left central base", left_blacklisted_out),
        ("right central base", right_blacklisted_out),
    ] {
        let blacklisted_arr: Array3<f32> =
            read_midpoint_zarr_counts(output_dir.join("sites.midpoint_profiles.zarr"))?;
        assert_eq!(blacklisted_arr.shape(), &[1, 1, 2]);
        assert_eq!(
            blacklisted_arr.sum(),
            0.0,
            "midpoint blacklist filtering should remove an even fragment when the {case_name} is blacklisted"
        );
    }

    Ok(())
}

#[test]
fn keep_blacklisted_intervals_keeps_sites_but_still_filters_fragments() -> Result<()> {
    // Arrange:
    // The 11 bp fragment spans [45,56), so its midpoint is 50 and would count at position 0 in
    // the output interval [50,57). The blacklist [45,46) overlaps the fragment span, and it also
    // lies inside the interval-level safety margin. Setting `keep_blacklisted_intervals` proves
    // that the site remains available while fragment-level blacklist filtering still removes the
    // fragment contribution.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 45, 11, 5)],
        Vec::new(),
        "midpoints_keep_blacklisted_intervals_fragment_filter",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    let blacklist_path = temp.path().join("blacklist.bed");
    write_bed(&bed_path, &[("chr1", 50, 57, "groupA")])?;
    std::fs::write(&blacklist_path, "chr1\t45\t46\n")?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![11, 12]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.blacklist = Some(vec![blacklist_path]);
    cfg.blacklist_strategy = BlacklistStrategy::Any;
    cfg.set_keep_blacklisted_intervals(true);

    // Act
    run(&cfg)?;

    // Assert
    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;
    assert_eq!(arr.shape(), &[1, 1, 7]);
    assert_eq!(
        arr.sum(),
        0.0,
        "fragment-level blacklist filtering should still remove the only fragment"
    );

    Ok(())
}

#[test]
fn midpoint_prefilter_fails_clearly_when_blacklist_drops_all_intervals() -> Result<()> {
    // Arrange:
    // With length bin [11,12), the interval blacklist margin is ceil(11 / 2) = 6 bp.
    // The output interval [50,57) plus that margin is [44,63), which overlaps [45,46).
    // Because `keep_blacklisted_intervals` is false, the only interval is prefiltered away.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 45, 11, 5)],
        Vec::new(),
        "midpoints_all_intervals_dropped_by_blacklist",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    let blacklist_path = temp.path().join("blacklist.bed");
    write_bed(&bed_path, &[("chr1", 50, 57, "groupA")])?;
    std::fs::write(&blacklist_path, "chr1\t45\t46\n")?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![11, 12]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.blacklist = Some(vec![blacklist_path]);
    cfg.blacklist_strategy = BlacklistStrategy::Any;

    // Act
    let error = run(&cfg).expect_err("blacklist prefiltering should drop the only interval");
    let message = error.to_string();

    // Assert
    assert!(
        message.contains("No midpoint intervals remain after filtering"),
        "unexpected all-dropped error: {message}"
    );
    assert!(
        message.contains("Blacklist prefiltering dropped 1 interval(s)"),
        "all-dropped error should name blacklist prefiltering: {message}"
    );

    Ok(())
}

#[test]
fn midpoint_command_smooths_full_resolution_profile_before_final_binning() -> Result<()> {
    // Arrange:
    // Build an expanded 13-position profile with counts 1..13 around output interval [50,57).
    // A 7 bp SavGol window has radius 3, so smoothing preserves the linear profile at retained
    // centers 4..10. Final 3 bp binning then averages [4,5,6] -> 5, [7,8,9] -> 8, [10] -> 10.
    let mut fragments = Vec::new();
    for midpoint in 47_i64..=59 {
        let count_at_position = (midpoint - 46) as usize;
        for _ in 0..count_at_position {
            fragments.push(paired_fragment_on_tid(0, midpoint - 5, 11, 5));
        }
    }
    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
        "midpoints_savgol_then_bin_command",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 50, 57, "groupA")])?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![11, 12]);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_smoothing(MidpointSmoothing::SavGol { window_bp: 7 });
    cfg.set_bin_size(3);

    // Act
    run(&cfg)?;

    // Assert
    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;
    assert_eq!(arr.shape(), &[1, 1, 3]);
    assert!((arr[[0, 0, 0]] - 5.0).abs() <= MIDPOINT_F32_TOL);
    assert!((arr[[0, 0, 1]] - 8.0).abs() <= MIDPOINT_F32_TOL);
    assert!((arr[[0, 0, 2]] - 10.0).abs() <= MIDPOINT_F32_TOL);

    let settings_path = temp.path().join("sites.midpoint_settings.json");
    let settings_text = std::fs::read_to_string(&settings_path)?;
    let settings: Value = serde_json::from_str(&settings_text)?;
    assert_eq!(
        settings["array_axes"],
        serde_json::json!(["group", "length_bin", "position"])
    );
    assert_eq!(settings["length_axis"]["min_fragment_length"], 11);
    assert_eq!(settings["length_axis"]["max_fragment_length"], 11);
    assert_eq!(settings["length_axis"]["n_bins"], 1);
    assert_eq!(settings["position_axis"]["output_interval_length_bp"], 7);
    assert_eq!(settings["position_axis"]["counted_interval_length_bp"], 13);
    assert_eq!(settings["position_axis"]["n_bins"], 3);
    assert_eq!(settings["position_axis"]["bin_size_bp"], 3);
    assert_eq!(settings["position_axis"]["bin_aggregation"], "mean");
    assert_eq!(settings["position_axis"]["last_bin_width_bp"], 1);
    assert_eq!(settings["smoothing"]["method"], "savitzky_golay");
    assert_eq!(settings["smoothing"]["polynomial_order"], 3);
    assert_eq!(settings["smoothing"]["window_bp"], 7);
    assert_eq!(settings["smoothing"]["computation_flank_bp"], 3);
    assert_eq!(settings["smoothing"]["applied_before_binning"], true);

    Ok(())
}

#[test]
fn length_bin_start_end_list_format_is_rejected() {
    // Arrange: This format was intentionally removed.
    let mut config = base_midpoints_config_for_length_bins();
    config.set_length_bins_spec("30-80,80-150");

    // Act
    let error = config
        .resolve_length_bins()
        .expect_err("start-end list format should fail");

    // Assert
    assert!(
        format!("{error}").contains("explicit start-end lists are not supported"),
        "Unexpected error message: {error}"
    );
}

#[test]
fn midpoint_profiles_written_with_group_index() -> Result<()> {
    // Arrange:
    // This fixture pins the public midpoint profile contract, not the current implementation
    // details. Positive fragments cover both length bins, both chromosomes, and one group whose
    // windows are split across chromosomes. The negative fragments prove that plausible
    // off-contract records do not leak into the output.
    let fixture = midpoint_axis_contract_fixture()?;
    let temp = TempDir::new()?;
    let cfg = midpoint_axis_contract_config(&fixture, temp.path(), 2, "sites");

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    assert!(counts_path.exists());
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    let expected = expected_midpoint_axis_contract_counts();
    assert_eq!(arr.shape(), &[3, 2, 40]);
    assert_eq!(
        arr, expected,
        "midpoint profile axes must be (group, length_bin, position), with excluded fragments contributing no mass"
    );
    assert_eq!(arr.sum(), 6.0);

    let map_path = temp.path().join("sites.group_index.tsv");
    let zarr_group_idx = read_midpoint_zarr_i32_1d(&counts_path, "/group")?;
    assert_eq!(
        zarr_group_idx,
        vec![0, 1, 2],
        "Zarr group coordinate must match group_idx order in group_index.tsv"
    );
    let group_to_idx = read_group_index_map(&map_path)?;
    let group_to_eligible_intervals = read_group_index_eligible_intervals(&map_path)?;
    assert_eq!(
        group_to_idx,
        HashMap::from([
            ("groupA".to_string(), 0usize),
            ("groupC".to_string(), 1usize),
            ("groupB".to_string(), 2usize)
        ])
    );
    assert_eq!(
        group_to_eligible_intervals,
        HashMap::from([
            ("groupA".to_string(), 2usize),
            ("groupC".to_string(), 2usize),
            ("groupB".to_string(), 2usize)
        ]),
        "eligible_intervals should count profile intervals, independent of fragment count"
    );

    let settings_path = temp.path().join("sites.midpoint_settings.json");
    assert!(settings_path.exists());

    Ok(())
}

#[test]
fn midpoint_profiles_zarr_directory_is_replaced_on_rerun_with_same_prefix() -> Result<()> {
    // Arrange:
    // Zarr outputs are directories. This test covers the final-output move path that removes the
    // previous complete Zarr directory before moving the new one into place.
    let fixture = midpoint_axis_contract_fixture()?;
    let temp = TempDir::new()?;
    let cfg = midpoint_axis_contract_config(&fixture, temp.path(), 2, "sites");
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");

    // Act
    run(&cfg)?;
    let stale_marker = counts_path.join("stale_file_from_previous_run.txt");
    std::fs::write(&stale_marker, "old run")?;
    assert!(stale_marker.exists());
    run(&cfg)?;

    // Assert
    assert!(
        !stale_marker.exists(),
        "rerunning midpoints with the same prefix should replace the old Zarr directory"
    );
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr, expected_midpoint_axis_contract_counts());

    Ok(())
}

#[test]
fn group_index_counts_eligible_intervals_after_prefilter() -> Result<()> {
    // Arrange:
    // Length bin [11,12) gives an interval blacklist margin of ceil(11 / 2) = 6 bp.
    //
    // groupA has two input intervals:
    // - [50,57) is dropped because [50,57) expanded by 6 bp -> [44,63), overlapping [44,46).
    // - [120,127) is retained.
    //
    // groupB has one input interval:
    // - [80,87) is dropped because [80,87) expanded by 6 bp -> [74,93), overlapping [74,75).
    //
    // The BAM contains no fragments. The expected counts therefore prove that `eligible_intervals`
    // describes the retained profile interval set, not observed coverage.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        Vec::new(),
        "midpoints_group_index_eligible_intervals",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    let blacklist_path = temp.path().join("blacklist.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 50, 57, "groupA"),
            ("chr1", 80, 87, "groupB"),
            ("chr1", 120, 127, "groupA"),
        ],
    )?;
    std::fs::write(&blacklist_path, "chr1\t44\t46\nchr1\t74\t75\n")?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![11, 12]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.blacklist = Some(vec![blacklist_path]);
    cfg.blacklist_strategy = BlacklistStrategy::Any;

    // Act
    run(&cfg)?;

    // Assert
    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;
    assert_eq!(arr.shape(), &[2, 1, 7]);
    assert_eq!(
        arr.sum(),
        0.0,
        "the fixture has no fragments, so eligible interval counts must not depend on observed coverage"
    );

    let group_to_eligible_intervals =
        read_group_index_eligible_intervals(&temp.path().join("sites.group_index.tsv"))?;
    assert_eq!(
        group_to_eligible_intervals,
        HashMap::from([
            ("groupA".to_string(), 1usize),
            ("groupB".to_string(), 0usize)
        ]),
        "eligible_intervals should reflect intervals left after interval-level blacklist prefiltering"
    );

    Ok(())
}

#[cfg(feature = "cli")]
#[test]
fn midpoint_profiles_are_identical_across_thread_counts() -> Result<()> {
    // Arrange:
    // Run the command in separate processes so each invocation owns a fresh Rayon global pool.
    // The intended user contract is deterministic output, independent of how tile work is
    // scheduled across workers. CI runners can expose only one core, so only request thread counts
    // that the current machine can actually run.
    let fixture = midpoint_axis_contract_fixture()?;
    let mut observed_outputs = Vec::new();
    let bam_path = path_text(&fixture.bam.bam);
    let intervals_path = path_text(&fixture.intervals_path);
    let available_parallelism = std::thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(1);
    let thread_counts: Vec<usize> = [1_usize, 2, 4, 8]
        .into_iter()
        .filter(|n_threads| *n_threads <= available_parallelism)
        .collect();
    if thread_counts.len() < 2 {
        return Ok(());
    }

    // Act
    for n_threads in thread_counts {
        let output_dir = TempDir::new()?;
        let output_dir_text = path_text(output_dir.path());
        let output_prefix = format!("sites_threads_{n_threads}");
        let n_threads_text = n_threads.to_string();
        let output = command_output(
            "midpoints",
            &[
                "--bam",
                bam_path.as_str(),
                "--output-dir",
                output_dir_text.as_str(),
                "--chromosomes",
                "chr1",
                "chr2",
                "--n-threads",
                n_threads_text.as_str(),
                "--intervals",
                intervals_path.as_str(),
                "--min-mapq",
                "30",
                "--length-bins",
                "20",
                "60",
                "120",
                "--smoothing",
                "none",
                "--tile-size",
                "1000000",
                "--output-prefix",
                output_prefix.as_str(),
            ],
        )?;
        assert_success_with_logs(
            &output,
            &format!("cfdna midpoints with {n_threads} thread(s)"),
        );

        let arr: Array3<f32> = read_midpoint_zarr_counts(
            output_dir
                .path()
                .join(format!("{output_prefix}.midpoint_profiles.zarr")),
        )?;
        let group_index_text = std::fs::read_to_string(
            output_dir
                .path()
                .join(format!("{output_prefix}.group_index.tsv")),
        )?;
        observed_outputs.push((n_threads, output_dir, arr, group_index_text));
    }

    // Assert
    let expected = expected_midpoint_axis_contract_counts();
    let reference_group_index_text = observed_outputs[0].3.clone();
    for (n_threads, _output_dir, arr, group_index_text) in observed_outputs {
        assert_eq!(
            arr, expected,
            "{n_threads} thread(s) should produce the hand-derived midpoint profile"
        );
        assert_eq!(
            group_index_text, reference_group_index_text,
            "{n_threads} thread(s) should produce the same group index artifact"
        );
    }

    Ok(())
}

#[test]
fn group_index_axis_matches_first_group_encounter_order_and_collapsed_counts() -> Result<()> {
    // Arrange:
    // BED rows are sorted by chromosome/start as required, but group names are intentionally
    // interleaved:
    //   chr1  [45,56)   groupB   -> first new group encountered, so index 0
    //   chr1  [65,76)   groupC   -> second new group encountered, so index 1
    //   chr2  [85,96)   groupA   -> third new group encountered, so index 2
    //   chr2  [105,116) groupA   -> same group, so it reuses index 2
    //
    // Fragments are chosen so every midpoint lands at position 5 inside its window:
    // - [20,81)  midpoint 50  -> groupB window [45,56)   -> position 5
    // - [40,101) midpoint 70  -> groupC window [65,76)   -> position 5
    // - [60,121) midpoint 90  -> groupA window [85,96)   -> position 5
    // - [80,141) midpoint 110 -> groupA window [105,116) -> position 5
    //
    // Therefore the collapsed profiles must be:
    // - groupB (axis 0): one count at position 5
    // - groupC (axis 1): one count at position 5
    // - groupA (axis 2): two counts at position 5
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 61, 20),
            paired_fragment_on_tid(0, 40, 61, 20),
            paired_fragment_on_tid(1, 60, 61, 20),
            paired_fragment_on_tid(1, 80, 61, 20),
        ],
        Vec::new(),
        "midpoints_group_axis_contract",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 45, 56, "groupB"),
            ("chr1", 65, 76, "groupC"),
            ("chr2", 85, 96, "groupA"),
            ("chr2", 105, 116, "groupA"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(40);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[3, 1, 11]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    assert_eq!(
        group_to_idx,
        HashMap::from([
            ("groupB".to_string(), 0usize),
            ("groupC".to_string(), 1usize),
            ("groupA".to_string(), 2usize),
        ])
    );

    let expected_rows = [
        (
            "groupB",
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
        (
            "groupC",
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
        (
            "groupA",
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
    ];
    for (group_name, expected_row) in expected_rows {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_eq!(
            row, expected_row,
            "array axis for {group_name} must match the written group index map"
        );
    }
    assert_eq!(arr.sum(), 4.0);

    Ok(())
}

#[test]
fn real_ref_gc_bias_then_gc_bias_package_is_neutral_in_single_bin_case_for_midpoints() -> Result<()>
{
    // Arrange:
    // Use one odd-length fragment so midpoint placement is deterministic rather than randomly
    // split across the two central bases of an even-length fragment.
    //
    // Fragment:
    // - span [20, 81), length 61
    // - midpoint = 20 + floor(61 / 2) = 50
    //
    // Window:
    // - [45, 56), length 11
    // - midpoint position inside the window = 50 - 45 = 5
    //
    // Real GC artifact derivation:
    // - `ref-gc-bias` is run for exactly one fragment length: 61 bp
    // - `gc-bias` is then run on exactly one 61 bp sample fragment over the same repeated reference
    // - all mass therefore lands in one GC-by-length cell on both sides
    // - after normalization and ratio, the produced correction is 1.0
    //
    // So the final midpoint profile must be exactly one count at position 5.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_real_gc_neutral",
    )?;
    let reference = simple_reference_twobit()?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;
    let gc_path = build_real_neutral_gc_package(&bam.bam, &reference.path, temp.path(), 61)?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);
    assert_eq!(arr.sum(), 1.0);
    assert_eq!(
        arr.slice(ndarray::s![0, 0, ..]).to_vec(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    assert_eq!(
        group_to_idx,
        HashMap::from([("groupA".to_string(), 0usize)])
    );

    Ok(())
}

#[test]
fn gc_file_late_tile_site_uses_reference_coordinates_after_fetch_narrowing() -> Result<()> {
    // Arrange:
    // - The only site is [925,936), well inside a tile whose fetch origin is still 0.
    // - The reference is shorter than the BAM chromosome, but long enough for the narrowed
    //   window-derived fetch span. Reading the full tile reference would overrun the reference.
    // - The fragment [900,961) has deterministic midpoint 930, so it lands at position 5.
    // - The fragment interval [900,961) is all C, so it lands in the high-GC correction bin with
    //   weight 7.0. Using prefix-local origin 0 would see A-only sequence instead.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1_500)],
        vec![paired_fragment_on_tid(0, 900, 61, 20)],
        Vec::new(),
        "midpoints_late_tile_gc_origin",
    )?;
    let reference = twobit_from_sequences(
        "midpoints_late_tile_gc_origin_ref",
        vec![("chr1".to_string(), late_origin_gc_reference_sequence())],
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("late_site.bed");
    let gc_path = temp.path().join("two_bin_gc_package.zarr");
    write_bed(&bed_path, &[("chr1", 925, 936, "late")])?;
    write_two_bin_gc_package(
        &gc_path,
        61,
        2.0,
        7.0,
        twobit_contig_footprint(&reference.path)?,
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    run(&cfg)?;
    let arr: Array3<f32> =
        read_midpoint_zarr_counts(temp.path().join("sites.midpoint_profiles.zarr"))?;

    // Assert
    assert_eq!(arr.shape(), &[1, 1, 11]);
    assert_eq!(
        arr.slice(ndarray::s![0, 0, ..]).to_vec(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(arr.sum(), 7.0);
    Ok(())
}

#[test]
fn real_ref_gc_bias_then_gc_bias_package_changes_midpoints_in_expected_direction() -> Result<()> {
    // Arrange:
    // Build a real non-neutral GC package, then consume it through `midpoints`.
    //
    // Reference genome:
    // - chr1[0,100)    = all A
    // - chr1[100,101)  = one-base spacer that is excluded from the reference BED windows
    // - chr1[101,201)  = all C
    //
    // Real producer setup:
    // - fragment length is fixed at 61, so midpoint placement is deterministic
    // - valid starts are 0..=140, because 201 - 61 + 1 = 141
    // - reference BED windows keep only starts that both fall inside the BED row and leave room
    //   for the full 61 bp fragment:
    //     [0,100)     -> starts 0..=39     -> GC%=0
    //     [101,201)   -> starts 101..=140  -> GC%=100
    // - starts 40..=100 would cross the spacer or boundary and are intentionally excluded
    //
    // So the reference-side counts are exactly balanced:
    // - 40 starts at GC%=0
    // - 40 starts at GC%=100
    //
    // Producer BAM:
    // - one A-only fragment [10,71)   -> GC%=0
    // - nine C-only fragments [110,171) -> GC%=100
    //
    // The real produced GC package is therefore the same two-bin non-neutral package as in the
    // corresponding `gc-bias`/`fcoverage`/`lengths` tests:
    // - GC%=0   -> weight 5.0
    // - GC%=100 -> weight 5/9
    //
    // Consumer BAM:
    // - one A-only fragment [10,71), midpoint 40
    // - one C-only fragment [110,171), midpoint 140
    //
    // Consumer windows:
    // - [35,46)   -> midpoint 40 lands at position 5, groupA
    // - [135,146) -> midpoint 140 lands at position 5, groupC
    //
    // No genomic scaling is applied, so the final midpoint profile must contain:
    // - groupA: 5.0   at position 5
    // - groupC: 5/9   at position 5
    let reference = twobit_from_sequences(
        "midpoints_real_non_neutral_reference",
        vec![(
            "chr1".to_string(),
            format!("{}N{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    // The producer deliberately stacks nine identical C-only fragments at one start. Use the
    // strict helper so those molecules keep distinct qnames in the BAM pairing layer.
    let producer_bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 201)],
        {
            let mut fragments = vec![paired_fragment_on_tid(0, 10, 61, 20)];
            for _ in 0..9 {
                fragments.push(paired_fragment_on_tid(0, 110, 61, 20));
            }
            fragments
        },
        Vec::new(),
        "midpoints_real_non_neutral_producer",
    )?;
    let consumer_bam = bam_from_specs(
        vec![("chr1".to_string(), 201)],
        vec![
            paired_fragment_on_tid(0, 10, 61, 20),
            paired_fragment_on_tid(0, 110, 61, 20),
        ],
        Vec::new(),
        "midpoints_real_non_neutral_consumer",
    )?;
    let temp = TempDir::new()?;
    let gc_path = build_real_non_neutral_gc_package(
        &producer_bam.bam,
        &reference.path,
        temp.path(),
        61,
        "chr1\t0\t100\nchr1\t101\t201\n",
        // Chromosome length 201 and fragment length 61 give 141 valid starts in total. With
        // `ref-gc-bias` BED mode, only starts that fit entirely inside the same BED row count, so
        // these windows contribute exactly:
        // - starts 0..=39    -> 40 pure-A starts
        // - starts 101..=140 -> 40 pure-C starts
        141,
    )?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[("chr1", 35, 46, "groupA"), ("chr1", 135, 146, "groupC")],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: consumer_bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[2, 1, 11]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    let expected = [("groupA", 5.0_f32), ("groupC", (5.0_f32 / 9.0_f32))];
    for (group_name, expected_weight) in expected {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_midpoint_profile_row_matches(&row, expected_weight, group_name);
    }
    assert!(
        (arr.sum() - (5.0_f32 + 5.0_f32 / 9.0_f32)).abs() <= MIDPOINT_F32_TOL,
        "unexpected total midpoint mass {}",
        arr.sum()
    );

    Ok(())
}

#[test]
fn midpoints_rejects_gc_package_when_length_bins_are_outside_supported_range() -> Result<()> {
    // Arrange:
    // The midpoint command resolves its fragment length range from the configured bin edges:
    //   [61, 62] -> counted fragment lengths are exactly 61 bp.
    //
    // We then hand-build the smallest valid GC package that only covers lengths 10..=60.
    // The shared GC loader should therefore reject the package before any per-tile counting:
    //   requested range [61,61] is outside package range [10,60].
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_gc_length_range_mismatch",
    )?;
    let reference = simple_reference_twobit()?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;
    let gc_path = temp.path().join("too_short_gc_package.zarr");
    write_minimal_gc_package_excluding_length_61(
        &gc_path,
        twobit_contig_footprint(&reference.path)?,
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    let err = run(&cfg).expect_err("out-of-range GC package should fail");

    // Assert
    let msg = err.to_string();
    assert!(
        msg.contains("fragment length range [61-61] is outside the range covered by the correction package [10-60]"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn midpoints_rejects_gc_package_with_schema_version_mismatch() -> Result<()> {
    // Arrange:
    // Build the smallest valid GC correction package shape, but with an intentionally
    // incompatible schema version. `midpoints` should fail while loading the package, before
    // reading any GC weights or accumulating profile mass.
    let bam = fixtures::simple_inward_bam()?;
    let reference = simple_reference_twobit()?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;
    let gc_path = temp.path().join("gc_pkg_bad_version.zarr");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION + 1,
        end_offset: 0,
        length_edges: vec![10, 200],
        gc_edges: vec![0, 101],
        correction_matrix: array![[1.0_f64]],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
    };
    package.write_zarr(&gc_path)?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![60, 61]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    let err = run(&cfg).expect_err("schema version mismatch should fail");

    // Assert
    let msg = err.to_string();
    assert!(
        msg.contains("GC correction package schema version mismatch"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn coverage_weights_tsv_changes_midpoints_by_full_fragment_average_not_window_overlap() -> Result<()>
{
    // Arrange:
    // Producer BAM:
    // - `simple_inward_bam()` has exactly one fragment [20, 80) on a 200 bp chromosome.
    // - We run `coverage-weights` with `bin_size = stride = 20`.
    // - In that identity case, `avg_overlapping_pos_cov == avg_pos_cov` for every stride bin.
    // - The producer therefore writes per-bin scaling factors:
    //     [0,20):  0   (no coverage)
    //     [20,40): 1   (covered at depth 1, global mean over non-zero bins is also 1)
    //     [40,60): 1
    //     [60,80): 1
    //     [80,200): 0
    //
    // Consumer BAM:
    // - One odd-length fragment [20, 81), length 61.
    // - Odd length makes midpoint deterministic:
    //     midpoint = 20 + floor(61 / 2) = 50.
    // - One window [45, 56), so the midpoint lands at window position:
    //     50 - 45 = 5.
    //
    // Crucial scaling derivation for `midpoints`:
    // - `midpoints` averages scaling over the full fragment span, not only over the midpoint
    //   window or over the fragment/window overlap.
    // - The consumer fragment overlaps scaling bins as:
    //     [20,40): 20 bp with factor 1
    //     [40,60): 20 bp with factor 1
    //     [60,80): 20 bp with factor 1
    //     [80,81):  1 bp with factor 0
    // - Average scaling over the fragment is therefore:
    //     (20*1 + 20*1 + 20*1 + 1*0) / 61 = 60 / 61.
    // - No GC weighting is applied, so the final midpoint profile mass must be exactly 60/61 at
    //   position 5.
    let producer_bam = fixtures::simple_inward_bam()?;
    let consumer_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_scaling_consumer",
    )?;
    let temp = TempDir::new()?;
    let weights_out_dir = temp.path().join("coverage_weights");
    std::fs::create_dir_all(&weights_out_dir)?;
    let scaling_cfg = make_simple_coverage_weights_config(&weights_out_dir, &producer_bam.bam);
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;

    // Act
    run_coverage_weights(&scaling_cfg)?;
    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

    let mut midpoints_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: consumer_bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    midpoints_cfg.set_output_prefix("sites");
    midpoints_cfg.set_length_bins(vec![61, 62]);
    midpoints_cfg.set_smoothing(MidpointSmoothing::None);
    midpoints_cfg.set_tile_size(1_000);
    midpoints_cfg.set_min_mapq(0);
    midpoints_cfg.set_require_proper_pair(false);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    midpoints_cfg.set_scale_genome(scale_genome);

    run(&midpoints_cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);
    assert_eq!(arr.slice(ndarray::s![0, 0, ..]).len(), 11);

    let expected_weight = 60.0_f32 / 61.0_f32;
    let row = arr.slice(ndarray::s![0, 0, ..]).to_vec();
    assert_midpoint_profile_row_matches(&row, expected_weight, "single-group scaling");
    assert!(
        (arr.sum() - expected_weight).abs() <= MIDPOINT_F32_TOL,
        "expected total midpoint mass {expected_weight}, got {}",
        arr.sum()
    );

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_multi_chromosome_coverage_weights_tsv_is_applied_per_chromosome_in_midpoints() -> Result<()>
{
    // Arrange:
    // Build a real multi-chromosome scaling artifact, then consume it through `midpoints`.
    //
    // Producer BAM:
    // - chr1 has one 61 bp fragment [20, 81)
    // - chr2 has two identical 61 bp fragments [20, 81)
    //
    // We use `coverage-weights` with `bin_size = stride = 20`, so each TSV row is just the
    // average positional coverage inside one 20 bp bin.
    //
    // Per-bin producer coverage is therefore:
    // - chr1:
    //     [20,40): 1
    //     [40,60): 1
    //     [60,80): 1
    //     [80,100): 1/20   (only the last base 80 is covered)
    // - chr2:
    //     [20,40): 2
    //     [40,60): 2
    //     [60,80): 2
    //     [80,100): 2/20 = 1/10
    //
    // Shared global mean over the 8 non-zero bins:
    //   ((3 * 1) + 1/20 + (3 * 2) + 1/10) / 8
    // = (3 + 1/20 + 6 + 1/10) / 8
    // = (61/20 + 61/10) / 8
    // = 183/160.
    //
    // The written scaling factors are mean / avg_pos_cov:
    // - chr1 full bins: (183/160) / 1    = 183/160
    // - chr1 tail bin:  (183/160) / 1/20 = 183/8
    // - chr2 full bins: (183/160) / 2    = 183/320
    // - chr2 tail bin:  (183/160) / 1/10 = 183/16
    //
    // Consumer BAM:
    // - one 61 bp fragment [20,81) on chr1
    // - one 61 bp fragment [20,81) on chr2
    // - odd length makes midpoint deterministic:
    //     20 + floor(61 / 2) = 50
    // - both windows are [45,56), so each midpoint lands at profile position:
    //     50 - 45 = 5
    //
    // `midpoints` averages scaling over the full fragment span:
    // - chr1 average scaling:
    //     (20*(183/160) + 20*(183/160) + 20*(183/160) + 1*(183/8)) / 61
    //   = (183/160) * (60 + 20) / 61
    //   = 183/122
    //   = 1.5
    // - chr2 average scaling:
    //     (20*(183/320) + 20*(183/320) + 20*(183/320) + 1*(183/16)) / 61
    //   = (183/320) * (60 + 20) / 61
    //   = 183/244
    //   = 0.75
    //
    // No GC weighting is applied, so the final midpoint profile must contain:
    // - group_chr1: 1.5 at position 5
    // - group_chr2: 0.75 at position 5
    // chr2 deliberately stacks two identical fragments at one start. Use strict identity so the
    // producer really contains three molecules and the derived scaling TSV is correct.
    let producer_bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 61, 20),
            paired_fragment_on_tid(1, 20, 61, 20),
            paired_fragment_on_tid(1, 20, 61, 20),
        ],
        Vec::new(),
        "midpoints_multichrom_scaling_producer",
    )?;
    let consumer_bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 61, 20),
            paired_fragment_on_tid(1, 20, 61, 20),
        ],
        Vec::new(),
        "midpoints_multichrom_scaling_consumer",
    )?;
    let temp = TempDir::new()?;
    let weights_out_dir = temp.path().join("coverage_weights");
    std::fs::create_dir_all(&weights_out_dir)?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 45, 56, "group_chr1"),
            ("chr2", 45, 56, "group_chr2"),
        ],
    )?;

    let mut scaling_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: producer_bam.bam.clone(),
            output_dir: weights_out_dir.clone(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2"]),
    );
    scaling_cfg.set_bin_size(20);
    scaling_cfg.set_stride(20);
    scaling_cfg.set_min_mapq(0);
    scaling_cfg.set_require_proper_pair(false);
    scaling_cfg.set_output_prefix("coverage".to_string());
    {
        let frag = scaling_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    run_coverage_weights(&scaling_cfg)?;
    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

    let mut midpoints_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: consumer_bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2"]),
        bed_path,
    );
    midpoints_cfg.set_output_prefix("sites");
    midpoints_cfg.set_length_bins(vec![61, 62]);
    midpoints_cfg.set_smoothing(MidpointSmoothing::None);
    midpoints_cfg.set_tile_size(1_000);
    midpoints_cfg.set_min_mapq(0);
    midpoints_cfg.set_require_proper_pair(false);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    midpoints_cfg.set_scale_genome(scale_genome);
    run(&midpoints_cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[2, 1, 11]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    let expected_total = 1.5_f32 + 0.75_f32;
    for (group_name, expected_weight) in [("group_chr1", 1.5_f32), ("group_chr2", 0.75_f32)] {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_midpoint_profile_row_matches(&row, expected_weight, group_name);
    }
    assert!(
        (arr.sum() - expected_total).abs() <= MIDPOINT_F32_TOL,
        "expected total midpoint mass {expected_total}, got {}",
        arr.sum()
    );

    Ok(())
}

#[test]
fn gc_tag_pair_average_sets_midpoint_profile_weight() -> Result<()> {
    // Arrange:
    // - One paired fragment spans [20, 81), length 61, so the midpoint is deterministic:
    //     20 + floor(61 / 2) = 50.
    // - One window [45, 56) therefore receives the fragment at position:
    //     50 - 45 = 5.
    // - Mate GC tags are 2.0 and 4.0.
    // - The shared fragment-level GC-tag rule is to average two valid mate weights:
    //     (2.0 + 4.0) / 2 = 3.0.
    // - No genomic scaling is applied, so the final midpoint profile must contain exactly 3.0 at
    //   position 5 and 0 elsewhere.
    let base_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_gc_tag_base",
    )?;
    let tagged_bam = bam_with_gc_tags(
        &base_bam.bam,
        "midpoints_gc_tag_paired_avg",
        &[Some(2.0), Some(4.0)],
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: tagged_bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);

    let row = arr.slice(ndarray::s![0, 0, ..]).to_vec();
    assert_midpoint_profile_row_matches(&row, 3.0, "gc-tag");
    assert!(
        (arr.sum() - 3.0).abs() <= MIDPOINT_F32_TOL,
        "expected total midpoint mass 3.0, got {}",
        arr.sum()
    );

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn gc_file_and_scaling_tsv_weights_multiply_in_midpoints() -> Result<()> {
    // Arrange:
    // Producer BAM:
    // - `simple_inward_bam()` contains one fragment [20, 80) on chr1.
    // - Run `coverage-weights` with a neutral real GC package over its full configured
    //   fragment length range so the written scaling TSV is GC-compatible with the consumer
    //   command, but the numerical scaling profile stays unchanged.
    // - With `bin_size = stride = 20`, the written scaling TSV is therefore still the identity
    //   profile:
    //     [20,40): 1
    //     [40,60): 1
    //     [60,80): 1
    //     everything else: 0
    //
    // Consumer BAM:
    // - One odd-length fragment [20, 81), length 61.
    // - The midpoint is therefore deterministic:
    //     20 + floor(61 / 2) = 50.
    // - One window [45, 56) receives that midpoint at profile position:
    //     50 - 45 = 5.
    //
    // Scaling derivation:
    // - `midpoints` averages scaling over the full fragment span [20, 81):
    //     [20,40): 20 bp at factor 1
    //     [40,60): 20 bp at factor 1
    //     [60,80): 20 bp at factor 1
    //     [80,81):  1 bp at factor 0
    // - Average scaling over the fragment is therefore:
    //     (20 + 20 + 20 + 0) / 61 = 60 / 61.
    //
    // GC derivation:
    // - Use the smallest valid GC package for the only supported fragment length 61:
    //     length_edges = [61, 62]
    //     gc_edges     = [0, 101]
    //     correction_matrix = [[3.0]]
    // - Every accepted fragment therefore gets GC weight 3.0.
    //
    // Combination contract from `midpoints.rs`:
    // - The command increments counts by `scaling_weight * gc_weight`.
    // - So the only non-zero output cell must be:
    //     3.0 * (60 / 61) = 180 / 61
    //   at group 0, length bin 0, position 5.
    let producer_bam = fixtures::simple_inward_bam()?;
    let consumer_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_gc_and_scaling_consumer",
    )?;
    let reference = simple_reference_twobit()?;
    let temp = TempDir::new()?;
    let weights_out_dir = temp.path().join("coverage_weights");
    std::fs::create_dir_all(&weights_out_dir)?;
    let mut scaling_cfg = make_simple_coverage_weights_config(&weights_out_dir, &producer_bam.bam);
    let weights_gc_path = build_real_neutral_gc_package_for_range(
        &producer_bam.bam,
        &reference.path,
        temp.path(),
        10,
        200,
    )?;
    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");
    let gc_path = temp.path().join("constant_gc_pkg.zarr");
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![61, 62],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
        correction_matrix: array![[3.0_f64]],
    };
    package.write_zarr(&gc_path)?;
    scaling_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(weights_gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    scaling_cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    run_coverage_weights(&scaling_cfg)?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: consumer_bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    cfg.set_scale_genome(scale_genome);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);

    let expected_weight = 180.0_f32 / 61.0_f32;
    let row = arr.slice(ndarray::s![0, 0, ..]).to_vec();
    assert_midpoint_profile_row_matches(&row, expected_weight, "gc x scaling");
    assert!(
        (arr.sum() - expected_weight).abs() <= MIDPOINT_F32_TOL,
        "expected total midpoint mass {expected_weight}, got {}",
        arr.sum()
    );

    Ok(())
}

#[cfg(feature = "cmd_bam_to_bam")]
#[test]
fn bam_to_bam_gc_file_output_drives_midpoints_gc_tag_same_as_original_gc_file() -> Result<()> {
    // Arrange:
    // One paired fragment spans [20, 81), length 61, so the midpoint is deterministic:
    //   20 + floor(61 / 2) = 50
    // One window [45, 56) therefore receives the fragment at profile position:
    //   50 - 45 = 5
    //
    // We use the smallest GC package that assigns a constant weight 3.0 to every 61 bp fragment:
    // - length_edges = [61, 62]
    // - gc_edges     = [0, 101]
    // - correction_matrix = [[3.0]]
    //
    // Then we compare two logically equivalent released workflows:
    // 1. original paired BAM -> `midpoints --gc-file <pkg>`
    // 2. original paired BAM -> `bam-to-bam --gc-file <pkg>` ->
    //    `midpoints --gc-tag GC`
    //
    // Because the package gives the only supported fragment a constant weight 3.0, both
    // workflows must produce the same midpoint profile:
    // - shape [1, 1, 11]
    // - exactly 3.0 at position 5
    // - total mass 3.0
    let source_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "midpoints_bam_to_bam_gc_source",
    )?;
    let reference = simple_reference_twobit()?;
    let temp = TempDir::new()?;
    let tagged_out_bam = temp.path().join("tagged_gc.bam");
    let gc_path = temp.path().join("constant_gc_pkg.zarr");
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![61, 62],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
        correction_matrix: array![[3.0_f64]],
    };
    package.write_zarr(&gc_path)?;

    let mut bam_to_bam_cfg = BamToBamConfig::new(
        source_bam.bam.clone(),
        tagged_out_bam.clone(),
        base_chromosomes(&["chr1"]),
    );
    bam_to_bam_cfg.min_mapq = 0;
    bam_to_bam_cfg.set_gc(cfdnalab::run_like_cli::common::ApplyGCArgFileOnly {
        gc_file: Some(gc_path.clone()),
        neutralize_invalid_gc: false,
    });
    bam_to_bam_cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = bam_to_bam_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 61;
        frag.max_fragment_length = 61;
    }

    let mut original_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: temp.path().join("orig_out"),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path.clone(),
    );
    original_cfg.set_output_prefix("origsites");
    original_cfg.set_length_bins(vec![61, 62]);
    original_cfg.set_smoothing(MidpointSmoothing::None);
    original_cfg.set_tile_size(1_000);
    original_cfg.set_min_mapq(0);
    original_cfg.set_require_proper_pair(false);
    original_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    original_cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act 1: write the tagged BAM from the real `bam-to-bam` producer and index it for fetch-based
    // downstream consumers.
    run_bam_to_bam(&bam_to_bam_cfg)?;
    build_bai_for_test_bam(&tagged_out_bam)?;

    // Act 2: compare original `--gc-file` consumption with downstream `--gc-tag GC`.
    run(&original_cfg)?;
    let mut tagged_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: tagged_out_bam,
            output_dir: temp.path().join("tagged_out"),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    tagged_cfg.set_output_prefix("taggedsites");
    tagged_cfg.set_length_bins(vec![61, 62]);
    tagged_cfg.set_smoothing(MidpointSmoothing::None);
    tagged_cfg.set_tile_size(1_000);
    tagged_cfg.set_min_mapq(0);
    tagged_cfg.set_require_proper_pair(false);
    tagged_cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });
    run(&tagged_cfg)?;

    // Assert
    let original_arr: Array3<f32> = read_midpoint_zarr_counts(
        &temp
            .path()
            .join("orig_out/origsites.midpoint_profiles.zarr"),
    )?;
    let tagged_arr: Array3<f32> = read_midpoint_zarr_counts(
        &temp
            .path()
            .join("tagged_out/taggedsites.midpoint_profiles.zarr"),
    )?;

    assert_eq!(original_arr, tagged_arr);
    assert_eq!(original_arr.shape(), &[1, 1, 11]);
    assert_eq!(original_arr[[0, 0, 5]], 3.0);
    assert_eq!(original_arr.sum(), 3.0);

    Ok(())
}

#[test]
fn scaling_tsv_must_cover_requested_chromosome_end_in_midpoints() -> Result<()> {
    // Arrange:
    // `simple_inward_bam()` uses chr1 length 200.
    // `midpoints` loads scaling factors through the same shared TSV contract as the other
    // released consumers. A TSV that stops at 100 is malformed even if the counted fragment
    // and interval both lie inside that prefix.
    //
    // We use one interval [45,56) that would otherwise count the fixture fragment midpoint, so a
    // successful run would produce a single non-zero profile cell. The correct behavior here is to
    // fail before any counting because the scaling artifact does not cover the full chromosome.
    let bam = fixtures::simple_inward_bam()?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;
    let scaling_path = temp.path().join("truncated_scaling.tsv");
    std::fs::write(
        &scaling_path,
        "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t100\t2.0\n",
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        bed_path,
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![60, 61]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    cfg.set_scale_genome(scale_genome);

    // Act
    let err = run(&cfg).expect_err("truncated scaling TSV should fail");

    // Assert:
    // `midpoints` also wraps the shared loader with `load scaling factors`, so the actionable
    // artifact-contract message is only visible in the full error chain.
    let msg = format!("{err:#}");
    assert!(
        msg.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn midpoint_fetch_narrowing_preserves_tile_halo_near_chromosome_end_on_three_chromosomes()
-> Result<()> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 95),
            ("chr2".to_string(), 95),
            ("chr3".to_string(), 95),
        ],
        vec![
            paired_fragment_on_tid(0, 84, 11, 3),
            paired_fragment_on_tid(1, 84, 11, 3),
            paired_fragment_on_tid(2, 84, 11, 3),
        ],
        Vec::new(),
        "midpoints_chrom_end_halo_three_chr",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows_three_chr_near_end.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 89, 95, "groupA"),
            ("chr2", 89, 95, "groupB"),
            ("chr3", 89, 95, "groupC"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        bed_path.clone(),
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![10, 15]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(40);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Manual expectations:
    // - Each chromosome ends with a 6 bp site [89,95), which falls in the last tile [80,95).
    // - The only fragment on each chromosome is [84,95), length 11, midpoint 89.
    // - The midpoint lies at window position 89 - 89 = 0, so each group gets one count at
    //   length-bin [10,15) and position 0.
    // - This command-level fixture checks that narrowing to the extreme midpoint sites does not
    //   discard the fetch halo already carried by the last tile near chromosome end.
    // - It does not isolate the separate `halo_bp` argument to the narrowing helper, because the
    //   tile fetch band was already built with the same maximum-fragment-length halo.
    run(&cfg)?;

    let counts_path = temp.path().join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[3, 1, 6]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;

    assert_eq!(group_to_idx.len(), 3);
    assert_eq!(arr.sum(), 3.0);
    for group_name in ["groupA", "groupB", "groupC"] {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_eq!(
            row,
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "{group_name} should have exactly one midpoint count at position 0"
        );
    }

    Ok(())
}

#[test]
fn midpoint_fetch_narrowing_reads_all_eligible_fragments_near_chromosome_end_on_three_chromosomes()
-> Result<()> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 95),
            ("chr2".to_string(), 95),
            ("chr3".to_string(), 95),
        ],
        vec![
            paired_fragment_on_tid(0, 79, 11, 3),
            paired_fragment_on_tid(0, 80, 11, 3),
            paired_fragment_on_tid(0, 82, 11, 3),
            paired_fragment_on_tid(0, 84, 11, 3),
            paired_fragment_on_tid(1, 79, 11, 3),
            paired_fragment_on_tid(1, 80, 11, 3),
            paired_fragment_on_tid(1, 82, 11, 3),
            paired_fragment_on_tid(1, 84, 11, 3),
            paired_fragment_on_tid(2, 79, 11, 3),
            paired_fragment_on_tid(2, 80, 11, 3),
            paired_fragment_on_tid(2, 82, 11, 3),
            paired_fragment_on_tid(2, 84, 11, 3),
        ],
        Vec::new(),
        "midpoints_chrom_end_fetch_reads_all_eligible",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp
        .path()
        .join("windows_three_chr_fetch_read_coverage.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 85, 95, "groupA"),
            ("chr2", 85, 95, "groupB"),
            ("chr3", 85, 95, "groupC"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        bed_path,
    );
    cfg.set_output_prefix("sites_fetch_reads_all");
    cfg.set_length_bins(vec![10, 15]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(40);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Manual expectations:
    // - Each chromosome has one site [85,95), which lies in the last tile [80,95).
    // - Four fragments are present per chromosome, all length 11:
    //     * [79,90) midpoint 84 -> outside the site, so it must not be counted
    //     * [80,91) midpoint 85 -> counted at site position 0
    //     * [82,93) midpoint 87 -> counted at site position 2
    //     * [84,95) midpoint 89 -> counted at site position 4
    // - The narrowing step therefore has to preserve enough of the tile fetch band to read all
    //   three eligible fragments, not just the one closest to chromosome end.
    // - Each group row must therefore be exactly [1,0,1,0,1,0,0,0,0,0].
    run(&cfg)?;

    let counts_path = temp
        .path()
        .join("sites_fetch_reads_all.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[3, 1, 10]);

    let map_path = temp.path().join("sites_fetch_reads_all.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;

    assert_eq!(group_to_idx.len(), 3);
    assert_eq!(arr.sum(), 9.0);
    for group_name in ["groupA", "groupB", "groupC"] {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_eq!(
            row,
            vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "{group_name} should count exactly the three eligible near-end fragments"
        );
    }

    Ok(())
}
