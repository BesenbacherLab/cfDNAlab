#![cfg(feature = "cmd_midpoints")]

mod fixtures;

use anyhow::Result;
#[cfg(feature = "cmd_bam_to_bam")]
use cfdnalab::commands::bam_to_bam::{
    bam_to_bam::run_inner as run_bam_to_bam, config::BamToBamConfig,
};
use cfdnalab::commands::cli_common::{
    ApplyGCArgs, ChromosomeArgs, IOCArgs, ScaleGenomeArgs,
};
#[cfg(feature = "cmd_coverage_weights")]
use cfdnalab::commands::coverage_weights::{
    config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
};
use cfdnalab::commands::gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage};
use cfdnalab::commands::midpoints::config::MidpointsConfig;
use cfdnalab::commands::midpoints::midpoints::run;
use fixtures::{
    FragmentSpec, ReadSpec, bam_from_specs, build_real_neutral_gc_package,
    build_real_non_neutral_gc_package, complex_bam_fixture, simple_reference_twobit,
    twobit_from_sequences, write_bed,
};
use ndarray::array;
use ndarray::Array3;
use ndarray_npy::read_npy;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read, Reader};
use std::collections::HashMap;
use std::path::PathBuf;
use tempfile::TempDir;

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
            mapq: 60,
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
            mapq: 60,
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

fn write_minimal_gc_package_excluding_length_61(path: &std::path::Path) -> Result<()> {
    // Smallest possible valid GC package that only covers fragment lengths 10..=60 and a single
    // GC bin spanning 0..=100.
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 60],
        gc_edges: vec![0, 101],
        correction_matrix: array![[1.0_f64]],
        length_bin_frequencies: array![1.0_f64],
    };
    package.write_npz(path)?;
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
fn default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero() -> Result<()> {
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
    let bam = bam_from_specs(
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
        let counts_path = dir
            .path()
            .join(format!("{prefix}.midpoint_profiles.npy"));
        read_npy(&counts_path).map_err(Into::into)
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
    let paired_arr: Array3<f32> = read_npy(paired_out.path().join("sites.midpoint_profiles.npy"))?;
    let unpaired_arr: Array3<f32> =
        read_npy(unpaired_out.path().join("sites.midpoint_profiles.npy"))?;

    assert_eq!(paired_arr, unpaired_arr);
    assert_eq!(paired_arr.shape(), &[1, 1, 11]);
    assert_eq!(paired_arr[[0, 0, 5]], 1.0);
    assert_eq!(paired_arr.sum(), 1.0);

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
    let bam = complex_bam_fixture()?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 40, 80, "groupA"),
            ("chr1", 180, 220, "groupA"),
            ("chr2", 20, 60, "groupB"),
            ("chr2", 60, 100, "groupB"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2"]),
        bed_path.clone(),
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![20, 60, 120]);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    run(&cfg)?;

    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    assert!(counts_path.exists());
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[2, 2, 40]); // groups, length bins, window size
    assert!(arr.sum() > 0.0);

    let map_path = temp.path().join("sites.group_index.tsv");
    let map_text = std::fs::read_to_string(&map_path)?;
    assert!(map_text.contains("groupA"));
    assert!(map_text.contains("groupB"));

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
    cfg.set_tile_size(40);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
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
        ("groupB", vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        ("groupC", vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        ("groupA", vec![0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
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
fn real_ref_gc_bias_then_gc_bias_package_is_neutral_in_single_bin_case_for_midpoints()
-> Result<()> {
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
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);
    assert_eq!(arr.sum(), 1.0);
    assert_eq!(arr.slice(ndarray::s![0, 0, ..]).to_vec(), vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    assert_eq!(group_to_idx, HashMap::from([("groupA".to_string(), 0usize)]));

    Ok(())
}

#[test]
fn real_ref_gc_bias_then_gc_bias_package_changes_midpoints_in_expected_direction() -> Result<()> {
    // Arrange:
    // Build a real non-neutral GC package, then consume it through `midpoints`.
    //
    // Reference genome:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    //
    // Real producer setup:
    // - fragment length is fixed at 61, so midpoint placement is deterministic
    // - valid starts are 0..=139, because 200 - 61 + 1 = 140
    // - reference BED windows keep only pure-start intervals:
    //     [0,40)   -> starts 0..=39   -> GC%=0
    //     [100,140) -> starts 100..=139 -> GC%=100
    // - starts 40..=99, which would cross the A/C boundary and create intermediate GC%, are
    //   intentionally excluded
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
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let producer_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
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
        vec![("chr1".to_string(), 200)],
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
        "chr1\t0\t40\nchr1\t100\t140\n",
        // Chromosome length 200 and fragment length 61 give:
        //   200 - 61 + 1 = 140 valid starts.
        140,
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
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[2, 1, 11]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    let expected = [("groupA", 5.0_f32), ("groupC", (5.0_f32 / 9.0_f32))];
    for (group_name, expected_weight) in expected {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        for (position, value) in row.iter().enumerate() {
            let expected_value = if position == 5 { expected_weight } else { 0.0 };
            assert!(
                (value - expected_value).abs() <= 1e-6,
                "unexpected midpoint weight for {group_name} at position {position}: expected {expected_value}, got {value}"
            );
        }
    }
    assert!(
        (arr.sum() - (5.0_f32 + 5.0_f32 / 9.0_f32)).abs() <= 1e-6,
        "unexpected total midpoint mass {}",
        arr.sum()
    );

    Ok(())
}

#[test]
fn midpoints_rejects_gc_package_when_length_bins_are_outside_supported_range() -> Result<()> {
    // Arrange:
    // The midpoint command resolves its fragment-length range from the configured bin edges:
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
    let gc_path = temp.path().join("too_short_gc_package.npz");
    write_minimal_gc_package_excluding_length_61(&gc_path)?;

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
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
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
    let gc_path = temp.path().join("gc_pkg_bad_version.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION + 1,
        end_offset: 0,
        length_edges: vec![10, 200],
        gc_edges: vec![0, 101],
        correction_matrix: array![[1.0_f64]],
        length_bin_frequencies: array![1.0_f64],
    };
    package.write_npz(&gc_path)?;

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
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
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
fn coverage_weights_tsv_changes_midpoints_by_full_fragment_average_not_window_overlap()
-> Result<()> {
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
    let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");

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
    midpoints_cfg.set_tile_size(1_000);
    midpoints_cfg.set_min_mapq(0);
    midpoints_cfg.set_require_proper_pair(false);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    midpoints_cfg.set_scale_genome(scale_genome);

    run(&midpoints_cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);
    assert_eq!(arr.slice(ndarray::s![0, 0, ..]).len(), 11);

    let expected_weight = 60.0_f32 / 61.0_f32;
    let row = arr.slice(ndarray::s![0, 0, ..]).to_vec();
    for (position, value) in row.iter().enumerate() {
        let expected = if position == 5 { expected_weight } else { 0.0 };
        assert!(
            (value - expected).abs() <= 1e-6,
            "unexpected midpoint weight at position {position}: expected {expected}, got {value}"
        );
    }
    assert!(
        (arr.sum() - expected_weight).abs() <= 1e-6,
        "expected total midpoint mass {expected_weight}, got {}",
        arr.sum()
    );

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_multi_chromosome_coverage_weights_tsv_is_applied_per_chromosome_in_midpoints()
-> Result<()> {
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
    let producer_bam = bam_from_specs(
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
    let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");

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
    midpoints_cfg.set_tile_size(1_000);
    midpoints_cfg.set_min_mapq(0);
    midpoints_cfg.set_require_proper_pair(false);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    midpoints_cfg.set_scale_genome(scale_genome);
    run(&midpoints_cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[2, 1, 11]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;
    let expected_total = 1.5_f32 + 0.75_f32;
    for (group_name, expected_weight) in [("group_chr1", 1.5_f32), ("group_chr2", 0.75_f32)] {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        for (position, value) in row.iter().enumerate() {
            let expected = if position == 5 { expected_weight } else { 0.0 };
            assert!(
                (value - expected).abs() <= 1e-6,
                "unexpected midpoint weight for {group_name} at position {position}: expected {expected}, got {value}"
            );
        }
    }
    assert!(
        (arr.sum() - expected_total).abs() <= 1e-6,
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
    let tagged_bam = bam_with_gc_tags(&base_bam.bam, "midpoints_gc_tag_paired_avg", &[Some(2.0), Some(4.0)])?;
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
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());
    cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        drop_invalid_gc: false,
    });

    // Act
    run(&cfg)?;

    // Assert
    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);

    let row = arr.slice(ndarray::s![0, 0, ..]).to_vec();
    for (position, value) in row.iter().enumerate() {
        let expected = if position == 5 { 3.0 } else { 0.0 };
        assert!(
            (value - expected).abs() <= 1e-6,
            "unexpected midpoint GC-tag weight at position {position}: expected {expected}, got {value}"
        );
    }
    assert!(
        (arr.sum() - 3.0).abs() <= 1e-6,
        "expected total midpoint mass 3.0, got {}",
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
    let gc_path = temp.path().join("constant_gc_pkg.npz");
    let bed_path = temp.path().join("windows.bed");
    write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![61, 62],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        correction_matrix: array![[3.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut bam_to_bam_cfg = BamToBamConfig::new(
        source_bam.bam.clone(),
        tagged_out_bam.clone(),
        base_chromosomes(&["chr1"]),
    );
    bam_to_bam_cfg.skip_chromosome_sort = true;
    bam_to_bam_cfg.min_mapq = 0;
    bam_to_bam_cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
        gc_file: Some(gc_path.clone()),
        drop_invalid_gc: false,
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
    original_cfg.set_tile_size(1_000);
    original_cfg.set_min_mapq(0);
    original_cfg.set_require_proper_pair(false);
    original_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
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
    tagged_cfg.set_tile_size(1_000);
    tagged_cfg.set_min_mapq(0);
    tagged_cfg.set_require_proper_pair(false);
    tagged_cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        drop_invalid_gc: false,
    });
    run(&tagged_cfg)?;

    // Assert
    let original_arr: Array3<f32> =
        read_npy(&temp.path().join("orig_out/origsites.midpoint_profiles.npy"))?;
    let tagged_arr: Array3<f32> =
        read_npy(&temp.path().join("tagged_out/taggedsites.midpoint_profiles.npy"))?;

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

    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
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
        .join("sites_fetch_reads_all.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
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
