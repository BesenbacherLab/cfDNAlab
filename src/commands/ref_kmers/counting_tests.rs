use super::*;
use crate::{
    commands::cli_common::DistributionWindowSpec,
    shared::{
        interval::IndexedInterval,
        kmers::kmer_codec::{build_kmer_specs, build_left_aligned_codes_for_spec, KmerSpec},
        windowing::DistributionWindowContext,
    },
};
use std::path::PathBuf;

struct CountFixture {
    kmer_size: u8,
    spec: KmerSpec,
    codes: crate::shared::kmers::kmer_codec::KmerCodes,
    none: u64,
    n: u64,
}

fn build_count_fixture(sequence: &[u8], kmer_size: u8) -> Result<CountFixture> {
    let mut specs = build_kmer_specs(&[kmer_size])?;
    let spec = specs
        .remove(&kmer_size)
        .expect("spec exists for requested k-mer size");
    let codes = build_left_aligned_codes_for_spec(sequence, &spec);
    let none = spec.sentinel_none();
    let n = spec.sentinel_n();
    Ok(CountFixture {
        kmer_size,
        spec,
        codes,
        none,
        n,
    })
}

fn run_count(
    fixture: &CountFixture,
    windows: &[IndexedInterval<u64>],
    window_spec: &DistributionWindowSpec,
    chr_idx_offset: u64,
    chrom_len: u64,
    assign_by: WindowAssigner,
) -> Result<KmerCountsByWindow> {
    let enc = Enc {
        k: fixture.kmer_size,
        codes: &fixture.codes,
        none: fixture.none,
        n: fixture.n,
    };
    let window_context = DistributionWindowContext {
        spec: window_spec,
        windows: (!windows.is_empty()).then_some(windows),
        chr_idx_offset,
    };
    let mut counts_by_window = KmerCountsByWindow::default();
    let mut selected_counts_by_window = SelectedKmerCountsByWindow::default();
    let mut window_pointer = 0usize;
    count_kmers_by_window(
        &mut counts_by_window,
        &mut selected_counts_by_window,
        &enc,
        &window_context,
        &mut window_pointer,
        0..chrom_len,
        0,
        chrom_len,
        assign_by,
        None,
    )?;
    assert!(
        selected_counts_by_window.is_empty(),
        "unselected count path should not fill selected counts"
    );
    Ok(counts_by_window)
}

fn fixture_kmer(fixture: &CountFixture, motif: &[u8]) -> Kmer {
    forward_kmer(fixture.kmer_size, fixture.spec.encode_kmer_bytes(motif))
}

fn row_counts_by_motif(
    counts: &KmerCountsByWindow,
    row_idx: u64,
    fixture: &CountFixture,
) -> FxHashMap<String, f64> {
    counts
        .get(&row_idx)
        .map(|row| {
            row.counts
                .iter()
                .map(|(kmer, weight)| (fixture.spec.decode_kmer(kmer.code), *weight))
                .collect()
        })
        .unwrap_or_default()
}

fn observed_weight(counts: &KmerCountsByWindow, row_idx: u64, kmer: Kmer) -> f64 {
    counts
        .get(&row_idx)
        .and_then(|row| row.counts.get(&kmer).copied())
        .unwrap_or_default()
}

fn total_row_weight(counts: &KmerCountsByWindow, row_idx: u64) -> f64 {
    counts
        .get(&row_idx)
        .map(|row| row.counts.values().sum())
        .unwrap_or_default()
}

fn assert_close(observed: f64, expected: f64) {
    assert!(
        (observed - expected).abs() < 1e-12,
        "observed {observed}, expected {expected}"
    );
}

#[test]
fn global_any_counts_each_valid_dinucleotide_start() -> Result<()> {
    // Derivation for ACGTAC with k = 2:
    //   0 AC, 1 CG, 2 GT, 3 TA, 4 AC.
    // Global `any` counts each valid k-mer start once into row 0.
    let fixture = build_count_fixture(b"ACGTAC", 2)?;
    let window_spec = DistributionWindowSpec::Global;

    let counts = run_count(&fixture, &[], &window_spec, 0, 6, WindowAssigner::Any)?;

    let row = row_counts_by_motif(&counts, 0, &fixture);
    assert_eq!(row.len(), 4);
    assert_close(row["AC"], 2.0);
    assert_close(row["CG"], 1.0);
    assert_close(row["GT"], 1.0);
    assert_close(row["TA"], 1.0);
    assert_close(total_row_weight(&counts, 0), 5.0);
    Ok(())
}

#[test]
fn count_overlap_adds_fractional_kmer_mass_at_both_window_edges() -> Result<()> {
    // Sequence starts for k = 4 in ACGTNNNNACGT:
    //   0: ACGT -> [0, 4)
    //   8: ACGT -> [8, 12)
    // Starts 1..7 all overlap N and are removed.
    //
    // Window row 7 is [3, 9). The first ACGT overlaps at base 3, so it contributes 1/4.
    // The second ACGT overlaps at base 8, so it contributes 1/4.
    let fixture = build_count_fixture(b"ACGTNNNNACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let windows = IndexedInterval::from_tuples(&[(3, 9, 7_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        12,
        WindowAssigner::CountOverlap,
    )?;

    assert_close(observed_weight(&counts, 7, kmer), 0.50);
    assert_close(total_row_weight(&counts, 7), 0.50);
    assert_eq!(counts.len(), 1);
    Ok(())
}

#[test]
fn count_overlap_accumulates_manual_bed_weights_across_multiple_windows() -> Result<()> {
    // AAAAAAAAAAAAAA has eleven 4-mer starts, all encoding AAAA.
    //
    // Window row 10, [2, 6):
    //   starts 0..5 contribute 2/4 + 3/4 + 4/4 + 3/4 + 2/4 + 1/4 = 3.75.
    // Window row 11, [6, 9):
    //   starts 3..8 contribute 1/4 + 2/4 + 3/4 + 3/4 + 2/4 + 1/4 = 3.00.
    // Window row 12, [8, 12):
    //   starts 5..10 contribute 1/4 + 2/4 + 3/4 + 4/4 + 3/4 + 2/4 = 3.75.
    let fixture = build_count_fixture(b"AAAAAAAAAAAAAA", 4)?;
    let kmer = fixture_kmer(&fixture, b"AAAA");
    let windows =
        IndexedInterval::from_tuples(&[(2, 6, 10_u64), (6, 9, 11_u64), (8, 12, 12_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        14,
        WindowAssigner::CountOverlap,
    )?;

    assert_close(observed_weight(&counts, 10, kmer), 3.75);
    assert_close(observed_weight(&counts, 11, kmer), 3.00);
    assert_close(observed_weight(&counts, 12, kmer), 3.75);
    assert_close(total_row_weight(&counts, 10), 3.75);
    assert_close(total_row_weight(&counts, 11), 3.00);
    assert_close(total_row_weight(&counts, 12), 3.75);
    assert_eq!(counts.len(), 3);
    Ok(())
}

#[test]
fn any_counts_one_per_touched_window() -> Result<()> {
    // Same sequence/window as the count-overlap test. Both surviving ACGT k-mers touch row 7, and `any`
    // converts each qualifying overlap to count 1.0.
    let fixture = build_count_fixture(b"ACGTNNNNACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let windows = IndexedInterval::from_tuples(&[(3, 9, 7_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        12,
        WindowAssigner::Any,
    )?;

    assert_close(observed_weight(&counts, 7, kmer), 2.0);
    assert_close(total_row_weight(&counts, 7), 2.0);
    assert_eq!(counts.len(), 1);
    Ok(())
}

#[test]
fn all_requires_the_full_kmer_span_inside_the_window() -> Result<()> {
    // Window row 7 is [0, 4). Only the first ACGT k-mer is fully inside it. The second ACGT
    // starts at 4, touches the boundary, and has zero overlap with this half-open window.
    let fixture = build_count_fixture(b"ACGTNNNNACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let windows = IndexedInterval::from_tuples(&[(0, 4, 7_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        12,
        WindowAssigner::All,
    )?;

    assert_close(observed_weight(&counts, 7, kmer), 1.0);
    assert_close(total_row_weight(&counts, 7), 1.0);
    assert_eq!(counts.len(), 1);
    Ok(())
}

#[test]
fn proportion_threshold_uses_kmer_base_fraction() -> Result<()> {
    // Window row 7 is [2, 10). Starts 1..7 are removed by N-containing k-mers.
    //   ACGT [0, 4) overlaps by 2/4 and passes proportion=0.5.
    //   ACGT [8, 12) overlaps by 2/4 and passes proportion=0.5.
    // Proportion assignment counts each accepted k-mer as 1.0, not by overlap fraction.
    let fixture = build_count_fixture(b"ACGTNNNNACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let windows = IndexedInterval::from_tuples(&[(2, 10, 7_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        12,
        WindowAssigner::Proportion(0.5),
    )?;

    assert_close(observed_weight(&counts, 7, kmer), 2.0);
    assert_close(total_row_weight(&counts, 7), 2.0);
    assert_eq!(counts.len(), 1);
    Ok(())
}

#[test]
fn proportion_threshold_counts_only_kmers_meeting_query_fraction() -> Result<()> {
    // AAAAAAAA has five 4-mer starts. Window row 7 is [2, 6).
    //   start 0 overlaps by 2/4 and fails threshold 0.51.
    //   starts 1, 2, and 3 overlap by 3/4, 4/4, and 3/4, so they pass.
    //   start 4 overlaps by 2/4 and fails.
    // Proportion assignment counts each passing k-mer as 1.0.
    let fixture = build_count_fixture(b"AAAAAAAA", 4)?;
    let kmer = fixture_kmer(&fixture, b"AAAA");
    let windows = IndexedInterval::from_tuples(&[(2, 6, 7_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        8,
        WindowAssigner::Proportion(0.51),
    )?;

    assert_close(observed_weight(&counts, 7, kmer), 3.0);
    assert_close(total_row_weight(&counts, 7), 3.0);
    assert_eq!(counts.len(), 1);
    Ok(())
}

#[test]
fn midpoint_requires_odd_kmer_size_and_counts_center_base_window() -> Result<()> {
    // For k = 3, sequence ACGTAC has one ACG k-mer at [0, 3). Its center base is coordinate 1.
    // Window row 11 is [1, 2), so midpoint assignment counts that k-mer once.
    let fixture = build_count_fixture(b"ACGTAC", 3)?;
    let mut specs = build_kmer_specs(&[3])?;
    let spec = specs.remove(&3).expect("spec exists for k = 3");
    let kmer = forward_kmer(3, spec.encode_kmer_bytes(b"ACG"));
    let enc = Enc {
        k: 3,
        codes: &fixture.codes,
        none: fixture.none,
        n: fixture.n,
    };
    let windows = IndexedInterval::from_tuples(&[(1, 2, 11_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));
    let window_context = DistributionWindowContext {
        spec: &window_spec,
        windows: Some(windows.as_slice()),
        chr_idx_offset: 0,
    };
    let mut counts_by_window = KmerCountsByWindow::default();
    let mut selected_counts_by_window = SelectedKmerCountsByWindow::default();
    let mut window_pointer = 0usize;

    count_kmers_by_window(
        &mut counts_by_window,
        &mut selected_counts_by_window,
        &enc,
        &window_context,
        &mut window_pointer,
        0..6,
        0,
        6,
        WindowAssigner::Midpoint,
        None,
    )?;

    assert_close(observed_weight(&counts_by_window, 11, kmer), 1.0);
    assert_eq!(counts_by_window.len(), 1);
    Ok(())
}

#[test]
fn midpoint_rejects_even_kmer_size() -> Result<()> {
    // k = 4 has no single center base. The counting path must reject midpoint mode rather than
    // silently picking one of the two middle bases.
    let fixture = build_count_fixture(b"ACGT", 4)?;
    let windows = IndexedInterval::from_tuples(&[(1, 3, 11_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));
    let mut counts_by_window = KmerCountsByWindow::default();
    let mut selected_counts_by_window = SelectedKmerCountsByWindow::default();
    let mut window_pointer = 0usize;
    let enc = Enc {
        k: 4,
        codes: &fixture.codes,
        none: fixture.none,
        n: fixture.n,
    };
    let window_context = DistributionWindowContext {
        spec: &window_spec,
        windows: Some(windows.as_slice()),
        chr_idx_offset: 0,
    };

    let error = count_kmers_by_window(
        &mut counts_by_window,
        &mut selected_counts_by_window,
        &enc,
        &window_context,
        &mut window_pointer,
        0..4,
        0,
        4,
        WindowAssigner::Midpoint,
        None,
    )
    .expect_err("even k-mer midpoint assignment should fail");

    assert!(
        error.to_string().contains("requires an odd `--kmer-size`"),
        "unexpected error: {error:#}"
    );
    assert!(counts_by_window.is_empty());
    assert!(selected_counts_by_window.is_empty());
    Ok(())
}

#[test]
fn grouped_bed_aggregates_by_stored_group_index() -> Result<()> {
    // Grouped BED windows keep group_idx in IndexedInterval.idx().
    //
    // Group 20:
    //   [0, 4) fully overlaps the first ACGT -> 1.0.
    //   [8, 12) fully overlaps the second ACGT -> 1.0.
    // Group 21:
    //   [3, 9) overlaps both ACGT k-mers by 1/4 each -> 0.5 under count-overlap.
    let fixture = build_count_fixture(b"ACGTNNNNACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let windows = IndexedInterval::from_tuples(&[(0, 4, 20_u64), (3, 9, 21_u64), (8, 12, 20_u64)])?;
    let window_spec = DistributionWindowSpec::GroupedBed(PathBuf::from("groups.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        12,
        WindowAssigner::CountOverlap,
    )?;

    assert_close(observed_weight(&counts, 20, kmer), 2.0);
    assert_close(observed_weight(&counts, 21, kmer), 0.50);
    assert_close(total_row_weight(&counts, 20), 2.0);
    assert_close(total_row_weight(&counts, 21), 0.50);
    assert_eq!(counts.len(), 2);
    Ok(())
}

#[test]
fn fixed_size_windows_apply_chromosome_row_offset() -> Result<()> {
    // Fixed-size windows of width 4 on chr-local coordinates:
    //   local row 0: [0, 4) gets the first ACGT.
    //   local row 1: [4, 8) gets the second ACGT.
    // With chr_idx_offset = 10, output rows are 10 and 11.
    let fixture = build_count_fixture(b"ACGTACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let window_spec = DistributionWindowSpec::Size(4);

    let counts = run_count(&fixture, &[], &window_spec, 10, 8, WindowAssigner::All)?;

    assert_close(observed_weight(&counts, 10, kmer), 1.0);
    assert_close(observed_weight(&counts, 11, kmer), 1.0);
    assert_eq!(counts.len(), 2);
    Ok(())
}

#[test]
fn n_containing_and_tail_sentinel_kmers_are_not_counted() -> Result<()> {
    // For ACGNACGT and k = 4:
    //   starts 0..3 contain N and use the N sentinel.
    //   start 4 is ACGT and counts once.
    //   starts 5..7 do not have a full k-mer and use the none sentinel.
    let fixture = build_count_fixture(b"ACGNACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let window_spec = DistributionWindowSpec::Global;

    let counts = run_count(&fixture, &[], &window_spec, 0, 8, WindowAssigner::Any)?;

    assert_close(observed_weight(&counts, 0, kmer), 1.0);
    assert_close(total_row_weight(&counts, 0), 1.0);
    assert_eq!(counts.len(), 1);
    Ok(())
}

#[test]
fn all_counts_nothing_when_window_is_shorter_than_kmer() -> Result<()> {
    // The only full k-mer is ACGTAA at [0, 6). Window row 9 is [0, 4), so `all` rejects it
    // because the complete k-mer span is not inside the window.
    let fixture = build_count_fixture(b"ACGTAA", 6)?;
    let windows = IndexedInterval::from_tuples(&[(0, 4, 9_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        6,
        WindowAssigner::All,
    )?;

    assert!(counts.is_empty());
    Ok(())
}

#[test]
fn chromosome_shorter_than_kmer_has_no_countable_starts() -> Result<()> {
    // No 3-mer can start in a 2 bp chromosome. The packed code vector contains only the
    // no-full-k-mer sentinel, so the global row remains empty.
    let fixture = build_count_fixture(b"AC", 3)?;
    let window_spec = DistributionWindowSpec::Global;

    let counts = run_count(&fixture, &[], &window_spec, 0, 2, WindowAssigner::Any)?;

    assert!(counts.is_empty());
    Ok(())
}

#[test]
fn kmer_crossing_chromosome_end_is_not_clipped_and_counted() -> Result<()> {
    // This fixture intentionally gives start 0 a real ACGT code while the declared chromosome
    // length is only 3 bp. Counting must require the complete [0, 4) k-mer span to fit. If the
    // overlap cursor were allowed to clip that span to [0, 3), the global row would get a count.
    let fixture = build_count_fixture(b"ACGT", 4)?;
    let window_spec = DistributionWindowSpec::Global;

    let counts = run_count(&fixture, &[], &window_spec, 0, 3, WindowAssigner::Any)?;

    assert!(counts.is_empty());
    Ok(())
}

#[test]
fn all_counts_one_kmer_when_window_length_equals_kmer_size() -> Result<()> {
    // ACGT has exactly one 4-mer start. Window row 4 is exactly [0, 4), so `all` accepts that
    // single span and no tail sentinels are counted.
    let fixture = build_count_fixture(b"ACGT", 4)?;
    let kmer = fixture_kmer(&fixture, b"ACGT");
    let windows = IndexedInterval::from_tuples(&[(0, 4, 4_u64)])?;
    let window_spec = DistributionWindowSpec::Bed(PathBuf::from("windows.bed"));

    let counts = run_count(
        &fixture,
        windows.as_slice(),
        &window_spec,
        0,
        4,
        WindowAssigner::All,
    )?;

    assert_close(observed_weight(&counts, 4, kmer), 1.0);
    assert_close(total_row_weight(&counts, 4), 1.0);
    assert_eq!(counts.len(), 1);
    Ok(())
}
