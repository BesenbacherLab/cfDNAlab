use super::*;

#[test]
fn outside_kmer_clip_strategy_warning_skips_zero_outside_bases() {
    // Arrange + Act + Assert: without any outside bases requested, the clip strategy does not
    // matter because there is no outside-reference lookup to become ambiguous.
    assert_eq!(
        outside_kmer_clip_strategy_warning(0, ClipStrategy::Aligned),
        None
    );
    assert_eq!(
        outside_kmer_clip_strategy_warning(0, ClipStrategy::RawShiftedBoundary),
        None
    );
}

#[test]
fn outside_kmer_clip_strategy_warning_skips_skip_strategy() {
    // Arrange + Act + Assert: `skip` already discards soft-clipped motifs, so it should not emit
    // the outside-position ambiguity warning.
    assert_eq!(outside_kmer_clip_strategy_warning(2, ClipStrategy::Skip), None);
}

#[test]
fn outside_kmer_clip_strategy_warning_mentions_each_non_skip_strategy() {
    // Arrange + Act: every non-skip strategy can retain soft-clipped motifs, so all of them
    // should emit the same warning with their CLI-visible strategy name.
    let aligned_warning = outside_kmer_clip_strategy_warning(2, ClipStrategy::Aligned)
        .expect("aligned should warn when outside bases are requested");
    let raw_aligned_warning =
        outside_kmer_clip_strategy_warning(2, ClipStrategy::RawAlignedBoundary)
            .expect("raw-aligned-boundary should warn when outside bases are requested");
    let raw_shifted_warning =
        outside_kmer_clip_strategy_warning(2, ClipStrategy::RawShiftedBoundary)
            .expect("raw-shifted-boundary should warn when outside bases are requested");

    // Assert
    assert!(aligned_warning.contains("`--clip-strategy aligned`"));
    assert!(raw_aligned_warning.contains("`--clip-strategy raw-aligned-boundary`"));
    assert!(raw_shifted_warning.contains("`--clip-strategy raw-shifted-boundary`"));
    assert!(aligned_warning.contains("more noise than signal"));
    assert!(raw_aligned_warning.contains("outside motif actually lies on the reference"));
    assert!(raw_shifted_warning.contains("Prefer `--clip-strategy skip`"));
}

#[test]
fn record_counted_fragment_stats_counts_each_fragment_once_and_each_end_once() {
    // Arrange: the same fragment was counted in multiple windows, but only its left and right
    // end motifs should contribute once each to the summary statistics.
    //
    // Mental derivation:
    // - `counted_fragments` is a per-fragment statistic, so any non-empty motif contribution
    //   should increment it by exactly 1
    // - `counted_motifs` is a per-end statistic, so two counted ends should contribute 2
    let mut counter = EndsCounters::default();
    let counted_end_flags = CountedEndFlags {
        left_counted: true,
        right_counted: true,
    };

    // Act
    record_counted_fragment_stats(&mut counter, counted_end_flags);

    // Assert
    assert_eq!(counter.base.counted_fragments, 1);
    assert_eq!(counter.counted_motifs, 2);
}

#[test]
fn record_counted_fragment_stats_skips_fragments_without_any_counted_motif() {
    // Arrange: if neither end produced a motif, then the fragment reached the counting stage
    // but contributed nothing to the output. Both public statistics should therefore stay at 0.
    let mut counter = EndsCounters::default();

    // Act
    record_counted_fragment_stats(&mut counter, CountedEndFlags::default());

    // Assert
    assert_eq!(counter.base.counted_fragments, 0);
    assert_eq!(counter.counted_motifs, 0);
}
