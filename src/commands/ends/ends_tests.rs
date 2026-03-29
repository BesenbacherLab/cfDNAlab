use super::*;

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
