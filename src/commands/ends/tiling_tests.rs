use super::*;
use crate::{
    commands::ends::counting::EndCountsByWindow,
    shared::kmers::motifs_file::EncodedMotifKey,
};
use fxhash::FxHashMap;
use tempfile::TempDir;

fn count_record_signature(
    count_records: &[TileWindowEndCounts],
) -> Vec<(u64, Vec<(u64, u64, bool, f64)>)> {
    count_records
        .iter()
        .map(|window| {
            (
                window.original_idx,
                window
                    .entries
                    .iter()
                    .map(|entry| {
                        (
                            entry.inside_code,
                            entry.outside_code,
                            entry.reverse_on_decode,
                            entry.value,
                        )
                    })
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn build_tile_count_records_sorts_windows_and_entries_deterministically() {
    // Arrange: insert windows and motif keys in hash-map order so the output must rely on the
    // explicit sort rather than insertion order.
    let mut counts_by_window: FxHashMap<u64, EndMotifCounts> = FxHashMap::default();
    counts_by_window.insert(
        9,
        EndMotifCounts {
            counts: FxHashMap::from_iter([
                (
                    EncodedMotifKey {
                        inside_code: 5,
                        outside_code: 1,
                        reverse_on_decode: true,
                    },
                    2.0,
                ),
                (
                    EncodedMotifKey {
                        inside_code: 2,
                        outside_code: 9,
                        reverse_on_decode: false,
                    },
                    1.0,
                ),
            ]),
        },
    );
    counts_by_window.insert(
        3,
        EndMotifCounts {
            counts: FxHashMap::from_iter([(
                EncodedMotifKey {
                    inside_code: 4,
                    outside_code: 4,
                    reverse_on_decode: false,
                },
                7.0,
            )]),
        },
    );

    // Act
    let count_records = build_tile_count_records(counts_by_window);

    // Assert
    assert_eq!(
        count_record_signature(&count_records),
        vec![
            (3, vec![(4, 4, false, 7.0)]),
            (9, vec![(2, 9, false, 1.0), (5, 1, true, 2.0)]),
        ]
    );
}

#[test]
fn merge_tile_count_records_merges_counts_by_window_and_key() {
    // Arrange
    let tile_count_records = [
        vec![TileWindowEndCounts {
            original_idx: 7,
            entries: vec![TileEndMotifCountEntry {
                inside_code: 1,
                outside_code: 2,
                reverse_on_decode: false,
                value: 1.5,
            }],
        }],
        vec![TileWindowEndCounts {
            original_idx: 7,
            entries: vec![TileEndMotifCountEntry {
                inside_code: 1,
                outside_code: 2,
                reverse_on_decode: false,
                value: 2.0,
            }],
        }],
    ];
    let mut reduced = EndCountsByWindow::default();

    // Act
    for count_records in tile_count_records {
        merge_tile_count_records(&mut reduced, count_records)
            .expect("count record merge should work");
    }

    // Assert
    let counts = reduced.get(&7).expect("window 7 should be present");
    let key = EncodedMotifKey {
        inside_code: 1,
        outside_code: 2,
        reverse_on_decode: false,
    };
    assert_eq!(counts.counts.get(&key), Some(&3.5));
}

#[test]
fn merge_tile_count_records_merges_multiple_windows_without_cross_talk() {
    // Arrange: each window should merge only with itself, even when count records arrive interleaved.
    let tile_count_records = [
        vec![
            TileWindowEndCounts {
                original_idx: 7,
                entries: vec![TileEndMotifCountEntry {
                    inside_code: 1,
                    outside_code: 2,
                    reverse_on_decode: false,
                    value: 1.5,
                }],
            },
            TileWindowEndCounts {
                original_idx: 9,
                entries: vec![TileEndMotifCountEntry {
                    inside_code: 3,
                    outside_code: 4,
                    reverse_on_decode: true,
                    value: 2.0,
                }],
            },
        ],
        vec![
            TileWindowEndCounts {
                original_idx: 7,
                entries: vec![TileEndMotifCountEntry {
                    inside_code: 1,
                    outside_code: 2,
                    reverse_on_decode: false,
                    value: 0.5,
                }],
            },
            TileWindowEndCounts {
                original_idx: 9,
                entries: vec![TileEndMotifCountEntry {
                    inside_code: 8,
                    outside_code: 1,
                    reverse_on_decode: false,
                    value: 3.0,
                }],
            },
        ],
    ];
    let mut reduced = EndCountsByWindow::default();

    // Act
    for count_records in tile_count_records {
        merge_tile_count_records(&mut reduced, count_records)
            .expect("count record merge should work");
    }

    // Assert
    assert_eq!(reduced.len(), 2);
    let window_7 = reduced.get(&7).expect("window 7 should be present");
    assert_eq!(
        window_7.counts.get(&EncodedMotifKey {
            inside_code: 1,
            outside_code: 2,
            reverse_on_decode: false,
        }),
        Some(&2.0)
    );

    let window_9 = reduced.get(&9).expect("window 9 should be present");
    assert_eq!(
        window_9.counts.get(&EncodedMotifKey {
            inside_code: 3,
            outside_code: 4,
            reverse_on_decode: true,
        }),
        Some(&2.0)
    );
    assert_eq!(
        window_9.counts.get(&EncodedMotifKey {
            inside_code: 8,
            outside_code: 1,
            reverse_on_decode: false,
        }),
        Some(&3.0)
    );
}

#[test]
fn serialize_tile_counts_round_trips_the_sorted_count_records() {
    // Arrange
    let out_dir = TempDir::new().expect("tempdir");
    let path = out_dir.path().join("tile.counts.bin");
    let count_records = vec![
        TileWindowEndCounts {
            original_idx: 3,
            entries: vec![TileEndMotifCountEntry {
                inside_code: 4,
                outside_code: 4,
                reverse_on_decode: false,
                value: 7.0,
            }],
        },
        TileWindowEndCounts {
            original_idx: 9,
            entries: vec![
                TileEndMotifCountEntry {
                    inside_code: 2,
                    outside_code: 9,
                    reverse_on_decode: false,
                    value: 1.0,
                },
                TileEndMotifCountEntry {
                    inside_code: 5,
                    outside_code: 1,
                    reverse_on_decode: true,
                    value: 2.0,
                },
            ],
        },
    ];

    // Act
    serialize_tile_counts(&path, &count_records).expect("serialisation should work");
    let restored = deserialize_tile_counts(&path).expect("deserialisation should work");

    // Assert
    assert_eq!(
        count_record_signature(&restored),
        count_record_signature(&count_records)
    );
}
