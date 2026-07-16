use super::*;
use crate::commands::ref_kmers::counting::SelectedKmerCounts;
use fxhash::FxHashMap;
use tempfile::TempDir;

fn kmer(k: u8, code: u64, orientation: KmerOrientation) -> Kmer {
    Kmer {
        k,
        code,
        orientation,
    }
}

fn count_record_signature(
    count_records: &[TileWindowKmerCounts],
) -> Vec<(u64, Vec<(u8, u64, KmerOrientation, f64)>)> {
    count_records
        .iter()
        .map(|window| {
            (
                window.original_idx,
                window
                    .entries
                    .iter()
                    .map(|entry| (entry.k, entry.code, entry.orientation, entry.value))
                    .collect(),
            )
        })
        .collect()
}

fn selected_record_signature(
    count_records: &[TileWindowSelectedKmerCounts],
) -> Vec<(u64, Vec<(u32, f64)>)> {
    count_records
        .iter()
        .map(|window| {
            (
                window.original_idx,
                window
                    .entries
                    .iter()
                    .map(|entry| (entry.target_idx, entry.value))
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn build_tile_count_records_sorts_windows_and_kmers_deterministically() {
    // Arrange: insert rows and k-mers in a hash map. The output order must come from explicit
    // sorting, not from insertion order.
    let mut counts_by_window = KmerCountsByWindow::default();
    counts_by_window.insert(
        9,
        KmerCounts {
            counts: FxHashMap::from_iter([
                (kmer(4, 5, KmerOrientation::Reverse), 2.0),
                (kmer(3, 9, KmerOrientation::Forward), 1.0),
                (kmer(4, 5, KmerOrientation::Forward), 3.0),
            ]),
        },
    );
    counts_by_window.insert(
        3,
        KmerCounts {
            counts: FxHashMap::from_iter([(kmer(4, 4, KmerOrientation::Forward), 7.0)]),
        },
    );

    // Act
    let count_records = build_tile_count_records(counts_by_window);

    // Assert
    assert_eq!(
        count_record_signature(&count_records),
        vec![
            (3, vec![(4, 4, KmerOrientation::Forward, 7.0)]),
            (
                9,
                vec![
                    (3, 9, KmerOrientation::Forward, 1.0),
                    (4, 5, KmerOrientation::Forward, 3.0),
                    (4, 5, KmerOrientation::Reverse, 2.0),
                ],
            ),
        ]
    );
}

#[test]
fn merge_tile_count_records_sums_by_row_and_kmer() {
    // Arrange: two tile records contribute the same k-mer to row 7. A third entry targets row 9
    // and must not affect row 7.
    let tile_count_records = [
        vec![TileWindowKmerCounts {
            original_idx: 7,
            entries: vec![TileKmerCountEntry {
                k: 4,
                code: 11,
                orientation: KmerOrientation::Forward,
                value: 1.5,
            }],
        }],
        vec![
            TileWindowKmerCounts {
                original_idx: 7,
                entries: vec![TileKmerCountEntry {
                    k: 4,
                    code: 11,
                    orientation: KmerOrientation::Forward,
                    value: 2.0,
                }],
            },
            TileWindowKmerCounts {
                original_idx: 9,
                entries: vec![TileKmerCountEntry {
                    k: 4,
                    code: 12,
                    orientation: KmerOrientation::Forward,
                    value: 4.0,
                }],
            },
        ],
    ];
    let mut reduced = KmerCountsByWindow::default();

    // Act
    for records in tile_count_records {
        merge_tile_count_records(&mut reduced, records).expect("merge should work");
    }

    // Assert
    assert_eq!(reduced.len(), 2);
    assert_eq!(
        reduced
            .get(&7)
            .and_then(|row| row.counts.get(&kmer(4, 11, KmerOrientation::Forward))),
        Some(&3.5)
    );
    assert_eq!(
        reduced
            .get(&9)
            .and_then(|row| row.counts.get(&kmer(4, 12, KmerOrientation::Forward))),
        Some(&4.0)
    );
}

#[test]
fn merge_tile_count_records_rejects_negative_weight() {
    // Arrange: negative counts cannot be represented in the sparse reference background.
    let records = vec![TileWindowKmerCounts {
        original_idx: 7,
        entries: vec![TileKmerCountEntry {
            k: 4,
            code: 11,
            orientation: KmerOrientation::Forward,
            value: -1.0,
        }],
    }];
    let mut reduced = KmerCountsByWindow::default();

    // Act
    let error = merge_tile_count_records(&mut reduced, records)
        .expect_err("negative weights should be rejected");

    // Assert
    assert!(
        error.to_string().contains("negative"),
        "unexpected error: {error:#}"
    );
    assert!(reduced.is_empty());
}

#[test]
fn selected_tile_records_sort_and_merge_by_target_index() {
    // Arrange: selected target rows use target_idx directly, so sorting and merging are independent
    // of concrete k-mer labels.
    let mut counts_by_window = SelectedKmerCountsByWindow::default();
    counts_by_window.insert(
        5,
        SelectedKmerCounts {
            counts: FxHashMap::from_iter([(3, 1.25), (1, 2.0)]),
        },
    );
    counts_by_window.insert(
        2,
        SelectedKmerCounts {
            counts: FxHashMap::from_iter([(4, 0.75)]),
        },
    );

    // Act
    let records = build_selected_tile_count_records(counts_by_window);
    let mut reduced = SelectedKmerCountsByWindow::default();
    merge_selected_tile_count_records(&mut reduced, records.clone()).expect("merge should work");

    // Assert
    assert_eq!(
        selected_record_signature(&records),
        vec![(2, vec![(4, 0.75)]), (5, vec![(1, 2.0), (3, 1.25)])]
    );
    assert_eq!(reduced.len(), 2);
    assert_eq!(
        reduced.get(&2).and_then(|row| row.counts.get(&4)),
        Some(&0.75)
    );
    assert_eq!(
        reduced.get(&5).and_then(|row| row.counts.get(&1)),
        Some(&2.0)
    );
    assert_eq!(
        reduced.get(&5).and_then(|row| row.counts.get(&3)),
        Some(&1.25)
    );
}

#[test]
fn serialize_tile_counts_round_trips_sorted_records() {
    // Arrange
    let temp_dir = TempDir::new().expect("tempdir");
    let path = temp_dir.path().join("tile.ref_kmer_counts.bin");
    let count_records = vec![
        TileWindowKmerCounts {
            original_idx: 3,
            entries: vec![TileKmerCountEntry {
                k: 4,
                code: 4,
                orientation: KmerOrientation::Forward,
                value: 7.0,
            }],
        },
        TileWindowKmerCounts {
            original_idx: 9,
            entries: vec![
                TileKmerCountEntry {
                    k: 3,
                    code: 2,
                    orientation: KmerOrientation::Forward,
                    value: 1.0,
                },
                TileKmerCountEntry {
                    k: 4,
                    code: 5,
                    orientation: KmerOrientation::Reverse,
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

#[test]
fn serialize_selected_tile_counts_round_trips_sorted_records() {
    // Arrange
    let temp_dir = TempDir::new().expect("tempdir");
    let path = temp_dir.path().join("tile.selected_ref_kmer_counts.bin");
    let count_records = vec![
        TileWindowSelectedKmerCounts {
            original_idx: 3,
            entries: vec![TileSelectedKmerCountEntry {
                target_idx: 2,
                value: 7.0,
            }],
        },
        TileWindowSelectedKmerCounts {
            original_idx: 9,
            entries: vec![
                TileSelectedKmerCountEntry {
                    target_idx: 1,
                    value: 1.0,
                },
                TileSelectedKmerCountEntry {
                    target_idx: 4,
                    value: 2.0,
                },
            ],
        },
    ];

    // Act
    serialize_selected_tile_counts(&path, &count_records).expect("serialisation should work");
    let restored = deserialize_selected_tile_counts(&path).expect("deserialisation should work");

    // Assert
    assert_eq!(
        selected_record_signature(&restored),
        selected_record_signature(&count_records)
    );
}
