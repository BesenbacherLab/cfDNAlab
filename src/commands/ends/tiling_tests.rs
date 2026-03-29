use super::*;
use crate::commands::ends::counting::{EncodedEndMotifKey, EndCountsByWindow};

#[test]
fn merge_tile_payload_merges_counts_by_window_and_key() {
    // Arrange
    let tile_payloads = [
        vec![TileWindowEndCounts {
            original_idx: 7,
            entries: vec![TileEndMotifCountEntry {
                within_code: 1,
                outside_code: 2,
                reverse_on_decode: false,
                value: 1.5,
            }],
        }],
        vec![TileWindowEndCounts {
            original_idx: 7,
            entries: vec![TileEndMotifCountEntry {
                within_code: 1,
                outside_code: 2,
                reverse_on_decode: false,
                value: 2.0,
            }],
        }],
    ];
    let mut reduced = EndCountsByWindow::default();

    // Act
    for payload in tile_payloads {
        merge_tile_payload(&mut reduced, payload).expect("payload merge should work");
    }

    // Assert
    let counts = reduced.get(&7).expect("window 7 should be present");
    let key = EncodedEndMotifKey {
        within_code: 1,
        outside_code: 2,
        reverse_on_decode: false,
    };
    assert_eq!(counts.counts.get(&key), Some(&3.5));
}
