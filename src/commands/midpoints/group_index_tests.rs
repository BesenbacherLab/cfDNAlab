use super::{eligible_interval_counts_by_group, write_midpoint_group_index_tsv};
use crate::shared::{bed::GroupedWindows, interval::IndexedInterval};
use fxhash::FxHashMap;
use tempfile::TempDir;

#[test]
fn eligible_interval_counts_include_zero_count_groups() {
    // groupA has two retained intervals, groupB has none, and groupC has one.
    // This is the metadata users need for mean profiles after interval prefiltering.
    let mut intervals_by_chromosome = FxHashMap::default();
    intervals_by_chromosome.insert(
        "chr1".to_string(),
        GroupedWindows::new(
            vec![
                IndexedInterval::new(10, 20, 0).expect("valid test interval"),
                IndexedInterval::new(30, 40, 0).expect("valid test interval"),
                IndexedInterval::new(50, 60, 2).expect("valid test interval"),
            ],
            None,
        ),
    );
    intervals_by_chromosome.insert(
        "chr2".to_string(),
        GroupedWindows::new(
            Vec::<IndexedInterval<u64>>::new(),
            None,
        ),
    );

    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(0, "groupA".to_string());
    group_idx_to_name.insert(1, "groupB".to_string());
    group_idx_to_name.insert(2, "groupC".to_string());

    let counts = eligible_interval_counts_by_group(&intervals_by_chromosome, &group_idx_to_name);

    assert_eq!(counts.get(&0), Some(&2));
    assert_eq!(counts.get(&1), Some(&0));
    assert_eq!(counts.get(&2), Some(&1));
}

#[test]
fn write_midpoint_group_index_sorts_sanitizes_and_defaults_missing_counts() {
    // The sidecar is parsed by name downstream, so the writer must keep the public column name,
    // sort rows by group index, and keep malformed group names from breaking TSV structure.
    let temp = TempDir::new().expect("temp dir should be created");
    let output_path = temp.path().join("sites.group_index.tsv");
    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(2, "group\tB\nname".to_string());
    group_idx_to_name.insert(0, "groupA".to_string());
    group_idx_to_name.insert(1, "groupWithoutIntervals".to_string());

    let mut eligible_interval_counts = FxHashMap::default();
    eligible_interval_counts.insert(0, 4);
    eligible_interval_counts.insert(2, 1);

    write_midpoint_group_index_tsv(&output_path, &group_idx_to_name, &eligible_interval_counts)
        .expect("group index should write");

    let text = std::fs::read_to_string(output_path).expect("group index should be readable");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "group_idx\tgroup_name\teligible_intervals");
    assert_eq!(lines[1], "0\tgroupA\t4");
    assert_eq!(lines[2], "1\tgroupWithoutIntervals\t0");
    assert_eq!(lines[3], "2\tgroup    B name\t1");
    assert_eq!(lines.len(), 4);
}
