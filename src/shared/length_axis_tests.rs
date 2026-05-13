use super::{LengthAxis, stepped_range_step};
use serde_json::json;

#[test]
fn stepped_range_step_accepts_equal_bins_and_one_short_final_bin() {
    assert_eq!(stepped_range_step(&[30, 40, 50]), Some(10));

    // The final edge is still explicit in the settings JSON, so a reader can reconstruct
    // [30,40), [40,50), and the shorter final [50,55) bin from start=30, end=55, step=10.
    assert_eq!(stepped_range_step(&[30, 40, 50, 55]), Some(10));
}

#[test]
fn stepped_range_step_rejects_zero_decreasing_and_nonuniform_bins() {
    assert_eq!(stepped_range_step(&[30]), None);
    assert_eq!(stepped_range_step(&[30, 30]), None);
    assert_eq!(stepped_range_step(&[40, 30]), None);

    // Interior bins must match the first width exactly.
    assert_eq!(stepped_range_step(&[30, 40, 55, 60]), None);

    // The final bin may be shorter, but not wider than the first bin.
    assert_eq!(stepped_range_step(&[30, 40, 52]), None);
}

#[test]
fn length_axis_settings_serializes_stepped_range_with_short_final_bin() {
    let length_axis =
        LengthAxis::new(vec![30, 40, 50, 55]).expect("valid length axis should resolve");

    let settings =
        serde_json::to_value(length_axis.settings()).expect("settings should serialize");

    assert_eq!(settings["column_intervals"], json!("half_open"));
    assert_eq!(settings["min_fragment_length"], json!(30));
    assert_eq!(settings["max_fragment_length"], json!(54));
    assert_eq!(settings["n_bins"], json!(3));
    assert_eq!(settings["single_bp_bins"], json!(false));
    assert_eq!(
        settings["bin_definition"],
        json!({"kind": "stepped_range", "start": 30, "end": 55, "step": 10})
    );
}

#[test]
fn length_axis_settings_serializes_nonuniform_edges_explicitly() {
    let length_axis =
        LengthAxis::new(vec![30, 45, 80]).expect("valid length axis should resolve");

    let settings =
        serde_json::to_value(length_axis.settings()).expect("settings should serialize");

    assert_eq!(settings["n_bins"], json!(2));
    assert_eq!(
        settings["bin_definition"],
        json!({"kind": "explicit_edges", "edges": [30, 45, 80]})
    );
}
