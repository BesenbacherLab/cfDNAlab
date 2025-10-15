use cfdnalab::commands::visualize_positions::config::VisualizeSelectedRegionConfig;
use cfdnalab::pos_kmer_viz::{BasesFrom, OverlapResolution, ReferenceFrame, Style};

#[test]
fn build_uses_expected_defaults() {
    let cfg = VisualizeSelectedRegionConfig {
        frame: ReferenceFrame::Left,
        positions: "1..5".to_string(),
        step: 1,
        bases_from: BasesFrom::PreferRead,
        overlap_resolution: OverlapResolution::NearestRead,
        lengths: Some(vec![120]),
        length_range: None,
        style: Style::Ascii,
        width: None,
        height: None,
        output: None,
        label: None,
        show_index: false,
        show_half: false,
        hide_mid: false,
    };

    let viz = cfg.build().expect("config should build");
    assert_eq!(viz.frame, ReferenceFrame::Left);
    assert_eq!(viz.positions_input, "1..5");
    assert_eq!(viz.step.get(), 1);
    assert_eq!(viz.bases, BasesFrom::PreferRead);
    assert_eq!(viz.overlap_resolution, OverlapResolution::NearestRead);
    assert!(viz.show_mid);
    assert_eq!(viz.fragment_lengths, vec![120]);
    assert_eq!(viz.style, Style::Ascii);
    assert_eq!(viz.width, 100);
    assert_eq!(viz.height, 120);
}

#[test]
fn build_applies_overrides() {
    let cfg = VisualizeSelectedRegionConfig {
        frame: ReferenceFrame::Right,
        positions: "..half".to_string(),
        step: 3,
        bases_from: BasesFrom::NearestRead,
        overlap_resolution: OverlapResolution::BaseQuality,
        lengths: Some(vec![90, 120]),
        length_range: None,
        style: Style::Svg,
        width: Some(140),
        height: Some(200),
        output: None,
        label: Some("test".to_string()),
        show_index: true,
        show_half: true,
        hide_mid: true,
    };

    let viz = cfg.build().expect("config should build");
    assert_eq!(viz.frame, ReferenceFrame::Right);
    assert_eq!(viz.step.get(), 3);
    assert_eq!(viz.bases, BasesFrom::NearestRead);
    assert_eq!(viz.overlap_resolution, OverlapResolution::BaseQuality);
    assert!(!viz.show_mid);
    assert_eq!(viz.fragment_lengths, vec![90, 120]);
    assert_eq!(viz.style, Style::Svg);
    assert_eq!(viz.width, 140);
    assert_eq!(viz.height, 200);
    assert_eq!(viz.label.as_deref(), Some("test"));
    assert!(viz.show_index);
    assert!(viz.show_half);
}

#[test]
fn build_rejects_zero_step() {
    let cfg = VisualizeSelectedRegionConfig {
        frame: ReferenceFrame::Left,
        positions: "1..5".to_string(),
        step: 0,
        bases_from: BasesFrom::PreferRead,
        overlap_resolution: OverlapResolution::NearestRead,
        lengths: Some(vec![100]),
        length_range: None,
        style: Style::Ascii,
        width: None,
        height: None,
        output: None,
        label: None,
        show_index: false,
        show_half: false,
        hide_mid: false,
    };

    let err = cfg.build().expect_err("step of zero must error");
    assert!(err.to_string().contains("--step must be at least 1"));
}

#[test]
fn build_rejects_mid_with_nearest_read() {
    let cfg = VisualizeSelectedRegionConfig {
        frame: ReferenceFrame::Mid,
        positions: "-5..5".to_string(),
        step: 1,
        bases_from: BasesFrom::NearestRead,
        overlap_resolution: OverlapResolution::NearestRead,
        lengths: Some(vec![100]),
        length_range: None,
        style: Style::Ascii,
        width: None,
        height: None,
        output: None,
        label: None,
        show_index: false,
        show_half: false,
        hide_mid: false,
    };

    let err = cfg
        .build()
        .expect_err("mid frame with nearest-read should fail");
    assert!(
        err.to_string()
            .contains("`--bases-from nearest-read` is incompatible")
    );
}
