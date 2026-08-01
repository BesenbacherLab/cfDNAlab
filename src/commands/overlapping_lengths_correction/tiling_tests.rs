use super::*;

#[test]
fn segment_adds_raw_depth_full_fragment_length_and_signal_separately() -> Result<()> {
    let core = Interval::new(100, 110)?;
    let mut depth = vec![0; 11];
    let mut lengths = vec![0; 11];
    let mut signal = vec![0.0; 11];
    assert!(add_segment_to_core_deltas(
        Interval::new(98, 105)?,
        core,
        120,
        Some(0.5),
        &mut depth,
        &mut lengths,
        Some(&mut signal),
    ));
    assert_eq!(depth[0], 1);
    assert_eq!(depth[5], -1);
    assert_eq!(lengths[0], 120);
    assert_eq!(lengths[5], -120);
    assert_eq!(signal[0], 0.5);
    assert_eq!(signal[5], -0.5);
    Ok(())
}
