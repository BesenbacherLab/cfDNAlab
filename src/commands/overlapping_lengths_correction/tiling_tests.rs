use super::*;
use anyhow::Result;

#[test]
fn segment_adds_raw_depth_and_full_fragment_length() -> Result<()> {
    let core = Interval::new(100, 110)?;
    let mut depth = vec![0; 11];
    let mut lengths = vec![0; 11];
    assert!(add_segment_to_core_deltas(
        Interval::new(98, 105)?,
        core,
        120,
        &mut depth,
        &mut lengths,
    ));
    assert_eq!(depth[0], 1);
    assert_eq!(depth[5], -1);
    assert_eq!(lengths[0], 120);
    assert_eq!(lengths[5], -120);
    Ok(())
}
