use super::is_midpoint;
use crate::Result;
use crate::shared::interval::Interval;

fn midpoint_blacklist_hit(
    blacklist_interval: (u64, u64),
    fragment_interval: (u64, u64),
) -> Result<bool> {
    let blacklist_intervals = Interval::from_tuples(&[blacklist_interval])?;
    let fragment_interval = Interval::new(fragment_interval.0, fragment_interval.1)?;
    let mut blacklist_ptr = 0;

    Ok(is_midpoint(
        &blacklist_intervals,
        fragment_interval,
        0,
        &mut blacklist_ptr,
    ))
}

#[test]
fn midpoint_blacklist_checks_both_central_bases_for_even_intervals() -> Result<()> {
    // The even fragment [40, 50) has no single discrete midpoint base.
    // Its two central bases are 44 and 45, so blacklist overlap at either base
    // should remove the fragment.
    let fragment_interval = (40, 50);

    assert!(midpoint_blacklist_hit((44, 45), fragment_interval)?);
    assert!(midpoint_blacklist_hit((45, 46), fragment_interval)?);
    assert!(!midpoint_blacklist_hit((43, 44), fragment_interval)?);
    assert!(!midpoint_blacklist_hit((46, 47), fragment_interval)?);

    Ok(())
}

#[test]
fn midpoint_blacklist_uses_single_central_base_for_odd_intervals() -> Result<()> {
    // The odd fragment [40, 51) has one central base at 45.
    let fragment_interval = (40, 51);

    assert!(midpoint_blacklist_hit((45, 46), fragment_interval)?);
    assert!(!midpoint_blacklist_hit((44, 45), fragment_interval)?);
    assert!(!midpoint_blacklist_hit((46, 47), fragment_interval)?);

    Ok(())
}
