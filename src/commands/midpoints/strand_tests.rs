use super::*;

fn test_window() -> Interval<u64> {
    Interval::new(100, 105).expect("test interval should be valid")
}

#[test]
fn forward_and_unstranded_positions_use_genomic_offsets() -> Result<()> {
    let window = test_window();

    // For [100,105), left-to-right genomic offsets are 0, 2, and 4.
    assert_eq!(stranded_window_position(window, 100, Strand::Forward)?, 0);
    assert_eq!(stranded_window_position(window, 102, Strand::Forward)?, 2);
    assert_eq!(stranded_window_position(window, 104, Strand::Forward)?, 4);

    assert_eq!(
        stranded_window_position(window, 100, Strand::Unstranded)?,
        0
    );
    assert_eq!(
        stranded_window_position(window, 102, Strand::Unstranded)?,
        2
    );
    assert_eq!(
        stranded_window_position(window, 104, Strand::Unstranded)?,
        4
    );

    Ok(())
}

#[test]
fn reverse_positions_mirror_the_half_open_window() -> Result<()> {
    let window = test_window();

    // [100,105) has five positions. Reverse orientation maps 100->4, 102->2, 104->0.
    assert_eq!(stranded_window_position(window, 100, Strand::Reverse)?, 4);
    assert_eq!(stranded_window_position(window, 102, Strand::Reverse)?, 2);
    assert_eq!(stranded_window_position(window, 104, Strand::Reverse)?, 0);

    Ok(())
}

#[test]
fn same_window_accumulates_forward_and_reverse_patterns_at_mirrored_ends() -> Result<()> {
    let window = test_window();
    let midpoint_positions = [100_u64, 101_u64];
    let mut forward_counts = [0_u32; 5];
    let mut reverse_counts = [0_u32; 5];

    for midpoint_position in midpoint_positions {
        let forward_position =
            stranded_window_position(window, midpoint_position, Strand::Forward)?;
        let reverse_position =
            stranded_window_position(window, midpoint_position, Strand::Reverse)?;

        forward_counts[forward_position] += 1;
        reverse_counts[reverse_position] += 1;
    }

    // The same two genomic midpoints occupy the left edge for + and the right edge for -.
    assert_eq!(forward_counts, [1, 1, 0, 0, 0]);
    assert_eq!(reverse_counts, [0, 0, 0, 1, 1]);

    Ok(())
}
