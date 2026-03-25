use super::*;
use crate::shared::tiled_run::Tile;

fn make_template() -> GCCounts {
    GCCounts::new(1, 1, 0, (0, 0)).expect("failed to build template")
}

#[test]
fn prepares_fixed_size_streaming_buffers_for_last_partial_window() -> Result<()> {
    // Human verification status: unverified
    // Chromosome length leaves a final partial 100 kb window
    // [114300000,114364328). The tile core lies in that last window, so there is no
    // valid window after it.
    let template = make_template();
    let tile = Tile::from_coords(
        "chr1".to_string(),
        0,
        11,
        114_300_000,
        114_364_328,
        114_299_000,
        114_364_328,
    )
    .expect("test tile should be valid");

    let prepared = prepare_fixed_size_streaming_buffers(
        100_000,
        114_364_328,
        tile.core.try_to_u64()?,
        &template,
    );

    assert!(
        prepared.is_ok(),
        "last partial window should not require an out-of-bounds next buffer: {prepared:?}"
    );
    let (current, next) = prepared?;
    assert_eq!(current.start(), 114_300_000);
    assert_eq!(current.end(), 114_364_328);
    assert!(
        next.is_none(),
        "the last partial window should not invent a following window"
    );
    Ok(())
}

#[test]
fn advances_fixed_size_streaming_buffers_into_last_partial_window() -> Result<()> {
    // Human verification status: unverified
    // Current window is the second-to-last 100 kb bin, next is the last partial bin,
    // and advancing once more must not try to build [114400000,114364328).
    let template = make_template();
    let chrom_len = 114_364_328_u64;
    let window_bp = 100_000_u64;
    let core_interval = Interval::new(110_000_000_u64, chrom_len)?;

    let current = window_state_from_idx(1142, window_bp, chrom_len, core_interval, &template)?;
    let next = window_state_from_idx(1143, window_bp, chrom_len, core_interval, &template)?;

    let advanced = advance_fixed_size_streaming_buffers(
        current,
        next,
        window_bp,
        chrom_len,
        core_interval,
        &template,
    );

    assert!(
        advanced.is_ok(),
        "advancing into the last partial window should not construct an invalid interval: {advanced:?}"
    );
    let (current, next) = advanced?;
    assert_eq!(current.start(), 114_300_000);
    assert_eq!(current.end(), 114_364_328);
    assert!(
        next.is_none(),
        "advancing past the last real window should leave no next window"
    );
    Ok(())
}

#[test]
fn rejects_fixed_size_window_index_past_chromosome_end() -> Result<()> {
    // Human verification status: unverified
    let err = fixed_size_window_interval(1144, 100_000, 114_364_328)
        .expect_err("out-of-range fixed window index should fail");

    assert!(
        format!("{err}").contains("beyond chromosome length"),
        "unexpected error message: {err}"
    );
    Ok(())
}
