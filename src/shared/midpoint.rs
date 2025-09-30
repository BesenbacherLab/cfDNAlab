use rand::Rng;

/// Compute midpoint with a random tie-break for even-length fragments.
///
/// - For odd `length`, returns `start + length/2`.
/// - For even `length`, randomly returns the left or right center:
///   either `start + (length/2 - 1)` or `start + (length/2)`.
#[inline]
pub fn midpoint_random_even<R: Rng + ?Sized>(start: u32, length: u32, rng: &mut R) -> u32 {
    debug_assert!(length > 0, "Zero-length fragment");
    if length == 0 {
        return start;
    }

    let half = length / 2;
    if (length % 2) == 1 {
        // Odd length -> unique midpoint
        start.saturating_add(half)
    } else {
        // Even length -> randomly pick left or right center
        let right = start.saturating_add(half);
        if rng.random_bool(0.5) {
            right.saturating_sub(1)
        } else {
            right
        }
    }
}

/// Convenience wrapper using thread-local RNG.
#[inline]
pub fn midpoint_random_even_with_thread_rng(start: u32, length: u32) -> u32 {
    let mut rng = rand::rng();
    midpoint_random_even(start, length, &mut rng)
}
