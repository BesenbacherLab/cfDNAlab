use fxhash::hash64;
use rand::{Rng, SeedableRng, rngs::StdRng};

/// Compute midpoint with a random tie-break for even-length fragments.
///
/// - For odd `length`, returns `start + length/2`.
/// - For even `length`, randomly returns the left or right center:
///   either `start + (length/2 - 1)` or `start + (length/2)`.
#[inline]
pub(crate) fn midpoint_random_even<R: Rng + ?Sized>(start: u32, length: u32, rng: &mut R) -> u32 {
    debug_assert!(length > 0, "Zero-length fragment");
    if length == 0 {
        return start;
    }

    let half = length / 2;
    if (length % 2) == 1 {
        // Odd length -> unique midpoint
        start.saturating_add(half)
    } else {
        // Even length -> pick left or right center through the supplied random generator
        let right = start.saturating_add(half);
        if rng.random_bool(0.5) {
            right.saturating_sub(1)
        } else {
            right
        }
    }
}

/// Compute midpoint with a deterministic random seed.
///
/// This keeps the same random tie-break behavior for even-length fragments, while making the
/// choice reproducible for a given seed.
#[inline]
pub(crate) fn midpoint_random_even_with_seed(start: u32, length: u32, seed: u64) -> u32 {
    let mut rng = StdRng::seed_from_u64(seed);
    midpoint_random_even(start, length, &mut rng)
}

/// Compute midpoint with a deterministic coordinate-derived seed.
///
/// The seed is derived from chromosome, start, and length, so duplicate fragments with the same
/// coordinates choose the same center for even-length fragments. The tie-break still goes through
/// a random generator so the left and right centers remain approximately balanced across many
/// distinct coordinates.
#[inline]
pub(crate) fn midpoint_random_even_for_fragment(chromosome: &str, start: u32, length: u32) -> u32 {
    let seed = fragment_midpoint_seed(chromosome, start, length);
    midpoint_random_even_with_seed(start, length, seed)
}

fn fragment_midpoint_seed(chromosome: &str, start: u32, length: u32) -> u64 {
    hash64(&(chromosome, start, length))
}

#[cfg(test)]
mod tests {
    include!("midpoint_tests.rs");
}
