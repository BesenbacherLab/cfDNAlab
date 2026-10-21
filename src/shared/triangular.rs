//! Integer triangular weights shared by scalar smoothing and coverage-histogram fitting.

use std::sync::{Arc, Mutex};

use anyhow::{Context, Result, anyhow};
use fxhash::FxHashMap;

/// Reuse immutable triangular profiles for the different radii needed during a command run.
///
/// For radius 2 the weights are `[1, 2, 3, 2, 1]`. Radius 0 gives `[1]`.
/// A profile is constructed once per radius, including when parallel callers request the same
/// radius together. Returned references share the stored weights without copying their contents.
/// Own a cache per run so unused profiles are released when that run finishes.
///
/// To clip the unavailable left side while keeping the peak on the requested center, start at
/// index `radius - available_left_neighbors`, then use the number of available positions.
/// Weights are unnormalized. Each consumer normalizes using its actual contributing support.
#[derive(Default)]
pub(crate) struct TriangularWeightCache {
    profiles: Mutex<FxHashMap<usize, Arc<[usize]>>>,
}

impl TriangularWeightCache {
    /// Get a cached profile, reporting size or allocation failures before storing it.
    pub(crate) fn get(&self, radius: usize) -> Result<Arc<[usize]>> {
        let mut profiles = self
            .profiles
            .lock()
            .map_err(|_| anyhow!("triangular weight cache lock was poisoned"))?;
        if let Some(weights) = profiles.get(&radius) {
            return Ok(Arc::clone(weights));
        }

        let peak = radius
            .checked_add(1)
            .context("triangular kernel peak exceeds usize")?;
        let length = radius
            .checked_mul(2)
            .and_then(|length| length.checked_add(1))
            .context("triangular kernel length exceeds usize")?;
        let mut weights = Vec::new();
        weights
            .try_reserve_exact(length)
            .context("allocating triangular weights")?;
        weights.extend((1..=peak).chain((1..=radius).rev()));
        let weights: Arc<[usize]> = weights.into();
        // Keep construction under the lock so simultaneous requests cannot build duplicates
        profiles.insert(radius, Arc::clone(&weights));
        Ok(weights)
    }
}

#[cfg(test)]
mod tests {
    include!("triangular_tests.rs");
}
