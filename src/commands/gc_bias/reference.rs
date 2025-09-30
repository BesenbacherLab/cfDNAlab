use ndarray::{Array1, Array2};

/// What expected/target distribution to correct **to**.
#[derive(Debug, Clone)]
pub enum ExpectedSpec {
    /// Uniform per length: each GC bin gets equal expected share of that row’s total.
    UniformPerLength,
    /// Reference across GC only (same for all lengths). Must be length n_gc and non-negative.
    Reference1D(Array1<f32>),
    /// Reference per (length, GC). Must be shape (n_len, n_gc) and non-negative.
    Reference2D(Array2<f32>),
}
