pub mod error;
#[cfg(feature = "cli")]
pub mod cli_app;
pub mod commands;
pub mod shared;

pub use error::{Error, Result};

// Curate the top-level API:
// pub use lengths::compute_lengths;
// pub use ends::compute_ends;
