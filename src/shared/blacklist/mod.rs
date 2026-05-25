pub(crate) mod apply_mask;
pub(crate) mod load;
pub(crate) mod overlaps;
pub(crate) mod strategy;

// Curate the top-level API:
pub use apply_mask::apply_blacklist_mask_to_seq;
pub use load::load_blacklists;
pub use overlaps::{compute_blacklist_overlap, is_blacklisted};
pub use strategy::BlacklistStrategy;
