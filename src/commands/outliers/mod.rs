pub(crate) mod config;
mod distribution;
mod failure;
mod inspection;
mod model;
pub(crate) mod outliers;
mod outputs;
#[cfg(feature = "plotters")]
mod plotting;
mod preservation;
#[cfg(feature = "plotters")]
mod preserved_plotting;
mod preserved_regions;
mod readme;
mod striding;
mod zinb;
