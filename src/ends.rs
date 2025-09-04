use std::path::PathBuf;

pub struct Config {
    pub bam: PathBuf,
    pub out: PathBuf,
    pub motifs: Vec<String>,
}

pub fn run(cfg: Config) -> anyhow::Result<()> {
    // ... your ends/cutting-site logic ...
    Ok(())
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
/// Compute fragment end motifs / cleavage profiles.
pub struct EndsConfig {
    #[cfg_attr(feature = "cli", command(flatten))]
    pub common: crate::common::Common,

    #[cfg_attr(feature = "cli", clap(long))]
    pub bam: PathBuf,

    #[cfg_attr(feature = "cli", clap(long))]
    pub out: PathBuf,

    /// Comma-separated motifs (e.g., CCCA,CCCT)
    #[cfg_attr(feature = "cli", clap(long, value_delimiter = ',', num_args=1..))]
    pub motifs: Vec<String>,
}
