use cfdnalab::lengths::LengthsConfig; // ends::EndsConfig, gc::GCConfig,
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "cfdna", version)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    // GC(GCConfig),
    // RefGC(RefGCConfig), // Extract reference GC counts
    Lengths(LengthsConfig),
    // Ends(EndsConfig),
}

fn main() {
    // Catch and handle errors
    // Ensures that tempfile has time to remove the tmp dir
    if let Err(e) = match Cli::parse().cmd {
        // Cmd::GC(cfg) => cfdnalab::gc::run(cfg)?,
        Cmd::Lengths(cfg) => cfdnalab::lengths::run(cfg)?,
        // Cmd::Ends(cfg) => cfdnalab::ends::run(cfg)?,
    } {
        eprintln!("{:?}", e);
        std::process::exit(1);
    }
    std::process::exit(0);
}
