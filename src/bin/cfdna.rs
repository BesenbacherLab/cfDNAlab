use cfdnalab::gc::GCConfig;
use cfdnalab::lengths::LengthsConfig; // ends::EndsConfig
use cfdnalab::refgc::RefGCConfig;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "cfdna", version)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    GC(GCConfig),
    RefGC(RefGCConfig), // Extract reference GC counts
    Lengths(LengthsConfig),
    // Ends(EndsConfig),
}

fn main() {
    // Run selected subcommand and capture its Result (no `?` in main).
    let res: anyhow::Result<()> = match Cli::parse().cmd {
        Cmd::GC(cfg) => cfdnalab::gc::run(cfg),
        Cmd::RefGC(cfg) => cfdnalab::refgc::run(cfg),
        Cmd::Lengths(cfg) => cfdnalab::lengths::run(cfg),
        // Cmd::Ends(cfg) => cfdnalab::ends::run(cfg),
    };

    if let Err(e) = res {
        eprintln!("{:#}", e);
        std::process::exit(1);
    }

    std::process::exit(0);
}
