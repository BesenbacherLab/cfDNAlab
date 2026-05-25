#[cfg(not(feature = "cli"))]
fn main() {
    eprintln!("This binary requires --features cli");
    std::process::exit(1);
}

#[cfg(feature = "cli")]
fn main() {
    cfdnalab::run_cli();
}
