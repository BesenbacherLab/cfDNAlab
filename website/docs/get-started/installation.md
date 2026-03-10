# Installation

## Build from GitHub

```bash
cargo install --git https://github.com/ludvigolsen/cfDNAlab --release --features cli,plotters
```

Verify installation:

```bash
cfdna --help
```

## Build from source

```bash
cargo build --release --features cli,plotters
```

The binary is available at `target/release/cfdna`.

Verify installation:

```bash
./target/release/cfdna --help
```
