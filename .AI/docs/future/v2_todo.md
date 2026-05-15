# v0.2 TODO

## Release Packaging

### Hide internal docs generator from crates.io install UI

The crates.io install panel for v0.1.0 shows both `cfdna` and `gen_cli_docs` as installed binaries.
That is misleading because `gen_cli_docs` is an internal documentation-generation helper, not a user-facing CLI tool.

Before v0.2, change the docs-generation workflow so `cargo install cfdnalab` advertises and installs only the public `cfdna` binary. Likely options:

- Move `gen_cli_docs` out of `src/bin/` into an `xtask`, `examples`, or another unpublished helper location.
- Keep docs generation available for maintainers without declaring it as an installable package binary.
- Recheck the crates.io package/install UI after publishing to confirm only `cfdna` is listed.

### Defer Bioconda packaging

Bioconda distribution is useful but not release-day work for v0.1.0. Revisit after the crates.io/GitHub release is stable.
