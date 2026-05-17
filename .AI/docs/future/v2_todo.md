# v0.2 TODO

## Release Packaging

### Hide internal docs generator from crates.io install UI

Status: implemented in the repository by moving the generator to
`examples/gen_cli_docs.rs` and removing the package `[[bin]]` declaration.
Before publishing, recheck the crates.io package/install UI to confirm only
`cfdna` is listed.

The crates.io install panel for v0.1.0 shows both `cfdna` and `gen_cli_docs` as installed binaries.
That is misleading because `gen_cli_docs` is an internal documentation-generation helper, not a user-facing CLI tool.

### Defer Bioconda packaging

Bioconda distribution is useful but not release-day work for v0.1.0. Revisit after the crates.io/GitHub release is stable.
