# Release Checklist

This checklist is for coordinated cfDNAlab releases across the Rust CLI and
future Python/R helper packages.

## Package Boundaries

Keep package roles separate. Use the project/runtime names when talking about
release responsibilities:

```text
Rust CLI        installs the `cfdna` executable
Python helper   imports as `cfdnalab`
R helper        imports as the future R package name
```

Specific distribution names are registry-specific:

```text
crates.io / Bioconda Rust package:  cfdnalab
PyPI Python package:                cfdnalab
future conda Python recipe:         py-cfdnalab
future conda R recipe:              r-cfdnalab
```

The Python and R packages read cfDNAlab output files. They should not wrap,
vendor, or install the Rust CLI. If users later want one install target for the
whole ecosystem, make a separate metapackage such as `cfdnalab-suite`.

## Versioning Policy

Do not synchronize Rust, Python, and R package versions by default. The Rust CLI
will change more often than the downstream loaders because many CLI releases add
commands, algorithms, performance improvements, or docs without changing public
output schemas.

Use package SemVer independently:

```text
Rust CLI:       tracks the command-line tool
Python helper:  tracks the Python reader/helper API
R helper:       tracks the R reader/helper API
```

Coordinate releases around schema compatibility, not matching version numbers.
If Rust `cfdnalab 0.5.0` still writes `midpoint_profiles` schema version 1, the
Python helper may reasonably remain at `0.2.3`. A helper release is needed only
when the helper API changes, a reader bug is fixed, or a new/changed output
schema needs support.

If a Rust release changes a public output schema, the release is not complete
until one of these is true:

- the relevant helper package supports the new schema
- the compatibility table says helper support is pending
- the user-facing docs clearly say the output should be rerun instead of read
  with older helpers

Helper packages must validate and document supported `cfdnalab_schema` and
`cfdnalab_schema_version` values.

### Citation file

Remember to update the citation file version.

## Schema Compatibility Table

Maintain an internal schema compatibility table whenever a public output schema
exists. This table is the source of truth for answering: "I already generated
files; which helper version can read them?"

Suggested location:

```text
.AI/docs/future/schema_compatibility.md
```

When stable enough for users, distill the same information into public docs.

Template:

```text
Schema                  Schema version  Written by CLI versions  Python helper       R helper            Recommendation
midpoint_profiles        1               0.2.0..0.5.x             cfdnalab >=0.2.0    pending             Read normally.
midpoint_profiles        2               0.6.0..                  cfdnalab >=0.3.0    r-cfdnalab >=0.3.0  Prefer new helpers; rerun if comparing with v1 output.
end_motif_counts         1               0.6.0..                  pending            pending             Rerun only if schema changes before helper support lands.
```

The table should record:

- schema name from `cfdnalab_schema`
- integer schema version from `cfdnalab_schema_version`
- CLI versions that wrote the schema
- minimum Python helper version that can read it
- minimum R helper version that can read it
- whether users should read old files, install an older helper, install a newer
  helper, or rerun the command

Do not use Rust/Python/R package version matching as a substitute for this
table. The schema version is the output contract; package versions are only
release vehicles.

## Git Tags

Use component-scoped tags rather than one ambiguous tag namespace:

```text
rust-v0.2.0
py-v0.2.0
r-v0.2.0
```

If all components are released together from the same commit, also add an
optional suite tag:

```text
v0.2.0
```

The suite tag should mean "the coordinated cfDNAlab release state", not "the
Rust crate only". If only the Rust CLI is released, tag only `rust-v...`.

This avoids the long-term confusion where `v0.2.1` might refer to a Rust patch,
a Python reader patch, or an R helper patch.

## Rust CLI Release Checks

Before publishing the Rust CLI:

```bash
cargo check
cargo check --tests
cargo check --all-features
cargo test
cargo package --list
cargo package
cargo publish --dry-run
```

Also verify:

- `Cargo.toml` version is correct.
- `CHANGELOG.md` includes user-facing changes.
- generated CLI docs are current if command help changed.
- website build passes if docs changed.
- `cargo package --list` does not include `examples/gen_cli_docs.rs`, and the
  crate exposes only the public `cfdna` binary.
- package whitelist does not include internal planning docs, website sources, or
  generated local artifacts.

## Python Helper Release Checks

Before publishing the Python helper:

```bash
cd py-cfdnalab
uv lock --check
uv run python -m py_compile src/cfdnalab/__init__.py src/cfdnalab/midpoints.py
uv run pytest
uv build
```

Also verify:

- `pyproject.toml` version is correct.
- `README.md` states that this package does not install the Rust CLI.
- runtime dependencies are minimal and justified.
- the package imports as `cfdnalab`.
- examples work against a cfDNAlab-generated Zarr store.

## R Helper Release Checks

Before publishing an R helper package:

```r
devtools::check()
```

Also verify:

- package docs list supported schema names and versions.
- examples load a cfDNAlab-generated Zarr store.
- dependency choices are explicit, especially around Zarr readers and
  compression support.
- package docs state that the R package reads outputs and does not install the
  Rust CLI.

## Downstream Compatibility Checks

Run downstream compatibility checks before any release that changes public
output files:

```bash
# local or CI-equivalent
cargo test --no-default-features --features cmd_midpoints,cmd_ends \
  --test generate_downstream_zarr_fixtures \
  -- --ignored
```

Then run the downstream Python and R reader checks against the generated stores.

The downstream fixture must be produced by cfDNAlab itself. A hand-authored Zarr
fixture is not enough for release validation because it cannot catch writer
schema drift.

## Documentation Checks

Before a coordinated release:

- update `CHANGELOG`
- update package READMEs
- update user guides with any output-format changes
- update the generated website loader-doc mapping in
  `website/scripts/generate_loader_docs.py` when public Python/R loader symbols,
  output groups, quick examples, or output filenames change
- keep Python and R examples concise and aligned with the helper APIs
- avoid documenting several downstream packages when one recommended path per
  language is enough

Run the website build when docs changed:

```bash
npm --prefix website ci
npm --prefix website run build
```

## Release Order

For output-format releases:

1. Release the Rust CLI after writer tests and downstream fixture generation
   pass.
2. Release Python and R helpers in the same major/minor line.
3. Run downstream compatibility tests against the released or release-candidate
   packages.
4. Tag released components with component-scoped tags.
5. Add an optional suite tag only if the component releases are coordinated from
   the same commit.
6. Update package badges and install instructions after publication is live.

For helper-only patch releases, release and tag only the affected helper.
