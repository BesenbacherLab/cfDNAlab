# Downstream Compatibility Testing

This note sketches CI-only tests for public output formats. These tests should
verify that files written by `cfDNAlab` can be loaded by the downstream tools
users are likely to use, without making those tools part of the normal local
Rust development environment.

## Why This Belongs In Separate CI

The normal development loop should stay Rust-focused:

```text
cargo check
cargo check --tests
```

R packages, Python scientific stacks, system libraries, and Bioconductor release
timing are a different dependency surface. They are still important because the
public output contract is only useful if real downstream readers can load it.
The right compromise is a separate GitHub Actions workflow that runs on demand,
on pushes to `main`, and on PRs that touch public output code or downstream
fixtures.

## Proposed Layout

```text
downstream_tests/
  zarr_midpoints/
    README.md
    make_fixture.py
    test_python_zarr.py
    test_python_xarray.py
    test_r_cran_zarr.R
    test_r_bioc_rarr.R
```

During the design spike, `make_fixture.py` can create tiny hand-authored Zarr
stores. Once the Rust midpoint writer exists, fixture generation should use the
actual `cfdna midpoints` command with tiny BAM and BED inputs.

## GitHub Actions Workflow

Use a separate workflow, for example:

```text
.github/workflows/downstream-zarr.yml
```

Suggested triggers:

```yaml
on:
  workflow_dispatch:
  push:
    branches:
      - main
  pull_request:
    paths:
      - "src/commands/midpoints/**"
      - "src/shared/**zarr**"
      - "downstream_tests/**"
      - ".github/workflows/downstream-zarr.yml"
```

Run the workflow on pushes to `main` and relevant pull requests. Do not schedule
it with cron. The suite should be small enough that running it on code changes is
the right signal, and scheduled package-drift failures are likely to be noise
unless the output code is being changed.

Keep the workflow optional at first. Once the midpoint Zarr writer is
implemented and the supported downstream package set is stable, it can become
required for PRs that touch the Zarr writer or output schema.

## Midpoint Zarr Fixture

The fixture should be tiny but structurally complete:

```text
counts[group, length_bin, position]
group_name[group]
group_name_utf8[group, max_name_bytes]
group_name_nbytes[group]
eligible_intervals[group]
length_start_bp[length_bin]
length_end_bp[length_bin]
position_bin_start_bp[position]
position_bin_end_bp[position]
```

Use at least:

- two groups
- two length bins
- three position bins
- non-integer count values
- one zero count
- one group name with ordinary ASCII
- one group name that tests the chosen string encoding boundary

The fixture should be small enough that every test can also materialize a
dataframe slice.

## Required Python Checks

Python `zarr` checks:

- Store opens as Zarr V3.
- All expected arrays are present.
- Shapes and data types match the schema.
- Sliced reads return the expected values.
- Compression and chunking do not block reads.

Python `xarray` checks:

- `xr.open_zarr()` succeeds.
- `counts` has dimensions `("group", "length_bin", "position")`.
- Coordinate arrays are visible as coordinates or loadable data variables.
- A small subset can be converted to a dataframe.
- The dataframe contains count, group, length-bin, and position metadata.

The xarray check is essential. A raw Zarr group that opens with `zarr` but fails
with `xarray` is not good enough for the midpoint use case.

## Required R Checks

CRAN `zarr` checks:

- Store opens as Zarr V3.
- Numeric arrays can be read.
- Sliced reads work.
- Group names can be read or reconstructed from the fallback encoding.
- A small R `data.frame` can be built from a slice.

Bioconductor `Rarr` or `ZarrArray` checks:

- Numeric arrays can be read.
- Dimension metadata does not prevent reading.
- The selected compression codec works.
- String arrays work, or the fallback byte encoding can be decoded.

The R checks should avoid assuming that every package exposes the same object
model. The contract is successful loading and reconstruction, not identical R
APIs.

## Zarr V3 Safe Subset To Test First

Start with a deliberately small Zarr V3 subset:

- regular chunk grid only
- no sharding
- no storage transformers
- no extension data types
- no structured arrays
- no datetime or timedelta arrays
- no complex arrays
- numeric arrays limited to `float32`, `int32`, and `uint32`
- avoid public `uint64` unless values can exceed `int32`
- V3 `dimension_names` metadata on every array
- coordinate arrays stored as arrays, not only as attrs

For codecs:

- First fixture: uncompressed or gzip.
- Preferred production candidate: zstd.
- CI matrix should include zstd before the format is finalized.

Zstd is the likely production choice if Python, CRAN `zarr`, and Bioconductor
readers all pass. It is fast and compresses numeric scientific arrays well. The
compatibility test still needs a boring baseline because it separates "the
schema is wrong" from "this reader cannot handle this codec".

## Group Name Encoding Matrix

Test two options before choosing the public schema.

Native string array:

```text
group_name              string[group]
```

Fallback byte encoding:

```text
group_name_utf8         uint8[group, max_name_bytes]
group_name_nbytes       int32[group]
```

Use native strings only if they pass in Python `zarr`, Python `xarray`, CRAN
`zarr`, and Bioconductor. If not, use the byte encoding. It is less pleasant,
but it is explicit and easy to reconstruct everywhere.

## Chunking Checks

The compatibility fixture should not settle performance chunking, but it should
exercise chunked reads. Use chunks that split at least one dimension:

```text
counts chunks: [1, 1, 2]
```

The production midpoint writer can choose larger chunks later, probably biased
toward reading a small number of groups across all length bins and positions.

## Pass Criteria

The downstream test suite passes when:

- Python `zarr` can read all arrays.
- Python `xarray` sees named dimensions on `counts`.
- R can read counts and coordinate arrays without Python.
- R can read or reconstruct group names.
- zstd-compressed fixtures pass in every supported reader.
- A small dataframe can be built in both Python and R.

## Non-Goals

- Do not test every R package that might read Zarr.
- Do not make this part of the normal local `cargo check` workflow.
- Do not benchmark compression or chunking here.
- Do not test large production-sized outputs in this workflow.
- Do not require local developers to install R, Bioconductor, or Python
  scientific stacks just to work on unrelated Rust code.
