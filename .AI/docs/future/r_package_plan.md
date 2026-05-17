# R cfDNAlab Loader Package Plan

## Goal

Create an R package, `cfdnalab`, that reads cfDNAlab public analysis outputs and returns standard R objects for downstream plotting and analysis.

The package should not wrap or install the Rust CLI. It should be clearly described as a loader and analysis-helper package for files already produced by `cfdna`.

## Design Principles

- Use an R-native interface: plain functions, S3 classes, and standard R objects.
- Avoid Python-style object methods in the public API.
- Return data-frame-compatible objects for tabular metadata and plotting.
- Return `Matrix::sparseMatrix` for sparse end-motif counts.
- Do not silently densify sparse Zarr stores.
- Make large materialization explicit in function names and documentation.
- Keep tidyverse-friendly behavior without requiring the tidyverse as the core dependency.
- Validate schema metadata before exposing objects so malformed stores fail early with package-specific messages.

## Object Model

Use S3 classes, not S4, for the first version.

S3 is enough for the current scope:

- loaders return small list-backed objects with a class vector
- generic functions can dispatch by output type and row mode
- users can inspect objects easily with `str()`
- documentation is simpler than S4 for a small loader package

Suggested class structure:

```r
class(x) <- c("cfdnalab_midpoint_profiles", "cfdnalab_zarr_store")

class(x) <- c(
  "cfdnalab_global_end_motif_counts",
  "cfdnalab_end_motif_counts",
  "cfdnalab_zarr_store"
)

class(x) <- c(
  "cfdnalab_windowed_end_motif_counts",
  "cfdnalab_end_motif_counts",
  "cfdnalab_zarr_store"
)

class(x) <- c(
  "cfdnalab_grouped_end_motif_counts",
  "cfdnalab_end_motif_counts",
  "cfdnalab_zarr_store"
)
```

The loader should choose the most specific end-motif subclass from root metadata:

- `row_mode = "global"` -> `cfdnalab_global_end_motif_counts`
- `row_mode = "size"` or `"bed"` -> `cfdnalab_windowed_end_motif_counts`
- `row_mode = "grouped_bed"` -> `cfdnalab_grouped_end_motif_counts`

## Package Dependencies

Keep the first dependency set small:

- `zarr`: CRAN Zarr v3 reader backend
- `Matrix`: sparse end-motif matrices
- `methods`: only if needed for `Matrix`
- `jsonlite`: JSON root and array attributes if the Zarr backend does not expose them cleanly

Avoid hard dependency on `tibble`, `dplyr`, or `data.table` in the first version.

Return base `data.frame` objects by default. They work naturally with:

```r
tibble::as_tibble(x)
data.table::as.data.table(x)
```

This avoids making tidyverse or data.table users pay for dependencies they do not use.

Current candidate Zarr readers:

- CRAN `zarr`: native R implementation with Zarr v3 support.
- Bioconductor `Rarr`: promising for `DelayedArray`, but `Rarr` 1.6.0 expects
  Zarr v2 `.zarray` metadata and cannot read the current Zarr v3 stores.
- CRAN/R-universe `pizzarr`: useful to test, but do not assume full compatibility until downstream tests cover our stores.

First implementation decision:

- Start with CRAN `zarr` as the internal backend because it has native Zarr v3 support and keeps installation simpler for ordinary R users.
- Keep all backend access in `R/zarr_helpers.R` so we can swap or add `Rarr` later if downstream tests reveal missing support.

Decision rule:

- Use one backend internally if it can read all required numeric arrays, attrs, dimensions, and chunks reliably.
- If no single backend is enough, keep the public API stable and isolate backend-specific code in a small internal layer.

Sources checked while drafting:

- [`zarr`](https://cran.r-universe.dev/zarr): native R Zarr implementation with required Zarr v3 support.
- [`Rarr`](https://bioconductor.org/packages/release/bioc/vignettes/Rarr/inst/doc/Rarr.html): Bioconductor reader with `DelayedArray` integration, but not currently compatible with the V3 stores written here.
- [`pizzarr`](https://cran.r-universe.dev/articles/pizzarr/v3-read.html): Zarr v3 reader support through the Zarr developers R-universe.

## Public Loaders

```r
read_midpoints(path)
read_end_motifs(path)
```

Both loaders should:

- require `path` to exist
- require `path` to be a directory
- require `path` to end in `.zarr`
- read root `zarr.json`
- validate `cfdnalab_schema`
- validate `cfdnalab_schema_version`
- give package-specific errors for wrong schema, missing arrays, or unsupported storage modes

## Shared Generics

Use generic functions where names make sense across outputs:

```r
schema_version(x)
storage_mode(x)
groups(x)
length_bins(x)
positions(x)
motifs(x)
```

Only implement methods where the concept exists. For example:

- `groups()` exists for midpoint profiles and grouped end motifs.
- `positions()` exists for midpoint profiles.
- `windows()` exists only for windowed end motifs.
- `motifs()` exists only for end motifs.

Do not invent generalized `rows()` methods for users. The schema has rows internally, but users think in groups, windows, motifs, and profiles.

## Indexing Policy

R-facing selectors and metadata must be one-based. Users should not need to
know that the Zarr schema stores zero-based axis coordinates, and the package
must not expose public zero-based selectors where an off-by-one mistake can
silently return the wrong data.

Rules:

- Public function arguments and return columns use `*_idx` names and one-based
  values.
- Internal zero-based values use `*_idx0` or `*_index0` names.
- Loaded object internals may keep zero-based arrays when they mirror the Zarr
  schema, but those fields must have a zero-based suffix.
- Public data frames should not include `*_idx0` or `*_index0` columns.
- Schema/debug helpers may expose raw schema coordinates only if the function
  name and column names make the low-level nature explicit.
- Direct `+ 1L` or `- 1L` conversions should be isolated in helper functions,
  not scattered across output-specific methods.

Conversion helpers:

```r
cf_validate_r_index(index, size, name)
cf_r_index_to_index0(index)
cf_index0_to_r_index(index0)
```

The output-specific code should validate public one-based input, convert once
to `index0`, then use `index0 + 1L` only at the final R array/data-frame
subsetting boundary through shared helpers. Tests must use asymmetric fixtures
where choosing index `1` versus zero-based index `1` returns different values.

## Midpoint Profiles API

Loader:

```r
midpoints <- read_midpoints("sample.midpoint_profiles.zarr")
```

Metadata:

```r
groups(midpoints)
length_bins(midpoints)
positions(midpoints)
group_idx(midpoints, group_name)
length_bin_idx(midpoints, length)
```

Expected metadata data frames:

```r
groups(midpoints)
# group_idx, group_name, eligible_intervals

length_bins(midpoints)
# length_bin_idx, length_start_bp, length_end_bp

positions(midpoints)
# position_idx, position_bin_start_bp, position_bin_end_bp
```

Profile extraction:

```r
profile_data_frame(midpoints, group, length_bin_idx)
profile_data_frame(midpoints, group_idx, length_bin_idx)
profile_array(midpoints, group_idx, length_bin_idx)
```

`profile_data_frame()` should return one row per position bin:

```r
position_idx
position_bin_start_bp
position_bin_end_bp
count
group_idx
group_name
length_bin_idx
length_start_bp
length_end_bp
```

Full-array extraction should be explicit:

```r
midpoint_array(midpoints)
```

Document that this loads the full 3D count tensor into memory.

## End-Motif Counts API

Loader:

```r
ends <- read_end_motifs("sample.end_motifs.zarr")
```

Shared metadata:

```r
storage_mode(ends)
motifs(ends)
motif_idx(ends, motif)
has_motif(ends, motif)
```

`motifs()` should return:

```r
motif_idx
motif
```

`has_motif()` should return a scalar logical without scanning count arrays.

### Global End Motifs

For `cfdnalab_global_end_motif_counts`:

```r
dense_counts_vector(ends)
dense_data_frame(ends)
sparse_counts_matrix(ends)
sparse_data_frame(ends)
```

`dense_counts_vector()` should return a named numeric vector when dense output is available or explicitly requested.

`dense_data_frame()` should return:

```r
motif_idx
motif
count
```

### Windowed End Motifs

For `cfdnalab_windowed_end_motif_counts`:

```r
windows(ends)
dense_counts_matrix(ends)
dense_data_frame_for_window(ends, window_idx)
dense_data_frame_for_motif(ends, motif)
sparse_counts_matrix(ends)
sparse_data_frame(ends)
sparse_data_frame_for_window(ends, window_idx)
sparse_data_frame_for_motif(ends, motif)
```

`windows()` should return:

```r
window_idx
chromosome
chromosome_name
window_start_bp
window_end_bp
blacklisted_fraction
```

Sparse helpers should avoid constructing the full dense matrix.

### Grouped End Motifs

For `cfdnalab_grouped_end_motif_counts`:

```r
groups(ends)
group_idx(ends, group_name)
dense_counts_matrix(ends)
dense_data_frame_for_group(ends, group)
dense_data_frame_for_motif(ends, motif)
sparse_counts_matrix(ends)
sparse_data_frame(ends)
sparse_data_frame_for_group(ends, group)
sparse_data_frame_for_motif(ends, motif)
```

`groups()` should return:

```r
group_idx
group_name
eligible_windows
blacklisted_fraction
```

## Dense and Sparse Behavior

Rules:

- If the store is dense, dense helpers read dense arrays.
- If the store is sparse, sparse helpers read sparse COO arrays and build `Matrix::sparseMatrix`.
- Dense helpers on sparse stores should either error by default or require `allow_densify = TRUE`.
- Sparse helpers on dense stores may be supported, but the documentation must say that the dense matrix is read first.
- Long data frames from large sparse stores should contain only non-zero rows unless the function name clearly says dense.

Recommended first behavior:

```r
dense_counts_matrix(x, allow_densify = FALSE)
```

If `x` is sparse and `allow_densify = FALSE`, error with a message telling the user to use `sparse_counts_matrix()` or set `allow_densify = TRUE`.

## Internal Layout

Suggested files:

```text
R/
  zarr_helpers.R
  schema.R
  midpoints.R
  end_motifs.R
  print.R
  utils.R
tests/
  testthat/
    test-midpoints.R
    test-end-motifs.R
    test-schema-validation.R
```

`zarr_helpers.R` should own:

- root metadata reads
- array metadata reads
- backend-specific array slicing
- JSON label extraction
- `.zarr` path validation

`schema.R` should own:

- supported schema versions
- schema compatibility table
- root schema validation

Output-specific files should own biological semantics and public helper functions.

Documentation convention:

- Roxygen comments in `R/*.R` are the source of truth for user-facing docs,
  internal helper docs, imports, exports, S3 registration, `NAMESPACE`, and
  generated `man/*.Rd` files.
- Do not hand-edit generated `.Rd` files.
- Do not hand-edit generated `NAMESPACE` entries except as a temporary repair
  that must be moved back into roxygen comments.
- Internal helpers should still have roxygen blocks with `@noRd` so the source
  documents intent without creating user help pages.

## Testing Plan

Use actual cfDNAlab-generated fixtures, not hand-written Zarr stores, for integration tests.

Minimum tests:

- midpoint loader validates schema and version
- midpoint `groups()`, `length_bins()`, `positions()`
- midpoint `profile_data_frame()` for one group and length bin
- midpoint `profile_array()` matches the corresponding data frame counts
- end-motif loader returns correct subclass by `row_mode`
- end-motif `motifs()`, `motif_idx()`, `has_motif()`
- dense global end-motif data frame
- dense windowed end-motif window slice
- sparse end-motif `Matrix::sparseMatrix` round-trip
- sparse motif slice does not densify
- invalid schema gives a useful error
- sparse densification requires explicit opt-in

Local package install from the repository root:

```bash
Rscript -e 'install.packages("pak")'
Rscript -e 'pak::pak("./r-cfdnalab")'
```

Downstream GitHub Action should install the R dependencies directly and run:

```bash
Rscript -e 'testthat::test_local("r-cfdnalab")'
```

The current downstream Zarr workflow runs `testthat::test_local("r-cfdnalab")`
after generating the shared midpoint and end-motif fixtures with cfDNAlab. It
also installs `r-cfdnalab` and runs
`downstream_tests/R/test_cfdnalab_r_package.R`, which exercises the package as
an external downstream user would.

## README Examples

Keep examples short and user-facing.

Midpoint example:

```r
library(cfdnalab)

midpoints <- read_midpoints("sample.midpoint_profiles.zarr")
groups(midpoints)

profile <- profile_data_frame(
  midpoints,
  group = "LYL1",
  length_bin_idx = 1
)
```

End-motif example:

```r
library(cfdnalab)

ends <- read_end_motifs("sample.end_motifs.zarr")
motifs(ends)

motif_counts <- sparse_data_frame_for_motif(ends, motif = "_AA")
```

Do not make the README an exhaustive API reference. The package help pages should carry details about dense vs sparse behavior.

## Open Decisions

- Should GC Zarr package loaders be included in the first R package release or deferred?
- Should the R package include plotting helpers now, or only return analysis-ready objects?

Resolved for the first implementation:

- Use CRAN `zarr` as the default backend.
- Return base `data.frame` objects.
- Dense helpers on sparse stores error by default and require `allow_densify = TRUE`.

## First Milestone

Implement a minimal package that can:

1. Load midpoint profile Zarr stores.
2. Load end-motif Zarr stores.
3. Return metadata data frames.
4. Extract one midpoint profile as a data frame.
5. Extract end-motif sparse counts as `Matrix::sparseMatrix`.
6. Pass downstream tests on actual cfDNAlab-generated fixtures.

Plot helpers and GC package loaders can wait until the core readers are stable.

Current implementation status:

- Initial S3 loader implementation exists under `r-cfdnalab/`.
- The first backend is CRAN `zarr`.
- `testthat` integration tests use cfDNAlab-generated downstream Zarr fixtures when present and skip when fixtures have not been generated.
- Additional validation tests cover path checks, schema attributes, dimension names, labels, zero-based axes, sparse indices, sparse counts, interval metadata, scalar index validation, and motif ASCII decoding.
- Remaining validation should come from running those tests in the downstream workflow with the generated fixtures.

Pre-release cleanup:

- Replace placeholder `Authors@R` and `LICENSE` holder values with final package maintainer details.
- Re-run roxygen after any public API or documentation changes and keep generated `NAMESPACE` / `man/*.Rd` in sync.
- Convert the R public API from zero-based `*_idx` selectors and columns to
  one-based `*_idx` selectors and columns. Keep zero-based schema values only
  in clearly named internal `*_idx0` or `*_index0` fields and helper-local
  variables.
