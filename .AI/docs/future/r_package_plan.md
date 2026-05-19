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

The current public loader API is specified in
[`../specs/package_loader_api.md`](../specs/package_loader_api.md). Keep this
plan focused on R package structure, backend choices, plotting ideas, and future
extensions rather than duplicating every selector mapping.

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
- `data.table`: fast TSV reader for length-count outputs

Avoid hard dependency on `tibble` or `dplyr` in the first version.

Return base `data.frame` objects from public data frame helpers by default. The
loader may keep efficient internal representations, such as the numeric count
matrix and parsed metadata, but public data frame helpers should not require
users to know `data.table` syntax. Base data frames work naturally with:

```r
tibble::as_tibble(x)
data.table::as.data.table(x)
```

For `read_lengths()`, keep compressed TSV reading behind an internal helper.
Use `data.table::fread()` for TSV parsing. For `.tsv.zst` files, use
`data.table::fread(cmd = "zstd -dc <path>")` with proper shell quoting and
validate that `Sys.which("zstd")` finds a system `zstd` binary before reading.
If `zstd` is unavailable, error clearly instead of falling back to a slower or
partial path. The public API should be `read_lengths(path)`, not a generic
decompression utility.

If R plotting helpers are added, keep plotting packages out of `Imports`:

- `ggplot2` can be a `Suggests` dependency because plotting is optional.
- `ggridges` should not be a hard dependency. Either use it only when installed
  for ridge-style plots or draw a simple internal polygon-based ridge view.

Current candidate Zarr readers:

- CRAN `zarr`: native R implementation with Zarr v3 support.
- Bioconductor `Rarr`: promising for `DelayedArray`, but `Rarr` 1.6.0 expects
  Zarr v2 `.zarray` metadata and cannot read the current Zarr v3 stores.
- CRAN/R-universe `pizzarr`: useful to test, but do not assume full compatibility until downstream tests cover our stores.

First implementation decision:

- Start with CRAN `zarr` as the internal backend because it has native Zarr v3 support and keeps installation simpler for ordinary R users.
- Keep all backend access in `R/helpers.R` so we can swap or add `Rarr` later if downstream tests reveal missing support.

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
read_lengths(path)
```

Zarr loaders should:

- require `path` to exist
- require `path` to be a directory
- require `path` to end in `.zarr`
- read root `zarr.json`
- validate `cfdnalab_schema`
- validate `cfdnalab_schema_version`
- give package-specific errors for wrong schema, missing arrays, or unsupported storage modes

`read_lengths()` should load `<prefix>.length_counts.tsv.zst` outputs from the
`cfdna lengths` command. It should infer its public metadata from the TSV
itself, not from `<prefix>.length_settings.json`:

- Count columns named `count_<length>` represent half-open bins
  `[length, length + 1)`.
- Count columns named `count_<start>_<end>` represent half-open bins
  `[start, end)`.
- Global output has only count columns.
- Windowed output has `chrom`, `start`, and `end`, plus optional
  `blacklisted_fraction`.
- Grouped output has `group_name` and `eligible_windows`, plus optional
  `blacklisted_fraction`.

The settings JSON is not part of the public `read_lengths()` contract. Users
who need command provenance can read the JSON directly with `jsonlite`. The
loader should only parse and validate the TSV shape needed for downstream R
workflows.

`read_lengths()` should return mode-specific S3 classes:

```r
class(x) <- c(
  "cfdnalab_global_length_counts",
  "cfdnalab_length_counts"
)

class(x) <- c(
  "cfdnalab_windowed_length_counts",
  "cfdnalab_length_counts"
)

class(x) <- c(
  "cfdnalab_grouped_length_counts",
  "cfdnalab_length_counts"
)
```

Do not expose a generic row-metadata concept for lengths. Use `window_metadata()` for
windowed outputs and `group_metadata()` for grouped outputs, matching the rest
of the R API.

## Shared Generics

Use generic functions where names make sense across outputs:

```r
schema_version(x)
storage_mode(x)
group_metadata(x)
length_bins(x)
positions(x)
motifs(x)
window_metadata(x)
```

Only implement methods where the concept exists. For example:

- `group_metadata()` exists for midpoint profiles and grouped end motifs.
- `group_metadata()` also exists for grouped length counts.
- `positions()` exists for midpoint profiles.
- `window_metadata()` exists for windowed end motifs and windowed length counts.
- `position_bins()` may be added later for range-based window lookup, but is
  not part of the first length-count implementation.
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
to `index0`, then use `index0 + 1L` only at the final R array/data frame
subsetting boundary through shared helpers. Tests must use asymmetric fixtures
where choosing index `1` versus zero-based index `1` returns different values.

## Midpoint Profiles API

Loader:

```r
midpoints <- read_midpoints("sample.midpoint_profiles.zarr")
```

Metadata:

```r
group_metadata(midpoints)
length_bins(midpoints)
positions(midpoints)
group_idx(midpoints, group_name)
length_bin_idx(midpoints, length)
```

Expected metadata data frames:

```r
group_metadata(midpoints)
# group_idx, group_name, eligible_intervals

length_bins(midpoints)
# length_bin_idx, length_start_bp, length_end_bp

positions(midpoints)
# position_idx, position_bin_start_bp, position_bin_end_bp
```

Tabular extraction:

```r
midpoint_data_frame(
  midpoints,
  groups = NULL,
  group_idxs = NULL,
  with_lengths = NULL,
  length_bin_idxs = NULL
)
profile_array(midpoints, group_idx, length_bin_idx)
```

`midpoint_data_frame()` should return one row per selected group, length bin,
and position bin:

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

Use `with_lengths` for user fragment lengths that should be resolved to
containing bins. Use `length_bin_idxs` when selecting bins directly.

Full-array extraction should be explicit:

```r
midpoint_array(midpoints)
```

Document that this loads the full 3D count tensor into memory.

## Length Counts API

Loader:

```r
lengths <- read_lengths("sample.length_counts.tsv.zst")
```

The loader should parse the wide TSV written by `cfdna lengths` and expose a
length-bin-aware object for reshaping and plotting. It should not require or
implicitly read the settings JSON.

Shared metadata and extraction:

```r
length_bins(lengths)
length_bin_idx(lengths, length)
length_counts_matrix(lengths)
```

`length_bins()` should return:

```r
length_bin_idx
length_start_bp
length_end_bp
length_midpoint_bp
length_width_bp
```

`length_counts_matrix()` should return the numeric count matrix only. Its rows
are the output units from the TSV in file order, and its columns are the parsed
length bins. This helper is allowed to use row order internally, but user-facing
metadata should stay expressed as windows, groups, or global output.

`length_data_frame()` adapts to the loaded output mode. Global outputs have no
row selector. Windowed outputs accept `window_idxs`. Grouped outputs accept
`groups` or `group_idxs`.

By default it should return one record per selected output unit and length bin:

```r
length_bin_idx
length_start_bp
length_end_bp
length_midpoint_bp
length_width_bp
count
```

For windowed length counts, include the window columns:

```r
window_idx
chrom
start
end
blacklisted_fraction
```

For grouped length counts, include the group columns:

```r
group_idx
group_name
eligible_windows
blacklisted_fraction
```

The `blacklisted_fraction` column is present only when the TSV includes it.

If `keep_wide = TRUE`, `length_data_frame()` should return one row per output
unit and keep one value column per length bin. `value = "count"` should keep
the original `count_*` column names. `value = "fraction"` should return
`fraction_*` columns, and `value = "density"` should return `density_*`
columns, using the same length-bin suffixes as the source count columns.

`max_blacklisted_fraction` should filter output units before reshaping. It must
be a single finite fraction in 0..1. The default `1.0` keeps all rows with valid
blacklist fractions and effectively disables filtering. Values below `1.0`
keep only windows or groups with
`blacklisted_fraction <= max_blacklisted_fraction`. If the loaded output does
not contain `blacklisted_fraction`, the default `1.0` should keep all rows, and
stricter cutoffs should error with a clear message because the package cannot
evaluate the requested filter.

`value = "fraction"` should compute the within-output-unit fraction:

```r
count / sum(count)
```

`value = "density"` should compute that fraction divided by
`length_width_bp`. This makes unequal length-bin widths comparable.

Do not use `eligible_windows` for normalization. It is useful for filtering and
quality control. The TSV does not contain eligible bases, so the R package
should not expose per-base or per-eligible-window normalization helpers for
length outputs.

### Global Length Counts

For `cfdnalab_global_length_counts`:

```r
length_counts_vector(lengths)
length_data_frame(
  lengths,
  value = c("count", "fraction", "density"),
  keep_wide = FALSE
)
```

`length_counts_vector()` should return a named numeric vector with one value per
length bin.

### Windowed Length Counts

For `cfdnalab_windowed_length_counts`:

```r
window_metadata(lengths)
length_data_frame(
  lengths,
  window_idxs = NULL,
  value = c("count", "fraction", "density"),
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0
)
```

`window_metadata()` should return:

```r
window_idx
chrom
start
end
blacklisted_fraction
```

`window_idx` is one-based and follows the output order in the TSV. Genomic
range lookup is a later development stage. In the first implementation,
windowed length-count selection should use `window_idxs`, and ad hoc filtering
should happen through ordinary R subsetting on `window_metadata(lengths)`.
`length_data_frame()` should accept optional `window_idxs` selection for the
common case where users want one or a few windows.

In a later stage, `position_bins()` may return rows from `window_metadata(lengths)`
that overlap a requested genomic range:

```r
position_bins(lengths, chrom = "chr1", start = 50000000, end = 75000000)
```

The returned data frame should include:

```r
position_bin_idx
window_idx
chrom
start
end
```

Include `blacklisted_fraction` when the loaded TSV has that column.

`position_bin_idx` is a one-based public index for selecting bins on the
windowed genomic track. It should follow the plotted bin order. For fixed-size
windows this is the same order as `window_metadata(lengths)`. The helper should use the
same `chrom`, `start`, and `end` coordinate system as `window_metadata(lengths)` only
for range lookup. It should not expose or convert internal zero-based schema
indices. Validate that `chrom`, `start`, and `end` are supplied together when
`ranges = NULL`, that `start < end`, and that at least one output window
overlaps the requested range. For multiple requested ranges, `ranges` should be
a data frame with `chrom`, `start`, and `end` columns. Use either `ranges` or
the scalar `chrom`/`start`/`end` arguments, not both. Return the union of
overlapping position bins in plotted order.

### Grouped Length Counts

For `cfdnalab_grouped_length_counts`:

```r
group_metadata(lengths)
group_idx(lengths, group_name)
length_data_frame(
  lengths,
  groups = NULL,
  group_idxs = NULL,
  value = c("count", "fraction", "density"),
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0
)
```

`group_metadata()` should return:

```r
group_idx
group_name
eligible_windows
blacklisted_fraction
```

`length_data_frame()` should accept optional `groups` or `group_idxs`
selection. Use either `groups` or `group_idxs`, not both.

### Length-Bin Ratio Helpers

Ratios between user-selected length bins are common downstream summaries.
Short/long ratios are one example, but the API should be general enough for
any ratio of two length bins.

Shared helper:

```r
length_ratio_data_frame(
  x,
  num_length_bin_idx,
  denom_length_bin_idx,
  denom_zero = c("NA", "error"),
  max_blacklisted_fraction = 1.0
)
```

`num_length_bin_idx` and `denom_length_bin_idx` should be one-based indices
from `length_bins(x)`. For example, a short/long-style ratio can be created by
running `cfdna lengths --length-bins 100 151 221` and often `--by-size 5000000` 
and then using length bin 1 over length bin 2. The R helper should not reconstruct broad ranges by summing many bins. Users should choose the ratio bins during the command call when they want this derived feature.

Validation rules:

- `num_length_bin_idx` and `denom_length_bin_idx` must each be a scalar
  one-based length-bin index.
- Both indices must exist in `length_bins(x)`.
- The numerator and denominator bin indices must be different.
- Denominator zero handling must be explicit. The default can return `NA` for
  undefined ratios, but the returned data frame should include
  `num_count` and `denom_count` so users can filter low-support windows or
  groups.
- `max_blacklisted_fraction` should filter windowed or grouped outputs before
  computing ratios, using the same rules as `length_data_frame()`.

Returned columns should include the output-unit metadata plus:

```r
num_length_bin_idx
num_length_start_bp
num_length_end_bp
denom_length_bin_idx
denom_length_start_bp
denom_length_end_bp
num_count
denom_count
ratio
```

For windowed length counts, include:

```r
position_bin_idx
window_idx
chrom
start
end
blacklisted_fraction
```

`position_bin_idx` is one-based and follows the plotted genomic-track bin
order. For fixed-size windows, it should match the row order used by
`window_metadata(lengths)`.

For grouped length counts, include:

```r
group_idx
group_name
eligible_windows
blacklisted_fraction
```

Do not hide low-support behavior. Documentation should tell users to inspect
`num_count`, `denom_count`, their sum, and blacklist fraction before treating a
ratio as reliable.

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

Tabular extraction uses one public function:

```r
end_motif_data_frame(
  ends,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL
)
```

Mode-specific methods add only the selectors that make sense for the loaded
mode:

```r
end_motif_data_frame(
  ends,
  window_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0
) # windowed

end_motif_data_frame(
  ends,
  groups = NULL,
  group_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0
) # grouped
```

For sparse stores, `densify = FALSE` returns stored non-zero motif-count rows.
`densify = TRUE` adds explicit zero-count rows for the selected rows and the
selected observed motifs. Dense stores ignore `densify`.

### Global End Motifs

For `cfdnalab_global_end_motif_counts`:

```r
dense_counts_vector(ends)
sparse_counts_matrix(ends)
end_motif_data_frame(ends, densify = FALSE)
end_motif_data_frame(ends, densify = TRUE)
```

`dense_counts_vector()` should return a named numeric vector when dense output is available or explicitly requested.

`end_motif_data_frame()` should return:

```r
motif_idx
motif
count
```

### Windowed End Motifs

For `cfdnalab_windowed_end_motif_counts`:

```r
window_metadata(ends)
dense_counts_matrix(ends)
sparse_counts_matrix(ends)
end_motif_data_frame(
  ends,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  densify = FALSE,
  max_blacklisted_fraction = 1.0
)
```

`window_metadata()` should return:

```r
window_idx
chrom
start
end
blacklisted_fraction
```

Use the same public genomic-window column names for windowed end motifs and
windowed length counts. `chrom`, `start`, and `end` match the standard
BED-like table shape. They are genomic coordinates, not R indices.

Sparse helpers should avoid constructing the full dense matrix.
Windowed end-motif data frame helpers should accept
`max_blacklisted_fraction = 1.0` with the same 0..1 validation and keep-all
default as `length_data_frame()`.

### Grouped End Motifs

For `cfdnalab_grouped_end_motif_counts`:

```r
group_metadata(ends)
group_idx(ends, group_name)
dense_counts_matrix(ends)
sparse_counts_matrix(ends)
end_motif_data_frame(
  ends,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  densify = FALSE,
  max_blacklisted_fraction = 1.0
)
```

`group_metadata()` should return:

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
- Dense helpers on sparse stores should make densification explicit in the
  helper name.
- Sparse helpers on dense stores may be supported, but the documentation must say that the dense matrix is read first.
- End-motif data frames from sparse stores contain only non-zero rows by
  default.
- `end_motif_data_frame(..., densify = TRUE)` adds explicit zero-count rows for
  selected rows and selected observed motifs.
- Dense and sparse end-motif data frames use the same columns for a given row
  mode; only the returned rows differ.

## Length Plotting API

Plotting should be domain-aware and built on `cfdnalab_length_counts` objects.
Do not add generic plotting wrappers for arbitrary user data frames.

Use one S3 generic:

```r
plot_lengths <- function(x, ...) {
  UseMethod("plot_lengths")
}
```

Methods can have different arguments by class. Keep separate public plotting
functions out of the first version unless the display families become too
different to document under one generic.

Shared value choices:

```r
value = c("fraction", "count", "density")
```

Definitions must match `length_data_frame()`:

- `count`: loaded count or weighted mass.
- `fraction`: within-output-unit fraction over the length axis.
- `density`: fraction divided by `length_width_bp`.

The plotting API should not expose per-base or per-window normalization for
length outputs because eligible bases are not available in the TSV.

For `cfdnalab_global_length_counts`:

```r
plot_lengths(
  x,
  value = c("fraction", "count", "density"),
  display = c("line", "bar"),
  ...
)
```

For `cfdnalab_grouped_length_counts`:

```r
plot_lengths(
  x,
  groups = NULL,
  value = c("fraction", "count", "density"),
  max_blacklisted_fraction = 1.0,
  display = c("line", "facet", "ridges", "tiles"),
  max_groups = 12,
  ...
)
```

For `cfdnalab_windowed_length_counts`:

```r
plot_lengths(
  x,
  window_idxs = NULL,
  value = c("fraction", "count", "density"),
  max_blacklisted_fraction = 1.0,
  display = c("line", "facet", "ridges", "tiles"),
  max_windows = 12,
  ...
)
```

Display modes:

- `line`: length on the x-axis and value on the y-axis, with groups or selected
  windows mapped to color.
- `facet`: one length-distribution panel per group or selected window.
- `ridges`: stacked vertical distributions built from the binned length values.
  This is not a kernel-density estimate over raw fragments.
- `tiles`: group or selected window on the y-axis, length on the x-axis, and
  value mapped to fill.

Default behavior:

- Global output plots the single distribution.
- Grouped output may plot all groups when the number of groups is at most
  `max_groups`. Otherwise require `groups`.
- Windowed output should require selected `window_idxs` unless the number of
  windows is at most `max_windows`.
- All plotting helpers return `ggplot` objects so users can adjust themes,
  labels, colors, and export settings.
- Built-in defaults should be polished enough for reports, including meaningful
  axis labels, legend labels, and sensible scales. This is part of the value of
  package plotting, not a replacement for user customization.

### Length-Bin Ratio Plotting

Use a separate plotting generic for ratios because these plots show derived
summaries across output units, not length distributions themselves:

```r
plot_length_ratio <- function(x, ...) {
  UseMethod("plot_length_ratio")
}
```

The first useful plotting method should be for windowed length counts:

```r
plot_length_ratio(
  x,
  num_length_bin_idx,
  denom_length_bin_idx,
  position_bin_idx = NULL,
  denom_zero = c("NA", "error"),
  min_num_count = 0,
  min_denom_count = 0,
  min_total_count = 0,
  max_blacklisted_fraction = 1.0,
  support_color = c("none", "total_count", "num_count", "denom_count", "blacklisted_fraction"),
  require_regular_windows = TRUE,
  display = c("genome", "chromosome_facets"),
  ...
)
```

`num_length_bin_idx` and `denom_length_bin_idx` use the same one-based
length-bin selectors as `length_ratio_data_frame()`.

`plot_length_ratio()` should be a convenience wrapper. It should call
`length_ratio_data_frame()` and then pass the result to the explicit data frame
plotter:

```r
plot_length_ratio_track(
  ratios,
  position_bin_idx = NULL,
  min_num_count = 0,
  min_denom_count = 0,
  min_total_count = 0,
  max_blacklisted_fraction = 1.0,
  support_color = c("none", "total_count", "num_count", "denom_count", "blacklisted_fraction"),
  require_regular_windows = TRUE,
  display = c("genome", "chromosome_facets"),
  ...
)
```

`plot_length_ratio_track()` should accept the output of
`length_ratio_data_frame()` for windowed length counts. This gives users a
direct path for computing ratios, filtering or annotating the data frame, and
then plotting the same genomic ratio track without re-reading or recomputing
the length-count object.

Ratio track plotting should not expose `window_idx` selection. Window indices
are arbitrary output-order identifiers and are not meaningful for selecting a
range on a genome-scale track. Use `position_bin_idx` for selecting plotted
position bins.

`position_bin_idx` should be `NULL` or a one-based integer vector of plotted
position-bin indices. Users can get these indices from
`position_bins(lengths, chrom, start, end)` for windowed length counts, then
pass `bins$position_bin_idx` to `plot_length_ratio()` or
`plot_length_ratio_track()`. The plotting helper should validate that requested
indices exist in the ratio data frame and should error clearly if no ratio rows
remain after subsetting. Users who need more specialized subsetting can filter
the ratio data frame before calling `plot_length_ratio_track()`.

Windowed ratio plotting should target fixed-size genomic windows, for example
5Mb bins. Because the TSV itself does not record `window_mode`, the plotting
method should infer whether the window table is regular enough to plot as a
genomic track:

- Windows must have `chrom`, `start`, and `end`.
- Within each chromosome, windows must be sorted, non-overlapping, and
  contiguous.
- Window widths should be consistent, allowing only the final window on a
  chromosome to be shorter.

If those checks fail and `require_regular_windows = TRUE`, error with a message
telling the user to use `length_ratio_data_frame()` and custom plotting for
irregular BED windows.

`display = "genome"` should plot windows in chromosome order as one concatenated
genome-like track. `display = "chromosome_facets"` should facet by chromosome
and use each window's genomic midpoint on the x-axis.

The returned plot should keep low-support behavior visible. The ratio data
should include `num_count`, `denom_count`, and `ratio`, and the plotting method
should document how `min_num_count`, `min_denom_count`, `min_total_count`, and
`max_blacklisted_fraction` filter unstable or heavily blacklisted windows before
plotting.

`support_color` should optionally color plotted windows by a support or QC
column using a continuous gradient scale and a clear legend. It should default
to `"none"` so the basic ratio track stays visually simple. Useful choices are:

- `"total_count"`: colors by `num_count + denom_count`, which is usually the
  best first support metric.
- `"num_count"` and `"denom_count"`: color by support in either selected
  length bin.
- `"blacklisted_fraction"`: colors by blacklist overlap when available, and
  errors clearly if the loaded output has no blacklist metadata.

## Internal Layout

Suggested files:

```text
R/
  helpers.R
  schema.R
  midpoints.R
  end_motifs.R
  lengths.R
  plotting.R
  print.R
  utils.R
tests/
  testthat/
    test-midpoints.R
    test-end-motifs.R
    test-lengths.R
    test-length-plots.R
    test-schema-validation.R
```

`helpers.R` should own:

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

`lengths.R` should own TSV parsing, length-bin extraction, length-count
objects, and length data frame helpers.

`plotting.R` should own optional plotting helpers and must not make plotting
packages required for non-plotting use.

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

Use actual cfDNAlab-generated fixtures, not hand-written Zarr stores, for
Zarr integration tests. Length TSV parser unit tests may use small hand-written
TSV fixtures, with downstream tests covering real `cfdna lengths` outputs.

Minimum tests:

- midpoint loader validates schema and version
- midpoint `group_metadata()`, `length_bins()`, `positions()`
- midpoint `midpoint_data_frame()` for selected groups and length bins
- midpoint `profile_array()` matches the corresponding data frame counts
- end-motif loader returns correct subclass by `row_mode`
- end-motif `motifs()`, `motif_idx()`, `has_motif()`
- global end-motif data frame
- windowed end-motif data frame with `window_idxs`
- sparse end-motif `Matrix::sparseMatrix` round-trip
- sparse motif slice does not densify
- invalid schema gives a useful error
- sparse densification requires explicit opt-in
- length TSV loader detects global, windowed, and grouped outputs
- length TSV loader parses `count_<length>` and `count_<start>_<end>` columns
  as half-open bins
- length `length_bins()`, `length_bin_idx()`, `window_metadata()`, and
  `group_metadata()` use one-based public indices
- length `length_data_frame()` computes `count`, `fraction`, and `density`
  without using `eligible_windows` as a normalization denominator
- length `length_data_frame(..., keep_wide = TRUE)` returns selected rows in
  wide shape for `count`, `fraction`, and `density`, with `fraction_*` and
  `density_*` columns when transformed values are requested
- length `length_data_frame()` and ratio helpers use
  `max_blacklisted_fraction = 1.0` as the keep-all default, filter by
  `blacklisted_fraction` for stricter cutoffs, and error for stricter cutoffs
  if the output has no blacklist metadata
- length ratio helpers use one-based numerator and denominator length-bin
  indices from `length_bins()`
- length ratio helpers expose `num_count` and `denom_count`, and handle
  denom-zero behavior explicitly
- length plotting methods return `ggplot` objects for global, grouped, and
  selected windowed outputs when `ggplot2` is installed
- grouped and windowed length plots require explicit selections when output
  size exceeds `max_groups` or `max_windows`
- windowed length-ratio plotting validates fixed-size-like genomic windows
  before drawing genome-track plots
- `plot_length_ratio_track()` plots already-computed windowed ratio data frames
  from `length_ratio_data_frame()`
- ratio track plotting supports one-based `position_bin_idx` selection and does
  not expose arbitrary `window_idx` selection
- windowed length-ratio plotting can add support/QC gradient coloring for
  total count, numerator count, denominator count, and blacklist fraction

Local package install from the repository root:

```bash
Rscript -e 'install.packages("pak")'
Rscript -e 'pak::pak("./r-cfdnalab")'
```

Downstream GitHub Action should install the R dependencies directly and run:

```bash
Rscript -e 'testthat::test_local("r-cfdnalab")'
```

The current downstream workflow runs `testthat::test_local("r-cfdnalab")`
after generating the shared midpoint, end-motif, and length fixtures with
cfDNAlab. It also installs `r-cfdnalab` and runs
`downstream_tests/R/test_cfdnalab_r_package.R`, which exercises the package as
an external downstream user would.

## README Examples

Keep examples short and user-facing.

Midpoint example:

```r
library(cfdnalab)

midpoints <- read_midpoints("sample.midpoint_profiles.zarr")
group_metadata(midpoints)

profile <- midpoint_data_frame(
  midpoints,
  groups = "LYL1",
  length_bin_idxs = 1
)
```

End-motif example:

```r
library(cfdnalab)

ends <- read_end_motifs("sample.end_motifs.zarr")
motifs(ends)

motif_counts <- end_motif_data_frame(ends, motifs = "_AA")
```

Length example:

```r
library(cfdnalab)

lengths <- read_lengths("sample.length_counts.tsv.zst")
length_bins(lengths)

length_distribution <- length_data_frame(
  lengths,
  value = "fraction",
  max_blacklisted_fraction = 0.1
)
plot_lengths(lengths)

ratios <- length_ratio_data_frame(
  lengths,
  num_length_bin_idx = 1L,
  denom_length_bin_idx = 2L,
  max_blacklisted_fraction = 0.1
)
```

Do not make the README an exhaustive API reference. The package help pages should carry details about dense vs sparse behavior.

## Open Decisions

- Should GC Zarr package loaders be included in the first R package release or deferred?

Resolved for the first implementation:

- Use CRAN `zarr` as the default backend.
- Return base `data.frame` objects.
- End-motif sparse data frames return stored non-zero rows by default and
  require `densify = TRUE` for explicit zero-count rows.
- Defer generic plotting helpers, but add a focused `plot_lengths()` generic
  for length-count outputs. Length plotting should return `ggplot`
  objects and stay optional through `Suggests`.

## First Milestone

Implement a minimal package that can:

1. Load midpoint profile Zarr stores.
2. Load end-motif Zarr stores.
3. Load length-count TSV outputs.
4. Return metadata data frames.
5. Extract midpoint, end-motif, and length-count data frames.
6. Extract midpoint arrays, end-motif matrices, and length-count matrices or vectors.
7. Pass downstream tests on actual cfDNAlab-generated fixtures.

Plot helpers and GC package loaders can wait until the core readers are stable.

## Next R Extension Milestone

After the core readers are stable, add length-specific derived summaries and
optional plotting:

1. Add length-bin ratio helpers for user-selected numerator and denominator
   length-bin indices.
2. Add `plot_lengths()` with line, facet, ridge-style, and tile displays for
   the output modes where those displays make sense.
3. Add `plot_length_ratio_track()` for already-computed windowed ratio data
   frames and `plot_length_ratio()` as the length-count object convenience
   wrapper.
4. Add genomic range lookup helpers for windowed outputs if repeated workflows
   show that `window_metadata()` plus ordinary data frame filtering is not
   enough.

Current implementation status:

- Initial S3 loader implementation exists under `r-cfdnalab/`.
- The first backend is CRAN `zarr`.
- `testthat` integration tests use cfDNAlab-generated downstream Zarr fixtures when present and skip when fixtures have not been generated.
- Additional validation tests cover path checks, schema attributes, dimension names, labels, zero-based axes, sparse indices, sparse counts, interval metadata, scalar index validation, and motif ASCII decoding.
- Remaining validation should come from running those tests in the downstream workflow with the generated fixtures.

Pre-release cleanup:

- Replace placeholder `Authors@R` and `LICENSE` holder values with final package maintainer details.
- Re-run roxygen after any public API or documentation changes and keep generated `NAMESPACE` / `man/*.Rd` in sync.
