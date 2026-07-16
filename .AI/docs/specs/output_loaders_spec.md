# Rust Output Loaders Spec

This spec describes the current public Rust loaders under
`cfdnalab::output_loaders`. These loaders read files written by cfDNAlab
commands and return typed Rust objects for downstream code.

The loaders are strict readers for cfDNAlab outputs. They are not generic TSV or
Zarr importers.

## Feature Gates

The `cfdnalab::output_loaders` module is compiled when at least one supported
loader feature is enabled. Each loader is compiled with the cargo feature for
the command that writes the matching output:

- `cmd_lengths` exposes `load_lengths_output`.
- `cmd_fcoverage` exposes `load_fcoverage_output` and
  `load_fcoverage_output_with_group_index`.
- `cmd_ends` exposes `load_ends_output`.
- `cmd_midpoints` exposes `load_midpoints_output`.
- `cmd_ref_kmers` exposes `load_ref_kmers_output`.

`DenseMatrix`, `LengthBin`, and `WindowRow` are exported whenever
`output_loaders` is compiled. `DenseArray3` is exported with `cmd_midpoints`.

## Error Contract

Public loader functions and public fallible loader methods return
`OutputLoaderResult<T>`.

`OutputLoaderError` is the stable public error type. It keeps the contextual
parser error internally and exposes it through `as_anyhow()` for callers that
want the full message chain.

Error messages should include enough path, line, Zarr-array, or selector
context for a user to fix the malformed output. Exact message text is not a
public compatibility promise unless a test pins a specific behavior.

## Indexing And Axes

All public loader indices are zero-based.

Text output rows are returned in file order. Zarr axes are returned in stored
axis order after validation that axis arrays are zero-based, contiguous, and
shape-compatible with the count arrays.

Ranges use half-open intervals `[start, end)`.

Selectors preserve requested order and reject duplicate explicit selectors. A
missing selector means "select all values on this axis". Explicit row selectors
are rejected for global lengths, global end-motif outputs, and global reference
k-mer outputs because those outputs do not have a meaningful selectable row axis.

## Shared Containers

`DenseMatrix<T>` stores row-major two-dimensional values with explicit row and
column counts. It provides shape-aware row iteration, value lookup, and access
to row-major values.

`DenseArray3<T>` stores row-major three-dimensional values with the last axis
varying fastest. Midpoint counts use `(group, length_bin, position)` order.

## Lengths Loader

`load_lengths_output(path)` reads plain, gzip, or zstd-compressed
`cfdna lengths` TSV files.

Supported row modes:

- Global: one count row and no row metadata.
- Windows: `chrom`, `start`, `end`, optional `blacklisted_fraction`, and count
  columns.
- Groups: `group_name`, `eligible_windows`, optional `blacklisted_fraction`, and
  count columns.

Count columns define the fragment length axis. Supported headers are
`count_<length>` for a single-bp bin and `count_<start>_<end>` for a half-open
bin. Length bins must be contiguous, start at or above the minimum supported
fragment length, and end at or below the maximum supported fragment length edge.

Counts must be finite and non-negative. Group names must be unique.

Selections return `LengthCountSelection` with selected row metadata, selected
length-bin metadata, and an owned dense count matrix.

## Fcoverage Loader

`load_fcoverage_output(path)` reads non-positional aggregate `cfdna fcoverage`
TSV files. Plain text, gzip, and zstd compression are supported.

`load_fcoverage_output_with_group_index(path, group_index_path)` additionally
loads the matching `group_index.tsv` file for grouped outputs, attaches group
names, and enables name-based group selection.

Supported row modes:

- Windows: `chromosome`, `start`, `end`, value or summary columns, and support
  metadata.
- Groups: `group_idx`, support metadata, and value or summary columns.

Supported data modes:

- Scalar average values from `average_<signal>`.
- Scalar total values from `total_<signal>`.
- Summary statistics from `total_<signal>`, `total_squared_<signal>`,
  `average_<signal>`, `variance_<signal>`, `sd_<signal>`, and
  `coefficient_of_variation_<signal>`.

Positional fcoverage bedGraph and per-window positional outputs are rejected by
this API.

Support metadata must match the current writer contract:

- `blacklisted_positions <= span_positions`.
- `eligible_positions = span_positions - blacklisted_positions` where both
  fields are present.
- `nonzero_positions <= eligible_positions`.
- Finite `covered_fraction` values must be in `[0, 1]`.

Scalar totals must be finite and non-negative. Scalar averages must be finite
and non-negative, except `NaN` is allowed for zero-support rows.

Summary totals, squared totals, averages, variances, and standard deviations
must be finite and non-negative, except zero-support rows may contain `NaN`.
Coefficient of variation may be `NaN`, a finite non-negative value, or a writer
display cap such as `>1e6`.

Grouped aggregate outputs must contain at most one row per `group_idx`.
Group-index files must use `group_idx` and `group_name` columns. They must not
contain duplicate `group_idx` values, duplicate `group_name` values, empty names,
or names for missing TSV rows.

`FCoverageFilenameMetadata` reports command-mode hints parsed from canonical
filenames:

- `average`, `total`, and `summary_stats` map to ordinary aggregation.
- `average_on_unique_bases`, `total_on_unique_bases`, and
  `summary_stats_on_unique_bases` map to unique-base aggregation.
- `length_normalized.fcoverage` maps to unit-mass length normalization.
- `length_normalized.restored_mean.fcoverage` maps to restored-mean length
  normalization.
- Renamed files report `Unknown` for fields that are not present in the
  filename.

Selections return `FCoverageSelection`, which is either scalar values or summary
statistics with selected row metadata.

## End-Motif Loader

`load_ends_output(path)` reads `cfdna ends` Zarr stores with the supported
end-motif schema version.

Supported row modes:

- Global.
- Windows from fixed-size or BED windowing.
- Groups from grouped BED inputs.

Supported motif axes:

- Concrete motifs from `motif_ascii`.
- Motif-group labels from JSON attributes on `motif_index`.

Supported storage modes:

- Dense count matrix in `/counts`.
- Sparse COO values in the `/sparse` group.

Dense and sparse counts must be finite and non-negative. Dense count shapes must
match row and motif axes. Sparse coordinates must be sorted, unique,
zero-based, and within the declared sparse shape. Missing in-bounds sparse
coordinates are implicit zero counts.

Group names and motif labels are public selectors and must be unique where they
are loaded. JSON labels must not contain control characters. Concrete labels in
`motif_ascii` must be ASCII.

Sparse data exposes:

- `entries()` for ordered COO iteration.
- `count(row_index, motif_index)` for occasional point lookup by binary search.
- `to_lookup_index()` for repeated point lookup with a RAM-heavier hash index.
- `to_dense_matrix()` for explicit densification.

Selections preserve dense or sparse storage mode. Sparse selections copy only
stored entries whose source row and motif are both selected. Missing selected
cells remain implicit zero counts.

## Reference K-mer Loader

`load_ref_kmers_output(path)` reads `cfdna ref-kmers` Zarr stores with the
supported reference k-mer schema version.

Supported row modes:

- Global.
- Windows from fixed-size or BED windowing.
- Groups from grouped BED inputs.

Supported motif axes:

- Concrete reference k-mers from `motif_ascii`.
- Motif-group labels from JSON attributes on `motif_index`.

Supported storage modes:

- Dense frequency matrix in `/frequencies`.
- Sparse COO values in the `/sparse` group.

Stored values are frequencies. Each row has a `row_scaling_factor`, and counts
are reconstructed as `frequency * row_scaling_factor[row]`.

Dense and sparse frequencies must be finite values in `[0, 1]`.
`row_scaling_factor` values must be finite and non-negative. Dense frequency
shapes must match row and motif axes. Sparse coordinates must be sorted, unique,
zero-based, and within the declared sparse shape. Missing in-bounds sparse
coordinates are implicit zero frequencies.

Concrete reference k-mer labels must match `kmer_size`, contain only A/C/G/T
bases, and match the canonical representation when `canonical = true`.
Motif-group labels are public selectors and must be unique. JSON labels must not
contain control characters.

`all_motifs` means zero-frequency targets are retained on the motif axis. Without
a motifs file, those targets are all A/C/G/T k-mers for the configured `k`. With
a motifs file, they are the motifs or motif groups listed in that file.

Sparse data exposes:

- `entries()` for ordered COO iteration.
- `frequency(row_index, motif_index)` for occasional point lookup by binary
  search.
- `to_lookup_index()` for repeated point lookup with a RAM-heavier hash index.
- `to_dense_matrix()` for explicit frequency densification.

Counts can be reconstructed with `count()`, `sparse_count_entries()`, and
`to_dense_count_matrix()`.

Selections preserve dense or sparse storage mode, copy row and motif metadata,
and preserve requested selector order.

## Midpoints Loader

`load_midpoints_output(path)` reads `cfdna midpoints` profile Zarr stores with
the supported midpoint schema version.

The loader reads metadata eagerly and reads count values only through
`read_all_counts()` or `select().read()`.

Axes:

- Group metadata from the `group` axis and `eligible_intervals`.
- Fragment length bins from `length_bin`, `length_start_bp`, and
  `length_end_bp`.
- Position bins from `position`, `position_bin_start_bp`, and
  `position_bin_end_bp`.

The count array has shape `(group, length_bin, position)` and stores `f32`
values. Count values must be finite. Finite negative values are allowed because
smoothed midpoint profiles can legitimately go below zero. Group, length-bin,
and position axes must be non-empty.

Group names must be unique and must not contain control characters. Length bins
and position bins must be zero-based, contiguous, and valid half-open intervals.

Selections return `MidpointCountSelection` with selected group metadata,
selected length bins, selected position bins, and an owned `DenseArray3<f32>`.
Explicit `positions(&[])` is rejected because midpoint selections need at least
one position bin.

## Out Of Scope

The output loaders do not run commands, infer missing command settings, or load
arbitrarily renamed Zarr directories. They read supported cfDNAlab output
schemas and fail when files do not match those schemas.

External Python and R package loaders are tracked separately. Current
downstream fixture checks live under `downstream_tests/`.
