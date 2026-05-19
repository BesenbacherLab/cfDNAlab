# Analysis Output Format Roadmap

This note tracks the public output-format cleanup across analysis commands. The
goal is to make outputs self-contained, directly loadable from R and Python, and
harder to misinterpret than raw arrays plus loosely related sidecars.

## Current Direction

- `lengths`: settled as a wide compressed TSV plus settings JSON.
- `midpoints`: first Zarr implementation is in place on the development branch.
  Downstream loading examples and CI-only compatibility tests now exist.
- `ends`: dense and sparse motif outputs now use one Zarr schema on the
  development branch.
- Commands without a human-readable settings JSON should get one for
  consistency.
- Internal temporary files can keep their current formats unless a separate
  performance or correctness reason appears.

## Command Settings JSON

Every public analysis command should write a plain JSON settings file. This file
is for human-readable command provenance and interpretation. Users should not
need to open a Zarr store, inspect a binary package, or parse a log file just to
understand the major settings that produced an output.

Recommended naming:

```text
<prefix>.<command>_settings.json
```

Examples:

```text
<prefix>.length_settings.json
<prefix>.midpoint_settings.json
<prefix>.end_settings.json
```

The JSON should record command-level settings, selected modes, key filters,
normalization flags, and enough schema information to understand the primary
output. It can duplicate axis summaries that also exist as TSV columns or Zarr
coordinate arrays. The machine-readable primary output remains authoritative for
the actual data and coordinates.

## Lengths Status

`cfdna lengths` is no longer the main design problem. Its main output is now:

```text
<prefix>.length_counts.tsv.zst
<prefix>.length_settings.json
```

The lasting contract belongs in `.AI/docs/specs/lengths_spec.md`. In short:

- Output rows are global, genomic windows, BED rows, or grouped BED groups.
- Non-grouped rows are keyed by `chrom`, `start`, and `end`.
- Grouped rows are keyed by `group_name` and include `eligible_windows`.
- Fragment length bins are count columns.
- Single-bp bins use `count_<length>`.
- Wider half-open bins use `count_<start>_<end>`.
- `--decimals` controls count-column formatting only.
- `blacklisted_fraction` is written only when a blacklist was used and is
  rounded to three decimals.

The former long-form TSV and Zarr alternatives are closed for `lengths` unless a
new use case changes the requirements.

## Midpoints Design Target

`cfdna midpoints` now writes on the development branch:

```text
<prefix>.midpoint_profiles.zarr/
<prefix>.group_index.tsv
<prefix>.midpoint_settings.json
```

The settings file remains a plain, human-readable JSON file. Users should be
able to inspect command settings and provenance without opening the Zarr store.

The main array has shape:

```text
counts[group, length_bin, position]
```

This is a better Zarr candidate than `lengths` because the output is genuinely
multidimensional and because axis metadata is required for correct analysis.

### Recommended Public Package

Prefer one primary Zarr output plus one human-readable settings file:

```text
<prefix>.midpoint_profiles.zarr/
<prefix>.midpoint_settings.json
```

Core arrays:

```text
counts                  float32[group, length_bin, position]
group                   int32[group]
eligible_intervals      int32[group]
length_start_bp         int32[length_bin]
length_end_bp           int32[length_bin]
position_bin_start_bp   int32[position]
position_bin_end_bp     int32[position]
```

Human-readable group names are stored as JSON attributes on the `group` array:

```json
{
  "label_field": "group_name",
  "labels": ["group_a", "group_b"]
}
```

This keeps the Zarr arrays numeric-only for the first public schema and avoids
depending on Zarr V3 string-array support in every R reader. Downstream helper
packages should expose these labels as normal group metadata.

`group` is the same zero-based `group_idx` used by `<prefix>.group_index.tsv`
and by axis 0 of `counts`. The writer must reject non-contiguous group indices
because otherwise metadata rows could describe the wrong count rows.

The Zarr package should be self-contained for downstream analysis. The
plain-text `group_index.tsv` can remain as an auxiliary human-readable summary,
but users should not need it for normal plotting once the Zarr reader examples
are settled.

### Root Metadata

Root attrs should record the interpretation contract:

```json
{
  "cfdnalab_schema": "midpoint_profiles",
  "cfdnalab_schema_version": 1,
  "primary_array": "counts",
  "count_units": "weighted_midpoint_count"
}
```

The count array itself carries native Zarr V3 `dimension_names`. The root attrs
should not duplicate them.

Do not move general command settings into Zarr attrs. Keep them in
`<prefix>.midpoint_settings.json` so the run can be inspected without any Zarr
reader. The Zarr store should carry the arrays and small attrs needed for
downstream analysis. The JSON file should carry human-readable command settings,
provenance, and interpretation notes. Some duplication is acceptable when it
improves human readability, but the Zarr arrays remain the machine-readable
source for axis coordinates.

The store should also be written with the dimension metadata expected by the
target Python reader. For xarray, a raw Zarr group is not enough. Each array must
carry valid dimension metadata so `xr.open_zarr()` can construct a labeled
dataset instead of failing or exposing anonymous axes.

### Position Axis

The position axis should use half-open interval-relative bins:

```text
position_bin_start_bp[position]
position_bin_end_bp[position]
```

For `--bin-size 1`, this is equivalent to base offsets:

```text
position_bin_start_bp = 0, 1, 2, ...
position_bin_end_bp   = 1, 2, 3, ...
```

For larger `--bin-size`, it records the actual output bins, including a shorter
last bin. This avoids hiding the current "last bin uses actual width" rule in a
settings object that users may not inspect.

### Length Axis

Use the same half-open semantics as `lengths`:

```text
length_start_bp[length_bin]
length_end_bp[length_bin]
```

Do not encode midpoint values for length bins. Plotting midpoints are derived
display values, not part of the count definition.

### Dense Or Sparse Public Output

The current final public midpoint array is dense after merging and smoothing.
Keep the first Zarr transition dense unless real outputs prove that the dense
array is too large.

Sparse tile partials are an implementation detail. Do not expose sparse public
midpoint output in v1 of this transition unless there is a concrete user need,
because smoothed profiles are generally dense anyway.

### Python Loading

The guide currently shows direct `zarr` loading because it is explicit and close
to the on-disk schema. The downstream workflow also tests `xarray` and
`dask.array`.

The best high-level Python experience should still be xarray-first once the
package behavior is confirmed in CI:

```python
import xarray as xr

profiles = xr.open_zarr("sample.midpoint_profiles.zarr")

table = (
    profiles["counts"]
    .sel(group=slice(None))
    .to_dataframe(name="count")
    .reset_index()
)
```

Expected analysis columns after data frame conversion:

```text
group
length_bin
position
count
group_name
eligible_intervals
length_start_bp
length_end_bp
position_bin_start_bp
position_bin_end_bp
```

Docs must tell users to subset before converting a full 3D array to a data frame.
A full data frame can be much larger than the Zarr store.

### Helper Packages

Raw Zarr access is a compatibility layer, not a good user-facing analysis API.
The Python and R packages now provide the user-facing analysis API for
midpoints, end motifs, and lengths. The current loader API belongs in
[`../specs/package_loader_api.md`](../specs/package_loader_api.md); this
roadmap should not duplicate every method signature.

Representative Python usage:

```python
import cfdnalab as cfl

profiles = cfl.read_midpoints("sample.midpoint_profiles.zarr")
profile = profiles.data_frame(group_idxs=0, length_bin_idxs=0)

ends = cfl.read_end_motifs("sample.end_motifs.zarr")
motif_counts = ends.data_frame(motifs="_AA")

lengths = cfl.read_lengths("sample.length_counts.tsv.zst")
length_distribution = lengths.data_frame(value="fraction")
```

Representative R usage:

```r
library(cfdnalab)

profiles <- read_midpoints("sample.midpoint_profiles.zarr")
profile <- midpoint_data_frame(profiles, group_idxs = 1L, length_bin_idxs = 1L)

ends <- read_end_motifs("sample.end_motifs.zarr")
motif_counts <- end_motif_data_frame(ends, motifs = "_AA")

lengths <- read_lengths("sample.length_counts.tsv.zst")
length_distribution <- length_data_frame(lengths, value = "fraction")
```

Required helper operations:

- validate schema names and supported schema versions
- expose domain metadata for groups, windows, motifs, length bins, and positions
- return arrays, vectors, and sparse or dense matrices with explicit
  materialization semantics
- build data frames with mode-specific selectors and no generic public `row`
  vocabulary
- handle sparse end-motif stores without silently densifying by default

The first package helpers should not include plotting. data frame builders are
the useful boundary. Plot helpers can be added later if repeated workflows
converge on the same plot shapes.

Keep downstream tests against cfDNAlab-generated output so the packages cannot
drift away from the Rust schema.

### R Loading

The guide currently shows `Rarr` because it has a direct array-level read path.
The downstream workflow also tests CRAN `zarr`.

The target user workflow should stay close to:

```r
library(Rarr)
library(jsonlite)

store <- "sample.midpoint_profiles.zarr"

counts <- read_zarr_array(file.path(store, "counts"), index = list(1:10, NULL, NULL))
group_labels <- fromJSON(file.path(store, "group", "zarr.json"))$attributes$labels
eligible_intervals <- read_zarr_array(file.path(store, "eligible_intervals"), index = list(1:10))
length_start_bp <- read_zarr_array(file.path(store, "length_start_bp"), index = list(NULL))
length_end_bp <- read_zarr_array(file.path(store, "length_end_bp"), index = list(NULL))
position_bin_start_bp <- read_zarr_array(file.path(store, "position_bin_start_bp"), index = list(NULL))
position_bin_end_bp <- read_zarr_array(file.path(store, "position_bin_end_bp"), index = list(NULL))
```

The design requirement is not a specific R package. The requirement is that an R
user can load counts and axis metadata without Python.

## Ends Status

`cfdna ends` has both dense and sparse public motif outputs. Both modes now use
Zarr so users do not have to learn `.npy`, SciPy-style `.npz`, motif text files,
`bins.tsv`, and `group_index.tsv` as separate contracts.

The target public package should be:

```text
<prefix>.end_motifs.zarr/
<prefix>.end_settings.json
```

The settings sidecar is `<prefix>.end_settings.json` for consistency with the
other command settings files. The settings JSON remains human-readable command
provenance. The Zarr store carries the machine-readable arrays needed for
downstream analysis.

Root attrs:

```json
{
  "cfdnalab_schema": "end_motif_counts",
  "cfdnalab_schema_version": 1,
  "storage_mode": "dense",
  "row_mode": "global|size|bed|grouped_bed",
  "primary_array": "counts",
  "count_units": "weighted_end_motif_count"
}
```

For sparse output, `storage_mode` should be `"sparse_coo"` and `primary_array`
should be omitted or replaced with:

```json
{
  "primary_group": "sparse",
  "sparse_format": "coo",
  "sparse_indices_base": 0
}
```

Do not store general command options only as Zarr attrs. Keep motif-definition,
filtering, correction flags, dense/sparse mode, row mode, and other provenance
in `<prefix>.end_settings.json`.

### Shared Axes

Both dense and sparse stores should carry the same row and motif metadata.

Motif axis:

```text
motif_index            int32[motif]
motif_byte             int32[motif_byte]
motif_ascii            uint8[motif, motif_byte]
```

Motif labels are fixed-width ASCII for a given run because every label has
`k_outside` bases, one underscore, and `k_inside` bases. Do not store large
motif label sets as JSON attrs. JSON attrs are loaded as metadata and become a
bad fit when dense outputs can enumerate thousands or millions of motifs.

`motif_ascii[motif, motif_byte]` stores one row per motif label. No
`motif_nbytes` array is needed because all motif labels have the same width in a
run. Downstream helpers decode each row as ASCII and expose normal string labels
through their public APIs. This keeps labels chunkable, compressible as array
payload, and independent of Zarr string dtype support.

Row axis, all modes:

```text
row                    int32[row]
```

`row` is a zero-based coordinate label for the matrix row. It is not the
domain-level key. The domain key depends on `row_mode`.

Global mode:

`row` has `label_field = "row_label"` and `labels = ["global"]`.

Fixed-size and BED modes:

```text
chromosome             int32[chromosome]
row_chromosome         int32[row]
row_start_bp           int64[row]
row_end_bp             int64[row]
blacklisted_fraction   float64[row]
```

Use a chromosome dictionary instead of repeating chromosome strings for every
row. `row_chromosome[row]` indexes `chromosome[chromosome]`, whose JSON attrs
store `label_field = "chromosome_name"` and the chromosome-name labels. For BED
mode, rows must preserve the original BED row order after chromosome filtering,
just as current `bins.tsv` is keyed by the preserved output index.

Grouped BED mode:

```text
group                  int32[row]
eligible_windows       int32[row]
blacklisted_fraction   float64[row]
```

`group[row]` must be the same zero-based `group_idx` used by the count rows.
Group names are stored as JSON attrs on `group` with `label_field =
"group_name"` and the label array. Reject non-contiguous group indices instead
of reordering counts blindly. `eligible_windows` is part of grouped row
metadata because grouped counts need the same basic support context as grouped
length outputs.

### Dense Output

Dense output:

```text
counts                  float64[row, motif]
```

Dense output corresponds to current `--all-motifs`. It should include all motif
labels in deterministic order, including motifs with zero counts. The dense size
guard should remain, but update the error message to mention Zarr rather than
`.npy`.

Chunking should favor common reads:

- small dense matrices: one chunk
- larger matrices: chunk by rows while keeping the motif axis whole when
  possible
- target chunk size should be in the same rough range as midpoint counts
  chunks, adjusted for `float64`

This keeps "read a subset of rows with all motifs" efficient, which is the
likely data frame and plotting workflow.

### Sparse Output

Sparse output should be explicit COO inside the Zarr package:

```text
sparse/row              int32[nnz]
sparse/motif            int32[nnz]
sparse/count            float64[nnz]
sparse/shape            int32[sparse_dimension]    # [n_rows, n_motifs]
sparse/sparse_dimension int32[sparse_dimension]    # labels attr: ["row", "motif"]
```

Use zero-based sparse coordinates on disk. R examples must add 1 when
constructing `Matrix::sparseMatrix`. Python examples can construct
`scipy.sparse.coo_matrix((count, (row, motif)), shape=shape)`.

COO entries should be sorted by `(row, motif)` before writing. The current sparse
`.npz` helper preserves the matrix mathematically, but the Zarr transition is an
opportunity to make sparse payloads deterministic and easy to diff. Duplicate
`(row, motif)` entries should be rejected in the final writer because reduction
should already have merged them.

Sparse output keeps the observed motif axis written by the command. Dense
`--all-motifs` remains the path for users who need the full motif universe.

### Downstream Loading Targets

Python dense:

```python
import pandas as pd
import zarr

store = zarr.open_group("sample.end_motifs.zarr", mode="r", zarr_format=3)
counts = store["counts"][:]
motifs = ["".join(map(chr, row)) for row in store["motif_ascii"][:]]
```

Python sparse:

```python
import scipy.sparse
import zarr

store = zarr.open_group("sample.end_motifs.zarr", mode="r", zarr_format=3)
sparse = store["sparse"]
mat = scipy.sparse.coo_matrix(
    (sparse["count"][:], (sparse["row"][:], sparse["motif"][:])),
    shape=tuple(sparse["shape"][:]),
)
```

R sparse:

```r
library(Matrix)
library(Rarr)

store <- "sample.end_motifs.zarr"
row <- read_zarr_array(file.path(store, "sparse", "row"), index = list(NULL))
motif <- read_zarr_array(file.path(store, "sparse", "motif"), index = list(NULL))
count <- read_zarr_array(file.path(store, "sparse", "count"), index = list(NULL))
shape <- read_zarr_array(file.path(store, "sparse", "shape"), index = list(NULL))

mat <- sparseMatrix(
  i = as.integer(row) + 1L,
  j = as.integer(motif) + 1L,
  x = count,
  dims = as.integer(shape)
)
```

Docs should show one dense and one sparse example only if both are common. If one
example is preferred, show sparse because it is the default output path and the
one most likely to trip users up.

### Ends Implementation Notes

Keep dense and sparse Zarr writers in the ends command, but extract small shared
Zarr helpers once the current schema has settled. Useful shared helpers:

- filesystem store creation
- root metadata writing
- generic Zarr V3 array creation with zstd and dimension names
- single-chunk coordinate/metadata array writing
- checked `usize`/`u64` to public `i32`, `int64`, and `uint32` conversions
- public label validation for labels stored in JSON attrs
- fixed-width ASCII label-array writing for high-cardinality motif labels
- 2D dense chunk-shape selection
- COO sparse matrix writing

Do not build a broad generic "Zarr dataset builder" yet. The midpoint and ends
schemas are different enough that command-specific writers should remain the
readable entry points. The shared code should only cover mechanical Zarr
serialization and low-level dtype validation.

Current reuse assessment after midpoint and ends:

- Extract now-ish: store creation, generic V3 array creation, single-chunk
  metadata arrays, JSON attribute validation, integer narrowing, and
  public label validation. These are duplicated and testable once.
- Extract with the motif-label change: fixed-width ASCII label-array writing and
  decoding tests. Motifs are high-cardinality enough that they should not share
  the JSON-attrs path used for group and chromosome labels.
- Keep command-specific: root schema attrs, axis names, row metadata
  construction, group ordering validation, chromosome dictionaries, and
  midpoint/ends chunk policy.
- Keep manual: sparse COO layout. Zarr itself does not define a sparse matrix
  schema, so cfDNAlab must own `sparse/row`, `sparse/motif`, `sparse/count`,
  `sparse/shape`, ordering, and duplicate rejection.
- Let `zarrs` own: Zarr V3 metadata format, dimension-name metadata,
  compression, chunk serialization, and boundary chunk handling.

Reuse existing row metadata logic instead of recomputing it inside the Zarr
writer. The current code writes `bins.tsv` from `WindowBinInfo` and grouped rows
through `write_group_index_with_blacklist_tsv`. For Zarr, split the data
construction from TSV writing so both TSV and Zarr paths can consume the same
row summaries. This avoids the grouped-row ordering and blacklist-fraction drift
risks already encountered in the midpoint/length output work.

The sparse COO writer can later be reused by fragment-kmer positional sparse
outputs and any remaining SciPy `.npz` sparse outputs, but do not migrate those
commands as part of the first ends Zarr transition.

## Zarr Schema Cleanup Before Release

Write a concise current schema spec before the next release, either as
`.AI/docs/specs/zarr_schema.md` or split by command specs if that stays clearer.
The spec must list, for each public Zarr store:

- root attributes, their types, and allowed values
- array names, shapes, dtypes, dimension names, and array attributes
- which arrays are coordinates and which arrays are payloads
- the R conversion strategy for integer arrays outside native R `integer`

R needs explicit dtype guidance because it has no native 64-bit integer vector.
The R helper should convert on purpose instead of inheriting whatever a reader
package does:

- `int64` genomic coordinates should remain exact, either as `bit64::integer64`
  or as checked numeric values.
- Sparse coordinates are `int32` because they must become ordinary signed
  integer vectors when constructing `Matrix::sparseMatrix`, which uses one-based
  signed indices.
- Public small non-negative metadata should use `int32` when it is intended to
  become an ordinary R integer.

Root metadata should not duplicate native Zarr V3 dimension metadata. Midpoint
previously had a root-level `dimension_names` copy while each array also
carried native V3 dimension names. The root copy has been dropped; treat the
array metadata as the source of truth.

Root metadata should use symmetric discovery keys. For ends, write both
`primary_array` and `primary_group`, with `null` for the inactive representation,
or choose one explicit object such as:

```json
{
  "primary_data": {
    "storage_mode": "sparse_coo",
    "group": "sparse"
  }
}
```

The chosen form should be used consistently across dense and sparse stores so
naive readers do not need special-case missing keys.

Python and R loaders should require both the public `.zarr` suffix and the
expected schema attrs. The suffix is part of the cfDNAlab output contract, while
the schema attrs prove that the directory is the expected cfDNAlab Zarr store.
We do not need to support users renaming output directories.

Use the same required-array probing strategy in every helper. Prefer the ends
style that checks `store["path/to/array"]` and catches lookup errors, because it
works for nested arrays such as `sparse/row`. Do not use `array_keys()` for
required arrays if nested paths are possible.

Sparse arrays should validate their own dimension names. Ends already validates
dense `counts` dimension names; add equivalent checks that `sparse/row`,
`sparse/motif`, and `sparse/count` use the `nnz` dimension and that
`sparse/shape` and `sparse/sparse_dimension` use the `sparse_dimension`
dimension.

Document schema-version compatibility before the first helper-package release.
Exact-match version checks are fine for schema version 1, but future loaders
should have an explicit `MIN_SUPPORTED_SCHEMA_VERSION` and
`MAX_SUPPORTED_SCHEMA_VERSION`, or a written decision that helper versions are
strictly tied to schema versions.

Python helper cleanup before release:

- `MidpointProfiles.length_bin_idx(length)` resolves a fragment length to a
  length-bin index without being confused with the stored `length_bin` axis.
  Done.
- Expose a chunked/lazy handle for dense ends output, for example
  `dense_counts_zarr_array`, so users who need chunked iteration do not have to
  reach into private fields. Done.
- Add docstring warnings to `dense_*` ends methods that they may load or
  reconstruct dense data. Done.
- Optimize sparse end-motif slices by filtering stored COO arrays directly
  instead of converting the full matrix through CSR. Done.
- Route end-motif tabular output through mode-specific `.data_frame(...)`
  methods, with sparse stores returning non-zero rows by default and
  `densify=True` adding explicit zero-count rows when requested. Done.
- Add one Python fixture for `storage_mode == "dense"` and `row_mode ==
  "global"` so the helper tests cover the storage-mode by row-mode cross product
  used by the public schema. Done.

## Shared Rules

- Simple two-dimensional tables should stay TSV when TSV is compact enough.
- Multidimensional outputs should carry named dimensions and coordinate arrays.
- Do not expose row indices as the main public key when domain keys exist.
- Use `group_name` as the normal grouped-output key.
- Use genomic coordinates as the normal non-grouped window key.
- Store half-open coordinate arrays explicitly in Zarr.
- Put settings and provenance in the output package, but do not make attrs the
  only place where important axis coordinates live.
- Prefer one primary machine-readable output per command.
- Keep temporary sparse formats separate from public output contracts.

## Open Decisions For Midpoints

- Whether downstream CI confirms the current Zarr V3 safe subset across the
  target Python and R package set.
- Whether zstd remains acceptable across the target R readers.
- Whether `counts` should stay `float32` on disk or move to `float64`.
- Whether JSON-attribute label arrays are acceptable for long-term R/Python
  downstream use, or whether a future schema needs separate label sidecars.
- Whether the package should include a small `README.txt` with the schema version
  and loading pointers.
- Whether any optional export should produce a long TSV for small selected
  profiles. This should not replace the primary Zarr output.

## Open Decisions For Ends

- Whether the output store should include only one representation per run, or
  allow both dense and sparse groups in the same store when `--all-motifs` is
  requested. Prefer one representation per run unless a user need appears.
- Whether chromosome row metadata should keep the current dictionary
  (`row_chromosome` plus labels on `chromosome`) or repeat `chrom` per row in a
  helper data frame. Prefer the dictionary on disk unless downstream examples
  become too clumsy.
- `blacklisted_fraction` may be absent when no blacklist was used. Helper
  packages should keep all rows at the default cutoff and error clearly for
  stricter cutoffs when the column is unavailable.
- Grouped ends include `eligible_windows` in row metadata.
- Whether sparse coordinate `int32` limits are acceptable for all practical
  public sparse stores. Current decision: use `int32` to avoid noisy and fragile
  64-bit handling in R sparse workflows.
- Exact dense and sparse chunk shapes.
- Whether to preserve `end_motifs` in the store name or rename to `end_counts`.
  Prefer `end_motifs.zarr` because it matches existing user-facing terminology.

## Python/R Loader API Harmonization

Current cross-language loader decisions:

- Use per-schema supported-version ranges in both helper packages. A schema bump
  in one output type must not make another loader accept a schema it has not
  been updated to read.
- Use mode-specific data frame methods. Python uses object methods named
  `.data_frame(...)`; R uses command-prefixed functions such as
  `midpoint_data_frame()`, `end_motif_data_frame()`, and
  `length_data_frame()`.
- Use plural selector names for vector-capable selectors.
- Use `densify` only for sparse end-motif data frames. Sparse data frames return
  stored non-zero rows by default; densified data frames add explicit
  zero-count rows for selected observed motifs.
- Add Python `repr()` summaries matching the intent of R `print()` methods so
  interactive users can see schema version, shape, storage mode, and row mode
  without inspecting internals.
- Tighten Python scalar input validation where R already rejects ambiguous
  values, especially `length_bin_idx(length=...)` with non-integer numeric
  input.

Public API decisions to make deliberately:

- Decide whether `motifs()` returns just labels in both languages or metadata
  tables in both languages. Current state: Python returns `list[str]` and has
  `motif_metadata()`, while R returns a metadata `data.frame`.
- Decide whether sparse matrix naming should be harmonized. Current state:
  Python exposes SciPy COO-specific helpers, while R exposes `Matrix` sparse
  matrices and sparse data frames without promising a COO class.
- Decide whether R should add selected dense/sparse slice convenience helpers
  that match R idioms. Python can keep its broader method surface because method
  discovery is more natural there.

## Implementation Order

Midpoints:

1. Implement midpoint Zarr writing behind the normal public output path. Done.
2. Add website loading examples for Python and R. Done.
3. Add CI-only downstream compatibility checks. Done as an initial workflow.
4. Validate the JSON-attribute label design in downstream tests before release.
5. Confirm Zarr V3 subset, codec, chunk layout, and label-attribute handling
   before release.
6. Add Python and R helper packages, then update the guide to use one
   helper-based example per language.
7. Apply the Zarr schema cleanup items: remove duplicate root dimension names,
   keep strict suffix and schema validation, unify required-array probing,
   and document dtype conversion for R. Done.

Ends:

1. Split ends row metadata construction from TSV writing so Zarr and optional
   TSV outputs cannot drift. Done.
2. Extract only the low-level Zarr helpers already duplicated by midpoint and
   needed by ends. Done.
3. Implement dense `counts[row, motif]` Zarr output for `--all-motifs`. Done.
4. Implement sparse COO Zarr output for the default path. Done.
5. Add ends fixtures to the downstream compatibility workflow for Python
   `zarr`, Python sparse reconstruction, R `Rarr`, and R `Matrix`.
6. Move motif labels from JSON attrs to fixed-width `motif_ascii` arrays and
   update Python/R helper decoding. Done for the Rust writer and Python helper.
7. Refactor the Python ends helper into mode-specific classes returned by the
   loader. Done.
8. Update the ends guide with the sparse loading example first, and dense only
   if it remains useful for users.
9. Apply the Zarr schema cleanup items: symmetric primary-data attrs, strict
   suffix and schema validation, nested required-array probing, sparse
   dimension-name validation, direct sparse motif slicing, and dense-method
   docstring warnings. Done.
10. Distill the final ends schema into `.AI/docs/specs/ends_spec.md` after the
   motif-label and helper API decisions are implemented. Done.

The downstream compatibility workflow is tracked separately in
`downstream_testing.md`.

## Acceptance Criteria

The midpoint transition is successful when:

- The output is one self-contained package.
- Counts open with named axes in Python.
- Counts and coordinate arrays open in R without Python.
- Users can construct a data frame with group names, eligible intervals, length
  bins, position bins, and counts.
- Users do not have to manually reconcile `.npy`, `group_index.tsv`, and settings
  JSON files to trust the result.

The ends transition is successful when:

- Dense and sparse outputs share one Zarr schema for row and motif metadata.
- Sparse COO reconstructs the expected count matrix in Python and R.
- Motif labels are readable from fixed-width byte arrays without Zarr string
  dtype support or huge JSON metadata.
- Chromosome, sparse-dimension, row-label, and group labels are readable from
  JSON attributes without Zarr string-array support.
- Users can build window-indexed, group-indexed, and motif-indexed data frames
  from the Python helper without manually joining sidecar files.
