# Analysis Output Format Roadmap

This note tracks the public output-format cleanup across analysis commands. The
goal is to make outputs self-contained, directly loadable from R and Python, and
harder to misinterpret than raw arrays plus loosely related sidecars.

## Current Direction

- `lengths`: settled as a wide compressed TSV plus settings JSON.
- `midpoints`: first Zarr implementation is in place on the development branch.
  Downstream loading examples and CI-only compatibility tests now exist.
- `ends`: later target. Dense motif arrays and sparse motif arrays should move
  together so users do not have to learn both Zarr and SciPy-style `.npz` for the
  same command family.
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
- Fragment-length bins are count columns.
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
group_name              string[group]
eligible_intervals      uint32[group]
length_start_bp         uint32[length_bin]
length_end_bp           uint32[length_bin]
position_bin_start_bp   int32[position]
position_bin_end_bp     int32[position]
```

The current implementation also writes explicit byte storage for group names as
a temporary compatibility fallback:

```text
group_name_utf8         uint8[group, max_name_bytes]
group_name_nbytes       uint32[group]
```

Release blocker: before the next release, decide whether `group_name[group]`
encoded as Zarr `string` plus `vlen-utf8` is supported well enough by the target
Python and R readers. If yes, remove the byte fallback and document
`group_name` as the only public group-name array. If no, keep the fallback and
provide copy-paste decoding helpers in the R/Python docs. Do not ship the final
format with this undecided.

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
  "dimension_names": ["group", "length_bin", "position"],
  "count_units": "weighted_midpoint_count"
}
```

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

Expected analysis columns after dataframe conversion:

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

Docs must tell users to subset before converting a full 3D array to a dataframe.
A full dataframe can be much larger than the Zarr store.

### Sourceable Helper Scripts

Raw Zarr access is a compatibility layer, not a good user-facing analysis API.
Before turning anything into real Python or R packages, add sourceable midpoint
helper scripts that can load a `.midpoint_profiles.zarr` path and expose the
operations users are most likely to need.

Initial Python helper:

```python
from cfdnalab_midpoints import read_midpoint_profiles

profiles = read_midpoint_profiles("sample.midpoint_profiles.zarr")
group_frame = profiles.dataframe_for_group("CTCF")
length_frame = profiles.dataframe_for_length_bin((120, 180))
group_matrix = profiles.profile_for_group("CTCF")
all_counts = profiles.counts_array()
```

Initial R helper:

```r
source("cfdnalab_midpoints.R")

profiles <- read_midpoint_profiles("sample.midpoint_profiles.zarr")
group_frame <- data_frame_for_group(profiles, "CTCF")
length_frame <- data_frame_for_length_bin(profiles, c(120, 180))
group_matrix <- profile_for_group(profiles, "CTCF")
all_counts <- counts_array(profiles)
```

Required helper operations:

- validate midpoint schema name and version
- expose group metadata, length-bin metadata, and position-bin metadata
- return raw arrays or matrices for group-indexed and length-bin-indexed slices
- build long dataframes for one group across all length/position bins
- build long dataframes for one length bin across all groups/positions
- resolve groups by `group_idx` or `group_name`
- resolve length bins by index or by `(length_start_bp, length_end_bp)`

The first scripts should not include plotting. Dataframe builders are the useful
boundary. Plot helpers can be added later if repeated workflows converge on the
same plot shapes.

These helper scripts should live with downstream tests first and be tested in
the downstream workflow against cfDNAlab-generated Zarr output. Promote them to
actual Python/R packages only after the API has survived real examples and ends
needs the same style of helper layer.

### R Loading

The guide currently shows `Rarr` because it has a direct array-level read path.
The downstream workflow also tests CRAN `zarr`.

The target user workflow should stay close to:

```r
library(Rarr)

store <- "sample.midpoint_profiles.zarr"

counts <- read_zarr_array(file.path(store, "counts"), index = list(1:10, NULL, NULL))
group_name <- read_zarr_array(file.path(store, "group_name"), index = list(1:10))
eligible_intervals <- read_zarr_array(file.path(store, "eligible_intervals"), index = list(1:10))
length_start_bp <- read_zarr_array(file.path(store, "length_start_bp"), index = list(NULL))
length_end_bp <- read_zarr_array(file.path(store, "length_end_bp"), index = list(NULL))
position_bin_start_bp <- read_zarr_array(file.path(store, "position_bin_start_bp"), index = list(NULL))
position_bin_end_bp <- read_zarr_array(file.path(store, "position_bin_end_bp"), index = list(NULL))
```

The design requirement is not a specific R package. The requirement is that an R
user can load counts and axis metadata without Python.

## Ends Later

`cfdna ends` has both dense and sparse public motif outputs. Move both modes to
Zarr together so users do not have to learn `.npy`, SciPy-style `.npz`, motif
text files, `bins.tsv`, and `group_index.tsv` as separate contracts.

The target public package should be:

```text
<prefix>.end_motifs.zarr/
<prefix>.end_settings.json
```

The current settings sidecar is `<prefix>.end_motif_settings.json`. Rename it to
`<prefix>.end_settings.json` during the Zarr transition for consistency with the
other command settings files. The settings JSON should remain human-readable
command provenance. The Zarr store should carry the machine-readable arrays
needed for downstream analysis.

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
motif                  string[motif]
motif_utf8             uint8[motif, motif_byte]
motif_nbytes           uint32[motif]
motif_byte             int32[motif_byte]
```

`motif` is the user-facing label in `<outside>_<inside>` form. Keep the byte
fallback at least until downstream tests decide whether Zarr strings are safe
for all target R/Python readers. Motif labels are ASCII and fixed-width for a
given run, so the byte fallback is cheap and easier to decode than arbitrary
variable-width strings.

Row axis, all modes:

```text
row                    int32[row]
```

`row` is a zero-based coordinate label for the matrix row. It is not the
domain-level key. The domain key depends on `row_mode`.

Global mode:

```text
row_label              string[row]       # value: "global"
```

Fixed-size and BED modes:

```text
chromosome             int32[chromosome]
chromosome_name        string[chromosome]
chromosome_name_utf8   uint8[chromosome, chromosome_name_byte]
chromosome_name_nbytes uint32[chromosome]
chromosome_name_byte   int32[chromosome_name_byte]
row_chromosome         int32[row]
row_start_bp           uint64[row]
row_end_bp             uint64[row]
blacklisted_fraction   float64[row]
```

Use a chromosome dictionary instead of repeating chromosome strings for every
row. `row_chromosome[row]` indexes `chromosome_name[chromosome]`. For BED mode,
rows must preserve the original BED row order after chromosome filtering, just
as current `bins.tsv` is keyed by the preserved output index.

Grouped BED mode:

```text
group                  int32[row]
group_name             string[row]
group_name_utf8        uint8[row, group_name_byte]
group_name_nbytes      uint32[row]
group_name_byte        int32[group_name_byte]
eligible_windows       uint32[row]
blacklisted_fraction   float64[row]
```

`group[row]` must be the same zero-based `group_idx` used by the count rows.
Reject non-contiguous group indices instead of reordering counts blindly. Add
`eligible_windows` during the transition; current ends group metadata does not
include it, but it is needed to interpret grouped counts and keeps the command
consistent with the new `lengths` grouped output.

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
likely dataframe/plotting workflow.

### Sparse Output

Sparse output should be explicit COO inside the Zarr package:

```text
sparse/row              uint64[nnz]
sparse/motif            uint64[nnz]
sparse/count            float64[nnz]
sparse/shape            uint64[sparse_dimension]   # [n_rows, n_motifs]
sparse/dimension        string[sparse_dimension]   # ["row", "motif"]
```

Use zero-based sparse coordinates on disk. R examples must add 1 when
constructing `Matrix::sparseMatrix`. Python examples can construct
`scipy.sparse.coo_matrix((count, (row, motif)), shape=shape)`.

COO entries should be sorted by `(row, motif)` before writing. The current sparse
`.npz` helper preserves the matrix mathematically, but the Zarr transition is an
opportunity to make sparse payloads deterministic and easy to diff. Duplicate
`(row, motif)` entries should be rejected in the final writer because reduction
should already have merged them.

Sparse output should keep only observed motifs, matching current default ends
behavior. If users need a full motif universe, they should use dense
`--all-motifs` unless there is a later concrete use case for "sparse with all
motif labels".

### Downstream Loading Targets

Python dense:

```python
import pandas as pd
import zarr

store = zarr.open_group("sample.end_motifs.zarr", mode="r", zarr_format=3)
counts = store["counts"][:]
motifs = store["motif"][:]
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
Zarr helpers once the second command needs them. Useful shared helpers:

- filesystem store creation
- root metadata writing
- generic Zarr V3 array creation with zstd and dimension names
- single-chunk coordinate/metadata array writing
- checked `usize`/`u64` to public `i32`, `u32`, and `uint64` conversions
- string plus UTF-8 byte fallback writing
- 2D dense chunk-shape selection
- COO sparse matrix writing

Do not build a broad generic "Zarr dataset builder" yet. The midpoint and ends
schemas are different enough that command-specific writers should remain the
readable entry points. The shared code should only cover mechanical Zarr
serialization and low-level dtype validation.

Reuse existing row metadata logic instead of recomputing it inside the Zarr
writer. The current code writes `bins.tsv` from `WindowBinInfo` and grouped rows
through `write_group_index_with_blacklist_tsv`. For Zarr, split the data
construction from TSV writing so both TSV and Zarr paths can consume the same
row summaries. This avoids the grouped-row ordering and blacklist-fraction drift
risks already encountered in the midpoint/length output work.

The sparse COO writer can later be reused by fragment-kmer positional sparse
outputs and any remaining SciPy `.npz` sparse outputs, but do not migrate those
commands as part of the first ends Zarr transition.

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
- Whether `group_name` should be a string array or explicit UTF-8 byte matrix.
- Whether the package should include a small `README.txt` with the schema version
  and loading pointers.
- Whether any optional export should produce a long TSV for small selected
  profiles. This should not replace the primary Zarr output.

## Open Decisions For Ends

- Whether the output store should include only one representation per run, or
  allow both dense and sparse groups in the same store when `--all-motifs` is
  requested. Prefer one representation per run unless a user need appears.
- Whether `motif` and row label strings can rely on Zarr string arrays, or need
  documented byte-fallback decoding just like midpoint group names.
- Whether chromosome row metadata should use a dictionary (`row_chromosome` plus
  `chromosome_name`) or repeat `chrom` as a string array per row. Prefer the
  dictionary unless downstream examples become too clumsy.
- Whether `blacklisted_fraction` should be present for every row mode with
  all-zero values when no blacklist was used, or only when blacklists were
  supplied. Prefer always present for ordinary window and grouped modes.
- Whether grouped ends should add `eligible_windows` now. The answer should
  probably be yes, but this is a small behavior addition compared with the
  current `group_index.tsv`.
- Whether sparse coordinates should be stored as `uint64` or `int64`. `uint64`
  is semantically correct, but some R sparse constructors want signed integer
  vectors after loading.
- Exact dense and sparse chunk shapes.
- Whether to preserve `end_motifs` in the store name or rename to `end_counts`.
  Prefer `end_motifs.zarr` because it matches existing user-facing terminology.

## Implementation Order

Midpoints:

1. Implement midpoint Zarr writing behind the normal public output path. Done.
2. Add website loading examples for Python and R. Done.
3. Add CI-only downstream compatibility checks. Done as an initial workflow.
4. Run the downstream workflow and decide whether native `group_name` is safe
   enough to keep without the byte fallback.
5. Confirm Zarr V3 subset, codec, chunk layout, and group-name encoding before
   release.
6. Add sourceable Python and R midpoint helper scripts, then update the guide to
   use one helper-based example per language.

Ends:

1. Split ends row metadata construction from TSV writing so Zarr and optional
   TSV outputs cannot drift.
2. Extract only the low-level Zarr helpers already duplicated by midpoint and
   needed by ends.
3. Implement dense `counts[row, motif]` Zarr output for `--all-motifs`.
4. Implement sparse COO Zarr output for the default path.
5. Add ends fixtures to the downstream compatibility workflow for Python
   `zarr`, Python sparse reconstruction, R `Rarr`, and R `Matrix`.
6. Update the ends guide with the sparse loading example first, and dense only
   if it remains useful for users.
7. Distill the final ends schema into `.AI/docs/specs/ends_spec.md`.

The downstream compatibility workflow is tracked separately in
`downstream_testing.md`.

## Acceptance Criteria

The midpoint transition is successful when:

- The output is one self-contained package.
- Counts open with named axes in Python.
- Counts and coordinate arrays open in R without Python.
- Users can construct a dataframe with group names, eligible intervals, length
  bins, position bins, and counts.
- Users do not have to manually reconcile `.npy`, `group_index.tsv`, and settings
  JSON files to trust the result.
