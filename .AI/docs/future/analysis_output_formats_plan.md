# Analysis Output Format Roadmap

This note tracks the public output-format cleanup across analysis commands. The
goal is to make outputs self-contained, directly loadable from R and Python, and
harder to misinterpret than raw arrays plus loosely related sidecars.

## Current Direction

- `lengths`: settled as a wide compressed TSV plus settings JSON.
- `midpoints`: first Zarr implementation is in place on the development branch.
  The next step is downstream loading examples in R and Python.
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
<prefix>.length_counts.tsv.gz
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

### Python Loading Sketch

The best Python experience should be xarray-first:

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

### R Loading Sketch

R support should be tested before finalizing the exact Zarr version, string
encoding, and codec. The target user workflow should look like:

```r
library(zarr)

profiles <- open_zarr("sample.midpoint_profiles.zarr", mode = "r")

counts <- profiles[["counts"]][1:10, , ]
group_name <- profiles[["group_name"]][1:10]
eligible_intervals <- profiles[["eligible_intervals"]][1:10]
length_start_bp <- profiles[["length_start_bp"]][]
length_end_bp <- profiles[["length_end_bp"]][]
position_bin_start_bp <- profiles[["position_bin_start_bp"]][]
position_bin_end_bp <- profiles[["position_bin_end_bp"]][]
```

If the final R reader is `Rarr` or `ZarrArray` instead of `zarr`, adjust the docs
to that package. The design requirement is not a specific R package. The
requirement is that an R user can load counts and axis metadata without Python.

## Ends Later

`cfdna ends` has both dense and sparse public motif outputs. If it moves to Zarr,
move both modes under one conceptual format.

Dense output:

```text
counts                  float64[row, motif]
motif                  string[motif]
```

Sparse output should be explicit COO inside the Zarr package:

```text
sparse/row              uint64[nnz]
sparse/col              uint64[nnz]
sparse/data             float64[nnz]
sparse/shape            uint64[2]
motif                   string[motif]
```

Use zero-based sparse coordinates on disk. R examples should add 1 when
constructing `Matrix::sparseMatrix`.

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

- Exact safe Zarr V3 subset.
- Codec and chunk layout.
- Exact dimension metadata needed for `xarray.open_zarr()`.
- Whether `counts` should stay `float32` on disk or move to `float64`.
- Whether `group_name` should be a string array or explicit UTF-8 byte matrix.
- Whether the settings should live entirely in root attrs or as a JSON object/file
  inside the package.
- Whether the package should include a small `README.txt` with the schema version
  and loading pointers.
- Whether any optional export should produce a long TSV for small selected
  profiles. This should not replace the primary Zarr output.

## Implementation Order

1. Run a compatibility spike with a tiny midpoint Zarr package.
2. Verify loading in Python with `zarr` and `xarray`.
3. Verify loading in R with the current best Zarr reader.
4. Confirm Zarr V3 subset, codec, chunk layout, and group-name encoding.
5. Implement midpoint Zarr writing behind the normal public output path.
6. Update command docs and website loading examples.
7. Only then revisit `ends`.

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
