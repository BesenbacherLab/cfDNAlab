# Scaling Factors Spec

## Purpose

Record the final design choices behind the scaling-factor commands.

## Command Split

We keep two explicit commands:

- `cfdna coverage-weights`
- `cfdna fragment-count-weights`

This makes the modeling choice visible in workflows instead of hiding it behind
a mode flag.

## Counting Model

Both commands reuse internal `fcoverage` as the raw counting step.

They run internal:

- `fcoverage --by-size <stride> --per-window average`

and then read the resulting final by-size TSV back in.

The only difference between the two commands is whether internal `fcoverage`
uses `--normalize-by-length`:

- `coverage-weights`: `normalize_by_length = false`
- `fragment-count-weights`: `normalize_by_length = true`

This keeps all difficult fragment-counting logic in one place:

- fragment reconstruction
- unpaired handling
- mate-gap handling
- deletion / ref-skip handling
- GC-file weighting
- GC-tag weighting
- blacklist masking

## Why Internal `fcoverage`

The main technical reason for this design is RAM use with GC correction.

`fcoverage` already uses tiled counting with:

- tile core / fetch spans
- tile-local reference loading
- tile-local GC work
- per-tile temp outputs and reducers

Reusing internal `fcoverage` therefore solves the large-memory chromosome-wide
counting problem without reimplementing the scientific logic.

The `*-weights` commands only add:

1. read one value per stride bin from the final by-size TSV
2. apply triangular smoothing across stride bins
3. normalize and invert into scaling factors
4. write the final scaling TSV

## Output Files

The final output files are:

- `coverage-weights`: `<prefix>.coverage.scaling_factors.tsv`
- `fragment-count-weights`: `<prefix>.fragment_counts.scaling_factors.tsv`

Intermediate `fcoverage` output is written to a temp directory under the
user-chosen output directory, not to a system temp location.

## Smoothing

The raw by-size values are interpreted as stride-bin values.

Triangular smoothing is then applied across neighboring stride bins at the
chromosome level, after reading the full by-size output back in.

This preserves the existing smoothing logic while keeping counting tiled.

## Fragment Counts

`fragment-count-weights` is intended to reflect local fragment counts rather
than base-weighted coverage.

Strictly speaking this is still an approximation because fragments that overlap
multiple stride bins contribute partly to each of them, but in sufficiently
large bins the approximation error is tiny.

## Shared Config

The two commands share one internal base config:

- IO args
- unpaired args
- output prefix
- chromosomes
- fragment length filters
- mapping-quality filter
- proper-pair requirement
- blacklist args
- GC args
- reference 2bit path
- bin size / stride

Each public command still has its own top-level config docstring so the CLI help
can explain the intended use clearly.

## GC Mode Metadata

Scaling-factor files write `gc_mode` metadata.

This metadata is only used for:

- compatibility checking
- warnings / errors
- clearer diagnostics

It does not change the numerical application of scaling factors.

Missing metadata is treated as unknown.
