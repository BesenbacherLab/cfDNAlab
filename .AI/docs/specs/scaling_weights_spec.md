# scaling_weights Spec

This spec covers `cfdna coverage-weights`, `cfdna fragment-count-weights`, and the shared scaling-factor TSV contract consumed by analysis commands.

## Purpose

Scaling weights are multiplicative genomic normalization factors. The commands estimate broad genomic coverage or fragment-mass trends in stride bins, smooth them with a triangular kernel, normalize the smoothed values to global mean 1.0, then write the inverse as `scaling_factor`.

Use:

- `coverage-weights` when downstream features are coverage-like and long fragments should contribute more mass because they cover more bases.
- `fragment-count-weights` when downstream features should make long and short fragments contribute more equally.

## Internal fcoverage Source

- Both commands call internal `fcoverage` rather than duplicating fragment iteration.
- `coverage-weights` runs by-size stride bins with `--per-window average`.
- `fragment-count-weights` runs by-size stride bins with `--normalize-by-length=unit-mass --per-window total`.
- Internal `fcoverage` output is written under a guarded temporary directory inside the selected output directory.
- Intermediate fcoverage values are written with 12 decimals and then read back for smoothing.
- The top-level command owns the CLI banner and statistics. Nested `fcoverage` logs should appear as phase logs only.

## Fragment Geometry And Filters

- Fragment span semantics follow `fcoverage`: paired `[forward.pos, reverse.reference_end)`, unpaired `[read.pos, read.reference_end)`.
- Shared filters include MAPQ, fragment length range, duplicate/QC/secondary/supplementary filtering, inward pairing, optional proper-pair filtering, optional chromosome selection, and optional blacklists.
- `coverage-weights --ignore-gap` forwards to internal `fcoverage --ignore-gap` and cannot be used with `--reads-are-fragments`.
- `fragment-count-weights` does not expose `--ignore-gap`.
- Optional GC correction can be used when the intended downstream run also uses GC correction, to avoid baking GC bias into genomic smoothing factors.

## Smoothing

- `bin_size` is the large smoothing window size.
- `stride` is the output-bin size and internal fcoverage fixed-window size.
- `stride` must be less than or equal to `bin_size`.
- `bin_size` must be divisible by `stride`.
- The triangular kernel radius is `(bin_size / stride) - 1` stride bins.
- Kernel weights are integer triangular weights centered on each stride bin.
- Chromosome edges truncate the kernel and renormalize by the actually used finite support.
- Last stride bins can be shorter than `stride`; smoothing weights those rows by `bin_length / stride`.
- Non-finite raw stride values do not contribute to smoothing neighbors. They may still receive finite smoothed values from neighboring support, but their final `scaling_factor` is 0.

## Global Normalization

- The global mean is computed over finite, non-zero smoothed rows whose raw stride value is finite.
- The global mean is length-weighted by stride-bin length.
- Each row's smoothed value is divided by the global mean.
- Public `scaling_factor` is the inverse normalized value.
- Rows without usable finite support get `scaling_factor = 0`.
- If no usable finite non-zero smoothed mass exists, the command fails with a user-facing error.

## Output TSV

`coverage-weights` writes:

```text
<prefix>.coverage.scaling_factors.tsv
```

`fragment-count-weights` writes:

```text
<prefix>.fragment_counts.scaling_factors.tsv
```

Leading metadata lines:

```text
# gc_mode=uncorrected|corrected_file|corrected_tag
# ignore_gap=true|false
```

`ignore_gap` is written only when the source command has that concept.

Header:

```text
chromosome  start  end  stride_average_coverage|stride_fragment_mass  smoothed_coverage|smoothed_fragment_mass  scaling_factor
```

Rows are written in resolved chromosome order and contiguous stride order from 0 to chromosome length.

## Consumer Contract

- Scaling TSVs may start with `# key=value` metadata comments before the header.
- Unknown metadata keys are ignored.
- Duplicate known metadata keys are errors.
- The header is required and column names are matched case-insensitively.
- Required columns for consumers are `chromosome`, `start`, `end`, and `scaling_factor`.
- Coordinates are 0-based half-open intervals.
- Scaling factors must be finite and non-negative.
- For every selected chromosome, rows must start at 0, be perfectly contiguous, and end exactly at the BAM chromosome length.
- Known raw-vs-GC-corrected mismatches fail. File-based and tag-based corrected modes are considered mutually compatible.
- Missing GC metadata is accepted with a warning.
- `ignore_gap` metadata mismatches warn, because old files and count-based files may not know the setting.
