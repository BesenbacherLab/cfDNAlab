# fcoverage Spec

`cfdna fcoverage` is the canonical fragment-coverage command. It owns the current coverage geometry, length-normalization semantics, positional output behavior, and aggregate reducers used by downstream commands such as `coverage-weights` and `fragment-count-weights`.

## Fragment Geometry

- Paired fragments span `[forward.pos, reverse.reference_end)`.
- Unpaired read-as-fragment mode spans `[read.pos, read.reference_end)`.
- Default coverage includes the inter-mate gap as part of the fragment span. `--ignore-gap` removes `[forward.reference_end, reverse.pos)` from counted coverage and cannot be combined with `--reads-are-fragments`.
- Fragment length filters are applied to the fragment definition used by the command. The shared minimum is 10 bp, and default CLI filtering is 30..1000 bp inclusive.
- Always-on read filters remove secondary, supplementary, duplicate, QC-failed, unmapped, cross-tid, and non-inward paired fragments.

## Weighting

- GC correction can come from a correction package (`--gc-file` with `--ref-2bit`) or from a two-byte BAM AUX tag (`--gc-tag`). The two sources are mutually exclusive.
- Invalid or missing GC weights skip fragments by default. `--neutralize-invalid-gc` keeps those fragments with weight 1.0 and records failures in run statistics.
- Scaling factors are multiplicative per-base weights loaded from TSV files with full contiguous chromosome coverage. Known GC-mode mismatches fail during load.
- Coverage-weight TSVs can carry `ignore_gap` metadata. A mismatch with the current `fcoverage --ignore-gap` setting warns, not fails.

## Length Normalization

- `off`: per-base coverage mass is unnormalized.
- `unit-mass`: each fragment contributes total mass 1.0 across its counted bases.
- `restore-mean`: count like `unit-mass`, then multiply the final output by the observed mean normalization length.
- The normalization denominator includes blacklisted positions in the fragment span. Blacklist masking affects reported coverage support, not the fragment's length-normalization mass.
- Output filenames indicate normalization: `length_normalized` for unit mass and `length_normalized.restored_mean` for restored mean.

## Windows And Actions

- Global mode produces one aggregate row unless positional output is requested.
- Fixed-size mode uses contiguous chromosome bins. The final bin is clipped to chromosome length.
- BED mode preserves original BED row indices through reduction. Output order follows the original BED index, not necessarily coordinate order unless the input was already ordered that way.
- Grouped BED mode collapses rows by `group_idx`. Plain grouped aggregate actions preserve the command's grouped span semantics, while `*-on-unique-bases` actions merge overlapping or touching bases within each group before computing support.
- Grouped aggregate outputs write one row per `group_idx`. Duplicate `group_idx` rows are invalid.
- `--per-window average` writes mean coverage over eligible bases.
- `--per-window total` writes coverage sum over eligible bases.
- `--per-window summary-stats` writes the summary schema described below and is handled through reducer paths, not direct positional window extraction.
- BED positional output has two modes: unique positions and indexed per-window positions. Indexed output can repeat a genomic position for multiple BED windows.
- Positional outputs can split at tile boundaries. Reducers must preserve chromosome and tile-index order when merging.

## Blacklists

- Blacklisted positions are excluded from positional and aggregate output support.
- Aggregate rows carry `blacklisted_positions` so downstream consumers can distinguish no support from masked support.
- Fully blacklisted average rows are represented with `NaN` support behavior through the writer/reducer path.
- Summary statistics report both `span_positions` and `eligible_positions`, where `eligible_positions = span_positions - blacklisted_positions`.

## Output Schema

Basic aggregate rows:

```text
chromosome  start  end  average_<signal>|total_<signal>  blacklisted_positions
```

Grouped aggregate rows:

```text
group_idx  span_positions  blacklisted_positions  eligible_positions  average_<signal>|total_<signal>
```

`<signal>` is `coverage` for unnormalized coverage and `fragment_mass` for length-normalized outputs.

Windowed summary rows:

```text
chromosome
start
end
span_positions
blacklisted_positions
eligible_positions
nonzero_positions
covered_fraction
total_<signal>
total_squared_<signal>
average_<signal>
variance_<signal>
sd_<signal>
coefficient_of_variation_<signal>
```

Grouped summary rows:

```text
group_idx
span_positions
blacklisted_positions
eligible_positions
nonzero_positions
covered_fraction
total_<signal>
total_squared_<signal>
average_<signal>
variance_<signal>
sd_<signal>
coefficient_of_variation_<signal>
```

Numeric policies:

- Variance is `E[x^2] - E[x]^2`.
- Tiny negative variance from floating-point error is snapped to zero. Larger negative variance is an error.
- Coefficient of variation is `NaN` when the mean is effectively zero.
- Finite CV values above `1e6` are printed as `>1e6`.
- Restore-mean scales `total_<signal>` by the multiplier and `total_squared_<signal>` by multiplier squared.

## Reducer Invariants

- Reducers consume the explicit tile paths returned by workers. They must not discover work by scanning temp directories.
- BED reducer keys are original BED indices. Duplicate original indices are an error.
- Fixed-size reducer keys are full bin starts before final-bin clipping.
- Cross-index files list only rows that crossed tile boundaries. Missing cross-index keys mean the row had exactly one contribution.
- Aligned tile/window boundaries can bypass fixed-size reduction by concatenating tile finals, but restore-mean always forces reducer aggregation.
- Missing tile output for a chromosome that has expected windows is an error.

## Open Notes

! Warning: Sparse BED runs still build file-based GC reference prefixes for the full tile fetch span before no-window pruning. This is a tracked performance issue, not a correctness mismatch.
