# fcoverage Package Loader Plan

Status: future idea, not an accepted implementation plan.

This note captures possible R and Python package support for `cfdna fcoverage`
outputs. The core decision is whether package support would make common
downstream workflows clearer, or whether it would encourage users to load files
that should stay as indexed genomic tracks.

## Current Outputs

`fcoverage` writes two different output families.

Positional outputs:

- `<prefix>.fcoverage.per_position.bedgraph.zst`
- `<prefix>.fcoverage.per_position_per_window.tsv.zst`

Aggregate outputs:

- `<prefix>.fcoverage.average.tsv.zst`
- `<prefix>.fcoverage.total.tsv.zst`
- `<prefix>.fcoverage.summary_stats.tsv.zst`
- grouped variants such as `average_on_unique_bases`,
  `total_on_unique_bases`, and `summary_stats_on_unique_bases`

These should not share one high-level loader abstraction. Positional tracks are
indexed genomic signal data after conversion to BigWig. Aggregate TSVs are
ordinary row metadata plus numeric summaries.

## Recommendation

Add support for aggregate `fcoverage` TSVs first, if this feature is pursued.

Do not add an eager `read_fcoverage()` path that loads large positional
bedGraph files into memory. Full bedGraph outputs can be many GB and are better
converted to BigWig for indexed range queries.

## Aggregate TSV Loader

Aggregate TSV support fits the existing package pattern:

- validate the TSV schema
- split row metadata from value columns
- expose `window_metadata()` for windowed rows
- expose `group_metadata()` for grouped rows
- support row selectors and blacklist filters
- return ordinary data frames for analysis and plotting

Candidate Python API:

```python
coverage = cfdnalab.read_fcoverage(
    "sample.fcoverage.summary_stats.tsv.zst",
    group_index_path="sample.group_index.tsv",
)

coverage.data_frame(groups=["open"], max_blacklisted_fraction=0.2)
coverage.group_metadata()
coverage.window_metadata()
```

Candidate R API:

```r
coverage <- read_fcoverage(
  "sample.fcoverage.summary_stats.tsv.zst",
  group_index_path = "sample.group_index.tsv"
)

fcoverage_data_frame(coverage, max_blacklisted_fraction = 0.2)
group_metadata(coverage)
window_metadata(coverage)
```

Use a coverage-specific data-frame helper name in R. Avoid `counts_array()` and
`length_counts_matrix()` style names because `fcoverage` aggregate values are
coverage, fragment mass, and summary statistics, not count matrices.

## Modes To Detect

The loader can infer row mode from leading columns.

Windowed non-summary rows:

```text
chromosome  start  end  average_coverage|total_coverage|average_fragment_mass|total_fragment_mass  blacklisted_positions
```

Windowed summary-stat rows:

```text
chromosome  start  end  span_positions  blacklisted_positions  eligible_positions  ...
```

Grouped non-summary rows:

```text
group_idx  span_positions  blacklisted_positions  eligible_positions  average_coverage|total_coverage|average_fragment_mass|total_fragment_mass
```

Grouped summary-stat rows:

```text
group_idx  span_positions  blacklisted_positions  eligible_positions  nonzero_positions  ...
```

The loader should expose the signal label as metadata:

- `coverage`
- `fragment_mass`

It should also expose the value mode when it is inferable:

- `average`
- `total`
- `summary_stats`

Plain grouped output and `*-on-unique-bases` output may have identical columns.
That semantic is not fully self-describing from the main TSV. If the package
needs to expose it, use one of these approaches:

- infer from the filename when it follows cfDNAlab naming
- accept an explicit `aggregation_basis` or `action` argument
- add a future Rust settings JSON for `fcoverage` and read it when present

Do not guess silently when the distinction matters.

## Group Metadata

Grouped `fcoverage` aggregate files contain `group_idx`, not `group_name`.
The Rust command writes a separate `<prefix>.group_index.tsv` group-index file with:

```text
group_idx  group_name
```

Package behavior should be explicit:

- If `group_index_path` is supplied, validate and join group names.
- If it is omitted, try deterministic group-index file discovery only when the main file
  name follows the cfDNAlab pattern.
- If group names are requested and no group-index file is available, error clearly.
- Keep numeric group-index selection usable without group names.

R must convert public group indices to one-based indices. Store the raw file
indices internally as `group_idx0` or equivalent and expose public
`group_idx = group_idx0 + 1`.

Python should keep zero-based public indices, matching the current Python
loader conventions.

## Filtering

Useful row filters:

- `groups` and `group_idxs` for grouped outputs
- `window_idxs` for windowed outputs
- `chromosomes` for windowed outputs
- one half-open genomic interval filter, for example `chrom`, `start`, `end`
- `max_blacklisted_fraction`
- finite-only filtering for value columns when users want to drop rows with
  `NaN` derived metrics

For `fcoverage`, `blacklisted_fraction` is not a stored column. It can be
derived when `span_positions` is available:

```text
blacklisted_fraction = blacklisted_positions / span_positions
```

For non-summary windowed rows, `span_positions` is equivalent to `end - start`.
For grouped rows, `span_positions` is written directly.

If `span_positions` is zero, the derived fraction should be `NaN`, not zero.

## Positional BedGraph Outputs

Do not make compressed bedGraph the main package-level track interface.

Reasons:

- `.bedgraph.zst` has no genomic index.
- Region queries require scanning unless an external index exists.
- Eager loading can be too large for ordinary R or Python sessions.
- Tile boundaries can split bedGraph rows even when the coverage value is
  unchanged, so row shape is not a stable biological unit.

If support is added, make it explicitly streaming and low-level:

```python
track = cfdnalab.open_fcoverage_bedgraph("sample.fcoverage.per_position.bedgraph.zst")
for chunk in track.iter_runs(chrom="chr1"):
    ...
```

Do not provide a convenient full-data-frame method unless it requires an
explicit opt-in such as `allow_full_read=True`.

## BigWig Interop

The preferred positional workflow is:

1. Convert `fcoverage` bedGraph to BigWig.
2. Use existing BigWig readers for indexed range queries.

Python option:

- `pyBigWig`
- expose optional helpers only when installed
- useful methods map naturally to `values()`, `stats()`, and `intervals()`

R options:

- Bioconductor `rtracklayer` for `BigWigFile`, `import.bw()`, and `summary()`
- Bioconductor `megadepth` for region coverage matrices from BigWig or BAM

Do not make these hard runtime dependencies in the first pass. They are better
as optional integrations or documented workflows.

R-specific caveat:

- `rtracklayer` BigWig support is not available on Windows according to its
  manual, so hard-depending on it would make the R package less portable.

## Indexed Positions Output

`per_position_per_window.tsv.zst` keeps duplicate overlapping positions and
appends the original BED window index.

This output is useful for specialized workflows but not self-contained enough
for a friendly loader:

- it has no header
- it needs the original BED or derived window metadata to interpret window ids
- overlapping windows intentionally duplicate positions

Treat it as a low-level run table unless a future use case justifies a richer
object. If a richer object is added, require the source BED or a window metadata
file rather than guessing.

## Implementation Order

TODO: before treating R/Python aggregate `fcoverage` loaders as supported, add
cfDNAlab-generated downstream fixtures for average, total, summary-stats,
grouped output with `group_index.tsv`, and unique-base grouped output. The
fixtures should cover both `coverage` and `fragment_mass` headers where the
loader behavior differs.

1. Add aggregate TSV schema validation in Python and R.
2. Add grouped group-index file loading and group-name selectors.
3. Add `window_metadata()`, `group_metadata()`, and data-frame helpers.
4. Add derived `blacklisted_fraction` filtering.
5. Document BigWig conversion and existing R/Python BigWig reader options.
6. Consider optional BigWig interop helpers.
7. Consider low-level streaming bedGraph support only after aggregate TSV
   support is stable.

## Open Questions

- Should `fcoverage` write a settings JSON so package loaders can distinguish
  grouped plain and grouped unique-base outputs without filename inference?
- Should grouped aggregate TSVs eventually include `group_name` directly, like
  grouped length outputs, or is the separate group-index file acceptable?
- Should package loaders accept `chromosome` as a public column name for
  `fcoverage`, or normalize it to `chrom` for consistency with other package
  loaders?
- Should R use `fcoverage_data_frame()` or a more general `coverage_data_frame()`
  generic?
- Should Python expose one class, `FCoverageAggregates`, or mode-specific
  subclasses for windowed and grouped outputs?
