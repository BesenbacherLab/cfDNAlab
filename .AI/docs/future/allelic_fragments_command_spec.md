# `cfdna allelic-fragments` future spec

## Scope

This is a future design spec for a command that writes one row per fragment
overlapping one or more prespecified A/B target sites.

The command should be generalized enough for tumor-informed, germline, trio,
and haplotype-oriented analyses, while staying easy to run and producing a
compact output table.

Working command name:

```text
cfdna allelic-fragments
```

This keeps `fragments` as the head noun, so users should expect fragment rows.
It also avoids using `cfdna fragments`, which should remain available for a
more general future fragment extraction command.

## Use Cases

The command should be designed from the use cases outward. The default output
should serve the common fragment-level workflows, while optional features should
cover heavier analyses without making every run write wide rows.

### Tumor-informed fragmentomics

Users provide tumor-informed sites and compare fragment features by whether the
fragment carries the A allele, B allele, another observed base, or no callable
base.

The command should support downstream grouping by allele class without requiring
users to parse target-level lists first.

Design consequences:

- `allele_class` should be a default column, because it is the primary grouping
  variable.
- A/B allele labels must stay generic. The tumor allele might be A or B
  depending on how the target file was prepared.
- `O` and `U` must be distinct. An observed non-A/non-B base is different from a
  site with no callable base.
- Control fragment sampling matters, because tumor-informed sites are sparse
  and target-overlapping fragments need a background comparison set.
- GC correction and genomic scaling weights are useful optional columns for
  downstream models that compare feature distributions across genomic regions.

Useful features:

- default `allele_class`
- default `target_calls`
- optional `target_ids` when target records have stable names
- optional `gc_weight`, `gc_fraction`, and scaling weights
- optional control sampling with deterministic seeds
- optional target-call diagnostics for auditing low-support calls

### Trio and haplotype studies

A/B labels are not ref/alt. They may represent maternal and paternal alleles,
haplotypes, or other study-specific allele labels.

The command must not assume one allele is the reference allele. It should avoid
columns such as `ref`, `alt`, or `is_alt`.

Design consequences:

- The target file should use `A_allele` and `B_allele`, not `ref` and `alt`.
- Multi-target fragments are meaningful, because a fragment can provide direct
  evidence about local haplotype consistency.
- Disagreement across targets should be visible in `target_calls`, not hidden
  behind a single hard class.
- Mixed classes can be useful QC signals for phasing, but only when the
  overlapping targets are expected to belong to the same phased context. If the
  fragment crosses targets from different phase blocks, `A+B` is not
  automatically evidence that phasing failed.
- Fragment-relative target positions are useful, because the distance between
  targets inside one fragment can matter for phasing and molecule-level
  interpretation.
- Requiring a separate per-fragment-site table would make the command less
  convenient for the main fragment-level workflow, so list columns are the
  better compromise here.

Useful features:

- default `target_calls`, including multi-target calls such as `A;B`
- default fragment-relative `target_positions`
- optional `target_ids`
- possible future target-file `phase_block` column for phasing QC
- optional read-depth or call-support diagnostics for evaluating ambiguous
  haplotype evidence
- possible future haplotype-consistency summary for fragments with multiple
  callable targets

### Prenatal and fetal-fraction studies

In maternal plasma, allele-informative fragments can help separate fetal and
maternal contributions when the target panel encodes informative parental or
fetal alleles.

This use case overlaps with trio and haplotype studies, but it is worth naming
separately because the analysis question is often mixture estimation rather
than haplotype reconstruction.

Design consequences:

- A/B labels must be able to stand for fetal, maternal, paternal, inherited, or
  non-inherited alleles without changing command semantics.
- Control sampling and genomic correction are important, because fetal-fraction
  estimates can be biased by regional coverage and fragment feature differences.
- Target-call diagnostics may matter more than in ordinary fragment feature
  extraction, because low fetal fraction makes false allele support costly.
- The command should not hard-code fetal-specific terminology. The same output
  shape should serve other mixture settings.

Useful features:

- default `allele_class`
- optional `target_ids` for informative-marker panels
- optional control sampling
- optional GC and scaling weights
- optional target read depths and base-quality summaries

### Donor-recipient and transplant cfDNA

In transplant or donor-derived cfDNA analyses, A/B labels can represent donor
and recipient alleles at informative variants.

This is another mixture use case, but the biological interpretation differs
from prenatal cfDNA and tumor-informed analyses.

Design consequences:

- The command should remain label-agnostic. Donor and recipient labels belong
  in the target file or downstream metadata, not in the core output schema.
- The default fragment table is still useful, because downstream users can
  group fragments by `allele_class` and compare fragment features or estimate
  mixture proportions.
- Optional target IDs are more important when panels are reused across samples
  or donor-recipient pairs.

Useful features:

- default `allele_class`
- default `target_calls`
- optional `target_ids`
- optional target-call diagnostics
- optional target-panel QC summary

### Germline heterozygous-site profiling

Users may pass known heterozygous sites and compare fragmentomic features by
allele or haplotype.

Most fragments will overlap zero or one target. The default output should be
optimized for that common case.

Design consequences:

- Default rows should be compact. Count columns such as `n_a`, `n_b`,
  `n_other`, and `n_uncalled` are not worth the file size when `target_calls`
  already contains the same information.
- `fragment_length` should not be written by default, because `end - start`
  provides the aligned fragment length.
- Target positions should be fragment-relative if written, because absolute
  positions are recoverable from the target file and are less useful for
  fragment feature modeling.
- Optional target IDs are reasonable when the target panel has meaningful marker
  names, but synthetic IDs should not be invented by default.

Useful features:

- lean default columns
- optional target IDs
- optional end motifs, because allele-specific end signatures may be relevant
- optional `gc_fraction` for bias checks
- optional target-base-quality filtering

### Fragment-level model input

Users may want one compact table for model fitting or exploratory analysis.

The output should include enough information to group, filter, and model
fragments without re-reading the BAM, but it should not duplicate values that
are directly inferable from other output columns and the target file.

Design consequences:

- The table should stay one row per fragment.
- Default columns should be stable and low-cardinality where possible.
- Optional features should be explicit, because every added column multiplies
  file size by the number of target-overlapping and control fragments.
- Features that require reference sequence, such as outside motifs and
  `gc_fraction`, should be opt-in.
- Columns that are direct arithmetic transforms of default columns should not be
  defaults.

Useful features:

- default `allele_class` for grouping
- default `target_calls` for auditing the class
- default fragment-relative `target_positions` for feature engineering
- optional end motifs
- optional GC and scaling weights
- optional mismatch count
- optional controls with a reproducible sampling scheme

### Control fragment sampling

Users may want a matched or global sample of fragments that do not overlap any
target site, so target-overlapping fragment features can be compared against a
background set without writing every fragment from the BAM.

Design consequences:

- Control output should be opt-in. Writing every non-target fragment would turn
  this command into a broad fragment extractor and make output size
  unpredictable.
- Control rows should use `allele_class = control` and `.` for target-specific
  list columns.
- Sampling must be deterministic for reproducibility.
- Global, bin-matched, and bin-uniform controls should be the main modes. They
  answer different questions without turning control sampling into a large menu
  of distance-based variants.

Useful features:

- `--controls-per-case`
- `--n-controls`
- `--control-mode global`
- `--control-mode bin-matched=<bp>`
- `--control-mode bin-uniform=<bp>`
- `--control-seed`
- optional GC and scaling columns for background-adjusted comparisons

### Target-panel QC

Users may also want to use the command to audit whether a target panel behaves
well in a BAM.

This is not the main output mode, but it influences optional diagnostics.

Design consequences:

- The default fragment table should not become a panel-QC report.
- The command can still expose optional diagnostics that make bad targets
  visible without requiring another BAM pass.
- A future summary sidecar could report target-level aggregate QC without
  changing the fragment table shape.

Useful features:

- optional target read depths
- optional target base-quality summaries
- optional call reasons such as deletion, refskip, no read coverage, or base
  quality filtered
- optional run summary with counts by `allele_class`
- possible future target-level QC summary file

## Inputs

Required inputs:

- coordinate-sorted, indexed BAM
- target TSV

The target TSV should require a header:

```text
chrom	pos	A_allele	B_allele
```

Optional target TSV column:

```text
id
```

`pos` is 1-based, matching VCF-style position semantics. Internally, the
command converts each target to a 0-based single-base interval `[pos - 1, pos)`.

The initial target-calling semantics should be durable for SNPs:

- `A_allele` and `B_allele` must each be one of `A`, `C`, `G`, or `T`
- lower-case input is normalized to upper-case
- equal A and B alleles are an error
- duplicate `(chrom, pos, A_allele, B_allele)` targets are an error

Each input row is a distinct target record. Duplicate `(chrom, pos)` rows are
allowed because split multiallelic VCF records can produce several A/B
comparisons at the same physical locus.

When duplicate `(chrom, pos)` rows are present:

- the target TSV must include unique `id` values
- `target_ids` should be written automatically, because `target_positions`
  alone cannot identify which target row was called
- same-position rows must have a collapsible allele-label structure for
  fragment-level `allele_class`

A same-position target group is collapsible when either all rows share the same
`A_allele` or all rows share the same `B_allele`. Collapsible groups define a
site-level allele map:

- if all rows share `A_allele`, the site-level A set is the shared A base and
  the site-level B set is every B base listed at that position
- if all rows share `B_allele`, the site-level A set is every A base listed at
  that position and the site-level B set is the shared B base

Same-position groups that cannot be collapsed should fail with a clear error
instead of producing misleading fragment-level classes.

This distinction matters for split multiallelic-style input. For example,
targets `(A=C, B=T)` and `(A=C, B=G)` at the same position share the same A base.
An observed `T` is row-level `B` for the first target and row-level `O` for the
second target, but it is site-level `B` for fragment classification because `T`
is one of the B options at that physical site.

The command should not require `--ref-2bit` for allele calling, because A/B
alleles are supplied directly by the target TSV.

## Fragment and Target Semantics

Paired-end fragments use the existing directional project definition:

```text
forward.pos to reverse.reference_end
```

The emitted `start` and `end` columns use that fragment span.

A target overlaps a fragment when the target's 0-based position lies inside the
fragment span `[start, end)`.

A target can overlap the fragment span without being callable. This happens when
the target lies in the inter-mate gap, in a deletion or refskip, or in a part of
the fragment not covered by a concrete read base.

A fragment is a target-overlapping fragment when its span overlaps at least one
target site. This is true whether or not any overlapping target is callable.
Control fragments must pass the ordinary fragment filters and overlap zero
target sites.

Callability is not missing at random for fragmentomics. Longer fragments are
more likely to place interior targets in the inter-mate gap, so allele-stratified
length analyses should treat uncalled fragments as a meaningful QC category
rather than silently dropping them.

## Allele Calling

For each overlapping target, inspect observed read bases at the target reference
position.

Per-read behavior:

- a concrete observed base equal to `A_allele` contributes `A`
- a concrete observed base equal to `B_allele` contributes `B`
- a concrete observed base equal to neither A nor B contributes `O`
- deletion, refskip, absent read coverage, soft clipping over the reference
  position, and failed base-quality filtering do not contribute a concrete base
- an insertion adjacent to the target does not erase the aligned reference base
  at the target position

Mate-overlap behavior should use consensus by default:

- if only one mate has a passing concrete base at the target, use that base
- if both mates have passing concrete bases and they agree, use that base
- if both mates have passing concrete bases and disagree, the target call is
  `U`
- if neither mate has a passing concrete base, the target call is `U`

Mate disagreement should not create `A+B` by default, because this makes a
likely base-level error look like biological allele mixture. Optional call
diagnostics can preserve `mate_disagreement` as a reason for `U`.

Per-target call behavior after mate consensus:

- if no concrete base remains, the target call is `U`
- if the concrete base equals `A_allele`, the target call is `A`
- if the concrete base equals `B_allele`, the target call is `B`
- if the concrete base equals neither A nor B, the target call is `O`

Fragment-level `allele_class` is derived from site-level evidence after
same-position target rows are collapsed for classification:

- `A`: at least one A observation, no B observations, no O observations
- `B`: at least one B observation, no A observations, no O observations
- `A+B`: at least one A observation and at least one B observation, no O
  observations
- `O`: at least one O observation, no A observations, no B observations
- `mixed_other`: at least one O observation plus at least one A or B observation
- `uncalled`: all overlapping targets are `U`
- `control`: sampled non-target-overlapping control fragments

`U` does not make a fragment mixed. For example, target calls `A;U` classify as
`A`.

With the default mate-consensus policy, `A+B` means disagreement across multiple
target loci, or across a collapsible same-position target group. It should not
mean that two mates disagreed at a single target.

For unique target positions, row-level `target_calls` and site-level evidence
are the same. For duplicate target positions, `target_calls` remains aligned to
input target rows, while `allele_class` uses the collapsed site-level allele map
for that physical position.

## Default Output

The default output should be a compressed fragment-like TSV and a companion
header file, following the `bam-to-frag` output style.

Default columns:

```text
chromosome
start
end
min_mapq
read1_strand
allele_class
target_calls
target_positions
```

When the target TSV contains duplicate `(chrom, pos)` rows, `target_ids` should
also be written by default. Without it, `target_positions` cannot identify the
target row that produced each call.

`target_calls` is a `;`-separated list in genomic target order.

`target_positions` is a `;`-separated list in the same order as
`target_calls`. These are fragment-relative 0-based positions measured from the
left fragment coordinate, `start`:

```text
target_position = target_0based_pos - fragment_start
```

For control fragments, target-specific list columns should be `.`.

Default output intentionally omits count columns such as `n_a`, `n_b`,
`n_other`, and `n_uncalled`. They are redundant with `target_calls` and add file
size for the common one-target-overlap case.

Default output also omits `fragment_length`, because it is inferable as
`end - start` for the default aligned fragment span.

### R and Python Helpers

The CLI output should stay compact by default. R and Python helpers should not
expand `target_calls`, `target_positions`, or `target_ids` automatically.

The packages should instead offer explicit helper functions or read options for
users who want to parse list-like columns or expand the fragment table into a
long per-fragment-target data frame. Those helpers must document coordinate
semantics clearly: target TSV positions are 1-based genomic positions, while
`target_positions` in the fragment output are 0-based positions relative to the
fragment `start`.

## Optional Output Columns

### Target IDs

If the target TSV has an `id` column, users may request:

```text
target_ids
```

`target_ids` is a `;`-separated list aligned with `target_calls` and
`target_positions`.

The command should not synthesize IDs by default for target files with unique
positions. If a user requests target IDs and the target file has no `id` column,
the command should fail with a clear error.

When the target file contains duplicate `(chrom, pos)` rows, `id` is required
and `target_ids` should be written automatically.

### End Motifs

Users may request fragment end motifs:

```text
left_end_motif
right_end_motif
```

These should reuse the semantics of `cfdna ends` where possible:

- configurable `k_inside`
- configurable `k_outside`
- right-end motifs oriented consistently with `ends`
- clip handling policy aligned with existing `ends` behavior
- reference-backed outside bases require `--ref-2bit`

This command should only write per-fragment motif annotations. It should not
inherit the matrix-counting or window-assignment output model from `ends`.

### GC Features

Optional GC-related columns:

```text
gc_weight
gc_fraction
```

`gc_weight` should follow existing GC correction package behavior when a
`--gc-file` is supplied.

`gc_fraction` requires reference access and should be explicit, because it adds
work and output size.

### Genomic Scaling Weights

Optional scaling columns may mirror `bam-to-frag`:

```text
coverage_scaling_weight
count_scaling_weight
```

These are useful when downstream models need to carry region-level correction
weights with each fragment row.

### Mismatch Features

Optional mismatch columns:

```text
mismatch_count
```

The command should not write `mismatch_count_excluding_targets` by default.
Given A/B input, the command does not know which allele is the reference allele.
Therefore the target contribution to reference mismatches is not inferable from
A/B calls alone.

If future behavior uses `--ref-2bit` or a target file with reference alleles,
the command can add an explicit target-mismatch component rather than writing
both a value and its subtraction.

### Call Support Diagnostics

Optional diagnostics may be useful, but should not be default columns:

```text
target_read_depths
target_base_qualities
target_call_reasons
```

These would increase file size substantially. They should be added only when
there is a concrete downstream need.

## Speculative Feature Ideas

These ideas are not main recommendations for the core command surface. They are
design space worth keeping in mind if users find the basic fragment table useful
and start asking for richer allele-aware fragmentomics.

### Haplotype Consistency

For fragments overlapping multiple callable targets, the command could add a
compact consistency label:

```text
haplotype_consistency
```

Possible values:

- `single_target`
- `consistent_A`
- `consistent_B`
- `switch`
- `mixed_other`
- `uncalled`

This would be useful for trio and haplotype studies, but it should wait until
real examples show what labels are actually useful.

This feature probably needs an optional target-file column such as:

```text
phase_block
```

Without a phase-block label, a mixed fragment call can mean several different
things:

- a phasing error
- an A/B assignment error in the target file
- a read or mapping error
- a fragment crossing targets that were never meant to share one haplotype
  interpretation

So a phasing-QC feature should be explicit about which targets are expected to
be mutually consistent.

### Target Distance Features

For multi-target fragments, the command could optionally write compact distance
features instead of full per-site detail:

```text
target_span
nearest_target_to_left_end
nearest_target_to_right_end
```

These are model-friendly and smaller than carrying many target annotations.
They are still derivable from `target_positions`, so they should not be default
columns.

### Allele-Aware End Features

The command could stratify end-related annotations by target position:

```text
target_to_left_end
target_to_right_end
target_nearest_end
target_nearest_end_motif
```

This is biologically plausible for questions about mutation-bearing fragments
having different cleavage patterns near the observed allele. It is also easy to
overfit, so it should be optional and probably experimental at first.

### Read-Overlap Diagnostics

The command could expose how each target was observed:

```text
target_observation_mode
```

Possible values:

- `forward`
- `reverse`
- `both_agree`
- `both_disagree`
- `gap`
- `deletion`
- `refskip`
- `base_quality_filtered`

This would help debug allele calls and overlapping-mate conflicts. It is too
wide for default output if repeated for every target.

### Duplex-Like Support

If input data or tags support strand-family interpretation, a future mode could
separate weak single-read evidence from stronger bidirectional evidence.

Possible columns:

```text
target_strand_support
allele_support_tier
```

This should not be in the base design. It depends too much on sequencing
protocol and tagging conventions.

### Target-Panel Summary Sidecar

The command could optionally write a small target-level QC summary:

```text
target_id
chromosome
pos
fragments_overlapping
fragments_callable
calls_A
calls_B
calls_O
calls_U
median_target_base_quality
```

This would not be a per-fragment-site table. It would be an aggregate QC file
for judging whether a target panel is usable in a sample.

### Allele-Class Balanced Sampling

For model-building workflows, the command could optionally cap output per
allele class:

```text
--max-fragments-per-class <N>
```

This could help fast exploratory runs. It should be treated carefully, because
it changes the observed allele-class proportions.

### Distance-Based Control Modes

Several local control definitions are plausible:

```text
halo-only
exclude-halo
distance-weighted
```

`halo-only` would sample non-target fragments near target sites. `exclude-halo`
would do the opposite and deliberately avoid target neighborhoods.
`distance-weighted` would sample with a probability based on distance to the
nearest target.

These modes answer different questions. None should be added to the main API
unless a clear workflow requires it. The command should first optimize
`global`, `bin-matched`, and `bin-uniform`, because those cover the broad
background choices users are most likely to need.

### Nearby Context Annotations

Users may eventually want context columns such as:

```text
target_ref_context
target_gc_window
target_mappability
target_blacklist_distance
```

These are useful for QC and modeling, but they require extra reference or
annotation inputs and would move the command toward being a feature factory.
They belong behind explicit options or in downstream tooling unless a strong
common use case appears.

## Control Sampling

Control fragments are fragments that pass the normal read and fragment filters
but do not overlap any target site in their fragment span.

Control output should be opt-in.

Possible options:

```text
--controls-per-case <N>
--n-controls <N>
--control-mode <global|bin-matched=<bp>|bin-uniform=<bp>>
--control-seed <integer>
```

If neither `--controls-per-case` nor `--n-controls` is supplied, the command
should not run a counting pass. It should go directly to the full fragment pass
and write only target-overlapping fragments.

When controls are requested, control sampling should be deterministic for a
given seed, BAM, target file, chromosome selection, and filter configuration.

`control-mode` defines both the candidate universe and the sampling allocation:

- `global` samples uniformly across all non-target-overlapping fragments on the
  selected chromosomes.
- `bin-matched=<bp>` samples within fixed genomic bins of the requested size,
  proportional to target-overlapping fragments in each bin.
- `bin-uniform=<bp>` samples within fixed genomic bins of the requested size,
  with control allocation proportional to genomic bin weight rather than local
  case count or local sequencing depth.

These are different statistical choices. `global` is the right control when
users want a broad background over the observed fragment distribution. It is
uniform over fragments, so high-depth, amplified, or otherwise overrepresented
regions contribute more controls. `bin-matched` is the right control when users
want local genomic context to track the case fragments more closely.
`bin-uniform` is the right control when users want a background that is closer
to position-uniform across the selected genome while still outputting observed
fragments.

Other local definitions, such as halo-only controls, halo-excluded controls, or
distance-weighted controls, are plausible. They should stay out of the main API
until there is a clear need, because each adds another reasonable but
different control definition.

### Sampling Method

The recommended exact design is two-pass rank sampling. It avoids writing
temporary control rows and keeps the statistical definition independent of
processing tiles.

The first pass should use a minimal fragment representation. It should pair
reads, apply the ordinary read and fragment filters, apply blacklist filtering,
compute the project fragment span, and test whether the span overlaps any target
site. It should not inspect target bases, build allele calls, compute motifs,
compute GC annotations, or build final output rows.

For binned control modes, the first pass must also record `(bin, partition)`
control-candidate counts. This table is what makes it possible to recover the
correct bin-local rank when a processing partition starts in the middle of a
bin.

The second pass should use the full with-records fragment path. It writes all
target-overlapping fragments and writes only sampled control fragments.

The first and second pass must enumerate eligible control candidates in the
same canonical order. The canonical order should be independent of thread count,
tile size, and scheduling. It should be the same order used for sorted fragment
output: selected chromosome order, then fragment `start`, `end`, `read1_strand`,
and a stable fragment or template identity tie-breaker.

Execution partitions are implementation chunks, not part of the statistical
definition. They should be contiguous ranges in the canonical order so prefix
sums can map sampled ranks back to workers. Fetch halos may be used internally,
but a fragment is counted only once, in the partition whose core owns the
fragment start.

### Global Rank Sampling

For `--control-mode global`, pass 1 counts eligible control candidates per
execution partition:

```text
control_count_for_partition[p]
total_control_candidates = sum(control_count_for_partition)
total_cases = target-overlapping fragments after ordinary filters
```

For `--controls-per-case`, compute:

```text
requested_controls = ceil(total_cases * controls_per_case)
```

For `--n-controls`, use:

```text
requested_controls = n_controls
```

The command then samples `requested_controls` unique integer ranks from:

```text
0 .. total_control_candidates
```

These ranks refer to the canonical ordered list of eligible control candidates,
not to processing tiles. The implementation may choose the most appropriate
exact algorithm for sampling unique ranks, but the sampled rank set must be
deterministic for the control seed.

If `requested_controls` is larger than `total_control_candidates`, the command
should deliver all eligible controls and report the shortfall clearly in the run
summary.

Prefix sums map global ranks to partition-local ranks:

```text
partition_prefix[p] = sum(control_count_for_partition[q] for q < p)
```

A global rank belongs to partition `p` when:

```text
partition_prefix[p] <= rank < partition_prefix[p] + control_count_for_partition[p]
```

The local rank is:

```text
local_rank = rank - partition_prefix[p]
```

During pass 2, each worker enumerates eligible controls in the same local order
used by pass 1. It writes a control only when the current local rank is in that
partition's selected local rank set.

### Binned Rank Sampling

For `--control-mode bin-matched=<bp>` and
`--control-mode bin-uniform=<bp>`, bins are fixed genomic bins, not processing
tiles:

```text
bin = floor(fragment_start / bin_size)
```

The bin coordinate system is 0-based from the chromosome start, and a fragment
belongs to exactly one bin by its `start` coordinate. Bins do not need to be
divisible by tile size. A bin may cross execution partition boundaries.

Pass 1 counts:

```text
case_count_for_bin[bin]
control_count_for_bin_and_partition[bin, partition]
control_count_for_bin[bin] = sum(control_count_for_bin_and_partition[bin, p])
```

The `(bin, partition)` count matrix is required because execution partitions do
not need to start on bin boundaries. If partition `p` starts in the middle of a
bin, the controls before that partition are represented by the per-bin prefix
sum:

```text
bin_partition_prefix[bin, p] =
    sum(control_count_for_bin_and_partition[bin, q] for q < p)
```

During pass 2, each worker keeps a local counter for each bin. The bin-local
rank of the current control is:

```text
rank_in_bin =
    bin_partition_prefix[bin, partition] + local_control_counter_for_bin
```

The worker writes the control when `rank_in_bin` is selected for that bin.

### Bin-Matched Allocation

For `--controls-per-case`, compute:

```text
total_cases = sum(case_count_for_bin)
requested_controls_total = ceil(total_cases * controls_per_case)
raw_controls_for_bin = case_count_for_bin * controls_per_case
```

The per-bin requested counts should be allocated with floors plus largest
remainders so the total equals `requested_controls_total` without inflating
sparse bins by independently applying `ceil` to every bin.

For `--n-controls`, allocate `n_controls` across bins in proportion to
`case_count_for_bin`, again using floors plus largest remainders. `bin-matched`
with no target-overlapping fragments should fail with a clear error, because
there is no case distribution to match.

### Bin-Uniform Allocation

For `bin-uniform`, allocate controls across bins by genomic bin weight rather
than by observed case fragments.

The default bin weight should be the number of selected genomic bases in the
bin after clipping to chromosome bounds. If blacklist filtering has a
well-defined base-level excluded fraction for the command configuration, the
weight should use non-blacklisted bases. Target sites do not define the bin
weight, because target-overlap exclusion is already applied at the fragment
candidate level.

For `--controls-per-case`, compute:

```text
total_cases = target-overlapping fragments after ordinary filters
requested_controls_total = ceil(total_cases * controls_per_case)
```

For `--n-controls`, use:

```text
requested_controls_total = n_controls
```

Allocate `requested_controls_total` across bins in proportion to
`bin_weight[bin]`, using floors plus largest remainders. Bins with zero weight
should receive zero requested controls.

`bin-uniform` still samples observed fragments within each bin. It reduces
overrepresentation of high-depth regions across bins, but it does not make
fragment positions perfectly uniform inside a large bin. Smaller bins make the
mode closer to position-uniform and increase the chance that some bins have too
few eligible control fragments.

### Per-Bin Rank Selection

For each binned mode, once the requested count for a bin is known, sample that
many unique ranks from:

```text
0 .. control_count_for_bin[bin]
```

If a bin has fewer eligible controls than requested, the command must not fail
silently. The run summary should report requested and delivered controls
genome-wide and per affected bin. The default should not silently refill missing
controls from other bins, because that changes the chosen binned sampling
definition.

## Filtering and Defaults

Read and fragment filters should follow the existing BAM fragment commands:

- secondary, supplementary, duplicate, and failed-QC reads are excluded
- paired-end mode requires mapped mates on the same tid and inward orientation
- `--reads-are-fragments` supports unpaired fragment-like input
- fragment length filters use the existing defaults and validation
- `--min-mapq` should default to `30`, matching the more interpretive fragment
  commands rather than the raw `bam-to-frag` exporter
- blacklist filtering should reuse existing blacklist semantics

Target-base filters should be explicit:

```text
--min-target-base-quality <integer>
```

If all observed bases at a target fail this filter, the target call is `U`.

The command should include `uncalled` fragments by default, because they passed
the target-overlap requirement and may be informative for QC. Users can request
dropping them with an explicit option.

## Implementation Fit

The implementation should be closer to `bam-to-frag` than to `ends`.

Reusable pieces:

- `bam-to-frag` orchestration for per-chromosome or partitioned parallel
  processing, sorted chunks, final gzip concat, output prefix handling, and
  header writing
- the minimal fragment path for the first counting pass when controls are
  requested
- `fragments_with_records_from_bam`, because allele calling needs BAM records,
  CIGARs, sequences, and base qualities in the full output pass
- `Interval` and `IndexedInterval` helpers for target overlap
- existing GC correction and scaling helpers for optional feature columns
- `ends` motif semantics for optional per-fragment end motif annotations

Avoid:

- using `FragFileFragment` for allele calling, because it discards the records
  needed to inspect target bases
- copying the `ends` matrix-counting output model
- adding a second per-fragment-site table
- making processing tiles part of the statistical definition for control
  sampling

## Open Decisions

- Exact command name, with `allelic-fragments` as the current preferred name
- Exact name for the position-oriented binned control mode, with
  `bin-uniform=<bp>` as the current working name
- Whether `target_positions` should always be written by default or be
  controlled by a feature flag if file size becomes a problem
- Whether indel targets should be supported from the first implementation or
  added later as an extension to the SNP target semantics
