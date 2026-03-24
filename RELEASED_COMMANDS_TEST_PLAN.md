# Released Commands Test Plan

This plan is for the intended public command set from `CHANGELOG`, not the experimental feature-gated commands.

Scope:

- `cfdna fcoverage`
- `cfdna lengths`
- `cfdna midpoints`
- `cfdna gc-bias`
- `cfdna coverage-weights`
- `cfdna bam-to-bam`
- `cfdna bam-to-frag`
- `cfdna frag-to-bam`
- `cfdna ref-gc-bias`

The goal is not "more tests" in the abstract. The goal is confidence that these commands behave as intended under the real scientific invariants they claim to preserve.

What this plan optimizes for:

- First, correctness of biological and mathematical intent
- Second, invariance across tiling, windowing, scaling, masking, and file conversion boundaries
- Third, interoperability between released commands
- Last, CLI polish and small validation cases

What release confidence still lacks:

- One workflow spine per real user pipeline, using actual generated artifacts rather than hand-built stand-ins
- One invariance harness per tiled command, proving that chunking changes runtime only, not science
- One artifact-contract suite for every file format that other released commands consume
- One cross-command semantics suite proving that "fragment", "window overlap", "blacklist overlap", and weighting mean the same thing everywhere they are exposed
- One default-semantics suite proving that silent defaults which materially change outputs are intentional, stable, and explicitly different where the commands are meant to differ

## Current read on the codebase

The repo is not untested. The problem is that the strongest risks are not uniformly covered.

- `lengths`, `fcoverage`, and `frag-to-bam` already have large active suites
- `bam-to-bam` has a moderate command suite
- `bam-to-frag` is much thinner than a released format-boundary command should be
- `midpoints` and `coverage-weights` have some command coverage but still look structurally thin
- `gc-bias` has a lot of helper-level tests, but little evidence of full `run()` coverage
- `ref-gc-bias` is the most obviously under-tested released command. It has only a very small test surface relative to complexity and no meaningful command-level run coverage

Concrete evidence from the current tree:

- `gc-bias` tests are mostly helper- and reducer-level. There is little active evidence that the full command pipeline is exercised end to end
- There is no active test coverage of `ref_gc_bias::run()` itself; the current `ref-gc-bias` test file only exercises lower-level counting helpers
- `ref-gc-bias` currently has only a tiny test file relative to the command complexity
- `bam-to-frag` appears to have only a few active command tests, mostly smoke/global/BED shape checks
- `midpoints` tests currently emphasize fetch narrowing and output layout more than full signal correctness under combined filters and weights
- `coverage-weights` does have command-level `run()` tests, but they are still dominated by simple one-chromosome derivations and do not yet pin the full structural space of chromosome edges, blacklist effects, and cross-chromosome normalization
- `fcoverage` has a real by-size aligned fast path that bypasses the normal reducer, but the active suite still looks much stronger on general/reduced paths than on proving fast-path parity under all modifiers
- `bam-to-bam` currently has active `FLEN` and `COV` evidence, but little visible command-level evidence for `GC` tags or combined `GC` + `COV` + `FLEN` behavior
- The codebase uses two midpoint semantics in released behavior: randomized tie-breaking for even-length window assignment, but deterministic floor midpoint for blacklist-midpoint filtering. That may be intended, but it is not yet locked down by explicit release-level tests
- The GC artifact loaders (`load_reference_gc_data`, `load_gc_corrector`, `load_length_agnostic_gc_corrector`) have real compatibility and shape checks, but there is almost no active test evidence around those guardrails at the released-command level
- GC artifact versioning is asymmetric: `gc-bias` correction packages carry an explicit schema version, but the `ref-gc-bias` package currently relies on exact field names and shapes without a version field. That makes producer-consumer contract tests even more important
- Current active GC-consumer tests are still dominated by hand-built mini packages. There is much less evidence that real `ref-gc-bias` output feeds `gc-bias`, and real `gc-bias` output feeds all downstream consumers, without semantic drift
- The scaling TSV is intentionally consumed with different math in different commands: per-base scaling in `fcoverage`, full-fragment averaging in `midpoints` and the converters, and assignment-dependent overlap-vs-fragment averaging in `lengths`. That design needs explicit release-level proofs, not just isolated one-command tests
- Several released commands intentionally use different spans for different tasks on the same fragment: for example `lengths` can bin by indel-adjusted length while still using the full reference span for blacklist, scaling, and GC correction; `fcoverage` can omit the mate gap from counting while still applying positional masking/scaling over the counted bases. Those span contracts need explicit release tests
- BED semantics are intentionally not uniform across released commands: `gc-bias` and `ref-gc-bias` flatten to unique positions, `fcoverage` can either flatten or preserve duplicates depending on output mode, `lengths` preserves original BED rows, and `bam-to-bam`/`bam-to-frag` use BEDs only as inclusion filters. That matrix is scientifically important and not yet locked down at release level
- `gc-bias` has a meaningful default windowing choice (`--by-size 100000`, not global), but there is little active evidence that the default behavior is pinned against the explicit equivalent or against the global alternative
- Shared config code exposes materially different defaults across released commands: `fcoverage` and `lengths` default to global/windowless counting, `gc-bias` defaults to `--by-size 100000`, released counting commands usually default to `min_mapq = 30`, while released transformer commands intentionally default to `min_mapq = 0`
- Chromosome semantics are also intentionally non-uniform: the shared default selection is `chr1..chr22`, `--chromosomes all` resolves from the BAM header when possible, `bam-to-bam` sorts chromosomes lexicographically by default, and `frag-to-bam` builds its BAM header in `--chrom-sizes` order
- Current audit work has likely surfaced two real `gc-bias` bugs rather than missing tests: the fixed-size path appears non-invariant between aligned and misaligned tiling, and equivalent `--by-size` versus `--by-bed` windows appear not to agree in some cross-tile cases. Those should be treated as code fixes, not test rewrites, once the surrounding test setup noise is removed

That means the highest-value work is not "add 100 more small tests". It is to close the structural gaps that can still let scientifically wrong outputs pass.

## Priority 0: Structural test campaigns

These are the first things to build. They cut across commands and would have caught the kinds of bugs we have already found.

### P0.1 Tiling invariance campaign

Commands in scope:

- `fcoverage`
- `lengths`
- `midpoints`
- `gc-bias`
- `ref-gc-bias`

What must be proven:

- Output is invariant to tile size when all logical inputs are unchanged
- Aligned and misaligned tile boundaries produce the same final scientific result
- Halo/fetch narrowing does not lose eligible fragments near chromosome ends
- Crossing-window reduction produces the same answer as a hypothetical single-pass whole-chromosome computation
- Multi-chromosome processing does not reorder, duplicate, or drop results

Core test families:

- Same BAM/reference, many tile sizes, same final output
- Boundary-touching fragments around every relevant edge:
  - tile boundaries
  - window boundaries
  - chromosome ends
  - blacklist boundaries
- Aligned vs misaligned `--by-size` runs
- BED windows that overlap, touch, or straddle tile boundaries
- Sparse chromosomes and chromosomes with no windows

Why this is first:

- We have already found real reducer/tile bugs
- These commands all rely on chunking and merge logic
- Silent tile bugs are high-impact and hard to detect from superficial tests
- See also: [SA.4](#sa4) for oversized BED windows spanning 3+ tiles, and [SA.6](#sa6) for the separate NPZ-based cross-tile merge contract in `gc-bias`

### P0.2 Windowing semantics campaign

Commands in scope:

- `fcoverage`
- `lengths`
- `gc-bias`
- `ref-gc-bias`
- `bam-to-bam`
- `bam-to-frag`

What must be proven:

- `global`, `by-size`, and `by-bed` agree when they describe the same logical regions
- Overlapping BED windows obey the documented duplicate/merge semantics
- Window ordering, original indices, and chromosome grouping are preserved correctly
- Window assignment rules match the intended fragment semantics, not just the current implementation

Core test families:

- Equivalent BED and by-size windows produce equivalent aggregates
- BED windows with overlap, touching edges, and containment
- BED semantics matrix across released commands:
  - flatten-to-unique semantics in `ref-gc-bias` and `gc-bias`
  - unique-vs-indexed positional semantics in `fcoverage`
  - original-row preservation in `lengths`
  - inclusion-filter semantics in `bam-to-bam` and `bam-to-frag`
- Empty-window chromosomes and late chromosomes after empty ones
- Global-vs-window equivalence when there is exactly one window spanning the chromosome
- Assignment-mode parity tests for:
  - `count-overlap`
  - `any`
  - `all`
  - `midpoint`
  - `proportion`

Why this is first:

- The commands are defined by their window semantics
- A command can "run fine" and still be scientifically wrong if assignment or aggregation is off by one rule

### P0.3 Fragment semantics contract campaign

Commands in scope:

- `fcoverage`
- `lengths`
- `midpoints`
- `bam-to-bam`
- `bam-to-frag`
- `frag-to-bam`
- `gc-bias`

What must be proven:

- All commands agree on what a fragment is
- Paired-end and unpaired modes are consistent with their documented span definitions
- Inter-mate gap handling is correct where relevant
- Indel handling is correct in all supported modes and does not leak into unrelated computations
- Blacklist strategies behave the same way across commands that expose the same concept
- Even-length midpoint semantics are explicit and stable across released commands, including places where midpoint assignment and blacklist-midpoint filtering intentionally differ
- Shared blacklist loading semantics are stable wherever commands expose blacklist files: multi-file merge, touching-interval merge, minimum-size filtering, and any halo-expansion behavior must all match the intended scientific contract
- When a command intentionally uses different spans for different purposes, that contract stays stable: count span, exclusion span, weighting span, and reported fragment length must not silently drift apart

Core test families:

- Same synthetic fragments consumed by multiple commands, compared against hand-derived expectations
- Paired vs unpaired equivalence when representing the same physical fragment
- Indel-mode parity:
  - `ignore`
  - `adjust`
  - `skip`
- Blacklist strategy parity:
  - `any`
  - `all`
  - `midpoint`
  - `proportion`
- Boundary fixtures where an even-length fragment center lands exactly on a window or blacklist edge
- Shared blacklist fixtures with multiple BED inputs, touching intervals, and filtered short intervals
- Mate-gap edge cases for converters and coverage-like commands
- Span-contract fixtures where the same fragment is subjected to indel adjustment, midpoint assignment, mate-gap omission, blacklist filtering, and scaling so the intended “count span vs filter span vs weighting span” behavior is pinned explicitly

Why this is first:

- If fragment semantics drift between commands, all downstream interoperability becomes suspect

### P0.4 Weight composition campaign

Commands in scope:

- `fcoverage`
- `lengths`
- `midpoints`
- `bam-to-bam`
- `bam-to-frag`
- `gc-bias`

What must be proven:

- Genomic scaling and GC correction compose correctly
- Missing or invalid correction values follow the documented fallback or drop behavior
- Weighting does not silently change denominator logic in windowed outputs
- Mass conservation still holds when that is the intended contract
- The same scaling or GC artifact is interpreted with the same mathematical meaning by every released consumer that accepts it

Core test families:

- Scaling only
- GC only
- Scaling then GC together
- Invalid GC tags or invalid GC-file lookups
- Blacklist plus GC plus scaling together
- Rounded text outputs vs raw numeric outputs
- One shared fixture where `coverage-weights` and GC packages are consumed by multiple released commands and produce mutually consistent changes
- A scaling-semantics matrix fixture proving the intended differences:
  - per-base scaling in `fcoverage`
  - full-fragment averaging in `midpoints`, `bam-to-bam`, and `bam-to-frag`
  - overlap-aware vs full-fragment scaling in `lengths` depending on assignment mode

Why this is first:

- These are multiplicative modifiers with high scientific impact
- Bugs here often produce plausible-looking but wrong outputs
- See also: [SA.3](#sa3) for sparse-bin scaling blowups, and [SA.5](#sa5) for the concrete GC-fallback semantics already present in `bam-to-frag`

### P0.5 Cross-command interoperability campaign

Commands in scope:

- `bam-to-frag`
- `frag-to-bam`
- `bam-to-bam`
- `coverage-weights`
- `ref-gc-bias`
- `gc-bias`
- `fcoverage`
- `lengths`
- `midpoints`

What must be proven:

- Released commands work correctly together, not just in isolation
- Derived artifacts are accepted and interpreted correctly by downstream commands
- Real producer outputs, not just handcrafted stand-ins, survive full released workflows with the intended chromosome order, weighting semantics, and optional metadata intact

Core test families:

- `bam-to-frag -> frag-to-bam` roundtrip on realistic mixed cases
- `bam-to-bam` outputs are consumable by downstream commands and preserve intended semantics
- `coverage-weights` output changes downstream weighting as expected in `fcoverage`, `lengths`, and `midpoints`
- `ref-gc-bias` output drives `gc-bias` correctly
- `gc-bias` output drives `fcoverage` and `lengths` correctly
- Shared real-artifact fixture builders for:
  - one neutral `ref-gc-bias -> gc-bias` package
  - one non-neutral `ref-gc-bias -> gc-bias` package
  so consumer tests reuse the same produced artifacts instead of rebuilding the full producer
  pipeline independently in every file
- Workflow spine:
  - `ref-gc-bias -> gc-bias -> lengths`
  - `ref-gc-bias -> gc-bias -> fcoverage`
  - `bam-to-bam -> lengths`
  - `bam-to-frag -> frag-to-bam -> bam-to-bam` parity checks where representable
- `coverage-weights -> bam-to-bam` and `coverage-weights -> bam-to-frag`, so real producer TSVs are exercised through the released transformers and not only through the counters
- `bam-to-frag -> frag-to-bam -> lengths` and `bam-to-frag -> frag-to-bam -> fcoverage`, so external frag ingestion is validated as part of the released counting workflows

Why this is first:

- Users will chain these commands
- Release quality is determined by workflow reliability, not isolated unit correctness

### P0.6 Artifact contract campaign

Commands in scope:

- `ref-gc-bias`
- `gc-bias`
- `coverage-weights`
- `bam-to-frag`
- `bam-to-bam`
- `frag-to-bam`
- `midpoints`

What must be proven:

- Every released command that writes a reusable artifact writes the exact schema and metadata its downstream consumers expect
- Schema fields reflect actual run settings, not stale defaults or partially applied options
- Downstream commands produce numerically different results when the input artifact meaningfully changes
- Artifact readers fail clearly on malformed or shape-incompatible inputs
- Consumers enforce compatibility constraints consistently when an artifact is version-mismatched, range-mismatched, or shape-incompatible
- Chromosome order, axis order, and companion metadata are stable and intentional across artifacts, because different released commands deliberately preserve requested order, sort explicitly, or follow chrom-size/header order

Core test families:

- `ref-gc-bias` package loaded by real `gc-bias`, not just by helper readers
- `gc-bias` package loaded by real downstream consumers using both file-based and tag-based weighting paths where relevant
- `coverage-weights` TSV consumed by `fcoverage`, `lengths`, `midpoints`, `bam-to-bam`, and `bam-to-frag`
- `bam-to-frag` companion header consumed by `frag-to-bam`
- `midpoints` group-index TSV checked against the actual 3D array axis contract
- Real producer helpers should live in shared test fixtures when the produced artifact semantics
  are already proven elsewhere; downstream consumer tests should focus on their distinct behavior,
  not duplicate producer setup
- Corrupt or shape-mismatched artifacts rejected with specific actionable failures
- GC package version, edge-vector, and length-range compatibility failures checked through the real released commands that load them
- Reference-GC metadata fields (`skip_smoothing`, `skip_interpolation`, support masks, width corrections) round-tripped through real producer and consumer commands
- Exact reference-GC package field and shape contract checked through the real producer and consumer path, because this artifact currently has no explicit schema version field
- Scaling TSV compatibility failures checked through the real commands that load them: contig mismatch, chromosome subset mismatch, uncovered requested spans, malformed factors, and order-sensitive multi-chromosome fixtures
- Chromosome and order contracts checked through real artifacts: `coverage-weights` row order, `frag-to-bam` BAM header order from `--chrom-sizes`, `midpoints` group-axis order, and released consumer behavior under `--chromosomes all` vs explicit chromosome lists

Why this is first:

- Several released commands are primarily valuable because they feed other released commands
- A file can be syntactically valid while still encoding the wrong semantics

### P0.7 Default semantics campaign

Commands in scope:

- `fcoverage`
- `lengths`
- `gc-bias`
- `coverage-weights`
- `bam-to-bam`
- `bam-to-frag`
- `frag-to-bam`
- `midpoints`
- `ref-gc-bias`

What must be proven:

- Silent defaults that materially affect scientific output are pinned by explicit tests, not left as accidental behavior
- Commands that expose the same concept but intentionally choose different defaults are tested for both equivalence and intentional difference
- Default chromosome selection, chromosome order, and MAPQ thresholds are treated as release behavior, because they change which fragments are included and in what order artifacts are emitted

Core test families:

- Window-default matrix:
  - `gc-bias` no window flags equals explicit `--by-size 100000`
  - `fcoverage` and `lengths` no window flags equal explicit global mode
  - `bam-to-bam` and `bam-to-frag` no `--by-bed` behave as full-chromosome/global selection
- Chromosome-default matrix:
  - omitted chromosome selection uses `chr1..chr22`
  - `--chromosomes all` uses BAM-header order where applicable
  - `frag-to-bam` uses `--chrom-sizes` order for the BAM header and output iteration
  - `bam-to-bam` sorts lexicographically by default while other commands preserve the requested order unless documented otherwise
- Filter-default matrix:
  - released counting commands default to `min_mapq = 30`
  - released transformer commands intentionally default to `min_mapq = 0`
  - default blacklist strategy and default paired-end assumptions produce the intended cross-command behavior on the same fixture

Why this is first:

- A release can look correct in explicit test configurations while still being wrong in the actual default user path
- These defaults directly change inclusion, normalization, and artifact order, not just CLI ergonomics

## Priority 1: Command-specific high-risk gaps

### `ref-gc-bias`

Status:

- Most under-tested released command by far
- Very small current suite compared with complexity
- No meaningful active evidence that `run()` itself is being validated end to end

Needs:

- Full command-level tests for sampling, seeding, blacklist masking, by-bed restriction, smoothing, support masking, interpolation, and output package integrity
- Determinism tests with fixed seeds
- Determinism tests across `n_threads` and tile sizes when the seed is fixed
- Tile-size invariance tests
- Tests that compare expected support-mask behavior on crafted masked references
- Tests that prove end offsets and effective lengths are respected end to end
- Tests for the `sampling_density > 1.0` rejection path and other configuration guardrails that affect scientific validity
- Tests proving BED windows are merged/touch-flattened exactly as intended before counting
- Tests that load the written NPZ package and verify metadata fields consumed by downstream `gc-bias`
- Tests that prove the per-Mb support threshold behaves as intended when total covered ACGT support is tiny vs ample
- Tests that the sampled start-position process is independent of windowing and blacklisting as intended, while the retained counts/support still change in the expected way after those masks are applied
- Real consumer tests where `gc-bias` ingests the package and the observed correction changes in the expected direction
- Producer/consumer contract tests where the written package is loaded by the real `gc-bias` path, not just by helper readers, and metadata toggles change downstream behavior exactly as intended
- Exact package field-presence and scalar-shape tests for the written reference-GC package, because unlike the downstream GC correction package this artifact currently has no schema-version guard

Why high priority:

- It feeds the GC correction pipeline
- If it is wrong, downstream corrections can all be wrong in a coordinated way
- See also: [SA.2](#sa2) for the support-threshold step-function that can make test-scale and production-scale behavior diverge

### `gc-bias`

Status:

- Good helper-level coverage
- Weaker end-to-end coverage than the command complexity justifies
- Little active evidence that `run()` is covered across its main branches

Needs:

- Full run tests covering:
  - global vs windowed bias estimation
  - by-size aligned vs misaligned tiles
  - fixed-size streaming windows vs general window handling when they describe the same logical regions
  - crossing-window spillover reduction vs contained windows
  - all outlier methods and scopes
  - window weighting modes
  - minimum ACGT support thresholds per window
  - support-mask and interpolation interactions
  - end-offset interactions
  - paired/unpaired behavior
  - output package structure and downstream loadability
- Tests that prove reference-package metadata (`skip_smoothing`, `skip_interpolation`, masks, width corrections) is respected by the command
- Tests for `save_intermediates` that verify intermediate artifacts are coherent and correspond to the final matrix
- End-to-end comparisons against hand-built tiny reference/bias scenarios
- Stronger tests around the generated interpolation logic, which is explicitly marked for validation in code
- Tests that malformed or shape-incompatible reference packages fail early with clear messages
- Tests that `by-bed` overlapping/touching windows are flattened exactly once before counting
- Real workflow tests where the produced correction package measurably changes downstream `fcoverage` and `lengths` outputs in the expected way
- Tests that the produced GC package is accepted consistently by all released consumers (`fcoverage`, `midpoints`, `bam-to-bam`, `bam-to-frag`, and `lengths` after length marginalization)
- Compatibility-failure tests for GC package schema version, fragment-length range, and end-offset constraints through the real commands that load the package
- Command-level tests that the default no-window-flag behavior equals explicit `--by-size 100000`, and differs from explicit global mode only when the intended window-scaling semantics make it differ
- Full command-level tests for outlier configuration (`none`, quantile, IQR, stddev, MAD), outlier scope, and `save_intermediates`, since these are exposed scientific knobs with very little active run-level evidence
- Command-level tests for `min_window_acgt_pct` thresholds and the failure mode where no usable windows remain after filtering/blacklisting
- Likely active bug to fix after test cleanup: fixed-size `gc-bias` windowing appears to disagree across aligned vs misaligned tile sizes and across logically equivalent `--by-size` and `--by-bed` window definitions when tiles split the windows

Why high priority:

- This is a mathematically dense command with many tuning knobs
- Local helper coverage does not substitute for command-level validation
- See also: [SA.6](#sa6) for the NPZ-based cross-tile merge shape contract that sits underneath the command-level tiling story

### `coverage-weights`

Status:

- Some command coverage
- Still thin relative to the edge behavior of triangular smoothing
- Current tests appear to emphasize simple one-chromosome derivations more than whole-output invariants

Needs:

- Chromosome-edge truncation tests for the triangular weighting kernel
- Blacklist-aware scaling tests
- Multi-chromosome output ordering tests
- Multi-chromosome global-mean normalization tests to prove per-chromosome coverage is normalized against one shared denominator as intended
- Multi-chromosome tests where one chromosome is mostly zero, to lock in the intended exclusion of zero-coverage bins from the denominator
- Last-bin truncation tests where the final stride bin is shorter than `stride`, so length weighting cannot silently skew the global mean
- Zero-coverage region behavior tests
- Large `bin_size/stride` combinations beyond the simplest hand-derived case
- Paired vs unpaired mode parity where the underlying physical fragments are the same
- Interoperability tests with downstream consumers
- Consumer-consistency tests proving the written scaling TSV has the same practical effect in `fcoverage`, `lengths`, `midpoints`, `bam-to-bam`, and `bam-to-frag` where mathematically comparable
- A command-matrix test that locks in the intentional differences in how those consumers average scaling across fragment spans, overlaps, and per-base coverage
- Real producer-consumer tests where the exact TSV written by `coverage-weights` is loaded by every released scaling consumer, rather than relying mostly on handcrafted scaling fixtures
- Compatibility-failure tests for scaling TSVs through the real loader path: missing chromosomes, contig-length mismatches, uncovered requested spans, and malformed or non-finite factors

Why high priority:

- This command defines downstream normalization factors
- Small edge mistakes propagate into all smoothed analyses
- See also: [SA.3](#sa3) for the specific near-zero non-zero coverage case that can explode scaling factors

### `midpoints`

Status:

- Some active tests around grouping and fetch narrowing
- Still structurally thinner than `lengths` and `fcoverage`
- Current tests do not yet strongly prove whole-profile correctness under combined scaling, GC weighting, and blacklisting

Needs:

- Full end-to-end profile correctness tests with:
  - blacklist filtering
  - genomic scaling
  - GC weighting from file
  - GC weighting from tag
  - paired/unpaired modes
  - multiple groups and length bins
- Tie-handling tests for even-length midpoints with deterministic expectations
- Tile invariance and chromosome-end halo tests beyond the currently narrow cases
- Tests for the 3D view/layout contract in `counting_by_group.rs`, which is explicitly marked `TODO: Test!!`
- Tests proving the merged temporary tile arrays are numerically identical to a single logical accumulation
- Tests that the written group-index TSV and the first axis of the output array stay in exact agreement under nontrivial group ordering
- Tests that group-axis order remains stable under nontrivial BED encounter order, chromosome order, thread counts, and tile sizes, because the group index file is the only downstream contract for interpreting the array
- Tests that GC-from-tag and GC-from-file paths each produce the expected weighted profile on the same fixture
- Explicit command-level tests for even-length midpoint ties at window edges, including the fact that blacklist-midpoint filtering uses a different midpoint rule than profile placement
- Real GC-package consumer tests using packages generated by `gc-bias`, not only handcrafted matrices

Note:

- Strong midpoint tie tests probably require an explicit deterministic seam around midpoint randomization. Without that, the tests will stay weaker than they should be

Why high priority:

- Group-collapsed profiles are easy to get mostly right while still mis-indexing length or position axes
- See also: [SA.1](#sa1) for the concrete ndarray stride/view parity test that proves the written array is interpreted correctly

### `bam-to-frag`

Status:

- Thin active suite compared with released status
- Current active tests seem concentrated on smoke/global/BED chromosome handling, not on the command's main scientific branches

Needs:

- Output-column semantics tests when GC and scaling are absent, present individually, and present together
- Paired vs unpaired output equivalence where possible
- Blacklist behavior tests for all supported blacklist strategies
- GC correction tests from file, including invalid-GC drop/fallback behavior
- Scaling-factor tests and GC+scaling combined tests
- Blacklist and window filtering parity with `bam-to-bam`
- Multi-chromosome ordering and sorting guarantees
- Header-file contract tests for the written `.frag.header.tsv`
- More downstream interoperability checks with `frag-to-bam`
- Tests that row arity always matches the companion header exactly
- Tests that produced rows stay coordinate-sorted within chromosome after the bounded sorter flushes
- Roundtrip tests where `bam-to-frag -> frag-to-bam` preserves the intended `GC`/`COV`/`FLEN` information where representable
- Parity tests against `bam-to-bam` using the same scaling TSV and GC package, so the emitted frag weights match the BAM tags on the same physical fragments
- Real GC-package consumer tests using packages generated by `gc-bias`, not only handcrafted matrices
- Real scaling-artifact consumer tests using TSVs generated by `coverage-weights`, not only handcrafted scaling fixtures

Why high priority:

- This is a format boundary command
- Format boundary commands need strong contract tests, not just smoke tests
- See also: [SA.5](#sa5) for the specific GC-file fallback behavior that must be pinned before release

### `bam-to-bam`

Status:

- Moderate suite, but still not broad enough for a release-quality transformer

Needs:

- GC-tag writing tests, not just scaling/COV tag tests
- Combined `GC` + `COV` + `FLEN` tag scenarios
- Combined filter/correction/tag scenarios
- Unpaired-mode behavior where supported
- Stronger parity tests against `bam-to-frag` and downstream counting commands
- Output sort-order and chromosome-order guarantees under mixed selection settings
- Invalid-GC drop/fallback semantics tested at command level, not just assumed from logging
- Tests that both mates of the same fragment always receive identical `GC`, `COV`, and `FLEN` tags
- Tests that downstream consumers (`lengths`, `fcoverage`) see the expected numerical change after BAM tagging workflows
- Parity tests against `bam-to-frag` using the same scaling TSV and GC package, so both released transformers encode the same weights for the same fragments
- Real GC-package consumer tests using packages generated by `gc-bias`, not only handcrafted matrices
- Real scaling-artifact consumer tests using TSVs generated by `coverage-weights`, not only handcrafted scaling fixtures

Why high priority:

- Wrong outputs here pollute downstream workflows while still looking like valid BAM files

## Priority 2: Still-important work on commands that already look stronger

These commands already have more active coverage, so the goal is not to flood them with shallow tests.

### `lengths`

Needs:

- Cross-command parity with `fcoverage` and `bam-to-bam` on shared fragment/filter semantics
- More combinatorial tests for scaling + GC + blacklist together
- More invariance tests across many tile sizes
- More seeded/deterministic tests around midpoint/proportion style assignment interactions where relevant
- Artifact-driven tests where real `coverage-weights` and `gc-bias` outputs are consumed rather than mocked
- Explicit tests around even-length midpoint window assignment near bin boundaries, and around how that differs from blacklist-midpoint filtering
- Explicit span-contract tests showing that `indel_mode=adjust` changes the counted length bin while blacklist, scaling, and GC correction still follow the full reference fragment span
- Stronger tests for the length-marginalization contract: the same GC package should drive different expected results under `equal`, `coverage`, and `max-coverage`, and those should remain compatible with full-matrix consumers on the same fixture
- Real GC-package consumer tests using packages generated by `gc-bias`, not only handcrafted matrices

### `fcoverage`

Needs:

- Cross-command parity with `lengths`, `coverage-weights`, and `gc-bias`
- Stronger tests around aligned size-mode finals vs reduced outputs
- More combined-scenario tests for blacklist + scaling + GC corrections
- Artifact-driven tests where real `coverage-weights` and `gc-bias` outputs change per-position and per-window results in the expected way
- Fast-path parity tests proving the aligned by-size final-writer path matches the general reducer path under blacklist masking, scaling, GC correction, and chromosome-end clipping
- Explicit span-contract tests for `--ignore-gap`, proving that mate-gap omission changes only the counted coverage support while fragment inclusion, positional masking, and scaling semantics remain the intended ones
- Real GC-package consumer tests using packages generated by `gc-bias`, not only handcrafted matrices

### `frag-to-bam`

Needs:

- More true workflow roundtrips, not just syntax/header permutations
- Tests that downstream released commands consume the produced BAM exactly as intended
- Roundtrip tests starting from real `bam-to-frag` outputs, including optional `GC`/`COV`/`FLEN` columns
- Stronger order/blacklist/filter workflow tests so the converter is validated as part of a release pipeline, not just as a parser
- Explicit downstream parity tests where `frag-to-bam` output feeds `bam-to-bam`, `lengths`, and `fcoverage`, proving that external frag ingestion preserves the intended chromosome order, blacklist behavior, and optional tag semantics through the released workflows

## Second-opinion additions: gaps found by direct source review

The following were identified by reading the actual algorithm implementations, not just
the architectural descriptions. They are ordered by scientific impact.

<a id="sa1"></a>
### SA.1 `midpoints`: `view_ndarray3_group_len_pos` vs `to_3d_group_len_pos` parity

**Source**: `counting_by_group.rs` lines 304 and 343. The code itself marks the ndarray3
view with `// TODO: Test!!`.

What must be proven:

- `view_ndarray3_group_len_pos()` and `to_3d_group_len_pos()` return identical values for
  the same data, using a non-trivial fixture where all three axes have more than one element,
  groups, length bins, and positions are all different sizes, and the counts are not symmetric
- The zero-copy stride math `(group*p*l, 1, l)` for shape `(g, l, p)` reinterprets the
  internal `(group, position, length_bin)` flat buffer correctly

Why this matters above what is already in the plan:

The plan mentions the 3D view contract but does not pin the specific test that proves the
zero-copy non-contiguous stride view and the allocating copy produce the same array. A
bug in the stride math would produce silently wrong midpoint profiles for every downstream
consumer of the written NPY file. Adding this parity test is the minimal proof of the
`TODO: Test!!` marker in the code.

<a id="sa2"></a>
### SA.2 `ref-gc-bias`: support-threshold integer-division plateau

**Source**: `ref_gc_bias.rs` line 238: `let threshold_per_mb = 1 + opt.n_positions / 100000000;`

What must be proven:

- For any `n_positions < 100_000_000`, integer division yields 0, so `threshold_per_mb`
  is always 1 regardless of how small `n_positions` is
- The "per-Mb" semantics are only activated once `n_positions >= 100_000_000`
- Tests should verify the threshold value at the exact crossover boundaries:
  `99_999_999`, `100_000_000`, and `200_000_000`
- Tests should verify that support-mask coverage changes in the expected direction as
  `n_positions` crosses those boundaries

Why this matters above what is already in the plan:

The plan mentions support-threshold behavior but treats it as a continuous parameter. The
integer division creates a step-function where small `n_positions` values (including those
used in tests) all produce the same threshold. This is an invisible calibration difference
between test-scale and production-scale runs that could mask support-mask bugs entirely
in the test suite.

<a id="sa3"></a>
### SA.3 `coverage-weights`: near-zero non-zero bins producing very large or infinite scaling factors

**Source**: `striding.rs` lines 274-285. The `inverter` closure does `1.0 / x` for any
non-zero `x`. The caller skips `NaN`, `inf`, and exactly zero, but not near-zero values.

What must be proven:

- When a chromosome end has a stride-bin with very low but non-zero average coverage
  (e.g., 1-3 fragments in a 500kb bin), the inverted scaling factor is large but finite
  in the intended range
- Downstream consumers (`fcoverage`, `lengths`, `midpoints`) handle very large scaling
  factors gracefully: either the scientific result is intentional, or there is a cap
- At the f32 boundary: `avg_overlap_coverage` values that round to the smallest
  non-zero f32 when cast from f64 produce an inverted f32 of `inf`, which propagates
  into weighted sums downstream silently
- Tests should cover: single-read stride bin at chromosome end, inverted factor written
  to TSV, downstream consumer reads that TSV, output is finite or the behavior is pinned

Why this matters above what is already in the plan:

The plan covers zero-coverage regions and last-bin truncation. It does not cover the
near-zero non-zero case where the inversion is mathematically valid but produces a
scaling factor of, say, 80,000x or worse infinity as f32. A fragment in that stride bin
will dominate any weighted output by an arbitrary multiple of the genome-wide mean.
This is a real scenario at chromosome ends with shallow coverage.

Priority note:

- The first thing to pin here is the realistic sparse-coverage regime: very low but
  non-zero stride bins that produce huge yet finite scaling factors and can dominate
  downstream weighted outputs
- The explicit smallest-non-zero-`f32` to `inf` boundary is still worth testing, but it
  is a second-order numerical guardrail after the scientifically realistic huge-but-finite
  amplification case

<a id="sa4"></a>
### SA.4 `fcoverage` reducer: BED windows wider than 2 × tile size (spanning 3+ tiles)

**Source**: `reducer.rs` lines 232-264. `expected_contribs` is built by counting how many
`cross_index` files list a given `orig_idx`. A window spanning N tiles crosses N-1
boundaries, so it appears in N-1 cross-index files and gets `expected_contribs = N-1`.
But the window produces partial rows in N tiles. The accumulator is emitted after
`N-1` contributions, then the N-th contribution re-inserts the index into `accum_by_idx`,
which is caught by the safety check at line 364 as an error rather than a silent undercount.

What must be proven:

- A BED window strictly wider than 2 × `tile_size` (e.g., a 30Mb window with 10Mb tiles)
  either (a) produces the correct aggregated value through whatever mechanism the
  tiling layer uses for this case, or (b) triggers a clear error that informs the user
- The `expected_contribs` count correctly matches the actual number of partial rows the
  window generates for windows spanning exactly 2 tiles, exactly 3 tiles, and 4+ tiles
- The safety check (`accum_by_idx.is_empty()`) fires with an actionable error message
  rather than silently passing or panicking if the count is wrong

Why this matters above what is already in the plan:

The tiling campaign in P0.1 covers aligned/misaligned `--by-size` and tile-boundary
straddling but does not explicitly exercise a single BED window wider than 2 tile sizes.
Real user BED files (large TADs, chromosomal arm-scale windows) frequently exceed 20Mb.
If the cross-index count logic is off-by-one for spanning windows, the safety check
catching the error is the only guard. Verifying this path is concrete, low-cost, and
covers a plausible real-world input shape.

<a id="sa5"></a>
### SA.5 `bam-to-frag`: GC correction fallback is weight=1.0, not a drop

**Source**: `bam_to_frag.rs` lines 380-388. When `gc_file` is set and `correct_fragment`
returns `None`, the code sets `gc_weight = Some(1.0)` and increments
`gc_failed_fragments`, not `filtered_fragments`. The fragment is still written.

What must be proven:

- A fragment whose GC weight cannot be computed (length outside the GC matrix, or
  reference base ambiguity) contributes weight 1.0 to the frag output, not zero,
  and not a missing row
- The `gc_failed_fragments` counter matches the number of rows written with `gc_weight=1.0`
- The companion header file still lists `gc_weight` as a column when `gc_file` is set,
  even when all GC lookups fail
- `frag-to-bam` and downstream consumers that read the frag file receive 1.0 weights for
  those rows without any further silent fallback
- Parity test: `bam-to-bam` and `bam-to-frag` should agree on whether a fragment with
  a missing GC weight is included or excluded, and at what weight

Why this matters above what is already in the plan:

The plan mentions GC fallback behavior under `bam-to-frag` and `bam-to-bam` but calls it
"drop/fallback behavior" generically. The current implementation silently upgrades failed
GC lookups to weight 1.0 rather than dropping. This is a meaningful scientific choice:
it means GC-corrected outputs include unweighted fragments instead of excluding them,
which inflates their contribution relative to fragments with valid weights. Whether this
is intentional must be pinned by an explicit test before release.

<a id="sa6"></a>
### SA.6 `gc-bias`: NPZ-based cross-tile merge shape contract

**Source**: `gc_bias/cross_tile_parts.rs`. The cross-tile GC accumulation serializes
`GCCounts` arrays into a per-tile NPZ file, then merges them. The merge assumes every
tile's `GCCounts` was initialized from the same template (same shape).

What must be proven:

- When tiles span different chromosome regions (some with dense coverage, some sparse
  or entirely empty), the per-tile `GCCounts` arrays are initialized with the same shape
  before merging
- A tile that processes zero fragments still produces a zero-filled `GCCounts` of the
  correct shape (not a missing or differently shaped array), so the NPZ merge does not
  silently drop its contribution
- Cross-tile NPZ files from misaligned tiles (windows straddle a tile boundary) merge
  into the same final GC matrix as a non-tiled run

Why this matters above what is already in the plan:

The plan covers gc-bias tiling invariance at the command level but does not specifically
address the NPZ-based cross-tile merge path. The fcoverage and gc-bias cross-tile
mechanisms are structurally different (text partials vs NPZ arrays). A shape mismatch
in the NPZ merge for an empty tile would produce a silent count underestimate for
windows near tile boundaries, which is exactly the kind of bug the prior tile-alignment
fix was needed to address.

## Priority 3: Small and easy checks to do later

These should be done, but only after the structural work above.

- CLI error message quality
- Help text coverage beyond smoke-level presence
- Prefix/output naming consistency
- Compression-format combinations
- Temp-directory cleanup and failure cleanup behavior
- Minor parser edge cases
- Header/autodetection polish checks

## Concrete execution order

1. Build the cross-command tiling/windowing invariance harness
2. Build the fragment semantics contract harness shared by multiple commands
3. Build the artifact-contract harness for reusable packages, headers, and index files
4. Add full end-to-end suites for `ref-gc-bias` and `gc-bias`
5. Deepen `coverage-weights` and `midpoints`
6. Strengthen converter interoperability: `bam-to-bam`, `bam-to-frag`, `frag-to-bam`
7. Backfill `lengths` and `fcoverage` with cross-command parity and composition cases
8. Only then spend time on CLI/polish cases

## What to avoid

- Do not start with clap type parsing or help formatting
- Do not measure progress by test count
- Do not add unit tests that just restate implementation details without checking command intent
- Do not trust helper-level tests as proof that a whole command is correct

## First concrete deliverables from this plan

- A reusable tile-invariance fixture harness for released counting commands
- A reusable fragment-semantics fixture set shared across converters and counters
- A reusable artifact-contract harness for generated NPZ, TSV, header, and BAM-tag outputs
- A full `ref-gc-bias` command-level suite
- A full `gc-bias` command-level suite
- A workflow test that chains:
  - `ref-gc-bias`
  - `gc-bias`
  - `coverage-weights`
  - `fcoverage`
  - `lengths`

That workflow should become the release-confidence spine for the public command set.
