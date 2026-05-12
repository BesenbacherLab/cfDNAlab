# `cfdna midpoints` review

Date: 2026-04-24

Scope: `src/commands/midpoints/*`, the CLI dispatch for `midpoints`, and directly used shared helpers for BED loading, tiling, overlap lookup, GC correction, scaling, blacklist checks, midpoint placement, and grouped count merging. I also read the existing midpoint-focused tests in `tests/test_profile_groups_command.rs`, `tests/test_profile_groups_counts.rs`, `tests/test_tiling.rs`, `tests/test_cli_smoke.rs`, and cross-command roundtrip/artifact tests. I did not run tests.

Shared findings that affect this command:

- None active.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release semantic decision:

- None active. M-005 was resolved by defining shared midpoint blacklisting as central-base support: odd fragments check the unique central base, and even fragments check either central base.

Operational guardrail, only if oversized dense outputs are in the three-day target:

- M-D-001: dense final shape needs fail-fast checked allocation instead of late panic/OOM behavior.

Post-release performance/scalability:

- G-006: sparse-window GC reference pruning.
- M-001B: optional sparse final output would help very large group, width, and length-bin shapes.

## Findings

### M-001B - Medium - Very large final profile shapes need an optional sparse output format

The final output is intentionally dense today: `<prefix>.midpoint_profiles.npy` stores counts shaped `(group, length_bin, position)` after collapsing all intervals in each group. Dense output is straightforward for downstream numpy workflows and should remain the default for ordinary profile sizes.

Impact: dense final output becomes impractical when users request high positional resolution, many groups, and near-base-pair fragment length bins. For example, `2001 * 3000 * 971` cells is about 5.83 billion `f32` values, or roughly 21.7 GiB for the final counts alone. Sparse internal tile partials solve per-tile memory and temp disk pressure, but they do not solve the final artifact size when the requested output tensor itself is too large.

Recommended tasks:

- Add an explicit sparse final output option rather than changing the default dense artifact.
- Choose and document a stable sparse layout. A practical first version is a SciPy-compatible sparse `.npz` with rows representing `group_idx * num_length_bins + length_bin_idx` and columns representing `position`, plus sidecar metadata for group names and length-bin edges.
- Preserve zero-count groups in the sparse matrix shape so group indices remain stable even when a group has no observed counts.
- Keep dense plotting tied to dense output, or require densifying only selected groups for plotting. Do not accidentally densify the full sparse output just to make plots.
- Make the CLI estimate compare dense final size with expected sparse payload size when possible. If expected density is unknown before counting, report the dense final size and explain that sparse output avoids storing zeros in the final artifact.
- Add roundtrip tests that load the sparse final output, densify it in test code, and compare it against the dense `.npy` output for the same small fixture.
- Document downstream loading in Python, including how to reconstruct `(group, length_bin, position)` from the sparse two-dimensional representation.

## Sparse final output design note

The urgency of M-001B depends heavily on the intended axis sizes. In one recent benchmark, `--length-bins 30:1001:50` produced a dense midpoint profile file of about 191 MB. That range resolves to 20 length bins. A full 1 bp length axis for 80-600 bp inclusive would use `--length-bins 80:601:1`, which resolves to 521 bins. Holding groups and profile width fixed, that scales the dense output by about `521 / 20 = 26.05x`, or roughly 5.0 GB for the same sample. That is large but still manageable on cluster storage and memory for many single-sample workflows.

Sparse final output is therefore less urgent for one sample at this shape, but it still matters for cohort workflows with thousands of samples side by side, higher group counts, wider windows, or full 30-1000 bp resolution. Dense should remain the default because it is easy to load and slice in NumPy. A sparse option should optimize downstream ergonomics, not only file size.

Recommended sparse-output shape:

- Add an explicit option such as `--save-as-sparse`.
- Write a SciPy-compatible CSR `.npz`, e.g. `<prefix>.midpoint_profiles.sparse.npz`.
- Store rows as `group_idx * num_length_bins + length_bin_idx` and columns as `position`.
- Preserve empty groups and empty length bins in the matrix shape so group and bin indices remain stable.
- Keep sidecar metadata for group names and length-bin edges, likely reusing `<prefix>.group_index.tsv` and adding a small length-bin sidecar if needed.
- Document Python loading with `scipy.sparse.load_npz`, plus reconstruction of `(group, length_bin, position)` from sparse row and column indices.
- Keep plotting dense-only by default, or densify only the selected rows needed for plotting. Do not densify the full sparse matrix just to make plots.

# Claude

Date: 2026-05-03

Scope: deep dive into `src/commands/midpoints/{midpoints,counting_by_group,config,windows,plotting,mod}.rs` plus the directly-used shared helpers (`shared/midpoint.rs`, `shared/overlaps.rs`, `shared/scale_genome.rs`, `shared/window_fetch.rs`, `shared/tiled_run.rs`, `shared/length_axis.rs`, `shared/bed.rs`, `commands/cli_common.rs`). For comparison I cross-checked the analogous fragment-streaming loops in `lengths.rs` and `ends.rs`. I did not run tests.

## Findings


## Items I checked and ruled out

These were investigated but I did not find evidence of an actual problem:

- **`scaling_with_bin_idx` placeholder index of 0 plus `OverlappingWindow.idx` reuse**: works correctly because `find_overlapping_windows` returns the scan position (not the carried `idx`) in BED mode ([overlaps.rs:269-300](../../src/shared/overlaps.rs#L269-L300)), and the scaling vector is built in the same order as `scaling_chr`. The contract is documented in the inline comment at [midpoints.rs:438-444](../../src/commands/midpoints/midpoints.rs#L438-L444) and matches the same pattern used in `ends.rs` and `lengths.rs`.
- **`wd_ptr` and midpoints going backward across consecutive fragments**: `find_overlapping_windows` uses `query.start.saturating_sub(look_back)` with `look_back = max_fragment_length`, which strictly upper-bounds the maximum decrease in midpoint between two start-sorted fragments (`(max_len - min_len)/2 ≤ max_len`). No off-by-one in the streaming pointer.
- **`compute_per_window_scaling_over_fragment` filter pruning midpoint-overlap windows**: cannot occur. Any window returned by the midpoint query satisfies `window.end > fragment_start` and `window.start < fragment_end` because the midpoint lies inside both intervals.
- **Sparse merge correctness**: `read_sparse_profile_partial_file` validates length agreement, strictly ascending unique indices, `usize` fit, in-bounds indices, and shape match before any merge work begins ([counting_by_group.rs:554-632](../../src/commands/midpoints/counting_by_group.rs#L554-L632)). The chunked-Mutex merge in `merge_sparse_profile_partial_file` correctly partitions by `flat_idx / chunk_size` and only holds one chunk lock at a time, with the wrap-around split via `partition_point` preserving full coverage of every entry.
- **`gc_corrector.clone()` per tile**: clones an `Array2<f64>`, which is non-trivial, but the same pattern is used in `ends`, `lengths`, and `fcoverage` (see explicit comment "Quite small memory footprint" at [fcoverage.rs:562](../../src/commands/fcoverage/fcoverage.rs#L562)). Not a midpoints-specific concern.
- **`incr_weighted` casting `f64 -> f32`**: realistic GC × scaling products stay far below `f32::MAX`; no overflow in practice and consistent with the dense `.npy` output type.
- **Even-fragment blacklist midpoint vs counted midpoint divergence**: tracked as M-005 above and now resolved with central-base support semantics.

## Codex comments on Claude findings

Overall: the retained Claude finding is a low-severity performance cleanup rather than an output-count correctness issue.

- **M-C-003**: Agree, but this is a small optimization. If changed, keep the helper readable because the fetch-narrowing semantics are more important than saving one short per-tile pass.
- **Items ruled out**: Mostly agree. The note about `wd_ptr` only rules out backward movement across fragment order inside one window slice. It does not cover the separate tile-local versus chromosome-wide pointer bug found during the tile-size comparison, where `core_overlapping_windows` is a compact vector and must use a tile-local pointer starting at zero.

## Existing coverage notes

The command already has broad integration coverage: length-bin parsing, default MAPQ, paired/unpaired parity, group index ordering, even-length midpoint edge placement, blacklist midpoint behavior, real GC packages, GC tags, GC/scaling multiplication, scaling TSV validation, tile-boundary behavior, chromosome-end fetch narrowing, CLI smoke output, and cross-command BAM/fragment roundtrips are all represented.

The most important midpoint-specific missing test from this review is sparse final-output roundtrips if M-001B is implemented. The chosen even-midpoint blacklist semantics now have command-level regression coverage. The deferred skipped-tile GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.

## Re-review additions (2026-05-04)

The shared raw-chromosome temporary filename issue (G-019) and overlong `--gc-tag` issue (G-021) originally noted here have since been implemented. `midpoints` now uses shared temporary chromosome tokens and shared AUX-tag validation.

### Release triage additions

Pre-release correctness/safety:

- None active from this re-review.

Post-release performance:

- G-006 remains the only midpoint-specific performance item from this pass.

### Command-specific findings

No new midpoints-only counting correctness finding was added in this pass. The shared input/output contract issues found during the re-review have since been implemented.

## Fundamental-only re-review (2026-05-10)

Scope: re-read `src/commands/midpoints/{midpoints,counting_by_group,config,windows,plotting}.rs`, the directly used interval, tiling, BED, blacklist, GC, scaling, fragment, BAM, and artifact helpers, the midpoint spec, and the existing midpoint integration/unit tests. I did not run tests.

Result: I did not find a new midpoint-local fundamental count-correctness bug in this pass. The earlier M-D backlog was too broad for this review target: G-025 was retired after `MAX_SUPPORTED_FRAGMENT_LENGTH` was lowered to 50,000 bp, and the retained M-D entries are hardening, performance, or artifact-consistency notes rather than release-blocking midpoint-count defects.

Core invariants rechecked:

- Tile ownership is by the sampled midpoint falling in the tile core, so the same fragment can be fetched by neighboring halo tiles without double-counting the midpoint.
- Fetch narrowing from core-overlapping windows remains safe because a fragment of length at most `max_fragment_length` whose midpoint can hit one of those windows must lie inside the narrowed fetch band.
- The streaming overlap pointers remain safe despite paired-fragment yield order because fragment starts and midpoints cannot move backward by more than `max_fragment_length` between records yielded from coordinate-sorted pairs.
- Fragment length indexing goes through the shared `LengthAxis`, and out-of-range lengths are dropped before counting.
- GC and scaling weights are applied once per selected midpoint-window count; full-fragment coverage scaling is intentionally averaged over the entire fragment span, not just the 1 bp midpoint.

Remaining fundamental decisions:

- M-D-001 is an operational guardrail, not a counting bug. It matters before release only if the expected target runs include dense final tensors large enough to risk allocator aborts or late panics.

The remaining M-D entries below are retained only as non-fundamental backlog notes from the same read-through. They should not be treated as the answer to "what can make midpoint counts wrong?"

### M-D-004 - Low - Even-sized windows use an uncentered QC plot x-axis

`plot_midpoint_profiles` centers x-values only when `window_size` is odd. For even windows it plots positions as raw `0..window_size` values and heatmap edges as `-0.5..window_size-0.5` ([plotting.rs](../../src/commands/midpoints/plotting.rs#L63-L78)).

Impact: the dense `.npy` artifact is still correct and zero-based relative to the BED window start, but the QC plots silently switch conventions based on window parity. For even fixed-size windows, a user looking for midpoint enrichment around the site center can read the peak position incorrectly.

Recommended fix: choose one plot contract. Either always plot zero-based position to match the array artifact, or always plot center-relative coordinates with an explicit even-window convention. Add a plotting helper test for even and odd window sizes.

# Test coverage analysis (2026-05-11)

Date: 2026-05-11

Scope: every test currently exercising any code path in [src/commands/midpoints/](../../src/commands/midpoints/) — both inline `#[cfg(test)] mod tests` blocks and integration tests under [tests/](../../tests/) — plus the directly used shared helpers (`shared/midpoint.rs`, `shared/blacklist/overlaps.rs`, `shared/scale_genome.rs`, `shared/tiled_run.rs`, `shared/window_fetch.rs`, `shared/length_axis.rs`, `shared/bed.rs`, `shared/read.rs`, `shared/gc_tag.rs`). Goal: judge whether tests are anchored to **what the command should do** (intention) rather than **how it currently does it** (implementation), and write down everything that needs to be added, fixed, or hardened before release.

The recurring failure mode for these tests is hidden in the assertions, not in the fixtures. Many tests construct a believable input, then assert against shape and a single non-zero cell, but never assert that the rest of the array is zero, never assert the total mass, and never assert what the *intent* of the configuration was supposed to filter out. A test that says "the only non-zero count is at position 5, value 1" is much stronger than "the value at position 5 is 1". The fixes below are mostly tightening assertions to capture intent.

A second recurring failure mode is using the same fragment length, same window size, same tile size, and same single-bin length axis in nearly every test, which means many tests do not actually exercise the dimensions they appear to exercise (the length-bin axis, in particular, is essentially a single integration test plus index-of unit tests).

I have grouped the findings by the concrete behavior being tested or missed, and inserted "need to check" wherever the documented or implemented intention is ambiguous to me.

## What is well covered today

These are intention-level behaviors that already have a regression test I trust. They should be left as-is unless the underlying contract changes:

- **Length-bin parsing**: range spec vs explicit edges, max-supported-length cap, brace-style range rejection — covered by [config_tests.rs](../../src/commands/midpoints/config_tests.rs) and `length_bin_range_spec_matches_brace_expansion_edges`, `length_bin_start_end_list_format_is_rejected`.
- **Default MAPQ behavior** vs explicit MAPQ values — `midpoints_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero`.
- **Paired/unpaired parity for the same physical fragment span** — `unpaired_single_read_matches_paired_midpoint_profile_for_same_span`.
- **Tile ownership by midpoint and not by start/end**: a fragment is counted by exactly one tile no matter how the BED row falls — `bed_sites_mixed_core_and_halo_rows_keep_only_the_core_midpoint_count_across_tile_sizes`, `core_overlap_bed_site_is_kept_for_midpoints`, `later_tile_site_keeps_midpoint_count_when_window_span_starts_after_zero`.
- **Even-fragment midpoint tie behavior at the command level** — `even_length_midpoint_tie_counts_exactly_one_of_two_adjacent_edge_windows`.
- **Even-fragment central-base blacklist contract** — `blacklist_midpoint_filtering_checks_both_centers_for_even_fragments`.
- **Group axis encounter order and collapsed counts** — `group_index_axis_matches_first_group_encounter_order_and_collapsed_counts`.
- **GC file pipeline neutrality and non-neutrality** — `real_ref_gc_bias_then_gc_bias_package_is_neutral_in_single_bin_case_for_midpoints`, `real_ref_gc_bias_then_gc_bias_package_changes_midpoints_in_expected_direction`.
- **GC package rejection paths**: length-range mismatch and schema-version mismatch — `midpoints_rejects_gc_package_when_length_bins_are_outside_supported_range`, `midpoints_rejects_gc_package_with_schema_version_mismatch`.
- **GC × scaling multiplication semantics** — `gc_file_and_scaling_tsv_weights_multiply_in_midpoints`.
- **GC-tag pair-average weight semantics** — `gc_tag_pair_average_sets_midpoint_profile_weight`.
- **bam-to-bam → midpoints (--gc-file vs --gc-tag) cross-command parity** — `bam_to_bam_gc_file_output_drives_midpoints_gc_tag_same_as_original_gc_file`.
- **Scaling TSV must cover the chromosome** — `scaling_tsv_must_cover_requested_chromosome_end_in_midpoints`.
- **GC file fetch-narrowed reference origin** — `gc_file_late_tile_site_uses_reference_coordinates_after_fetch_narrowing`.
- **Chromosome-end halo preservation, single fragment per chromosome** — `midpoint_fetch_narrowing_preserves_tile_halo_near_chromosome_end_on_three_chromosomes`.
- **Chromosome-end halo preservation, multiple fragments per chromosome** — `midpoint_fetch_narrowing_reads_all_eligible_fragments_near_chromosome_end_on_three_chromosomes`.
- **`get_overlapping_sites_and_adapt_fetch_to_extremes` unit-level halo behavior at the chrom end** — `midpoint_fetch_span_preserves_tile_carried_halo_near_chromosome_end`, `midpoint_fetch_span_keeps_fragment_start_that_old_symmetric_halo_would_drop`.
- **Sparse partial-file format roundtrip, malformed-file rejection, parallel merge across chunks** — `counting_by_group_tests.rs`.
- **Cross-command real-artifact midpoint mass** — `midpoints_consumes_shared_real_artifacts_with_expected_profile_mass`.
- **CLI smoke** — `midpoints_cli_minimal_invocation_writes_profiles_and_group_index`.
- **Shared midpoint randomization is reproducible, balanced, and seeded by chr/start/length** — `midpoint_tests.rs`.

## Retained coverage bullets after filtering

These are Claude bullets kept as actionable midpoint coverage work. The original wording is preserved below each triage note.

### T-FIX-001 - `midpoint_profiles_written_with_group_index` only asserts shape and "sum > 0"

Codex triage: retain high. This is the best place to pin multi-group, multi-length-bin axis correctness.

[test_profile_groups_command.rs:769](../../tests/test_profile_groups_command.rs#L769) asserts shape `[2, 2, 40]` and `arr.sum() > 0`, plus the presence of the two group names in the TSV. This is a presence test, not a correctness test.

Why this matters: if the merge dropped half the entries, or if length-bin assignment was off by one, or if the second group accidentally received all counts, this test would still pass. The fixture is the `complex_bam_fixture()` so the expected count distribution is derivable by hand for the given length bins `[20, 60, 120]` and window sizes.

Recommended fix: derive the exact expected count per `(group, length_bin, position)` cell by hand for this fixture, then assert the full array equals the derived expectation. At a minimum, assert per-group totals and per-length-bin totals so a regression on either axis fails the test.

### T-FIX-006 - Multi-thread parity is not asserted at all

Codex triage: retain high. This is the highest-value guard for parallel merge and scheduling regressions.

Every `set_n_threads` call uses 1 or 2, and there is no test that runs the same fixture at `n_threads=1` and `n_threads=8` and asserts identical output. The current parallel merge contains a `start_chunk` jitter (`sparse_merge_start_chunk`) that is *intended* to be a perf-only optimization, but without a parity test, no automated check guards against a regression in that path that changes results.

Recommended fix: pick one mid-sized fixture (e.g. the `group_index_axis_matches_first_group_encounter_order_and_collapsed_counts` fixture or `complex_bam_fixture`) and run it with `n_threads ∈ {1, 2, 4, 8}`. Assert all four outputs are bit-identical, both for the `.npy` and the `group_index.tsv`. This is the single highest-value test to add for catching parallel-merge regressions, because the sparse-merge code path is otherwise only exercised by unit tests on synthetic partial files.

### T-FIX-016 - `neutralize_invalid_gc=true` is not tested for midpoints

Codex triage: retain high. `neutralize_invalid_gc = true` is an untested command branch.

For both the GC-file path and the GC-tag path, `--neutralize-invalid-gc` swaps invalid-GC behavior from "drop the fragment" to "count with weight 1.0". This option is *never* `true` in any midpoints test.

Why this matters: the neutralization path explicitly increments `gc_failed_fragments` but does not increment `counted_fragments` unless the fragment is also counted. The two counters interact in a way that is easy to break.

Recommended fix:
- One test: GC-tag mode, fragment with missing GC tag, `neutralize_invalid_gc=true`. Assert the fragment is counted with weight 1.0, `gc_missing_tags == 1`, `gc_failed_fragments == 1`, and the profile mass is 1.0 at the expected cell.
- One test: same fragment, same fixture, `neutralize_invalid_gc=false`. Assert the fragment is dropped, `gc_missing_tags == 1`, profile mass is 0.0.
- Same pair of tests for the GC-file path with a fragment whose corrected weight is `None` (e.g. an all-N reference slice — *need to check whether you want to expose a test path that produces a `None` weight from the file corrector; possibly easier via the GC-tag path*).

### T-NEW-002 - Fragments whose midpoint lies exactly on a tile core boundary

Codex triage: retain high. This is sharper than a broad tile-size sweep because it pins the half-open core-boundary contract.

`process_tile` keeps a fragment when `midpoint >= tile.core_start() && midpoint < tile.core_end()`. There is no test that constructs a fragment whose midpoint is exactly `tile.core_start()` of tile *N+1* (and exactly `tile.core_end()` of tile *N*) to prove the half-open ownership.

Add: tile size 10, midpoint at position 30, BED windows at both [29,30) (in tile N) and [30,31) (in tile N+1). Assert the count lands in the [30,31) window because of half-open ownership. Then change midpoint to 29 and assert the count lands in [29,30) of tile N. This proves the boundary contract.

### T-FIX-008 - Sparse-merge tests do not cover the wrap-around boundary

Codex triage: retain medium. Keep this near sparse merge unit tests, not as a broad command fixture.

[counting_by_group_tests.rs:167](../../src/commands/midpoints/counting_by_group_tests.rs#L167) exercises `add_from_sparse_npz_files_parallel_with_chunk_size` with two files and chunk size 3. The two files share both endpoints (first_idx and last_idx). The wrap-around split in `merge_sparse_profile_partial_file` (start chunk N → end, then 0 → N-1) is therefore not stressed: the test does not construct a `start_chunk` that lands in the middle of an entry distribution.

Why this matters: the wrap-around logic is the most subtle part of the sparse-merge path. The test for `sparse_merge_start_chunk` covers the chunk-picking function in isolation, but there is no end-to-end assertion that "starting at chunk N still visits every entry exactly once".

Recommended fix: extend the parallel-merge test to use enough chunks (e.g. 8 chunks) and a path index that forces a non-zero `start_chunk`. Construct one sparse partial file whose entries are spread across every chunk, and assert that all entries are summed exactly once regardless of `start_chunk`. A property-style version: loop over `start_chunk` values from 0..num_chunks and assert the dense output is identical to a single-threaded reference merge.

### T-FIX-009 - Blacklist tests only cover the `Midpoint` strategy

Codex triage: retain medium. This should stay command-level because blacklist strategy selection is user-facing.

The integration tests only exercise `BlacklistStrategy::Midpoint`. The default is `BlacklistStrategy::Any` ([config.rs:178](../../src/commands/midpoints/config.rs#L178)), but the default strategy is never explicitly asserted against an interval that overlaps a fragment without overlapping the midpoint. `All` and `Proportion=...` strategies have no command-level coverage.

Why this matters: documentation calls out that `midpoint` checks *both* central bases for even fragments. By contrast, what does `Any` actually do for even-length fragments with the midpoint outside the blacklist but a tail base inside? *Need to check whether the user-facing intention is "fragment touches blacklist anywhere drops it" or "midpoint touches blacklist drops it" for the default `Any` strategy.* The current documentation only commits to behavior for the `midpoint` strategy.

Recommended fix:
- Add one test per strategy (`Any`, `All`, `Midpoint` already exists, `Proportion=0.5`) that uses identical inputs and a blacklist tuned to drop or keep the fragment under each strategy. Assert that the same fragment produces 0 mass for the dropping strategy and the full mass for the keeping strategy.
- Add an explicit assertion that `--blacklist-strategy` is `Any` by default in `MidpointsConfig::new()`, and that running with `--blacklist-strategy any` vs. running with no `--blacklist-strategy` flag produces identical output.
- Add a `blacklist_min_size` test: build a blacklist with two intervals, one of size 1 bp and one of size 100 bp, run with `blacklist_min_size=10`, and assert the 1 bp interval is ignored. Today this path is silent.

### T-FIX-011 - No test asserts that fragments outside `--length-bins` are dropped

Codex triage: retain medium. This pins the agreement between fragment filtering and `LengthAxis` lookup.

The default fragment-length cap is `[30, 1001)`. The integration tests set custom narrow length bins (e.g. `[61, 62]`) and *implicitly* rely on the iterator filter to drop everything else. There is no test that says "given a fragment of length 1500 and length bins [30, 1001), the fragment is dropped, the run still succeeds, and the run statistics record one filtered fragment".

Why this matters: the `LengthAxis` `bin_index` lookup uses a `vec![usize::MAX; max_edge]` sentinel and depends on the fragment filter dropping out-of-range lengths *before* the count call. If the filter and axis disagree, the call to `incr_weighted` would silently `bail!` on the worker, which `?`-propagates and aborts the whole run.

Recommended fix: write one test with two fragments — one length 61 (in range), one length 200 (out of range with `--length-bins 30 70`). Assert the in-range fragment produces 1 count at the expected cell, the out-of-range fragment produces 0 counts (`arr.sum() == 1.0`), and the run does not error.

### T-FIX-014 - The `--require-proper-pair` conflict with `--reads-are-fragments` is not tested

Codex triage: retain medium. Cheap user-facing guard test.

[midpoints.rs:91-93](../../src/commands/midpoints/midpoints.rs#L91-L93) bails on the combination. `ends.rs`, `fcoverage.rs`, `fragment_kmers.rs`, and `lengths.rs` all have the same check and *do* have command-level tests asserting the bail. `midpoints.rs` does not.

Recommended fix: add one test that sets both flags and asserts `run(&cfg)` returns an error containing `--require-proper-pair cannot be used with --reads-are-fragments`.

### T-FIX-017 - `ensure_uniform_window_len` failure path is not asserted via midpoints

Codex triage: retain medium. Cheap user-facing guard test for the fixed-width window contract.

`ensure_uniform_window_len` in [windows.rs](../../src/commands/midpoints/windows.rs) bails if any two windows have different lengths. This is core to the midpoint contract (every window must have the same number of output positions). There is no test, integration or unit, asserting that mixed-length windows fail loudly.

Recommended fix: add one test using a BED with two windows of different lengths (e.g. `[10,20)` and `[30,50)`) and assert the run bails with the expected message ("Non-uniform window length detected..."). Optionally add a second test asserting that an empty BED produces "No windows found...".

### T-FIX-019 - `group_idx_to_name` stability across no-count groups is not asserted

Codex triage: retain medium. This matters for downstream group-index stability.

[load_grouped_windows_from_bed](../../src/shared/bed.rs#L458) enumerates groups in BED-encounter order. The midpoints output preserves this enumeration in the group axis. There is no test that asserts: "a group with zero counted fragments still appears as the correct row of zeros in the output, with the correct name in `group_index.tsv`".

Recommended fix: write one test where two groups are declared in the BED, fragments only hit the second group's windows, and the assertion checks `arr[group_idx["empty_group"], 0, ..] == zeros` and `arr[group_idx["hit_group"], 0, expected_position] == expected_value`. This proves group stability for downstream callers that index by name.

### T-NEW-003 - BED rows where windows in different groups overlap each other

Codex triage: retain medium, pending semantics decision. Different-group overlap behavior should be explicit.

A midpoint that falls inside two overlapping BED rows from different groups should add 1 count to each group's collapsed profile. There is no test for this. The current group-axis test uses non-overlapping BED rows.

Add: two BED rows `chr1 45 56 groupA` and `chr1 50 61 groupB` (overlapping), one fragment with midpoint 52. Assert `groupA[..., 52-45=7] == 1.0` and `groupB[..., 52-50=2] == 1.0`. Total mass should be 2.0. *Need to check that this is the intended semantics rather than "midpoint counts in at most one group"; my reading of [find_overlapping_windows](../../src/shared/overlaps.rs) and the loop at [midpoints.rs:742-764](../../src/commands/midpoints/midpoints.rs#L742-L764) is that all overlapping windows get incremented, but a regression test makes that explicit.*

### T-NEW-004 - BED rows where two rows of the same group overlap each other

Codex triage: retain medium, pending semantics decision. Same-group overlap behavior should be explicit or rejected.

The grouped BED loader uses `IndexedInterval::new` and pushes both rows. Both rows will then be returned by `find_overlapping_windows` for a midpoint that is inside both, and both will increment the *same* group axis row. This means a single midpoint can produce *two* counts in the same group's collapsed profile if the user supplied overlapping same-group BED rows.

Why this matters: *need to check whether this is the intended behavior or a documentation gap*. Either:
- (a) Overlapping same-group rows are intentional and add two counts (i.e. the user is allowed to "weight" sites by listing them multiple times).
- (b) Overlapping same-group rows are accidental and should be normalized first.

Add: a test that pins down whichever answer you choose, and document the answer in [config.rs](../../src/commands/midpoints/config.rs).

### T-NEW-007 - Multiple `--blacklist` paths concatenate and merge

Codex triage: retain medium. Multiple blacklist paths are a real CLI surface.

[config.rs:150-153](../../src/commands/midpoints/config.rs#L150-L153) accepts `Option<Vec<PathBuf>>` for multiple blacklist files. There is no test asserting that two blacklist files are loaded, merged, and applied as the union.

Add: two blacklist files, each containing one disjoint interval that overlaps a different fragment. Assert both fragments are dropped, i.e. the union semantics holds.

### T-FIX-018 - Plotting helper is not unit-tested

Codex triage: retain low. Wait until M-D-004 chooses the plot x-axis contract.

`plot_midpoint_profiles` is only exercised when the `plotters` feature is enabled, and there is no test of it in either form (unit or integration). M-D-004 already calls this out for the even-vs-odd x-axis issue. Coverage for the helper as a whole — plot file creation, bad `plot_groups` index bail path, empty `plot_groups` no-op path — is missing.

Recommended fix:
- One unit test for empty `plot_groups`: returns `Ok(())` without writing any plot files.
- One unit test for out-of-bounds `plot_groups`: bails with the expected error message.
- One unit test for `num_length_bins == 1`: only the line plot is produced, no heatmap.
- One unit test asserting that the x-axis values for both odd and even window sizes follow the chosen convention (this will be the M-D-004 regression test after the contract is fixed).

### T-NEW-013 - `MidpointsConfig::new()` defaults

Codex triage: retain low. Useful only if config defaults are worth snapshotting.

There is no test asserting the documented defaults of `MidpointsConfig::new()`: `length_bins == ["30", "1001"]`, `tile_size == 20_000_000`, `min_mapq == 30`, `require_proper_pair == false`, `blacklist_strategy == Any`, `blacklist_min_size == 1`, `plot_groups == [0]`.

Add: one snapshot-style test that asserts every default field. This catches regressions where a default is silently changed in a CLI refactor.

## Folded coverage bullets

These Claude bullets are not wrong, but they should be handled through one of the retained items rather than tracked separately.

### T-FIX-010 - Length-bin axis is barely exercised at the command level

Codex triage: fold into T-FIX-001. Full tensor assertions on `midpoint_profiles_written_with_group_index` cover the important length-axis integration risk.

Almost every command-level test uses `set_length_bins(vec![61, 62])` or a similar single-bin axis. `midpoint_profiles_written_with_group_index` is the only multi-bin test, and (see T-FIX-001) it only asserts `arr.sum() > 0`.

Why this matters: the length-bin axis is one of the three output dimensions. A bug that placed length 50 fragments in bin index 0 instead of bin index 1, but only in `process_tile`, would slip past every other test.

Recommended fix: add at least one focused test where the input is a small set of paired fragments at known lengths (e.g. 30, 40, 70) with `--length-bins 30 50 80 110`, every fragment placed so the midpoint lands at position 5 of a single window. Assert the full `[1, 4, 11]` array equals an explicit expected tensor with exactly one count at each `(0, expected_bin_idx, 5)` cell. Repeat with `--length-bins 30:80:10` and assert the same outcome.

### T-NEW-001 - End-to-end determinism: same input, multiple runs, byte-identical output

Codex triage: fold into T-FIX-006, but do not require raw byte identity unless byte-identical `.npy` output is a public contract.

There is no test that says "run the same midpoints command twice; the `.npy` files are bit-identical". The midpoint random tie-break is documented as reproducible per-coordinate, but the only proof of end-to-end determinism is the helper-level `midpoint_tests.rs`. The merge step combines tile partial files in a parallel order that depends on Rayon scheduling, so this is the right place to confirm determinism survives parallelism.

Add: run the same `MidpointsConfig` twice over the same fixture (paired and unpaired, GC and non-GC), compare the resulting `.npy` files byte-by-byte. *Need to check whether you accept bit-equality as the contract or whether `f32` summation order tolerances may be allowed; if so, narrow this to `assert_close` with a tight tolerance.*

### T-NEW-015 - The `view_ndarray3_group_len_pos` consumer contract

Codex triage: fold into T-FIX-001. Reading the written `.npy` and asserting the full tensor pins axis order.

`ProfileGroupsCounts::view_ndarray3_group_len_pos` exposes a `(group, length_bin, position)` view that is passed straight into both `write_npy` and the plotting code. There is one unit test in [test_profile_groups_counts.rs](../../tests/test_profile_groups_counts.rs) that walks every cell — that's good. There is no integration test that says "when I read the written `.npy` back with `ndarray-npy`, the axis order matches `(group, length_bin, position)` for a known multi-group, multi-length-bin, multi-position run".

The whole T-FIX-001 fix family addresses this implicitly but it deserves a dedicated test that focuses on axis order.

Add: one fixture with three groups, three length bins, five positions, and *distinct* nonzero counts at every cell. Assert that `read_npy::<Array3<f32>>(...)` returns the array shaped `[3, 3, 5]` and that every cell equals the expected value. This is the load-bearing test for downstream NumPy code.

### B-006 - `complex_bam_fixture` is reused without intent-anchored assertions

Codex triage: fold into T-FIX-001. This is the same fixture/assertion weakness.

The fixture is used in `midpoint_profiles_written_with_group_index`. That fixture has a complex, hand-coded set of fragments. Reusing it without a derived expected count tensor (see T-FIX-001) means each test relying on it can pass for the wrong reason if the fixture is silently mutated later. *Need to check whether you want to add a closed-form "the expected count per cell for `complex_bam_fixture` with `--length-bins 20 60 120` and BED `[…]` is ` …`" fixture-side helper, so that every test re-using the fixture asserts the same intent.*

## Filtered priority

If there is time for a small coverage pass, prioritize:

1. **T-FIX-006** plus folded **T-NEW-001**: multi-thread and rerun parity, using numeric equality unless byte-identical output is explicitly required.
2. **T-FIX-001** plus folded **T-FIX-010**, **T-NEW-015**, and **B-006**: exact full-tensor assertion for the broad multi-group, multi-length-bin fixture.
3. **T-FIX-016**: `neutralize_invalid_gc=true`.
4. **T-NEW-002**: half-open tile core-boundary ownership.
5. **T-FIX-008**: sparse merge wrap-around.
6. **T-FIX-014** and **T-FIX-017**: cheap user-facing guard tests.
7. **T-FIX-009** plus **T-NEW-007**: blacklist strategies and multiple blacklist paths.
8. **T-FIX-019**, **T-NEW-003**, and **T-NEW-004**: zero-count group stability and overlapping BED-row semantics.
