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

### M-005 - Implemented - `--blacklist-strategy midpoint` now checks both central bases for even fragments

Blacklist filtering happens before the command samples the counted midpoint ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L495-L510)). The old shared blacklist midpoint strategy used `start + (end - start) / 2`, which is the right-center base for even-length half-open fragments, while count placement randomly selected left or right center for the same even-length fragment.

Status: implemented by making the shared `midpoint` blacklist strategy check central-base support. Odd-length fragments check the unique central base. Even-length fragments check both central bases, so either central base overlapping a blacklist interval excludes the fragment. The command-level regression now proves both left-center and right-center blacklist hits remove the even-length fragment.

Impact: this keeps `midpoint` as a conservative coordinate-mask strategy without coupling shared blacklist filtering to the command-specific randomized midpoint tie break.

No remaining action for this finding.

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

### M-D-007 - Low - Dense output is written before the group-index sidecar

`run` writes `<prefix>.midpoint_profiles.npy` before it writes `<prefix>.group_index.tsv` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L295-L303)). If the sidecar write fails, the user can be left with a fresh primary count artifact that has no matching group-name map.

Impact: the count tensor is not self-describing. A missing or stale sidecar can make downstream group interpretation ambiguous, especially when group order depends on the first occurrence order in the BED.

Recommended fix: write both public artifacts to temporary paths and rename them into place only after both succeed. At minimum, write the sidecar first so a sidecar failure does not leave a new unlabeled count tensor.
