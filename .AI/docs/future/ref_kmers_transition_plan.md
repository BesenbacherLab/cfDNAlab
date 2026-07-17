# Ref-kmers Transition Plan

## Purpose

Move the standalone `../reference` reference k-mer counter into cfDNAlab as a native command:

```text
cfdna ref-kmers
```

The command should not be a blind port. It should use the normal cfDNAlab command API, row/window semantics, logging, output handling, and feature gates. It should also be designed as a reusable reference background for downstream correction or subtraction in:

- `cfdna ends`
- future `cfdna fragment-kmers`

The current `fragment-kmers` implementation is not the target contract. Future fragment-kmers work should follow `src/commands/fragment_kmers/positional_selection_logic.md`: k-independent eligibility, intersection, runs, final step, and only then k-dependent k-mer fit.

## Existing Standalone Reference Command

The current `../reference` code provides a useful core:

- reads a 2bit reference by chromosome
- supports `--global`, `--by-size`, and `--by-bed`
- masks blacklist intervals into the reference sequence
- builds left-aligned k-mer code arrays per requested k
- counts k-mer starts that fit completely inside each window
- writes dense `.npy` or sparse SciPy `.npz` count matrices plus motif lists
- writes a `bins.bed` sidecar for non-global rows

The core algorithm is small and worth reusing conceptually. The integration work is larger than the counting work.

Important differences from cfDNAlab expectations:

- It is a separate CLI with its own clap surface, output naming, progress, and thread-pool handling.
- Its BED and blacklist parsers are more permissive than cfDNAlab loaders.
- It silently clips some intervals and can produce bad metadata for windows that start beyond the chromosome.
- It reads and encodes whole chromosomes, which is simple but can become RAM-heavy for large k, many k values, or many worker threads.
- It writes raw array sidecars, not a self-contained cfDNAlab package.
- It has no grouped BED mode.
- Its canonical odd-k behavior differs from the current shared fragment-kmers collapse behavior.

## Target Command Shape

Add a new feature and module:

```text
cmd_ref_kmers
src/commands/ref_kmers/
```

Expected public CLI and Rust API:

```text
cfdna ref-kmers
```

```rust
RefKmersConfig
run_ref_kmers(&RefKmersConfig, RunOptions) -> Result<RefKmersRunResult>
```

The command should follow the established pattern:

- config struct with clap docstrings
- `new` constructor and explicit setters
- `ToCliCommand`
- `RunOptions`
- `CommandRunResult`
- `run_like_cli::ref_kmers` export
- CLI enum wiring in `src/cli_app.rs`
- feature wiring in `Cargo.toml` and `build.rs`

The command should be named `ref-kmers`, not `reference`, because `cfdnalab::reference` already means reference-genome utilities at the crate root.

The current `src/commands/ref_kmers/ref_kmers.rs` should be cleaned up incrementally. It already has pieces of the desired cfDNAlab command shape, but also contains copied `ends`, `ref-gc-bias`, and standalone `reference` fragments. The transition should fix those pieces in place instead of deleting the file and starting over.

## Config Surface

Initial config should be close to:

- `Ref2BitRequiredArgs`
- `output_dir`
- `output_prefix`
- `n_threads`
- `kmer_size`
- `DistributionWindowsArgs`
- `AssignToWindowArgs`
- `ChromosomeArgs`
- `blacklist`
- `blacklist_min_size`
- `canonical`
- `all_motifs`
- `motifs_file`
- `tile_size`
- `LoggingArgs`

`kmer_size` is singular on purpose. Supporting multiple k values while also supporting motif-file grouping makes the motif contract and output shape much less clear.

`all_motifs` controls whether the complete motif axis is requested. Dense output should be used when the output axis is known to contain the complete motif set for the configured k-mer contract. Otherwise sparse output is preferred. This is a mathematical/output-contract decision, not just a question of whether the command generated all motifs itself.

`motifs_file` is part of the target command. It should support selected motifs and motif groups. The docs should describe the accepted k-mer motif format directly, not inherit `ends` wording about inside and outside motif halves.

`canonical` is a boolean. It means final reverse-complement collapse of k-mer labels and counts. It should not become a larger internal mode enum unless a concrete second collapse rule is added.

`tile_size` is part of the config. The command should be tiled in the first full cfDNAlab implementation because the output is meant to be useful for large reference backgrounds.

## Window And Row Semantics

Use cfDNAlab windowing semantics, not the standalone reference parser:

- Global mode writes one output row.
- Fixed-size mode writes generated windows per selected chromosome, with final windows clipped to chromosome length.
- BED mode uses original BED row indices as row identity.
- Grouped BED mode uses `group_idx` as row identity.
- Selected BED and grouped BED input should fail early when no intervals remain after chromosome filtering.

Grouped BED should initially aggregate intervals by `group_idx` with interval-mass semantics:

- Count each retained grouped BED interval as its own reference interval.
- Add its k-mer counts into the group row.
- Overlapping same-group intervals contribute separately.
- Write group metadata with `group_idx`, `group_name`, `eligible_windows`, and interval-width weighted `blacklisted_fraction`.

Do not implement unique-base grouped semantics unless a separate flag is explicitly designed. `fcoverage` has unique-base grouped modes because coverage has a separate support interpretation. That should not be inherited silently here.

BED coordinates need an explicit chromosome-boundary policy before implementation. Recommended behavior:

- reject selected intervals with `start >= chromosome_length`
- clip `end` to chromosome length only if the metadata records the clipped interval consistently
- never divide by an empty clipped span
- make the error message include chromosome and interval coordinates

## K-mer Window Assignment

Do not add a ref-kmers-specific overlap enum when the package already has `WindowAssigner`.

The counted object for `ref-kmers` is the reference k-mer span `[p, p + k)`. Apply the existing assignment semantics to that span:

- `count-overlap`: count `overlap_bp / kmer_size` for each overlapped row.
- `any`: count `1.0` when any k-mer base overlaps the row.
- `all`: count `1.0` only when the full k-mer span is contained in the row.
- `proportion=<threshold>`: count `1.0` when at least that fraction of k-mer bases overlaps the row.
- `midpoint`: count `1.0` for the row containing the k-mer midpoint, and only accept this mode when `kmer_size` is odd.

Reject `midpoint` for even `kmer_size` with a direct error message. For even k-mers, there is no single center base and silently choosing one side would create an orientation bias.

For `ends`-style reference backgrounds, `count-overlap` should be documented as the mode that preserves equivalence between counting left-to-right and reverse-complementing the final output. For ordinary fully-contained k-mer opportunity tables, `all` is the stricter mode.

## Counting Invariants

`ref-kmers` counts reference k-mer background mass, then stores row-wise frequencies.

For ordinary k-mer mode:

- A k-mer start at reference coordinate `p` owns the reference span `[p, p + k)`.
- Tile ownership is by start coordinate. A k-mer is processed by the tile whose core contains `p`.
- Window assignment is independent of tile ownership and follows `AssignToWindowArgs`.
- K-mers containing `N`, masked blacklist bases, or sentinel-invalid positions do not count.
- Count mass can be fractional when `--assign-by count-overlap` is used.
- Count mass should be accumulated before converting to frequencies.
- Output frequencies are `count_mass[row, motif] / scaling_factor[row]`.
- The scaling factor is the row total count mass needed to reconstruct counts downstream.
- Blacklist masking removes every k-mer that overlaps a blacklisted base.
- K-dependent fit and overlap belong at the final counting step.

These invariants match the useful part of `../reference`, but should be implemented with cfDNAlab interval helpers and explicit validation.

## Background Compatibility

The command has to be useful for background correction in more than one downstream command. That requires recording the geometry used to define the motif axis, not just writing motif strings.

### Fragment-kmers

For future `fragment-kmers`, the first useful reference background is row-wise reference k-mer availability:

```text
row x motif
```

This can be joined to observed fragment-kmers by row and motif, provided both outputs use the same:

- k-mer size
- canonical flag
- motif ordering
- window assignment mode
- window/group row identity
- blacklist and chromosome selection assumptions

Do not overfit this to current fragment-kmers positional counting. If a later positional fragment-kmers output needs a position-specific reference background, it should reuse the future k-independent selection/run model from `positional_selection_logic.md`, not the current implementation.

### Ends

`ends` is not automatically compatible with ordinary left-aligned k-mer starts. End motifs have endpoint geometry:

- `k_outside`
- `k_inside`
- left and right endpoint orientation
- final inward-facing motif labels
- optional motif-file grouping
- optional complement/collapse behavior

If `ref-kmers` is expected to produce an accurate `ends` background, it needs either:

- an end-context mode that counts all possible endpoint contexts with the same motif contract as `ends`, or
- a documented statement that ordinary `ref-kmers` output is only a generic sequence-composition background and not an exact `ends` opportunity table.

The cleaner design is to make the counting layer generic over reference motif geometry:

```text
ReferenceMotifGeometry::KmerStart { k }
ReferenceMotifGeometry::EndContext { k_outside, k_inside, endpoint_side_policy }
```

The command can still be named `ref-kmers`; the important part is that the output metadata says which geometry produced the motif axis.

## Canonicalization

`canonical` is a final output collapse step.

The command should first count oriented reference k-mer labels according to the selected window assignment. After tile reduction, if `canonical` is enabled, reverse-complement-equivalent labels are collapsed and their count masses are summed. Frequencies and scaling factors are then computed from the collapsed count matrix.

The settings output must record whether canonical collapse was used. Downstream correction must reject mismatched canonical settings instead of silently joining incompatible motif axes.

## Shared K-mer Refactor

Most reusable code already lives under `src/shared/kmers`, but it is gated around `cmd_ends` and `cmd_fragment_kmers`.

Needed cleanup:

- Add `cmd_ref_kmers` to the relevant gates or introduce an internal `uses_kmers` cfg alias in `build.rs`.
- Avoid tying the core reference k-mer key to fragment orientation. A reference row count only needs `{ k, code }`; fragment-kmers may need orientation separately.
- Expose sentinel accessors and left-aligned code builders to `cmd_ref_kmers`.
- Keep reference count mass separate from weighted fragment counts. Fractional window assignment means ref-kmers needs `f64` count mass even though `any`, `all`, and `midpoint` produce integer-like increments.
- Factor motif preparation so canonicalization, motif grouping, and motif ordering can be shared across ordinary and selected-motif outputs.

The goal is a small shared k-mer core, not a large abstraction layer around every command.

## Output Contract

Because this command is intended to feed downstream correction, raw matrices plus loose sidecars are risky. The target public output should be self-describing.

Preferred target:

```text
<prefix>.ref_kmers.zarr/
<prefix>.ref_kmer_settings.json
```

Recommended Zarr structure:

```text
frequencies[row, motif]            float64 dense or sparse COO
scaling_factors[row]               float64
motif_index[motif]                 integer coordinate
motif_ascii[motif, motif_byte]     motif labels
row metadata arrays
```

Root metadata should include:

- `cfdnalab_schema = "reference_kmers"`
- schema version
- row mode: `global`, `by-size`, `by-bed`, or `by-grouped-bed`
- motif geometry: `kmer-start` or `end-context`
- value units: `reference_kmer_frequency`
- scaling factor units: `reference_kmer_count_mass`
- `kmer_size`
- window assignment mode
- canonical flag
- blacklist handling
- dense or sparse storage mode

Row metadata should include:

- global: one row label
- fixed-size and BED: chromosome axis, row chromosome, start, end, and optional `blacklisted_fraction`
- grouped BED: group axis, group labels, `eligible_windows`, and optional `blacklisted_fraction`

The loader API should expose frequencies by default and also support reconstructing count mass:

```text
count_mass[row, motif] = frequency[row, motif] * scaling_factors[row]
```

Rows with zero support have scaling factor `0`. Loaders must handle those rows explicitly and must not divide by zero silently.

Settings JSON should record human-readable command settings and enough provenance to decide whether a sample output can be corrected with the package:

- k-mer size
- motif geometry and geometry-specific parameters
- canonical flag
- window assignment mode
- selected chromosomes
- reference contig footprint
- window mode
- blacklist paths or at least whether blacklisting was used
- dense/sparse storage mode
- whether the motif axis is complete
- motif-file path or motif-file summary when selected motifs or groups are used

If raw `.npy`/`.npz` outputs are kept for convenience, they should be treated as secondary exports. They should not be the only public package contract for a reusable background.

## Memory And Tiling

The standalone implementation reads a whole chromosome and builds code arrays for every requested k. This is simple and can be a useful baseline, but it is not the ideal public implementation for high k or many k values.

The scalable design is tiled reference counting:

- tile cores own candidate k-mer starts
- sequence reads include left and right halos of `kmer_size`
- a k-mer is counted only when its start coordinate is inside the tile core
- the k-mer span is assigned to rows using `AssignToWindowArgs`
- blacklisted bases in the halo still invalidate the motif

This guarantees each candidate reference opportunity is counted once while avoiding full-chromosome code arrays per worker.

Tiling should preserve the same row identities as the non-tiled algorithm. Temporary row indices, tile-local ordering, and chromosome-local offsets must not leak into the public output.

## Implementation Stages

1. Define the public command contract.

   Confirm the frequency output contract, reconstruction scaling factors, window assignment semantics, motif-file grouping behavior, and chromosome-boundary validation.

2. Add the native command skeleton.

   Add `RefKmersConfig`, `run_ref_kmers`, `RefKmersRunResult`, CLI wiring, `ToCliCommand`, feature gates, and `run_like_cli` exports.

3. Move reusable reference counting logic into cfDNAlab style.

   Use shared reference readers, blacklist loading/masking, interval helpers, window loaders, progress, and thread-pool handling. Accumulate count mass before converting to frequencies.

4. Add row metadata and grouped BED aggregation.

   Build row summaries in final count-row order. Reject non-contiguous grouped row indices when the output axis requires contiguous rows.

5. Add output writing.

   Prefer a Zarr package plus settings JSON. Store frequencies and scaling factors. If raw arrays are written too, make them secondary and keep their row/motif axes tied to the same metadata.

6. Add downstream compatibility checks.

   Provide helper metadata or loader code so `ends` and future `fragment-kmers` can reject incompatible reference packages instead of applying the wrong background. Loaders should optionally reconstruct count mass from frequencies and scaling factors.

7. Add focused tests.

   Cover counting semantics, blacklist masking, BED row order, grouped BED aggregation, boundary validation, canonical mode differences, dense/sparse output shape, settings metadata, and CLI roundtrips.

## Tests To Add

Core counting tests:

- `count-overlap` gives fractional mass at window edges
- `any` counts a k-mer once per touched row
- `all` requires the full k-mer span inside the row
- `proportion=<threshold>` uses the fraction of k-mer bases inside the row
- `midpoint` works for odd `kmer_size`
- `midpoint` rejects even `kmer_size`
- k-mer starts are owned by tile core starts, not by window overlap
- final partial fixed-size windows are clipped correctly
- masked blacklist bases remove every crossing k-mer
- `N` bases remove every crossing k-mer

Window tests:

- BED row order follows original BED row indices
- grouped BED rows follow contiguous `group_idx`
- grouped intervals with the same group aggregate into one row
- overlapping same-group intervals contribute separately
- empty selected BED/grouped BED input fails

Output tests:

- motif labels match count columns
- frequencies and scaling factors reconstruct count mass
- sparse and dense outputs have the same logical matrix shape
- row metadata matches count row order
- settings JSON records canonical flag, motif geometry, k-mer size, window assignment mode, and motif axis completeness

Compatibility tests:

- future fragment-kmers correction rejects mismatched k-mer size, canonical flag, or window assignment mode
- ends correction rejects ordinary `kmer-start` packages when exact endpoint geometry is required
- end-context packages, if implemented, match `ends` motif labels for small synthetic references

## Not In Scope For The Transition

- Fixing the current `fragment-kmers` positional implementation.
- Implementing sample-level background subtraction in `ends` or `fragment-kmers`.
- Supporting backwards compatibility with the standalone `../reference` CLI.
- Adding unique-base grouped semantics without a separate design.
- Treating raw `.npy`/`.npz` sidecars as the long-term background package contract.

## Open Decisions

- Should the first version write only Zarr, or Zarr plus raw dense/sparse array exports?
- Is exact `ends` endpoint-context background part of the first `ref-kmers` command, or a follow-up mode?
- What exact policy should selected BED intervals beyond chromosome bounds use?
- Should dense output be allowed for a motif-file group axis only when every possible motif maps to exactly one output column, or should motif-file outputs always stay sparse unless `all_motifs` is set?
