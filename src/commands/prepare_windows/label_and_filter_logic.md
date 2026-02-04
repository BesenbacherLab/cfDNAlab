# prep-windows — labeling, compositions, and filtering

This document explains how `prep-windows` builds label columns, defines named compositions, and applies filtering rules. It also defines how merging, coordinate choices, cluster tagging, and label exclusions work together. It is written to be clear and exact so you can set things up once and rely on it.

---

## What you can control

1. **Which label columns are written**
   You can write any combination of the atomic parts and your own named compositions

2. **How composed labels are built**
   You can define any number of named compositions, each as an ordered list of parts

3. **How merge grouping is defined**
   You can choose which label key defines "within-group" for merging

4. **Which coordinates define merging**
   You can merge on original coordinates or resized coordinates

5. **Which coordinates are used for distances**
   You can choose whether binning uses original or resized coordinates (resized by default)

6. **How close windows are allowed to be**
   You can drop windows that are too close within a group

7. **Which windows are tagged as clusters**
   You can mark windows as clusters based on an overlap threshold

8. **Which minimum-per rules are enforced**
   You can set several min-per constraints at once, targeting atomic parts or any named composition

9. **Which terms are excluded**
   You can drop windows whose labels include specific terms

All arguments that accept multiple items are **space separated**.
The parts inside `--compose` are comma separated.

---

## Atomic parts

The base pieces you can reference anywhere

* `input` — original group from `--group-cols`.
  * When windows carry multiple groups, values are joined with `__` in stable order.
  * Missing group values are written as `[NA]` so each label keeps the same number of segments.
* `win-direction` — one of `- + =` describing where the window sits relative to the near interval **and its strand**:
  * `-` window upstream of the near interval
  * `+` window downstream of the near interval
  * `=` window overlaps the near interval
  Distances are oriented the same way - strand flips upstream/downstream for `-` strand intervals.

  * When a chromosome has no near intervals, this is `[NONE]`.
* `near-name` — the group name from the near file.
  * Missing group values are written as `[NA]` so each label keeps the same number of segments.
  * When a chromosome has no near intervals, this is `[NONE]`.
* `bin` — the distance bin label.
  * When a chromosome has no near intervals, this is `[NO-NEAR]` if `--distance-bins` is set.
* `cluster` — set to `cluster` for windows that meet the overlap threshold, otherwise `none`.

If a value is missing, it is empty.

---

## Label tuples and compact output

Each window carries one or more label tuples. A tuple is the set of atomic parts that belong together for one original window or near tie. When windows are merged, their tuples are combined into a list that is used for composition and filtering.

This matters for compositions and filtering:

* Compositions are built per tuple, not by cross product across parts.
* Multi-label windows are written in a compact form when only the `input` value differs.
* Lists are used only when other parts differ and we need to preserve the pairings.
* With merging completed before near lookup, bin values are single-valued. The only expected multi-bin case is annotated near ties with signed bins.
* Duplicate tuples created during merges or tie expansion are removed without changing their order.
* For atomic parts other than `input`, **all-equal values collapse to a single value** (e.g., two `win-direction` values `+,+` write as `+`). Pairings to `near-name` are preserved only via compositions such as `near=win-direction,near-name`. Use a composition when you need side–name association in the output or downstream filters.

Compact example:

* Tuple list: `[input=A, bin=prox]` and `[input=B, bin=prox]`
* `input` column: `A__B`
* `bin` column: `prox`
* `core = input bin` column: `A__B.prox`

List example (near ties with signed bins):

* Tuple list: `[input=A, win-direction=-, bin=up]` and `[input=A, win-direction=+, bin=down]`
* `input` column: `A`
* `bin` column: `up,down`
* `core = input bin` column: `A.up,A.down`

---

## Named compositions

A named composition is a label you define by listing parts in order. The tool joins those parts with a dot.

```bash
--compose <NAME=PART1,PART2,...>
```

* `<PARTS...>` can include atomic parts and names of **earlier** compositions.
* Compositions form a directed acyclic graph.
  * Cycles and references to unknown names are rejected with a clear error
* The join character is a dot between parts. We **do not drop empty parts**.
  * If a part is missing, we leave an empty segment so splitting later is consistent.

### Compositions when a window has multiple tuples

When a window has more than one label tuple (after a merge or a near tie), a composition is built for each tuple. The output is compacted when possible.

**Example setup**

```bash
--compose core=input,bin
```

**Window values**

* tuples: `[input=A, bin=prox]` and `[input=B, bin=prox]`

**Composed values**

* `core` -> `A.prox,B.prox`

**What gets written**

* If only `input` varies, we compact the label.
  * `core` column -> `A__B.prox`
* If other parts differ, we write a list to preserve pairings.
  * Example: near ties with signed bins -> `core` column -> `A.up,A.down`

**What gets counted**

* If you set `--min-per core=N`, counting uses the **expanded values**
  * The example window contributes `+1` to bucket `A.prox` and `+1` to bucket `B.prox`
* The window **passes** the `core=N` rule if **any** of its buckets reaches at least `N` in the final counts

**Recursive interaction with `--min-per`**

Min-per rules describe the desired sizes **in the final output**. We therefore evaluate them **recursively** until the result is stable.

Algorithm sketch:

1. Start with all windows and their current label tuples. This includes atomic parts and any compositions.

2. For each `--min-per key=K`:

   a. Compute counts for `key` over the **current** windows.

   b. Determine which label buckets for `key` meet the threshold `K`.

3. For each window, test each `key`:

   * If `key` is single-valued for the window, the window survives that key only if its bucket is in the kept set.
   * If `key` has multiple values for the window, the window survives that key if **any** of its values is in the kept set.

4. Remove anything that no longer qualifies:

   * If a label value for a window fails the rule for a key, **remove that value** from the window's assignments for that key.
   * If a composition uses that key, its values shrink accordingly.
   * If a window loses **all** values for any required `key`, drop the window entirely.

5. If anything changed in step 4, go back to step 2.
   Repeat until no counts, kept sets, or window assignments change.

6. Write the output from this final stable state:

   * Multi-label compositions use the compact form when possible (for example `A__B.prox`).
   * Otherwise they are serialized as lists (for example `A.up,A.down`).

Example outcome:

* With `--compose core=input,bin` and rules `--min-per input=1000 core=250`,
  a window with `core = A__B.prox` contributes to `A.prox` and `B.prox` while both labels are still in play.
  If `A` later fails the `input` rule during recursion, the tuple with `input=A` is removed, `core` for that window becomes just `B.prox`, counts are recomputed, and recursion continues until nothing else changes.

---

## Choosing which columns to write

```bash
--out-labels <TOKENS...>
```

Each token is either an atomic part or a composition name you defined with `--compose`. Columns are written after the coordinates in the order you specify.

* If you list a composition that was never defined, the command fails early.
* `input` uses `A__B__C` when only the input value varies across tuples.
  When other parts differ, it is written as a list to preserve pairings.
* Other atomic parts are written as a single value when all tuples agree, or as a `,`-separated list when they differ.
* Compositions use the compact form when only `input` varies, otherwise a `,`-separated list.

Examples

```bash
--out-labels input bin
--out-labels input win-direction near-name
--out-labels input core report
  --compose core=input,bin
  --compose report=core,win-direction
```

---

## Merge coordinates

Use `--merge-on` to choose which coordinates define merging.

* **Resized**:
  * Merge uses resized coordinates.
* **Original**:
  * Merge uses original coordinates.
  * If resizing is configured, the merged window is resized after merging.

When no resize or flank is configured, resized coordinates match the originals.

Blacklist checks are applied before merging using the **resized** coordinates. When merging uses original coordinates, we still compute a resized span for blacklist filtering and then recompute resized coordinates from the merged original span.
Near lookup and binning happen after merge and resize using the selected distance choice.

---

## Distance from for binning

Distance binning uses a single coordinate choice. Resized coordinates are the default.

Use `--distance-from` to select the coordinate set.

* **Resized**:
  * Distances and bins are computed on the resized coordinates.
* **Original**:
  * Distances and bins are computed on the original coordinates.

When no resize or flank is configured, resized coordinates match the originals.

Distances are computed after merging using the merged window coordinates.
Only the selected coordinates are carried forward for binning.
When a chromosome has no near intervals, `--distance-max` drops its windows.

---

## Merge grouping key

Within-group merging is controlled by a label key. The key can be:

* an atomic part, such as `input` or `bin`
* a named composition

Windows only merge within groups that share the same key value. Across-group merging ignores the key and merges purely by coordinates, but the merged window keeps the full list of label tuples.
The merge key only affects merging, not minimum-distance filtering or output ordering.

---

## Minimum distance within group

You can drop windows that are too close within an input group.

* The comparison uses the same coordinate choice as `--distance-from`.
* Grouping always uses the `input` label, independent of `--merge-key`.
* Spacing is applied after merging so the distance is measured on the merged window.

---

## Cluster labeling

Windows can be marked as clusters based on average position-wise overlap across groups on the same chromosome after merging.

* The overlap threshold is configurable.
* The average overlap is computed using `--cluster-on` coordinates.
* A window is tagged with `cluster=cluster` when its average overlap meets or exceeds the threshold.
  The average is the total overlap depth across the window divided by its length.
* Otherwise `cluster=none`.
* `--cluster-before-min-distance` controls whether clustering happens before or after `--min-distance-within-group`.

This label can be used in `--out-labels`, `--compose`, `--min-per`, or exclusion rules.

---

## Label exclusion rules

You can drop windows whose labels include specific terms before any min-per filtering. These rules apply to atomic parts or compositions.

* For single-valued labels, the window is dropped when the value matches a term.
* For multi-valued labels, the window is dropped when **any** value matches a term.

Examples:

* Drop all proximal bins: `bin=prox`
* Drop cluster windows: `cluster=cluster`
* Keep only cluster windows: `cluster=none`

## Minimum-per filtering

You can set several min-per rules. Each rule is a key and a threshold.

```bash
--min-per input=1000 core=200 near-name=30
```

* Keys can be atomic parts or composition names you defined.
* You can pass several pairs in one `--min-per` flag or repeat the flag.
* All rules are enforced together at the end. A window is only kept if it satisfies **every** rule.
* If a key is unknown or references an undefined composition the command fails early.

**Counting rules by key**

* `input`:
  * **Any-member** counting: A window with `input = A__B` adds one to `A` and one to `B`.
  * The window satisfies `input=N` if **any** of its input groups ends up with a count at least `N`.

* Atomic parts (`win-direction`, `near-name`, `bin`):
  * Each window adds one to every distinct value it carries for that part.
  * This can be more than one bucket when near ties create multiple tuples.
  * The window satisfies the rule if any of its buckets reaches the threshold.

* Named compositions:
  * If the composition produces multiple labels for a window, it uses the same any-member logic.
  * If the composition produces a single label, it behaves like a single-valued atomic key.

**Pruning label assignments during recursive filtering**

When any `--min-per` rule excludes label values for a key, those values are removed from every affected window's assignment for that key during the recursive evaluation. This applies to **all** keys. After each removal, any composition that uses that key is rebuilt from the updated values before the next counting pass. If a window loses all values for any required key, the window is dropped. The process repeats until no counts, kept sets, or window assignments change. The final output is written from this stable state, using the compact form when possible and falling back to lists when needed (for example `A.up,A.down`).

---

## Evaluation order

1. Parse input windows and build initial label tuples.
2. Compute resized coordinates as configured.
3. Apply blacklist checks using the resized coordinates.
4. Deduplicate identical windows by `(chrom,start,end,input)` before merging, using resized coordinates when resizing or flanking is enabled.
5. Merge windows within groups when `--merge-scope within` is selected.
6. If `--merge-on` is `original` and resizing is enabled, recompute resized coordinates from the merged original span.
7. If `--cluster-before-min-distance` is set, compute cluster labels across groups using `--cluster-on` coordinates.
8. Apply minimum-distance filtering within input groups.
9. If `--cluster-before-min-distance` is not set, compute cluster labels across groups using `--cluster-on` coordinates.
10. Merge windows across groups when `--merge-scope across` is selected.
11. If `--merge-on` is `original` and resizing is enabled, recompute resized coordinates from the merged original span.
12. Apply blacklist checks again on the output coordinate set after merging.
13. Compute near distances and bins on the selected distance choice.
14. Build named compositions from the current atomic parts for each tuple.
15. Apply label exclusion rules and drop any windows that match.
16. **Recursive filtering loop**

   Repeat the following until nothing changes:

* Compute counts for every key referenced by `--min-per` on the **current** window assignments.
* For each key, determine which label values meet its threshold.
* For every window and for every key:
  * If the window's assignment for that key is single-valued and that value does not meet the threshold, drop the window.
  * If the window's assignment for that key has multiple values, remove the values that do not meet the threshold. If no values remain for that key, drop the window.
* Rebuild any compositions affected by these removals (for example, compositions that include `input`).

17. Write the requested columns from the final stable state in the order given by `--out-labels`.

---

## Separators and formatting

* Parts inside a composition are joined with a dot.
  Example: `input.bin.win-direction`.
* Multiple `input` values are joined with `__`.
  Example: `A__B.prox`.
* Multiple distinct labels inside a single column are joined with `,`.
  Example: `A.up,A.down`.
* Missing atomic values are written as empty fields.
  * Missing parts in a composition appear as empty segments between dots so every composition keeps the same number of fields after splitting.
  * Example: `A..prox` means the middle part is empty.

---

## Validation and clear failures

* Undefined composition in `--out-labels` -> error that names the missing composition.
* Undefined or cyclic composition in `--compose` -> error that points to the first cycle or unknown.
* Unknown key in `--min-per` -> error that lists the allowed keys and known composition names.
* `--min-per` that targets a composition that expands only to empty labels for all windows -> results in zero kept windows by design.
* User-defined label values are restricted to ASCII alphanumerics. `none` is reserved and cannot be used for input groups, near names, bin labels, or composition names.

---

## Practical recipes

Enforce minimum counts for original groups and for a stratified composition:

```bash
--compose core=input,bin
--out-labels input core
--min-per input=2000 core=250
```

Write atomic parts only and filter by *near names*:

```bash
--out-labels input near-name bin
--min-per near-name=300
```

Two compositions and filtering on both:

```bash
--compose core=input,bin
--compose report=core,win-direction
--out-labels input report
--min-per core=150 report=80
```

Filter on input but output composition:

```bash
--compose core=input,bin
--out-labels core
--min-per input=1000
```

---

## Notes on merging

* **Within-group merges** use the `--merge-key` label value (default `input`) to decide which windows belong together.

* **Across-group merges** ignore `--merge-key` and merge by coordinates only.

* `--merge-label join` keeps all label tuples, so `input` can compact to `A__B` when only `input` differs.

* `--merge-label first` keeps only the first window's labels.

* Across-group merging runs after clustering and min-distance filtering, so cluster labels reflect pre-merge overlaps.

* When both merging and deduplication are enabled, deduplication happens first.

* Deduplication uses the input label and resized coordinates when resizing or flanking is enabled, otherwise original.

---

## Defaults

* No compositions are defined by default.
* `--out-labels` defaults to `input`.
* No min-per rules are applied by default.

---
