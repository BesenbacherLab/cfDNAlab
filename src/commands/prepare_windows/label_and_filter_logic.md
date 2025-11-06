# prep-windows — labeling, compositions, and filtering

This document explains how `prep-windows` builds label columns, defines named compositions, and applies independent minimum-per-group filters. It is written to be clear and exact so you can set things up once and rely on it.

---

## What you can control

1. **Which label columns are written**
   You can write any combination of the atomic parts and your own named compositions

2. **How composed labels are built**
   You can define any number of named compositions, each as an ordered list of parts

3. **Which minimum-per rules are enforced**
   You can set several min-per constraints at once, targeting atomic parts or any named composition

All arguments that accept multiple items are **space separated**.

---

## Atomic parts

The base pieces you can reference anywhere

* `input` — original group from `--group-cols`.
  - When windows are merged across groups this becomes a set of groups for that window.
* `near-side` — one of `- + =` relative to the near interval strand.
* `near-name` — the group name from the near file.
* `bin` — the distance bin label.

If a value is missing, the value is `.`

---

## Named compositions

A named composition is a label you define by listing parts in order. The tool joins those parts with a dot.

```bash
--compose <NAME> <PARTS...>
```

* `<PARTS...>` can include atomic parts and names of **earlier** compositions.
* Compositions form a directed acyclic graph.
  - Cycles and references to unknown names are rejected with a clear error
* The join character is a dot between parts. We **do not drop empty parts**.
  - If a part is missing, we place `.` in that position so splitting later is consistent.

### Expansion when `input` has multiple groups

When a window has more than one original `input` group (after an across-group merge) **and** a composition includes `input`, the composition expands **per input group** for logic, but is written in a compact way.

**Example setup**

```bash
--compose core input bin
```

**Window values**

* `input = A__B`
* `bin = prox`

**Expanded (logical) values**

* `core` expands to `{ A.prox, B.prox }`

**What gets written**

* For output, we serialize multi-`input` once, then append the other parts.
  - `core` column -> `A__B.prox`
  
* General rule: join all `input` values with `__`, then append the remaining parts joined by `.`

  * `input bin` -> `A__B.prox`
  * `input near-side bin` -> `A__B.-.prox`
  * `input near-name bin` -> `A__B.GENE1.prox`

**What gets counted**

* If you set `--min-per core=N`, counting uses the **expanded values**
  - The example window contributes `+1` to bucket `A.prox` and `+1` to bucket `B.prox`
* The window **passes** the `core=N` rule if **any** of its expanded buckets reaches at least `N` in the final counts

**Recursive interaction with `--min-per`**

Min-per rules describe the desired sizes **in the final output**. We therefore evaluate them **recursively** until the result is stable.

Algorithm sketch:

1. Start with all windows and their current labels. This includes `input` sets and any compositions.

2. For each `--min-per key=K`:

   a. Compute counts for `key` over the **current** windows.

   b. Determine which label buckets for `key` meet the threshold `K`.

3. For each window, test each `key`:

   * If `key` is single-valued for the window, the window survives that key only if its bucket is in the kept set.
   * If `key` expands to multiple values for the window (because it includes `input`), the window survives that key if **any** of its values is in the kept set.

4. Remove anything that no longer qualifies:

   * If an `input` group for a window fails the `input` rule, **remove that group from the window's `input` set**.
   * If a composition includes `input`, its expanded values for that window shrink accordingly.
   * If a window loses **all** values for any required `key`, drop the window entirely.

5. If anything changed in step 4, go back to step 2.
   Repeat until no counts, kept sets, or window assignments change.

6. Write the output from this final stable state:

   * Multi-`input` compositions are serialized in the compact form (for example `A__B.prox`).

Example outcome:

* With `--compose core input bin` and rules `--min-per input=1000 core=250`,
  a window with `input = A__B` and `bin = prox` contributes to `A.prox` and `B.prox` while both `A` and `B` are still in play.
  If `A` later fails the `input` rule during recursion, `A` is removed from that window's `input`, `core` for that window becomes just `B.prox`, counts are recomputed, and recursion continues until nothing else changes.

---

## Choosing which columns to write

```bash
--out-labels <TOKENS...>
```

Each token is either an atomic part or a composition name you defined with `--compose`. Columns are written after the coordinates in the order you specify.

* If you list a composition that was never defined, the command fails early.
* The `input` column for an across-merged window is written as `A__B__C` in stable order.
* A multi-valued composition is written the same way.

Examples

```bash
--out-labels input bin
--out-labels input near-side near-name
--out-labels input core report
  --compose core   input bin
  --compose report core near-side
```

---

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
  - **Any-member** counting: A window with `input = A__B` adds one to `A` and one to `B`.
  - The window satisfies `input=N` if **any** of its input groups ends up with a count at least `N`.

* Single-valued atomic parts (`near-side`, `near-name`, `bin`):
  - Each window adds one to exactly one bucket in that dimension.
  - The window satisfies the rule if the bucket it contributes to reaches the threshold.

* Named compositions:
  - If the composition includes `input`, it expands per window into multiple values (when present) and uses the same any-member logic.
  - If the composition does not include `input`, it behaves like a single-valued atomic key.

**Pruning label assignments during recursive filtering**

When any `--min-per` rule excludes label values for a key, those values are removed from every affected window's assignment for that key during the recursive evaluation. This applies to **all** keys, not just `input`. After each removal, any composition that uses that key is rebuilt from the updated values before the next counting pass. If a window loses all values for any required key, the window is dropped. The process repeats until no counts, kept sets, or per-window assignments change. The final output is written from this stable state, with multi-value compositions serialized in the compact form (for example `A__B.prox`).

---

## Evaluation order

1. Build atomic parts for each window.
2. Build named compositions from the current atomic parts.
   If a composition includes `input`, it can produce multiple values for a window.
3. **Recursive filtering loop**
   Repeat the following until nothing changes:

   * Compute counts for every key referenced by `--min-per` on the **current** window assignments.
   * For each key, determine which label values meet its threshold.
   * For every window and for every key:
     - If the window’s assignment for that key is single-valued and that value does not meet the threshold, drop the window.
     - If the window’s assignment for that key has multiple values (only possible via `input` or a composition that includes `input`), remove the values that do not meet the threshold. If no values remain for that key, drop the window.
   * Rebuild any compositions affected by these removals (for example, compositions that include `input`).
4. Write the requested columns from the final stable state in the order given by `--out-labels`.

---

## Separators and formatting

* Parts inside a composition are joined with a dot.
  Example: `input.bin.near-side`.
* Multiple values inside a single column are joined with `__`.
  Example: `A__B.prox`.
* Missing atomic values are written as `.`.
  - Missing parts in a composition appear as `.` between dots so every composition keeps the same number of fields after splitting.

---

## Validation and clear failures

* Undefined composition in `--out-labels` -> error that names the missing composition.
* Undefined or cyclic composition in `--compose` -> error that points to the first cycle or unknown.
* Unknown key in `--min-per` -> error that lists the allowed keys and known composition names.
* `--min-per` that targets a composition that expands only to `.` for all windows -> results in zero kept windows by design.

---

## Practical recipes

Enforce minimum counts for original groups and for a stratified composition:

```bash
--compose core input bin
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
--compose core input bin
--compose report core near-side
--out-labels input report
--min-per core=150 report=80
```

Filter on input but output composition:

```bash
--compose core input bin
--out-labels core
--min-per input=1000
```

---

## Notes on merging

* **Within-group merges** use the original `input` group as the merge key. The merged window keeps a single `input` group.

* **Across-group merges** union the window's input groups. The `input` column becomes `A__B` and any composition that includes `input` expands accordingly.

---

## Defaults

* No compositions are defined by default.
* `--out-labels` defaults to `input`.
* No min-per rules are applied by default.

---
