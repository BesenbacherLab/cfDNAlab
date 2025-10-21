# Position Selection: How It Works

This document explains how we turn the selection specs into the exact positions we count.

## The pipeline (in one glance)

1. **Per-spec eligible bases (k-independent):**
   Each spec (a frame + a positions string) marks bases on a fragment as **eligible**.
   *Eligible* means: "this base may participate in a k-mer for this frame."
   We do **not** check k-mer length yet.

2. **Intersect specs (k-independent):**
   If multiple specs are given, keep only bases that are eligible in **all** specs.
   The **first** spec is **primary**; later specs are filters only.

3. **Runs (k-independent):**
   Turn the surviving bases into **contiguous windows** ("runs"), e.g. `[start ... end]`.

4. **Final step (k-independent):**
   Apply a single `--step` **after** intersection, using the **first spec’s frame** as the stepping origin. This does not change the runs, it only marks every Nth eligible base as candidate k-mer start/anchor positions. Still no k-mer length involved.

5. **k-mer fit (k-dependent):**
   For each requested k, count only if the full k-mer fits **entirely inside a run** and **inside the fragment’s current segment and tile**. This is the only k-dependent stage.

This "filter -> make runs -> final step -> check k-mer fit" model keeps behavior predictable.

---

## Fragment definition

A fragment contains two inward-directed reads (one forward, one reverse) and starts at the left 5' start (forward.pos) and ends at the right 5' start (reverse.end).

```text
Reference 5' >>>>>>>>>>>>>>> '3
Fragment     |-------------|
Forward   5' |>>>>>>>| 3'     
Reverse        3' |<<<<<<<<| '5 
```

---

## Fragment indexing

* We index the fragment from **left to right**, starting at **0** and ending at **fragment length − 1**.
* All sets and runs are expressed in these left-based indices.

---

## Frames: how a spec marks **eligible bases**

Each frame turns its positions string into **eligible bases** (k-independent). Eligible bases are those that can be part of a counted k-mer.

* **Left**
  Interprets the specified positions from left 5' start (start of fragment) to the right '5 start (end of fragment).
  <br />**Stepping origin** when `Left` is first frame: index `0` at the left 5' end.

* **Right**
  Interprets the specified positions from right 5' start (end of fragment) to the left '5 start (start of fragment).
  <br />**Stepping origin** when `Right` is first frame: index `0` measured from the right end (we map this consistently to left-based indices).

* **Per-end**
  Eligible bases are the **union** eligible bases from `Left` and `Right`.
  <br />**Stepping origin** when `Per-end` is first frame: reflects left and right origins separately.

* **Nearest**
  Works in "distance from the nearest end". For each distance, include **both** fragment indices (one on each side).
  For odd fragment lengths, the exact middle base is **excluded** (both sides count up to, but not including, the center).
  <br />**Stepping origin** when `Nearest` is first frame: distance `0`, mirrored to both sides.
  
  Fragment:
  Origins:   .           .
  Positions: |>>>>>*<<<<<|
  Folded:
  Origin:    .
  Positions: |>>>>>|

* **Mid**
  Works in signed offsets around the middle. Counting is performed in the forward direction (left -> right).
  Let **center** be `floor(fragment length / 2)`.
  For **odd** lengths, index `center` is the physical middle.
  For **even** lengths (L = 2m), there is no single middle base; the center lies between indices `m-1` and `m`. We define the `Mid` frame so that:
    * Offset `0` maps to index `m` (the base right of center),
    * Offset `-1` maps to index `m-1`,
    * Offset `+1` maps to index `m+1`,
    * In general, offset `d` maps to index `m + d`.
  <br />**Stepping origin** when `Mid` is first frame: offset `0`; stepping is symmetric (..., −step, 0, +step, ...).

---

## Multiple specs

* The **first** spec is the **primary** spec. It decides how results are expressed:
  - Which frame the reported positions refer to (`Left`, `Right`, `Mid`, `Nearest`, `Per-End`).
  - Which strand/orientation is used when counting (forward vs reverse; reverse means bases are complemented).
  - The output structure:
    - One or two tracks (`Per-End`).
    - Folded around the middle (`Nearest`) or unfolded.
  - Which origin to `--step` from.
* Later specs **only filter**. They can only remove eligible bases.
* We intersect **eligible bases** across all specs before any stepping or k-mer checks.

---

## Runs (windows)

Intersecting several specs can produce several **disjoint regions**. We compress the surviving bases into inclusive **runs**:

```
[ start0 ... end0 ], [ start1 ... end1 ], ...
```

Runs encode all frame constraints and are respected by the final step and by the k-mer fit.

---

## Final `--step` (after intersection)

`--step N` is applied **once**, after runs are formed, using the **first spec’s frame** as the origin. We:

* Determine the stepping origin from the first frame (`Left`: index 0, `Right`’s right-origin, `Mid`’s offset 0 with symmetry, `Nearest`’s zero-distance mirrored, `Per-end`’s per-side origin).
* Do not alter the runs. Instead, within those runs, treat only the bases that land on the step lattice as candidate k-mer start/anchor positions.

Why one final step instead of per-spec steps?

* Matches expectations: **filter first, then downsample**.
* Avoids confusing phase interactions between multiple per-spec steps.

---

## K-mer fit

For each k-mer size:

* **Forward**
  A forward k-mer starting at index `s` spans `[s ... s + (k − 1)]` (both ends inclusive).
  It is valid only if this entire span stays **inside one run** and **inside the current fragment segment and tile**.

* **Reverse**
  A reverse k-mer ending at index `a` spans `[a − (k − 1) ... a]` (both ends inclusive).
  It is valid only if this entire span stays **inside one run** and **inside the current fragment segment and tile**.

A convenient way to check this in code:

* For **forward**, temporarily shrink each run’s **right edge** by `(k − 1)`; then a start is valid if it lies inside the shrunk run.
* For **reverse**, temporarily shrink each run’s **left edge** by `(k − 1)`; then an anchor is valid if it lies inside the shrunk run.

This guarantees k-mers never cross **breakpoints** between disjoint regions and never cross segment/tile boundaries.

---

## Example

In this example, we get the ±50bp around the middle of the fragment but never closer than 10bp away from the fragment ends.

Command:

```bash
--frame mid --positions "-50..50" \
--frame left --positions "10..-10" \
--step 1
```

1. **Eligible bases per spec:**

   * `Mid`: bases within ±50bp of the middle (per even/odd convention).
   * `Left`: bases at least 10bp from both ends.
2. **Intersection:** keep only bases present in **both** sets; typically two lobes near the middle.
3. **Runs:** compress those lobes into `[start ... end]` windows.
4. **Final step:** first frame is `Mid`, mark every Nth base (using `Mid`’s origin, symmetrically) as count site candidates.
5. **k-mer fit:** for k=5, we trim 4 bases from forward run ends (since `Mid` counts in the forward direction). The **runs themselves do not change** with k.

---

## Invariants

* Per-spec eligibility, intersection, runs, and final step are **k-independent**.
* Only the k-mer fit step depends on k.
* The first spec controls labeling and stepping origin; later specs only filter.
* If no bases remain after intersection and stepping, nothing is counted for that fragment length and k.
* Runs are never modified by --step; step only selects candidate positions within runs.

---

## FAQ

**Why not step per spec?**
Intersecting several periodic grids creates surprising sparsity and phase effects. A single final step is predictable and easy to explain.

**Why do sites near boundaries disappear for larger k?**
Because the entire k-mer must fit inside a run and inside the current segment/tile. Larger k leaves fewer valid positions near run boundaries. The runs themselves are unchanged.

**How is the exact middle handled?**
For odd lengths, the exact middle base may be excluded by the `Nearest` rules (both sides count up to it). For even lengths, we adopt the mapping above so the center has a clear origin for the `Mid` frame.

---

## Implementation checklist (for contributors)

1. For each fragment length:
   a. For each spec, compute **eligible bases** (set of left-based indices).
   b. Intersect across specs.
   c. Compress to **runs**.
   d. Apply the **final** `--step` using the first frame’s origin; keep the stepped bases.
   e. Keep the **first spec’s** `PositionSelection`s filtered to those bases (so orientation/group come from the first spec). Optionally store the runs for fast k-mer fit checks.

2. At count time (per fragment, per k):
   a. Clip to segment/tile.
   b. Enforce **k-mer fit** (shrink runs by `k − 1` on the appropriate side).
   c. Emit counts.

---
