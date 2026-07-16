# Ends Reference Correction Modes Plan

## Purpose

Define end-motif reference-correction modes for outputs whose motif labels
contain both outside and inside bases, such as `AC_TG`.

This feature is only a mode choice when both motif sides are present. One-sided
outputs already have only one motif axis, so they do not accept an
explicit mode value.

The public argument name is `two_sided_correction`.

```text
two_sided_correction = None
two_sided_correction = "joint"
two_sided_correction = "split"
two_sided_correction = "outside"
two_sided_correction = "inside"
```

`None` is not a correction mode. It means no two-sided mode was supplied.
Non-empty two-sided corrected outputs have no default correction mode. Users
must choose whether to correct full motifs, keep full labels while correcting
outside and inside labels separately, or aggregate to one side.

## Mode Availability

Mode availability depends on the loaded motif axis.

When both inferred side widths are positive:

- `two_sided_correction = None` is an error.
- `"joint"`, `"split"`, `"outside"`, and `"inside"` are valid.

When only the outside width is positive:

- `two_sided_correction = None` is required.
- Any explicit value is an error.
- Correction divides each `outside_` label by its matching reference
  denominator.

When only the inside width is positive:

- `two_sided_correction = None` is required.
- Any explicit value is an error.
- Correction divides each `_inside` label by its matching reference
  denominator.

Do not allow `joint`, `split`, `outside`, or `inside` on one-sided outputs.
Those values describe alternative two-sided behavior, and there is no two-sided
choice to make when the output has only one motif side.

Empty motif axes may be rejected by a loader before correction. If a loader
does expose such an output, return an empty result without inferring side
widths, because there are no corrected motif units. Do not require at least one
loaded motif label just to prove that an empty output is two-sided.
Since side widths cannot be inferred, do not validate whether an explicit
`two_sided_correction` value would have matched the missing sides. The chosen
value does not change the empty result. Motif selectors still target the mode
motif axis and must error if they request labels that are not present.

## Mode Summary

`joint` keeps full `outside_inside` motifs and corrects each full motif with
the matching reference k-mer denominator.

`split` keeps full `outside_inside` motifs, but corrects each full motif by
combining the separately support-normalized outside and inside reference
denominators. For `AC_TG`, the observed `AC_TG` count stays on the output axis,
while the correction uses the outside support for `AC` and the inside support
for `TG`. This mode is useful when exact joint reference support is too sparse.
It does not require the exact joint reference motif to have positive frequency
when both side denominators are positive.

`outside` returns outside motifs only. It sums joint end-motif counts over all
inside motifs with the same outside side, then corrects each outside motif with
the outside reference denominator.

`inside` returns inside motifs only. It sums joint end-motif counts over all
outside motifs with the same inside side, then corrects each inside motif with
the inside reference denominator.

## Invariants

- Keep one-sided correction behavior unchanged and require no explicit mode.
- Require an explicit two-sided mode when both sides are present.
- Keep selected row identity unchanged for every correction mode.
- Keep joint motif labels unchanged for `joint` and `split`.
- Return side-specific motif labels for `outside` and `inside`.
- Treat missing sparse reference entries as zero frequency.
- Keep `use_global_bias` semantics unchanged. A global reference can correct
  non-global end rows only with explicit opt-in.
- Compute reference support from the matched reference row and correction
  universe, not from motifs observed in the sample or selected by the loader.
- Reject positive corrected units with no positive reference denominator under
  the unsupported-reference policy.
- Motif-group end-motif outputs accept only `two_sided_correction = None`.
- Reject canonical reference k-mer outputs for all correction modes.

## Definitions

For an end-motif label:

```text
label = outside_inside
outside = bases before "_"
inside = bases after "_"
```

Side labels preserve the separator:

```text
outside label for "AC_TG" = "AC_"
inside label for "AC_TG" = "_TG"
```

This keeps side-specific output labels compatible with one-sided end-motif labels
when one side has length zero, and it prevents `AC` from being ambiguous between
outside and inside output.

Both `cfdna ends` and `cfdna ref-kmers` expose forward-oriented public motif
labels. Right-end motifs have already been reverse-complemented by `cfdna ends`.
Reference-correction loaders split and match public labels directly. They must
not apply another reverse-complement step.

Per matched row, define:

```text
count[outside, inside]
joint_reference_frequency[outside, inside]
```

Side reference frequencies are derived from the loaded joint reference row:

```text
outside_reference_frequency[outside] =
    sum over inside joint_reference_frequency[outside, inside]

inside_reference_frequency[inside] =
    sum over outside joint_reference_frequency[outside, inside]
```

Reference support counts are:

```text
joint_support_count =
    number of joint motifs with positive joint_reference_frequency

outside_support_count =
    number of outside motifs with positive outside_reference_frequency

inside_support_count =
    number of inside motifs with positive inside_reference_frequency
```

Reference denominators are:

```text
joint_reference_denominator[outside, inside] =
    joint_reference_frequency[outside, inside] * joint_support_count

outside_reference_denominator[outside] =
    outside_reference_frequency[outside] * outside_support_count

inside_reference_denominator[inside] =
    inside_reference_frequency[inside] * inside_support_count
```

The word "denominator" is intentional for internal discussion. Do not introduce
"scale" as the public name for the new two-sided modes.

## Correction Formulas

### Joint

Output axis:

```text
outside_inside
```

Corrected count:

```text
corrected[outside, inside] =
    count[outside, inside] /
    joint_reference_denominator[outside, inside]
```

### Split

Output axis:

```text
outside_inside
```

Corrected count:

```text
corrected[outside, inside] =
    count[outside, inside] /
    (outside_reference_denominator[outside] * inside_reference_denominator[inside])
```

The output still contains joint motifs. Only the reference correction
denominator changes.

### Outside

Output axis:

```text
outside_
```

Observed count:

```text
outside_count[outside] =
    sum over inside count[outside, inside]
```

Corrected count:

```text
corrected_outside[outside] =
    outside_count[outside] / outside_reference_denominator[outside]
```

### Inside

Output axis:

```text
_inside
```

Observed count:

```text
inside_count[inside] =
    sum over outside count[outside, inside]
```

Corrected count:

```text
corrected_inside[inside] =
    inside_count[inside] / inside_reference_denominator[inside]
```

## Reference Conditioning

`split`, `outside`, and `inside` derive side frequencies from the loaded joint
reference k-mer output. They do not require separate `ref-kmers` runs for
`k_outside` and `k_inside`.

This intentionally conditions side frequencies on the joint reference
opportunities represented by the loaded reference row. It can differ from
separate shorter-k reference runs near row edges or masked bases. That is
acceptable because correction stays tied to the same eligible joint motif space
as the end-motif output and avoids extra reference passes.

The correction universe comes from the reference package, not from sample
sparsity and not from loader motif selectors. If the reference output was built
from a motifs file, side frequencies are conditioned on that selected
full-motif set. Unlisted full motifs are not part of the denominator.

The mode motif axis comes from the loaded end-motif output, not from
reference-only labels. For `joint` and `split`, the mode motif axis is the
loaded joint end-motif axis. For `outside` and `inside`, the mode motif axis is
the side-specific labels represented by loaded joint end-motif labels.
Reference labels that are absent from the end-motif output can affect
correction denominators, but they do not create mode-axis rows or matrix
columns.

For `outside` and `inside`, dedupe side labels in the loaded end-motif axis
order. This keeps deterministic or motifs-file-defined source order while
making side-mode matrix columns inspectable after `read()`. Label selectors
return the requested side labels in requested selector order.

For large k-mer sizes, use a motifs file consistently for `ends` and
`ref-kmers` to define a tractable correction universe. Do not add a
selector-dependent denominator mode. The motifs-file route is the documented
way to define a large or custom correction universe.

## Validation

Apply these checks before correction when they are relevant to a non-empty mode
motif axis:

- Reference output must not be canonical.
- End-motif labels must contain exactly one `_`.
- All end-motif labels must agree on outside and inside widths.
- End-motif combined width must equal `ref_kmers.kmer_size`.
- Reference labels must split cleanly using the inferred outside and inside
  widths.
- Two-sided correction modes require positive inferred outside and inside
  widths.
- One-sided outputs require `two_sided_correction = None`.
- Motif-group end-motif outputs require `two_sided_correction = None`.
- End-motif and reference outputs must both use motif labels or both use
  motif-group labels.

Do not treat a missing side as a single empty motif for public `split`
correction. That would make `split` a redundant alias for one-sided correction.

## Unsupported Reference Policy

For `joint`, unsupported-reference handling is evaluated per selected joint
motif cell.

For `split`, unsupported-reference handling is evaluated on the product of the
outside and inside denominators. A positive joint count with zero or missing
exact joint reference frequency is still supported when both side denominators
are positive.

For `outside`, unsupported-reference handling is evaluated after summing counts
to the outside motif axis.

For `inside`, unsupported-reference handling is evaluated after summing counts
to the inside motif axis.

Positive counts with denominator zero follow the unsupported-reference policy:

```text
error   -> fail with the unsupported motif labels
keep_na -> return NA or NaN corrected counts for those positive units
drop    -> drop denominator-zero units from data-frame output
```

For `drop`, denominator-zero units are removed regardless of observed count.
For `error` and `keep_na`, zero counts with denominator zero remain corrected
zero.

`drop` is an R and Python data-frame policy. Rust corrected selections keep a
fixed shape and use only `Error` or `KeepNaN` unless a separate Rust drop API is
added later.

## Selector Semantics

Row selectors behave the same in all correction modes.

Motif label selectors always select the mode motif axis:

- `joint`: accepts only joint labels such as `AC_TG`.
- `split`: accepts only joint labels such as `AC_TG`.
- `outside`: accepts only outside labels such as `AC_`.
- `inside`: accepts only inside labels such as `_TG`.

Selector labels from the wrong axis are errors. For example, `AC_TG` is invalid
for `outside`, and `AC_` is invalid for `joint` or `split`.

For `outside`, selecting `AC_` returns only that outside motif. Internally, its
count is still:

```text
sum over all loaded joint motifs AC_inside
```

For `inside`, selecting `_TG` returns only that inside motif. Internally, its
count is still:

```text
sum over all loaded joint motifs outside_TG
```

Motif selectors must not change reference support. Selecting one mode-axis motif
must match running the full correction and filtering the result afterward.

Motif index selectors are valid for `joint`, `split`, and one-sided correction.
For `outside` and `inside`, reject motif index selectors. Side-mode label
selectors are clear, but there is no public method to inspect side-mode indices
before `read()`.

For `joint`, `split`, and one-sided correction, implementations may
apply motif selection before correction as an optimization, as long as reference
support and denominators are still computed from the full correction universe.
This keeps the same values as full correction followed by filtering.

Pre-aggregation full-motif filtering for `outside` and `inside` is not part of
this API. It changes the meaning of the side count and would need a separate
explicit API.

## Result Values

Corrected loaders return corrected values, not formula internals.

Matrix helpers return corrected counts. They have one matrix row per selected
end-motif output row and one matrix column per motif on the selected motif axis.

Corrected data-frame helpers keep selected row metadata and a `motif` column on
the selected motif axis. They keep the `count` column. For `outside` and
`inside`, `count` is the aggregated side count. Motif index metadata, where
returned, describes the selected mode axis rather than the source joint motif
axis. R keeps its public one-based index convention. They add:

```text
corrected_count
corrected_frequency
```

`corrected_frequency` is computed over the full mode motif axis after applying
the correction mode and unsupported-reference policy, but before motif
selection. Motif selection filters the corrected data frame afterward and does
not renormalize frequencies.

Per output row:

```text
corrected_frequency[motif] =
    corrected_count[motif] /
    sum over full mode motif axis corrected_count
```

If the finite full-axis corrected count total is zero, finite
`corrected_frequency` values for that row are zero. With `keep_na`, if any
full-axis corrected count needed for the row total is `NA` or `NaN`, all
`corrected_frequency` values for that row are `NA` or `NaN`. With `drop`,
denominator-zero units are removed before row totals are computed.
Implicit zero sparse cells do not need to be materialized for this row total,
because they do not change the denominator.

`joint` and `split` use the same result columns, because only the correction
denominator differs. `outside` and `inside` also use the same column shape, with
side-specific motif labels.

Do not include mode-specific formula columns such as `reference_frequency`,
`reference_denominator`, or correction motif counts in the default output. Do
not add a public diagnostic helper for these internals as part of this plan.
This intentionally replaces the unreleased reference-correction data-frame
columns that exposed formula internals.

## R Loader API

Add `two_sided_correction` to corrected R helpers.

Affected public methods:

- `end_motif_data_frame(..., ref_kmers = ref_kmers)`
- `dense_corrected_counts_matrix()`
- `sparse_corrected_counts_matrix()`

Validation:

- `two_sided_correction` with `ref_kmers = NULL` is an error.
- One-sided loaded outputs require the argument to be `NULL` or omitted.
- Non-empty two-sided loaded outputs require one of `"joint"`, `"split"`,
  `"outside"`, or `"inside"`.
- Empty loaded outputs that are exposed by the loader return empty without
  requiring a motif label only to prove the two-sided shape.
- Motif label selectors use the mode motif axis for the chosen mode.
- Motif index selectors are rejected for `"outside"` and `"inside"`.

## Python Loader API

Add the same `two_sided_correction` argument to Python corrected helpers.

Affected public methods:

- `data_frame(..., ref_kmers=ref_kmers)` paths that call reference correction
- `corrected_counts_array()`
- `sparse_corrected_counts_matrix()`

Validation and result values must match R. Passing `two_sided_correction`
without `ref_kmers` is an error.

## Rust Loader API

Add:

```rust
pub enum TwoSidedCorrectionMode {
    Joint,
    Split,
    Outside,
    Inside,
}
```

Add a builder method for selecting the two-sided correction model:

```rust
pub fn two_sided_correction(mut self, mode: TwoSidedCorrectionMode) -> Self
```

The method name is `two_sided_correction`.

Keep returning `EndMotifCountSelection`. Add an internal constructor path for
corrected selections that can replace the motif labels, motif indices, and count
shape together. For corrected side modes, motif labels and motif indices
describe the selected side axis in the corrected result, not the loaded
full-motif axis.

Validation:

- One-sided loaded outputs reject any explicit `TwoSidedCorrectionMode`.
- Non-empty two-sided loaded outputs require an explicit
  `TwoSidedCorrectionMode`.
- Empty loaded outputs return empty without requiring a motif label only to
  prove the two-sided shape.
- Motif-group outputs reject any explicit `TwoSidedCorrectionMode`.
- Motif label selectors use the mode motif axis for the chosen model.
- Motif index selectors are rejected for `Outside` and `Inside`.

## Implementation Notes

Build a per-reference-row correction cache:

```text
joint frequencies by joint motif
outside frequencies by outside motif
inside frequencies by inside motif
joint support count
outside support count
inside support count
```

Dense reference rows can be scanned directly. Sparse reference rows should be
processed from stored entries without expanding to a dense reference matrix.

For dense end-motif counts, side modes aggregate from all loaded joint motif
columns that contribute to the requested side-mode labels.

For sparse end-motif counts, side modes aggregate stored entries into
side-specific sparse rows. Missing sparse entries are implicit zero and do not
need materialization.

Row selection may happen before correction. Motif selection for `outside` and
`inside` is selection on the side mode axis, so it must not pre-filter joint
motif columns before side aggregation.

For side modes, build the full side-mode axis from loaded joint motif labels,
aggregate all loaded joint motif columns into that axis, then apply side-mode
motif label selection. This keeps side counts independent of selector
pre-filtering.

For `joint`, `split`, and one-sided correction, motif selection may
happen before correction for speed, provided denominators and reference support
do not use the selected motif set.

## Tests

Add public loader tests, not command reruns.

Rust tests:

- One-sided outputs accept no explicit mode and reject all explicit modes.
- Non-empty two-sided outputs require an explicit mode.
- `joint` keeps full `outside_inside` motif labels and uses the matching
  full-reference denominator.
- `split` keeps full `outside_inside` motif labels and combines the outside and
  inside side denominators.
- `split` supports a positive joint count when the exact joint reference
  frequency is zero but both side denominators are positive.
- `outside` returns outside labels and summed outside counts.
- `inside` returns inside labels and summed inside counts.
- Side-mode axes dedupe labels in loaded end-motif axis order.
- Sparse and dense fixtures produce matching corrected values.
- Motif label selectors select the mode motif axis and reject wrong-axis labels.
- Side-mode motif index selectors are rejected.
- Side-mode motif selectors match full-output filtering on the mode motif axis.
- Positive counts with denominator zero follow the unsupported-reference policy.
- Global reference opt-in works for every correction mode.
- Motif-group inputs accept no explicit mode and reject every explicit mode.
- Canonical reference inputs are rejected.
- Empty motif axes return empty outputs without side-width inference.

R tests:

- `two_sided_correction` without `ref_kmers` is an error.
- One-sided outputs accept no explicit mode and reject all explicit modes.
- Non-empty two-sided outputs require an explicit mode.
- Data-frame output exposes `corrected_count` and `corrected_frequency` for
  each correction mode.
- Side-mode data frames keep `count` as the aggregated side count.
- Dense and sparse matrix helpers keep the selected shape for each correction
  mode.
- `unsupported_motifs = "drop"` remains unsupported for fixed-shape matrix
  helpers.
- Motif label selectors select the mode motif axis and reject wrong-axis labels.
- Side-mode motif index selectors are rejected.
- Selector behavior matches the full-output-then-filter invariant.

Python tests:

- Match the R behavior for mode validation, returned columns, matrix shape,
  sparse output, unsupported motifs, motif selectors, and explicit
  `two_sided_correction` without `ref_kmers`.

Do not run full project tests for this plan. When implementing code later,
follow the project rule: run `cargo check` after code changes and
`cargo check --tests --features testing` after test code changes. Do not run
tests unless explicitly requested.
