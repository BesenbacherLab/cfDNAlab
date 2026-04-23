## `cfdna outliers` spec — review

Date: 2026-04-22

Reviewing `CFDNA_OUTLIERS_SPEC_2026-04-21.md` after reading the
LIONHEART code in `plans_and_specs/lionheart-main/` (Nature Communications
2025, `BesenbacherLab/lionheart`). This complements the spec's web-only
reading with what the code actually does.

The review covers the hard choices the spec identifies:

- what LIONHEART actually does vs what the spec says it does
- local vs global
- discreteness and the ZIP model
- architectural fit between a "clip" semantic and a "multiplicative weight"
  semantic
- open questions the spec leaves dangling

Where this review disagrees with the spec I state the disagreement plainly;
where it only tightens the spec I say so.

## Orientation: LIONHEART is inspiration, not a target

`cfdna outliers` is not obliged to reproduce LIONHEART. LIONHEART is a
convenient reference point because it is published, recent, and has code
we can read, but if `cfdna outliers` diverges that is a design choice,
not a failure. Read every "LIONHEART does X, we don't" passage below as
a data point about what a competent published pipeline chose, not as a
requirement for this command.

In particular:

- Where I say "pick either a LIONHEART-faithful port or a clean local
  method, don't hybrid", the force of the argument comes from
  "implementing two methodologies at once is costly and neither gets
  properly validated in v1". That argument still holds without
  LIONHEART as the anchor — it would just be "pick a local method
  carefully" vs "be free to hybrid".
- The technical facts about LIONHEART's ZIP class, order of operations,
  thresholds, and aggregation logic stay true regardless of whether
  `cfdna outliers` imitates them.
- Where the spec proposes features LIONHEART does not have (local MAD,
  fit-space discretization branches, interval-width-aware weighting),
  my critique is on the internal consistency of those features, not on
  their divergence from LIONHEART.

With that framing set, the rest of the review is unchanged.

## TL;DR

- LIONHEART's per-sample outlier step is **pure global chromosome-wise ZIP
  tail clipping on raw integer counts**. It is not a hybrid and contains no
  local-MAD / local-residual logic. The spec implies more nuance than the
  actual pipeline has.
- Their broad baseline (5 Mb window / 500 kb stride, division by mean, no
  centering) is applied **after** the clip and after GC/insert-size
  correction. So for them the clip space is raw counts, not residuals.
- Their ZIP class is not a standard MLE ZIP. `mu` is the overall mean
  (zeros included), and `p_nz = n_non_zero / n`. Upper-tail probability
  collapses to `p_nz * (1 - Poisson.cdf(k; overall_mean))`. That is a
  deliberate simplification; it biases the clip threshold downward
  (more conservative) relative to a textbook ZIP MLE.
- At **inference time**, the ZIP is fit **per chromosome** with a very strict
  tail threshold (`1 / 263_000_000` ≈ 3.8e-9, roughly `1 / N_bins`). At
  **resource creation** (building the recurrent blacklist), the ZIP is fit
  **genome-wide**, per sample, with a much looser tail threshold (`1e-4`),
  and candidates are then intersected-by-recurrence across samples. The spec
  collapses these into one sentence and loses that distinction.
- For a soft-weight v1 the honest recommendation is narrower than the spec's
  "hybrid": commit to one detector family and ship it well. A global ZIP
  tail on raw counts, an empirical per-chromosome quantile cap, or a local
  rolling-median/MAD residual are all defensible alone. The hybrid of
  "global tail + local MAD + broad-baseline normalization" is the most
  ambitious option and, in cfDNA specifically, the least validated.
  The point is commit-and-ship, not "match LIONHEART".
- There is a real architectural misalignment the spec flags but does not
  resolve: LIONHEART clips **counts post-correction**, which maps to a
  multiplicative weight of `min(1, cap / observed)` **per bin** — not per
  fragment. That clashes with cfDNAlab's fragment-averaged weight semantics
  for `midpoints`, `ends`, `lengths`, `bam_to_frag`. This should be called
  out in v1 docs, not papered over.

The rest of this document expands each bullet.

## 1. LIONHEART's actual outlier pipeline (from code, not paper text)

### 1.1 Inputs at inference time

- BAM file filtered to autosomes, MAPQ ≥ 20, fragment length 100–220 bp,
  fragment mode.
- `mosdepth --by 10` → integer fragment count per 10 bp bin.
- Sparse chromosome arrays are built from the nonzero bins.
- At load time, sparse arrays are rounded to `decimals=2` because
  downstream averaging can introduce small float error; the underlying
  counts are still integers.

### 1.2 Order of operations at inference

From `lionheart/features/create_dataset_inference.py`:

1. Load `include_indices` per chromosome. This already excludes
   ENCODE blacklist v2 + Umap k100 non-mappable bins + the recurrent
   outlier resource + the all-zero recurrent resource — they are all
   unioned into a single `exclude_paths` set passed in from the CLI.
2. Load sparse coverage, round to 2 decimals.
3. `poiss.reset().partial_fit(np.round(sample_cov).astype(np.int64))` —
   **per chromosome**, on already-filtered bins.
4. Iterate the fitted ZIP upward starting at `floor(mean)` until
   `1 - cdf(k) < 1 / 263_000_000`. Call that value `clipping_val`.
5. `sample_cov[sample_cov > clipping_val] = clipping_val` — **clip, not
   drop, not NaN.**
6. GC correction (chromosome-level).
7. Insert-size correction (three separate chromosome-level factors).
8. Megabin normalization: 5 Mb windows with 500 kb stride (i.e. 10
   overlapping windows per bin, averaged), divide coverage by that mean.
   No centering. `clip_above_quantile` is **None** in this call — the
   `<= 0.99` quantile clip mentioned in the megabin docstring is a
   capability that inference never exercises.
9. Remove consensus-site bin indices; compute Pearson R against each
   cell-type mask.

Important consequences:

- The blacklist / mappability / recurrent-outlier exclusions happen
  **first**, before the ZIP fit even looks at the data.
- The ZIP fit therefore operates on the already-masked counts, which makes
  its upper-tail estimate stricter in practice — the most extreme
  recurrent nuisance loci are already gone.
- The megabin baseline is computed **after** the clip, so broad baseline
  normalization does not see the original extreme values at all.

The paper's phrasing "fits chromosome-wise zero-inflated Poisson
distributions per sample, derives per-sample upper-tail thresholds and
truncates sample-specific extreme coverages, separately removes recurrent
empirical extreme bins across datasets, then divides coverage by a broad
5 Mb / 500 kb overlapping-average baseline" is **order-preserving but
misleading** — at inference the recurrent-extreme removal has already
happened as part of step 1, before the chromosome-wise ZIP; and the ZIP
clip already happens before the broad baseline divide.

### 1.3 Resource creation — candidate detection per sample

From `resource_creation/detect_outlier_candidates.sh` +
`calculate_tail_cdf.py`:

- mosdepth bin size passed as a workflow argument (10 bp per the paper).
- Coverage counts rounded to integers.
- `mean`, `n`, `n_nonzero`, `max_count`, `p_nonzero = n_nonzero / n`
  computed **genome-wide**, over all autosomes in one pass.
- Single tail-probability lookup table computed once using overall
  `mean` and `p_nonzero`.
- Threshold for "candidate" is `1e-4` (the shell script passes
  `1e-4` as the threshold argument).
- Only bins where `count > mean` are written out.
- Candidates: `chrom, chrom_idx, count, tail_prob`.
- All-zero bins: `chrom, chrom_idx`.

This is **not** chromosome-wise. The per-chromosome ZIP fit only exists at
inference time. The candidate-detection step is single-pass, genome-wide,
and uses a much looser threshold because the subsequent step requires
recurrence in ≥25% of samples to promote a candidate to an actual outlier.

### 1.4 Resource creation — cohort aggregation

From `collect_outliers.py`:

- Default `--thresholds 1e-4`, `--out_ofs 0.25`: a bin is a dataset-level
  outlier if it appears as a candidate in ≥ 25% of samples.
- Multiple `(threshold, out_of)` pairs can be passed; the union of the
  surviving bins is kept. The parameterization allows "strict threshold +
  low recurrence" combined with "loose threshold + high recurrence".
- All-zero bins are processed differently: kept only if **every** sample
  in the dataset has them as zero (intersection across samples).

From `collect_outliers_across_datasets.py`:

- Outliers across datasets: **union** by default.
- All-zero bins across datasets: **intersection** by default.

This asymmetry is deliberate and conservative in opposite directions:

- Outliers: broad net — any dataset can add a recurrent-outlier bin.
- All-zero bins: narrow net — only bins that are always zero everywhere
  are excluded as presumably unmappable.

The spec lumps these into one "recurrent empirical bin resource" and does
not note the asymmetry between outlier-union and zero-intersection, which
is relevant if the package ever builds a cohort resource of its own.

## 2. The ZIP model — what the code actually does

From `lionheart/features/correction/poisson.py`, class `ZIPoisson`:

- `mu = mean(x)` computed over the whole input including zeros.
- `p_nz = n_non_zero / n`, where `n_non_zero` is the number of
  strictly-positive entries.
- `pmf(k) = p_nz * Poisson(k; mu) + (1 - p_nz) * [k == 0]`
- `cdf(k) = (1 - p_nz) + p_nz * Poisson.cdf(k; mu)`
- Upper tail: `1 - cdf(k) = p_nz * (1 - Poisson.cdf(k; mu))`.

This is **not the MLE of a standard ZIP model**. A standard ZIP fits
`pi` (zero-inflation weight) and `lambda` (Poisson rate of the
non-degenerate component) jointly from the likelihood; `lambda` of the
MLE is typically larger than the sample mean because some of the zeros
are absorbed into the point-mass at zero. LIONHEART sidesteps this:
`lambda` is fixed at the overall sample mean and `1 - pi` is set
directly from the observed non-zero fraction.

Consequences:

- The Poisson component is **narrower** than the MLE would pick (smaller
  lambda). That makes any given `k` look farther into the tail than a
  proper ZIP-MLE would say. With a strict threshold like 3.8e-9 this
  just means slightly more bins get clipped — conservative in the
  "clip more, keep less" direction.
- At very strong zero inflation the simplification matters more; at
  moderate zero inflation (typical cfDNA at 10 bp fragment-mode) it
  is a small bias.
- The tail `p_nz * (1 - Poisson.cdf(k; mu))` is dominated by the
  Poisson-CDF term in the relevant regime (`k >> mu`). `p_nz` is a
  scalar shrinkage factor.

### 2.1 Discreteness — where the spec is right, and where it glosses

The spec's three-way fit-space distinction (true-count / pseudo-count /
continuous-residual) is genuinely useful. LIONHEART sits squarely in
"true-count": their ZIP sees integer fragment counts in 10 bp bins.

But LIONHEART also uses `np.round(sample_cov).astype(np.int64)`
before fitting at inference time. At that point `sample_cov` has already
been loaded from a sparse array with `decimals=2`. So they **do** round
a float back to int before the ZIP — the code treats this as a belt-and-
suspenders guarantee for upstream rounding error, not as a real
discretization of a continuous signal.

The spec's "pseudo-count branch" scenario — fitting ZIP on GC-corrected
or length-normalized values — is a **different** heuristic that
LIONHEART does not use and does not validate. Two reasons to be more
skeptical of it than the spec currently is:

- GC-corrected values are multiplicative rescalings of counts by
  GC-bin-specific factors, then clipped or corrected further. The
  resulting distribution is no longer Poisson-like even approximately:
  it is a mixture of rescaled Poissons with multiple rescaling factors.
  Fitting a single ZIP gives you a stand-in whose tail depends heavily
  on the correction factor distribution.
- For `--normalize-by-length` outputs, "rescale by representative
  fragment length" only restores a count-like scale in the mean. The
  variance structure is still ratio-like and the zeros are no longer
  hard zeros in general.

Recommendation that tightens the spec:

- v1 should operate **only** on true-count inputs for the parametric
  tail-score branch. Raw fragment counts, or raw base-coverage counts,
  in fixed bins, pre-GC and pre-length-normalization. This is what
  LIONHEART validated.
- The pseudo-count and continuous branches should be explicitly listed
  as "not supported in v1" rather than framed as coequal options.

### 2.2 A subtler point the spec misses

LIONHEART starts the iterator at `floor(mean)` and walks up until the
tail goes below threshold. It does **not** check that the resulting
`clipping_val` actually corresponds to a count present in the sample.
On shallow samples where `max(sample_cov)` is small, `clipping_val`
can exceed `max(sample_cov)`, in which case no bins are clipped at
all. That is fine mathematically but worth knowing — the effective
clip fraction varies with depth.

## 3. Local vs global — the spec's central design tension

The spec frames v1 around a "hybrid: broad baseline normalization →
local MAD residual + global tail companion → soft cap". Reading the
LIONHEART code, that framing is stronger than what LIONHEART itself did.

What LIONHEART actually does:

- **Global chromosome-wise ZIP clip** on raw integer counts.
- **No local MAD, no rolling median, no Hampel** anywhere in the
  pipeline.
- Broad baseline normalization (5 Mb / 500 kb) is a **separate
  step after the clip**, and it is a centering/scaling tool for
  CNA-like shifts, not an outlier detector.

So the spec's current framing "LIONHEART's method conceptually mixes
separate layers: ... single-sample upper-tail clipping + broad baseline
normalization" is factually right but gives equal weight to two things
that do very different work:

- the clip is the outlier step;
- the megabin normalization is the broad-baseline step;
- there is no local-residual step at all.

Three defensible v1 shapes, roughly in increasing ambition:

- **A. Global chromosome-wise ZIP tail clip on raw counts.** Fit ZIP
  per chromosome, cap at a strict tail threshold, emit
  `w_i = min(1, cap / c_i)`. Minimal scope. Validated analogue in
  LIONHEART, so at least one cfDNA pipeline has shown this shape works
  downstream. Good if you want a v1 that is obviously correct and easy
  to defend.

- **B. Empirical per-chromosome quantile cap on raw counts.** No
  parametric model. Pick a strict upper quantile (e.g. `1 - 1/N_bins`)
  and clip above it. Same emission shape as A. More robust to
  misspecification, but no reference pipeline for cfDNA. Good if you
  want to avoid parametric-model arguments entirely.

- **C. Local rolling-median / MAD residual on broad-baseline-normalized
  support.** The spec's preferred direction. The most ambitious option:
  two independent methodological choices (baseline + detector) that
  interact. Defensible in CNV literature (`copynumber`, LMADCNV) but
  not validated in cfDNA. Good if you want the command to detect
  **locally** weird bins rather than globally extreme ones, and if you
  are willing to accept that v1 is a novel methodological claim.

**My recommendation: pick one, don't hybrid.** The spec's current answer
— a hybrid — is the worst v1 shape because it commits to two
methodological risk surfaces at once and neither gets validated. That
recommendation stands independently of LIONHEART; it is about scoping.

Which of A / B / C to choose depends on the intended scientific
positioning:

- If `cfdna outliers` is scoped as "a residual cap after standard
  filtering", B is the cleanest v1 and has the fewest arguments to
  defend.
- If `cfdna outliers` is scoped as "detect what LIONHEART would call
  extreme and cap it", A is the natural starting point.
- If `cfdna outliers` is scoped as "detect locally extreme residuals
  and soft-cap them, which is different from global tail extremes",
  C is the right starting point — but then v1 should say explicitly
  that this is the methodological claim being made.

A and B can be ported to C later without contract-breaking changes
(same output format, same application path). C cannot be simplified
back to A/B without retracting a published claim. That asymmetry also
argues for starting with A or B.

### 3.1 If you really want local context, use coarse-bin normalization

The spec's Step 3b — coarse-bin normalization with breakpoint-bin
handling — is actually the better alternative to "local MAD" for this
codebase. Reasons:

- It keeps the scoring in a count-like space (no residual / ratio
  distribution to reason about), so the tail-probability fit stays
  honest.
- It is essentially LIONHEART's megabin step applied earlier, with
  an optional safety valve at breakpoints.
- It avoids inventing a new detector family.

The spec already articulates this well. The review position is that
if any "local" element is kept in v1, it should be this one, not
local-MAD. Local-MAD should be deferred.

## 4. The architectural misalignment the spec flags but does not resolve

LIONHEART clips **counts post-correction** on a per-10-bp-bin basis,
then uses the clipped count array directly. There is no notion of a
"weight track" in their pipeline.

cfDNAlab's scaling-factor contract multiplies a per-bin or per-fragment
weight into downstream accumulation. That translates a clip into one of
the following semantics:

- **Per-bin multiplicative weight** (coverage-style): `w_i = min(1,
  cap_i / c_i)`. Identical to LIONHEART in the positional-coverage
  case. Clean.
- **Per-fragment multiplicative weight** (fragment-averaged commands
  like `midpoints`, `ends`, `lengths`, `bam_to_frag`): a fragment
  spanning bins `i..j` sees the mean of `w_i..w_j`. A narrow clip is
  diluted across ~200 bp. Very different from what LIONHEART does to
  bin counts.

The spec is aware of this (Family A vs Family B vs Family C). Where it
stops short is:

- LIONHEART-style clipping on integer counts was validated on Family-A-
  like usage (Pearson R of bin coverage against chromatin masks). It
  was **not** validated on fragment-averaged workflows.
- Porting the same weights into a fragment-averaged command produces a
  different numerical effect — not incorrect, but different — and the
  downstream implication is that a LIONHEART-compatible weight track
  does **not** give LIONHEART-compatible results in Family B commands
  unless the clip resolution is coarse enough to not dilute.

Concrete consequence for `cfdna outliers`:

- v1 docs should state the clip granularity explicitly and say "this
  matches LIONHEART when applied to positional coverage; for
  fragment-averaged downstream commands, the effective clip weakens
  with fragment length / bin size".
- v1 should resist the urge to auto-tune for fragment-averaged
  commands. That would create a second output track with different
  science.

This aligns with the spec's own "coverage vs fragment-count mode"
Tier 1 ranking, but the framing I would use in v1 docs is tighter:

- `--counting-model coverage` produces weights that track LIONHEART's
  clip semantics.
- `--counting-model fragment-count` produces weights that approximate
  the clip for fragment-count-based support, at a **different** effective
  resolution. It is not "the same clip, applied differently"; it is
  "a different clip informed by the same detector".

Both are defensible. They should not share defaults.

## 5. Answers to the spec's open questions

I'll give an answer for each, with the understanding that these are
opinions informed by the code review, not final decisions.

- **Default scoring bin size?**
  10 bp if matching LIONHEART; otherwise 100–1000 bp for a more neutral
  general-purpose tool. Below 100 bp, fragment-averaged commands barely
  feel the weights.

- **Broad baseline computed internally or supplied externally?**
  Internally for v1. A broad baseline is a property of the detector,
  not of the user contract. Accepting external baselines is future
  scope.

- **Broad smoother or coarse-bin normalization?**
  Coarse-bin normalization (5 Mb / 500 kb, division by mean). That is
  what LIONHEART does, and it has one less knob than a general
  smoother. Breakpoint-bin handling is a Tier 2 feature.

- **ZIP or overdispersion-aware default?**
  ZIP with LIONHEART's simplified estimator. It matches published
  behavior. Adding negative-binomial adds a dispersion parameter that
  is hard to estimate stably under shallow cfDNA without cohort
  borrowing, which v1 does not do.

- **Breakpoint threshold for adjacent coarse-bin mean shift?**
  Defer. Not needed for v1. Most cfDNA samples relevant to this
  command will not have visible CNA-level shifts.

- **Interval width influence on score-to-weight mapping?**
  Defer. LIONHEART does not do this. It is a nice idea but adds a
  tuning axis with no reference point.

- **Positional vs fragment-averaged defaults?**
  Two explicit counting modes. Same detector, different output
  resolution / weights. Docs state they are **not** interchangeable.

- **Continuous scaling only, or greylist sidecar too?**
  Continuous only for v1. Emit the cap metadata in the header so a
  downstream user can threshold themselves if they want a BED.

- **Coverage-only or feature-aware for v1?**
  Coverage-only. Feature-aware scoring multiplies the scope and has
  no cfDNA validation parallel to LIONHEART.

- **Hard mask vs soft cap?**
  Soft cap by default. Hard mask only via explicit flag. LIONHEART's
  clip-not-drop choice is a hard cap at the boundary (`w = cap/c` at
  the boundary, `w = 1` everywhere else). That maps naturally to a
  soft floor on `w` in `[cap/max_count, 1]`.

- **One score mode internally in v1?**
  Yes. Start with the global chromosome-wise tail model; do not
  premature-hybridize.

- **Separate helper for pre-multiplying scaling tracks?**
  Not needed in v1. Pre-multiplication is one file-I/O command that
  users can script. A helper is easy to add later if demand shows.

- **Future `make-blacklist` / `recurrent-outliers` command?**
  Yes, as a separate command, but defer the decision on whether it
  lives here or as a sibling repo. LIONHEART built theirs via shell
  scripts + small Python aggregators; it would fit as a cfDNAlab
  subcommand, but it needs many samples and is a distinct product.

## 6. Smaller notes on the spec

- Section "What existing methods actually do" → "1b. Alignment-guided
  alternatives" is framed as orthogonal to `cfdna outliers`. Agreed, but
  worth adding that T2T + sponges **reduces** the load on this command
  in the first place. The spec says "do not oversell" — one sentence
  saying "and the residual nuisance signal this command targets should
  shrink as references improve" would make the long-term scoping more
  honest.

- Section "LIONHEART" → the ordered list is correct but the order of
  operations at inference differs (exclusions first, then ZIP clip,
  then corrections, then megabin normalize). The spec can state the
  two orderings separately (resource-build vs inference-apply) without
  changing its conclusions.

- Section "Step 4. Compute two scores, not ten" → "a global tail score
  on the broad-normalized support and a local robust residual score on
  the same broad-normalized support" is exactly what LIONHEART does
  **not** do. LIONHEART clips on raw-count support, then broad-
  normalizes. If the spec wants a LIONHEART-aligned v1, this step
  should score on raw-count support and reserve broad-normalized
  scoring for a secondary mode.

- Section "Step 4b. When Poisson or ZIP is scientifically defensible"
  is the best-written section in the spec. My only critique is
  demotional: the three-way distinction is genuinely useful, but v1
  should commit to true-count and explicitly reject the other two for
  this release rather than leaving all three as options.

- Recommended output contract → the `outlier_class` column is a nice
  touch. I would reserve values like `"global_tail"`, `"recurrent"`,
  `"clipped_raw_count"` even if v1 only emits one. Makes later columns
  additive.

- CLI direction → agree with "keep the public surface small". Of the
  listed flags, the v1 minimum is `--counting-model`, `--score-bin-size`,
  `--cap-policy`, and a single tail-threshold knob. Everything else can
  wait.

- Risks and anti-goals → the "do not present the output as truth"
  framing is good. Worth adding one more: **do not imply that matching
  LIONHEART's clip reproduces LIONHEART's downstream performance when
  the downstream command is not positional coverage**. The clip
  semantic does not travel across Family A → Family B unchanged.

## 7. What I would change in the v1 recommendation

Concretely, as a suggested rewrite of the "Best current recommendation"
section, assuming option A or B from Section 3 (the two least-ambitious
shapes). If you pick C, the structure below still applies; only the
detector step changes.

- one sample-specific command
- one explicit counting mode
- single detector step on raw integer-count support:
  - option A: chromosome-wise ZIP fit + strict tail cap (e.g. 1 / N_bins)
  - option B: chromosome-wise empirical upper-quantile cap
- emit weights as `min(1, cap / count)`; merge adjacent identical-weight
  bins
- broad baseline normalization is **not** part of this command; it is a
  downstream smoothing-weight job
- soft by default, hard mask only via explicit flag
- full-coverage scaling TSV output
- easy composition with existing smoothing factors by multiplication
- explicit docs that the clip semantic travels cleanly to positional
  commands and is diluted in fragment-averaged commands

That is a v1 we can defend on its own merits — not because LIONHEART
does the same thing, but because it is one clearly-scoped detector
instead of two entangled ones.

The hybrid-with-local-MAD version (spec's current preferred shape) is
a reasonable v2 or sidecar, **after** v1 has been deployed and we have
a sample of cases where the simpler method visibly under-catches real
residual spikes. Until then, it is a research direction, not a v1
default.

If you prefer to go directly to option C (local MAD residual with broad
baseline), the same "ship one detector, not two" logic applies — in
that case the global-tail score becomes the v2 companion instead of the
other way around.

## 9. Optimal outlier detection for cfDNAlab — reasoning from scratch

Sections 1–7 above take LIONHEART as a reference point and ask
"where should `cfdna outliers` fall relative to that?". That framing
was wrong for this project. This section drops LIONHEART, drops the
"v1" scoping talk, and reasons from the structure of the actual
signal cfDNAlab produces and the questions its downstream commands
actually ask.

### 9.1 Scope: cap extreme upper-tail support, nothing else

This command exists for one job: **detect bins with abnormally high
support and reduce their influence by emitting a multiplicative
weight below 1**. That is the whole product.

What that rules out, and keeps ruled out for the rest of this
section:

- It is not a general anomaly detector. Lower-tail weirdness,
  distributional shifts, pattern breaks — not in scope.
- It is not a CNV caller, a blacklist generator, or a mappability
  tool.
- It does not attempt symmetric handling. Low-support bins get
  weight 1 and move on. A multiplicative weight in `[0, 1]` is
  structurally asymmetric, and the command's job is to match that
  structure, not fight it.

A bin is a candidate for capping if its observed signal is **higher
than what would be predicted** for this locus from known biological
and technical effects, by a margin larger than sampling noise
explains. Two implications:

- The score is on a **residual**, not on the raw value. A principled
  method makes "predicted" explicit (broad baseline, GC effects,
  mappability).
- A cap threshold needs a notion of "sampling noise at the
  predicted level". Without one, "too high" is scale-arbitrary.
  That notion is what the detector choice in Section 9.5 supplies.

Everything downstream follows from this scope. If the scope changes
later — e.g. to include lower-tail capping — the design has to be
reopened, because the current structure is shaped by upper-tail-only.

### 9.2 The "no true count" constraint — and why rescaling solves it

cfDNAlab's signal at the point where this command runs is almost never
literally count-like. By the time a user is thinking about outliers:

- coverage may be GC-corrected
- length normalization may have been applied
- smoothing weights may have been applied
- broad-baseline normalization may have been applied
- different counting models (coverage / fragment-count /
  overlap-weighted) produce different variance structures

My first pass at this section concluded from that list that
parametric count models are out entirely. That was wrong. The spec
already anticipated the correct route (Step 4b, "pseudo-count space
after rescaling and discretization") and I under-weighted it. The
correction matters, because parametric count models give an
interpretable per-bin tail probability in a meaningful unit and a
principled threshold, which robust-scale methods do not.

Most of cfDNAlab's normalizations are **multiplicative rescalings with
known factors**. That fact is what rescues the count interpretation:

- GC correction multiplies each bin by a per-GC-bin factor `f_GC`.
  If the sample and its corrected coverages are bookkept together,
  dividing by `f_GC` returns a count-like quantity per bin.
- `--normalize-by-length` divides by fragment length. Multiplying
  back by a representative (sample-mean) length returns to a
  count-like scale.
- Smoothing with window `W` and mean aggregation multiplies by
  `1/W`. Multiplying back by `W` returns sums which are Poisson(Wλ)
  under the Poisson assumption.
- Broad-baseline division by `m_i` is a per-bin rescaling. Multiplying
  back by the sample mean of `m_i` returns to the original scale.

None of these are lossy in the sense that matters for tail scoring.
They are bookkept rescalings. If `cfdna outliers` is given (or can
reconstruct) the applied rescaling factor per bin, it can invert it
in place, round to integer, and fit a parametric count model honestly.
This is what the spec's "rescale to a count-like mean scale, round,
then fit the global tail model" sentence already says. It is a valid
pipeline, not just a heuristic.

The honest limits of the rescale-back approach:

- The **relative** variance structure within the rescaled signal must
  still be approximately Poisson-like (or NB-like with estimable
  dispersion). For GC correction this is fine when GC factors are
  close to 1; for heavy length-normalization it is fine because the
  representative-length rescale is a global scalar; for smoothing it
  is fine for sums but introduces correlation between neighbors that
  inflates tail rates unless accounted for.
- Correlation from smoothing is the only real subtlety. A fitted
  ZIP on a smoothed-then-rescaled signal will be anticonservative
  (more apparent outliers than true) because neighboring bins are no
  longer independent. The fix is to fit the model on a subsampled set
  of effectively-independent bins (stride ≥ smoothing window), then
  score all bins against the fitted parameters.
- For heavily transformed signals where the rescale factor is not
  tracked (e.g. imported from external tools), rescale-back is not
  available. This is the genuine continuous-residual case.

So the correct taxonomy is:

- **Count-space detector** (Poisson / ZIP / NB, fit on
  rescale-inverted counts): the default when cfDNAlab has the
  rescale bookkeeping to offer. Principled and interpretable.
- **Mean-variance-curve detector** (non-parametric variance
  regression, standardized residuals): the fallback when the signal
  has been through transforms whose rescaling cannot be inverted,
  but non-outlier bins are plentiful enough to estimate `σ²(μ)`.
- **Rank / robust-scale detector** (quantile cutoff, MAD): the
  fallback for small or unusually heterogeneous samples where even
  the variance curve is unstable.

These are ordered by how much structure the detector exploits. Pick
the highest level the signal supports; fall back only when it
doesn't.

This updates my earlier rejection of parametric count models. The
rejection was too strong. The correct position is that the **raw
signal at the scoring step** may not be count-like, but the
**rescale-inverted signal** usually is, and that's where the
parametric model should be fit.

### 9.3 Decompose, score, cap — non-negotiable order

The spec gets this right and then immediately muddies it by letting
the detector re-estimate broad baseline. An optimal design keeps the
three steps clean:

1. **Decompose**: split the signal into a broad component
   `m_i` (CNA / library-depth / arm-level structure, scale ≥ several
   hundred kb) and a local component `r_i = x_i / m_i` (or
   `x_i - m_i`, depending on what stabilizes variance). The broad
   component is estimated internally with a robust method
   (rolling median, optionally with breakpoint-aware handling).
2. **Score**: produce one or more standardized extremeness scores
   on `r_i`. Scores must be scale-equivariant so they survive
   whatever transform space the user passed in.
3. **Cap**: convert scores to multiplicative weights in `[0, 1]`
   via a smooth monotone saturation. Hard zeros only on explicit
   request.

The decomposition step is **part of the detector**, not part of the
downstream smoothing contract. The user's choice of downstream
smoothing weights is independent. That separation is clean only if
the detector commits to owning its broad-baseline estimate.

The log-ratio form `r_i = log((x_i + ε) / (m_i + ε))` has two
properties worth noting:

- It is symmetric around zero, so the same detector machinery
  applies to upper- and lower-tail anomalies without re-derivation.
- It is (approximately) variance-stabilizing when the mean-variance
  relationship is multiplicative, which is the expected regime after
  correction-based rescaling.

The additive form `r_i = x_i - m_i` is better when the residual
variance is close to constant across `m_i`, which is rare for raw
coverage but not uncommon after heavy normalization.

Both forms should be available. The choice is data-dependent and
worth letting the user decide or auto-pick based on a simple
mean-variance diagnostic.

### 9.4 Single-scale is enough; the scale is the user's bin size

An earlier version of this section argued for multi-scale detection
as a requirement. That was overreach. The scope is "cap bins that
are individually extreme enough to dominate downstream statistics",
and at a fixed bin size:

- a single-bin spike at that size crosses threshold and is capped,
- a plateau of moderately-high bins, where each is within normal
  sampling, emits `w = 1` everywhere in the plateau — correctly,
  because capping bins that sit within normal variation for their
  predicted mean is the same error as capping low-count bins
  (Section 9.5b): biasing the signal where the observation is
  statistically normal.

"Plateau of moderate bins that collectively look artifactual" is
not something the per-sample single-BAM detector can call from the
data alone. If you can't distinguish "moderate plateau is artifact"
from "moderate plateau is real biology" at the relevant scale, you
shouldn't cap it. That kind of plateau-level artifact belongs to
static blacklist / mappability resources or to cohort-level
recurrence (Section 9.8), not to the single-sample detector.

The other reason single-scale suffices: **the user picks the bin
size**. It is set by the downstream analysis, not by the detector.
Scoring at the user's bin size matches what the weights will be
applied to. Multi-scale scoring is a tool for problems where the
feature scale is unknown in advance; here it's given.

**Position:** run the detector at the user's chosen bin size.
Only one scoring scale, single score per bin, single cap function.
Multi-scale can be an optional flag for specific workflows (e.g.
a downstream consumer that integrates over several bin widths),
but it is not required by the optimal design.

Practical consequence: if the user's bin size is too fine to
clear the count gate (Section 9.5b) over most of the genome, the
command emits mostly `w = 1`. The correct response for the user
is to pick a coarser bin, not to layer on multi-scale sophistication.

### 9.5 Detector hierarchy: count-space first, then variance-curve, then MAD

Given Section 9.2, the detector at the scoring step has three
options in decreasing order of structural assumption:

**Count-space detector (preferred when applicable).**

1. Invert known multiplicative rescalings to get back to count-like
   scale (per bin, using tracked rescale factors).
2. Round to integer.
3. Fit a count model with zero handling — ZIP with the LIONHEART-
   style simplified estimator is fine; NB with estimated
   overdispersion is stricter and almost always worth the
   additional fit cost.
4. Compute per-bin upper-tail probability `P(X > k)` under the
   fitted model.
5. Threshold: bins below a target tail-probability are outliers.
6. If the input signal was smoothed, fit on a stride-subsampled
   bin set and score on all bins. This avoids the correlation-
   induced anticonservatism noted in Section 9.2.

The per-bin tail probability is directly interpretable as a
per-bin false-positive rate and gives a principled threshold
(e.g. `1/N_bins` for family-wise control, or a looser threshold
for discovery).

Choice between ZIP and NB: ZIP is natural when zero-inflation is
the main departure from Poisson (typical for fine-bin count data
where many bins carry no fragment). NB is natural when
overdispersion is the main departure (typical for coarser bins
or after corrections that inflate residual variance). A detector
that can fit both and picks by a simple likelihood criterion is
not hard to build.

**Mean-variance-curve detector (fallback when rescaling is
unavailable).**

1. Broad-baseline-detrend as in Section 9.3.
2. Estimate `σ²(μ)` non-parametrically from non-outlier residuals.
3. Standardize `z_i = (x_i - m_i) / σ(m_i)`.
4. Threshold against Gaussian or Student-t tails.

This handles heteroskedasticity cleanly but gives less
interpretable thresholds than the count model.

**Robust-scale detector (fallback when bin counts are low).**

1. Rolling MAD of residuals.
2. Standardize and threshold at MAD-units.

Default behavior: count-space when the input signal has tracked
rescaling and enough counts per bin to fit; variance-curve when
rescaling isn't trackable; MAD when the variance curve is too
noisy to estimate. The command can pick automatically based on
sample-level diagnostics, or the user can force one.

### 9.5b Low-depth behavior — handled by the detector, not a gate

Earlier drafts of this section argued for a minimum-expected-count
gate as a hard correctness rule. That was wrong. A properly
specified count-space detector (ZIP / NB / Poisson) handles
low-λ bins correctly via the tail test itself:

- The threshold is in probability units, not count units.
- At λ = 0.3, observing k = 3 is not extreme (tail ≈ 3e-3) and
  does not cross a strict threshold like 1/N_bins. No cap.
- At λ = 0.3, observing k = 8 is extreme (tail ≈ 1e-9) and does
  cross. Cap is warranted — the observation is genuinely
  inconsistent with the fitted model at this locus.

So low expected count does not imply "cannot detect". It implies
"the discrete tail has fewer relevant values, and the model sets
the bar where it should be".

Where low-depth still matters — and how to address it honestly:

- **MAD and variance-curve detectors do degenerate at low counts.**
  This is a property of those specific fallbacks. The fix is to
  prefer the count-space detector in that regime (Section 9.5),
  not to refuse to score.
- **Shallow whole-sample behavior.** The count model fits fine on
  shallow samples; the tail threshold is still meaningful. The
  command does not need to refuse shallow samples.
- **Essentially-zero regions** (λ close to 0 everywhere). Any
  non-zero observation crosses any strict threshold. This is an
  input-filtering concern — these regions are almost always either
  blacklist / unmappable / always-zero bins that the user's
  prefilter should have removed — not a detector concern.

Per-sample diagnostic log should still report sample-wide depth,
per-chromosome fit quality (e.g. residual variance vs Poisson
variance), and number of bins capped. That lets a user catch
model-misspecification regimes without the detector silently
refusing to score.

### 9.6 One counting-model mode, not two bundled tracks

An earlier draft recommended emitting two weight tracks per run —
one for coverage-geometry commands and one for fragment-geometry
commands — so each downstream family would get appropriately
shaped weights.

That was wrong. The two tracks require two input signals (per-base
coverage vs fragment-count-per-bin) and two full detector runs, so
emitting both from one invocation is ~2× the work of running the
command once. The user paying for both makes sense only if they use
both, and for most analyses they don't. Users who do need both can
run the command twice at the same total cost, with control over
when they pay.

Correct position: one `--counting-model` flag, one detection run
per invocation, one weight track out. The track is explicitly
tagged with the counting model it was built for. Users pick the
counting model that matches their downstream command family.

The cap shape is the same regardless of counting model — a smooth
saturation, not a hard clip:

```
w_i = 1                                 if Z_i <= τ
w_i = (τ / Z_i)^γ                       if Z_i > τ
```

- `τ` sets the detection threshold.
- `γ` sets how quickly the weight decays; `γ = 1` is a cap-at-τ
  semantic, `γ > 1` suppresses more aggressively, `γ < 1` more
  gently.

The floor is `τ / max_Z` in the sample, automatically adapting to
how extreme the worst observation is.

### 9.7 Interval width asymmetry and the score-to-weight map

Narrow and wide detected intervals should not get the same cap
shape. A 1-bin `z = 10` spike and a 100-bin plateau at `z = 10`
carry different evidence — the plateau is almost certainly
systematic, the spike might be a sampling accident.

Concrete mechanism:

- After scoring, identify contiguous runs of bins above threshold.
- For each run of width `L` bins, compute an evidence score
  `E_L = Z_mean * sqrt(L)` (the "sqrt(L)" reflects that independent
  samples of Z accumulate evidence like a mean, up to the effective
  independence ratio, which for smoothed signal is below the naive
  sqrt but still monotone in L).
- Map `E_L` to the cap:
  - Low `E_L` (narrow, barely-over-threshold) → gentle cap.
  - High `E_L` (wide, strongly-over-threshold) → aggressive cap.

Result: the command naturally does what an experienced analyst
would do — cautious with single-bin anomalies, firm with plateaus.

### 9.8 Cohort information is a prior, not a filter — but one cap, not two

The spec separates "recurrent blacklist" and "per-sample detector"
into two products. That is correct. What's missing is how they
**interact** at detection time, and specifically how to avoid
double-downweighting.

The central rule: **there must be exactly one cap function per
emitted weight track, and it consumes both sample and cohort
evidence as inputs to a single score.** Stacking caps is what
produces double-downweighting.

Wrong architecture (double-downweights):

- Cap the sample signal to get `w_sample_i`.
- Separately cap based on cohort recurrence to get `w_cohort_i`.
- Emit `w_combined_i = w_sample_i * w_cohort_i`.

This double-counts because each cap function fires independently,
and a bin that looks mildly extreme from both sources ends up
harshly capped even when neither source alone would cap it much.
Multiplication of weights is fine for composing **orthogonal**
corrections (broad-baseline smoothing × outlier capping), but
cohort and sample evidence about the same question are not
orthogonal — they are two pieces of evidence for the same claim.

Right architecture (single cap, dual input):

- Compute the sample-specific score `Z_i` per Sections 9.4–9.5.
- Compute a cohort-derived per-bin prior `π_i ∈ [0, ∞)`, built
  from held-out cohort data, representing how frequently bins like
  this were flagged in prior samples (with the count-gate from
  Section 9.5b applied during construction).
- Combine **in score space, before capping**:
  `Z'_i = combine(Z_i, π_i)`.
- Apply exactly one cap function to `Z'_i` → `w_i`.
- Emit `w_i`. There is no separate cohort weight track.

What `combine(·)` should be:

- Threshold-adjustment form: keep `Z'_i = Z_i` but shift the
  per-bin threshold `τ_i = τ / π_i`. A recurrent-looking bin
  (`π_i > 1`) gets a lower bar for suppression; a
  always-clean-looking bin (`π_i < 1`) needs stronger sample
  evidence.
- Score-addition form: `Z'_i = Z_i + α * log(π_i)` with `α`
  controlling how much weight the prior carries relative to the
  sample score. Mathematically close to the threshold-adjustment
  form for upper-tail thresholding.
- What `combine` should **not** be: `w_combined = w_sample *
  w_cohort`. That is the double-downweighting case.

Two additional protections against double-counting:

- The cohort prior must be built from samples **not including the
  sample currently being scored**, and ideally not from close
  biological or technical replicates of it. Otherwise the "prior"
  carries evidence that is not independent of the sample's own
  score.
- The prior must apply the same count-gate (Section 9.5b): only
  bins where the relevant prior sample cleared the gate contribute
  to `π_i`. A bin that was "never flagged as outlier" across the
  cohort because it was below the count gate in every sample is
  not a clean bin; it is an undetectable bin. Its `π_i` should be
  marked as uninformative (e.g. `π_i = 1` with a flag), not 0.

Consequences for the cohort-resource workflow:

- The cohort resource is a distribution of per-sample scores per
  bin (not a binary mask). That distribution lets downstream
  users derive either a prior `π_i` or a simpler blacklist
  cutoff, as needed.
- When the cohort resource is applied to a new sample, it enters
  only as input to the single cap function on `Z'_i`. Users who
  want to also apply a separately-maintained hard blacklist can
  do so, but the result is an exclusion (weight = 0 or = 1 per
  user choice), not another multiplicative soft cap layered on
  top.

The principled test for whether two downweighting sources are
independent enough to multiply: ask whether they are answering
different questions. Broad-baseline smoothing answers "what is
the normal expected signal at this locus?" and outlier capping
answers "does this sample's observation deviate from expected?"
Those are different questions, and the answers multiply cleanly.
Cohort recurrence and sample extremeness are both answers to
"is this bin an outlier?", and those should not multiply — they
should combine as co-evidence in a single score.

### 9.9 Upper-tail only — no lower-tail handling

Since the command's job is capping high values, the detector and
the scoring are upper-tail only:

- The residual score is signed, but only positive residuals
  (observed > predicted) can produce `w < 1`.
- Negative residuals (observed < predicted) emit `w = 1`. No
  sidecar. No annotation. They are simply out of scope.

This is the minimum surface that matches the output contract.
Expanding to symmetric handling would require a different output
representation (e.g. a `[0, 2]` weight or a separate depletion
track) and is not this command.

### 9.10 The concrete optimal pipeline

Putting the above together:

1. **Input**: per-chromosome signal track (any counting model, any
   prior normalization), plus an optional cohort prior `π`, plus
   optional mappability/blacklist exclusion indices.
2. **Prefilter**: mask excluded indices; they do not enter
   statistics and are emitted as `w = 1` (passthrough) or `w = 0`
   (hard mask) per user preference.
3. **Decompose**: robust broad baseline `m_i` via rolling median
   over a broad window (scale chosen to exceed expected CNA block
   size but stay well below chromosome length). Breakpoint-bin
   handling activates only when adjacent median shifts exceed a
   calibrated noise threshold.
4. **Residualize**: log-ratio `r_i = log((x_i + ε)/(m_i + ε))`,
   with ε chosen relative to sample mean.
5. **Score**: at the user's chosen bin size, pick the detector per
   Section 9.5 — count-space if rescalings are trackable
   (rescale-inverted counts fed to ZIP / NB / Poisson; the tail
   test naturally handles low-λ bins), mean-variance curve
   otherwise, MAD as last resort. Produce per-bin upper-tail
   statistic `Z_i` (or equivalently, tail probability).
6. **Cohort combination — in score space, before any cap**
   (Section 9.8): if a cohort prior `π_i` is supplied, combine via
   `Z'_i = Z_i + α log(π_i)` (or equivalently via a threshold
   adjustment `τ_i = τ / π_i`). Do not produce a separate cohort
   weight track.
7. **Run detection**: contiguous runs of `Z'_i > τ` → interval
   list with `(width L, mean Z within run)` per interval.
8. **Single evidence-weighted cap**: per run, compute `E_L`; a
   single cap function converts `(Z', L)` → `w`. This is the only
   place a weight below 1 is emitted. No cap is applied elsewhere
   in the pipeline.
9. **Weight emission**: full-coverage TSV with run-length merging
   of equal weights. Columns include `chromosome, start, end,
   scaling_factor, class, Z, L, m, x, fit_space`.
10. **Counting-model tag**: the emitted track is tagged with the
    counting model it was built from (coverage or fragment-count).
    Users who need a track in the other counting model run the
    command again with that flag.
11. **Upper-tail only**: only positive residuals produce `w < 1`.
    Negative residuals emit `w = 1`.
12. **Diagnostics**: per-sample log of (a) detector family used,
    (b) fit-quality summary (residual vs modeled variance),
    (c) number of runs detected, (d) weight distribution summary.

Defaults for the named parameters can be locked later by
calibration on a panel of clean and contaminated cfDNA samples.
The shape of the pipeline is what matters.

### 9.11 What this design explicitly rejects

- Parametric count models fit on the **raw post-correction signal
  without inverting rescalings**. The input at that stage is not
  count-like and a ZIP or Poisson fit there is a real
  misspecification. The fix is to invert tracked rescalings and fit
  on the count-like inverted signal, per Section 9.2 — not to
  abandon parametric models.
- Implicit or shared broad-baseline. The detector owns its
  broad-baseline estimate; downstream smoothing weights remain
  independent.
- Hard thresholds as the primary cap mechanism. Smooth saturation
  preserves information about how extreme each outlier is.
- Cohort information as exclusion, **or as a second multiplicative
  weight track layered on the sample cap**. Cohort information
  combines with sample evidence in score space under a single cap
  function; it must not produce its own independent weight factor
  (Section 9.8).
- Lower-tail handling of any kind. The command's job is capping
  high values; low-support bins get `w = 1` and are not the
  command's concern.

### 9.12 Relationship to the existing spec

The spec anticipated some of this — the decomposition principle, the
multi-scale intuition in the coarse-bin-normalization section, the
score-space discussion, the asymmetry note. Where this section goes
further:

- It keeps parametric count models as the preferred detector when
  rescalings can be inverted, but routes through explicit
  bookkeeping rather than the spec's looser "discretization" framing.
  Falls back to variance-curve or MAD when rescale inversion is
  unavailable.
- It introduces the evidence-weighted cap as a principled way to
  handle interval-width asymmetry, which the spec flags but does
  not resolve.
- It introduces cohort-prior-as-multiplier rather than
  cohort-as-exclusion, which is structurally different from the
  spec's "two separate products" framing.

The design still composes with the existing scaling-factor
infrastructure unchanged: the output is one (or two) full-coverage
TSVs of multiplicative weights in `[0, 1]`, most runs at exactly 1,
short runs softly capped below 1.

## 10. References actually consulted (beyond the spec's list)

- `plans_and_specs/lionheart-main/lionheart/features/correction/poisson.py`
- `plans_and_specs/lionheart-main/lionheart/features/correction/normalize_megabins.py`
- `plans_and_specs/lionheart-main/lionheart/features/create_dataset_inference.py`
- `plans_and_specs/lionheart-main/lionheart/commands/extract_features.py`
- `plans_and_specs/lionheart-main/resource_creation/detect_outlier_candidates.sh`
- `plans_and_specs/lionheart-main/resource_creation/calculate_tail_cdf.py`
- `plans_and_specs/lionheart-main/resource_creation/collect_outliers.py`
- `plans_and_specs/lionheart-main/resource_creation/collect_outliers_across_datasets.py`
- `plans_and_specs/lionheart-main/tests/test_poisson.py`
