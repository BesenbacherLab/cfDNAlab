# Outlier correction: what should count as safe enough?

2026-10-19. Proposal for discussion with the supervisor. This is not a new real-data result.

## What we want it to do

A small region with an extreme number of DNA fragments can dominate a coverage profile. We want
to reduce that technical excess while retaining useful biological differences. The command finds
unusually high raw coverage and assigns weights to overlapping fragments. Its preservation rule
leaves broad or clustered calls unchanged, but it cannot determine their cause from coverage alone.

The decision is therefore a tradeoff. Correcting more can remove useful signal. Preserving more
can leave technical artifacts. An acceptable tradeoff depends on the analysis using the coverage.

## Three possible harms

**Ordinary coverage beside a call can be reduced.** A fragment receives the full weight even if
only a base overlaps a call. Its contribution outside the call also falls. This follows directly
from the implemented rule. We have a hand-derived regression fixture for it. We have not measured
the resulting distortion of real biological profiles. The removed coverage need not form a simple
V shape, because fragment lengths, alignment gaps and overlapping calls affect the result.

**A genuine increase in DNA copies can be flattened.** Extra copies of a genomic region can
produce higher coverage. A short event, or scattered high bases inside a wider event, may be
corrected even when the event is biological. The preservation percentage measures called bases,
not the event's full width. Calling an event "preserved" also does not protect every fragment
inside it from a neighboring actionable call.

**GC correction may interact with the result.** Sequencing can represent DNA with different
proportions of G and C unevenly. GC weights and outlier weights multiply during counting. However,
the files used to estimate GC bias are not yet produced from the same outlier-weighted fragments.
Whether residual pileups materially distort those estimates is an open question, not an observed
failure of the current workflow.

## What evidence exists

The code establishes the overlap and preservation rules. New regression tests have been authored
and compiled, but have not been run by Codex. This update includes no new real-BAM measurements.

Historical simulations compared ZIP and ZINB on the same generated fragments. ZINB allows more
ordinary coverage variation than ZIP. It avoided ordinary-background correction in those tested
cases and retained several broad gains. It also
removed about 58% of biological signal in a tested 10 kb gain. That is a serious counterexample
to assuming biological safety. Those simulations used an older 10% preservation cutoff and a
different ZINB optimizer. They do not measure the current 5% rule or validate the Rust fitter.
See the [recorded methods and outcomes](outliers_historical_simulation_evidence.md).

## Measurements to agree on

Use matched raw and weighted coverage with identical filters and blacklists. Start without other
corrections, then repeat with GC correction. Include ordinary regions, independently supported
copy gains, likely technical pileups, and spikes inside gains. Coverage alone cannot supply the
truth labels. When independent evidence is absent, label the cause unknown.

For each group, measure technical excess removed, biological signal retained, coverage removed
outside exact calls, and changes in the profiles or features actually used downstream. Report
difficult individual examples as well as averages. Repeat across libraries and the intended depth
range. The proposed 0.1x to 200x range remains unvalidated. Record nominal sequencing depth and
counted fragment coverage separately.

Candidate acceptance criteria should cover:

- A named minimum event width and strength that must remain detectable.
- An agreed maximum biological signal loss, including outside calls and in difficult examples.
- An agreed maximum change in relevant downstream profiles or measurements.
- Enough technical-excess removal to justify those losses, including moderate artifacts.
- Explicit sample types, depths and settings for which the evidence supports use.

There is no defensible universal percentage to choose from the current evidence. Fitting a
distribution well, producing few called bases, or having a clean global plot does not establish
that the correction preserves the measurements we care about.

## Questions for the discussion

1. Which downstream measurement must this correction improve, and what would harmful distortion look like?
2. What is the smallest biological gain we must retain, and what independent evidence can identify it?
3. How much loss is acceptable inside that event and in neighboring ordinary regions?
4. Which sample types and depths must be supported before using this in the intended analysis?
