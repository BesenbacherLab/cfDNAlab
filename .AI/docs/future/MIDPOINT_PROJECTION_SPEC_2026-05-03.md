# Midpoint projection spec

Date: 2026-05-03

Scope: This note describes a possible opt-in midpoint definition for `cfdna
midpoints` that accounts for soft clipping and indels while still counting only
known aligned reference positions.

This is a design note, not an implementation decision. The main purpose is to
make the biological and coordinate assumptions explicit before discussing the
idea further.

## Goal

The current midpoint concept is reference-span based:

```text
paired fragment span = [forward.pos, reverse.reference_end)
midpoint = center of that aligned reference span
```

That is stable and easy to explain, but it ignores the fact that the observed
molecule can contain clipped bases, inserted bases, and deleted reference bases.

The adjusted concept considered here is:

```text
Adjusted midpoint = the midpoint of the observed molecule, projected back onto
the aligned reference span.
```

The counted coordinate must remain inside the aligned fragment coordinates. For
paired-end fragments, that means `forward.pos <= midpoint < reverse.reference_end`.
For unpaired read-as-fragment mode, that means `read.pos <= midpoint < read.reference_end`.

The adjusted midpoint is not allowed to become an outside-reference coordinate
from soft clipping. Soft clips can move the molecule midpoint, but the output
coordinate must still be a known aligned reference position. If projection
cannot produce a coordinate inside the aligned fragment span, skip the fragment
and count it in a separate projection-skipped statistics counter.

## Opt-in controls

Soft clipping and indels should be separate opt-ins, matching the direction of
the other commands.

Proposed behavior:

- `--clip-mode aligned`: ignore soft clips for midpoint placement.

- `--clip-mode skip`: skip fragments where the relevant fragment end is
  soft-clipped.

- `--clip-mode adjust`: include 5-prime soft clips as molecule bases that can
  move the projected midpoint.

- `--indel-mode ignore`: ignore indels for midpoint placement.

- `--indel-mode skip`: skip fragments with insertions or deletions.

- `--indel-mode adjust`: include indels when constructing the molecule axis and
  projecting the midpoint to reference coordinates.

These modes should combine independently. For example, a user should be able to
adjust clips while ignoring indels, adjust indels while using aligned clip
boundaries, or skip either class of fragment evidence.

## Mental model

Think of the fragment as two coordinate systems:

- Reference axis: positions in the genome.

- Molecule axis: bases in the observed fragment molecule.

Different CIGAR events consume these axes differently:

| Event | Consumes reference | Consumes molecule | Effect on projection |
| --- | --- | --- | --- |
| Match/equal/diff | yes | yes | ordinary aligned bases |
| Insertion | no | yes | molecule bases with no unique reference coordinate |
| Deletion/ref-skip | yes | no | reference bases absent from the molecule |
| Soft clip | no aligned coordinate | yes | molecule bases outside the aligned span |

The adjusted midpoint should be found on the molecule axis, then projected to an
aligned reference coordinate.

This is better than calculating an aligned midpoint and applying a single net
shift. A net shift loses where each event occurred. An insertion before the
midpoint and an insertion after the midpoint have opposite effects, and a total
insertion count alone cannot tell those cases apart.

## Directional intuition

The projection model gives the expected directional behavior:

- Left soft clipping tends to move the projected midpoint left.

- Right soft clipping tends to move the projected midpoint right.

- Insertions before the molecule midpoint tend to move the projected reference
  coordinate left, because extra molecule bases have already been consumed.

- Insertions after the molecule midpoint tend to move the projected reference
  coordinate right less often or not at all, depending on where the midpoint
  falls.

- Deletions before the molecule midpoint tend to move the projected reference
  coordinate right, because reference positions were crossed without consuming
  molecule bases.

- Deletions after the molecule midpoint usually do not affect the projected
  coordinate.

This intuition is useful for checking examples, but the implementation should be
a coordinate projection, not a set of ad hoc left and right shifts.

## Projection logic

For one fragment, construct a reference-ordered list of projection events over
the aligned fragment span. The walker tracks:

- current reference position

- current molecule offset

- target molecule midpoint offset

Soft clips contribute molecule bases before the first aligned reference position
or after the last aligned reference position, depending on which fragment end is
clipped and whether `--clip-mode adjust` is active.

Indels contribute according to `--indel-mode`:

- Insertions advance molecule offset but not reference position.

- Deletions advance reference position but not molecule offset.

- Matches/equal/diff advance both.

The target midpoint is selected from the adjusted molecule length. Even-length
fragments should keep the existing deterministic random tie-break idea, but the
tie-break should be based on the adjusted molecule coordinates rather than the
unadjusted reference span.

After selecting the target molecule offset, walk the projection until the target
is reached:

- If the target lands on a reference-consuming aligned base, count that reference
  position.

- If the target lands inside a soft-clipped prefix or suffix, count the nearest
  aligned reference coordinate inside the original aligned span.

- If the target lands inside an insertion, count the insertion anchor by a
  documented rule.

- If a deletion crosses the aligned midpoint, do not treat that alone as a
  reason to shift from the aligned midpoint. Deletions have no molecule bases.
  The walker should cross deleted reference positions and place the midpoint at
  the projected reference coordinate implied by the molecule midpoint. If that
  projection is ambiguous at a deletion boundary, use the documented deletion
  projection policy.

The final projected coordinate must be validated:

```text
fragment_start <= projected_midpoint < fragment_end
```

If the projected coordinate falls outside this interval, skip the fragment and
increment a separate counter for fragments whose adjusted midpoint projected
outside the aligned span. Do not clamp these fragments back into the aligned
span.

## Insertion fallback

An insertion has molecule bases but no unique reference coordinate. If the
adjusted midpoint lands inside an insertion, there are three plausible policies:

- Count the insertion anchor.

- Count the nearest flanking aligned reference base.

- Skip the fragment as ambiguous.

The most usable first policy is probably insertion-anchor projection, because it
keeps the fragment countable and makes the behavior deterministic. The drawback
is local pileup risk at insertion anchors in indel-rich regions.

Skipping is more conservative for positional profiles, but it changes the
fragment inclusion set and can make `--indel-mode adjust` behave closer to
`--indel-mode skip` around exactly the cases it was meant to handle.

Falling back to the original aligned midpoint is stable but less principled. If
earlier clips or indels already moved the molecule midpoint, returning to the
unadjusted midpoint can discard most of the opt-in adjustment.

## Deletion behavior

Deletions consume reference but not molecule bases. Because of that, the molecule
midpoint cannot truly land on a deleted base in molecule coordinates. There is
no observed molecule base there. But the projection back to reference can still
be ambiguous, because a deletion creates a jump across reference positions
without advancing the molecule offset.

Example:

```text
reference:  ... 98 99 [100 101 102 deleted] 103 104 ...
molecule:   ... 98 99                      103 104 ...
```

If the molecule midpoint is before the deletion, the projected reference
coordinate is before the deletion. If it is after the deletion, the projected
reference coordinate is after the deletion. The ambiguous case is when the
midpoint falls exactly at the jump. At that point, several reference-coordinate
policies are possible:

- Count the left flanking aligned base.

- Count the right flanking aligned base.

- Count the flank closest to the original aligned midpoint.

- Count a deleted reference coordinate.

- Skip the fragment as ambiguous.

Counting a deleted reference coordinate is allowed in a purely reference-space
sense, because the coordinate is inside known reference territory. It is weaker
biologically, because the molecule did not contain that base. This is different
from insertion-anchor projection, where the molecule contains bases but the
reference has no unique coordinate for them. Both cases are projections, but the
missing object is opposite.

This becomes less clear when multiple events offset each other. For example, a
large insertion before the midpoint can pull the projected reference coordinate
left, while a deletion elsewhere can push it right. At that point, the adjusted
midpoint is no longer "the original aligned midpoint plus a local correction".
It is a coordinate chosen from a warped molecule-to-reference projection. In
that setting, forbidding deleted reference coordinates is not a mathematical
requirement. It is a biological and reporting policy: do we require the counted
midpoint base to be molecule-backed, or is it acceptable to count a known
reference coordinate that the molecule skipped over during projection?

The most defensible first policy is probably to count the nearest flanking
aligned base, with a deterministic tie-break if both flanks are equally plausible.
That keeps the count on a molecule-backed aligned base and avoids placing
midpoint mass on reference bases absent from the molecule. If the intended
interpretation is stricter, skipping deletion-boundary ambiguities is cleaner
than counting inside the deletion.

If a deletion overlaps the old aligned midpoint, the adjustment should not
special-case that as "shift left" or "shift right". The correct question is
where the molecule midpoint projects after accounting for all events before it.

For documentation, this should be described as:

```text
Deleted reference bases are crossed during projection because they are not
present in the molecule. If the midpoint projection is ambiguous at the deletion
jump, use the documented deletion-boundary policy.
```

## Overlapping reads with differing CIGAR strings

Paired-end fragments can have an aligned mate-overlap. In that overlap, the two
reads may disagree about insertions or deletions. This matters more for midpoint
projection than for total adjusted length, because the event position can change
the projected coordinate.

The `lengths` command uses a conservative molecule-leaning policy:

- In non-overlap, count indels observed by the single covering read.

- In mate-overlap, count deletions only where the deletion intervals are
  supported by both reads.

- In mate-overlap, count insertions only where both reads have an insertion at
  the same reference anchor, using the shorter insertion length if lengths
  differ.

The midpoint projection should likely reuse the same support policy, but it
needs position-preserving events rather than only total insertion and deletion
counts.

Important implication:

```text
An adjusted midpoint implementation cannot reuse only the existing total indel
counts. It needs the resolved fragment-level insertion anchors and deletion
intervals after the overlap policy has been applied.
```

Differing CIGAR strings should be reported through counters, at minimum:

- fragments with any indel evidence

- fragments with overlap-disagreed indel evidence ignored by the support policy

- fragments where the adjusted midpoint landed inside an insertion fallback

- fragments skipped because projection could not produce an aligned coordinate
  inside the fragment's aligned span

These counters are needed because local profile changes could otherwise be
mistaken for biology when they are driven by alignment disagreement.

## Biological interpretation

The adjusted midpoint should be described as a reference-projected molecule
midpoint. That phrase matters.

It is not simply:

- the aligned reference midpoint

- the midpoint of an adjusted scalar length

- a midpoint outside the reference span

It is an attempt to answer:

```text
Where would the physical molecule midpoint land if we map the observed molecule
back onto known aligned reference positions?
```

This can better reflect molecule geometry when soft clips or indels represent
real molecule sequence. It can also make profiles worse when these CIGAR events
mostly reflect mapping artifacts, local repeats, poor alignment, or sample and
aligner-specific behavior.

## Expected profile effects

Most fragments will likely keep the same midpoint if they have no relevant
clipping or indels.

Adjusted fragments can move by one or more bases. The largest effects should
come from:

- asymmetric 5-prime soft clipping

- insertions or deletions before the molecule midpoint

- indel-heavy regions

- repetitive or low-mappability regions where CIGAR strings are less stable

Genome-wide average profiles may change little if adjusted fragments are rare.
Local profiles around indels, repeats, or problematic alignments can change
noticeably. That is exactly why the feature should remain opt-in and report
diagnostic counters.

## Drawbacks

- The adjusted midpoint is more complex to explain than the current aligned
  midpoint.

- Projection requires position-preserving indel events, not only adjusted length
  totals.

- Insertions have no unique reference coordinate, so any counted coordinate is a
  policy choice.

- Soft-clipped bases may be real molecule sequence, but they may also be
  alignment noise. Adjusting for them can move midpoint mass for the wrong
  reason.

- Mate-overlap CIGAR disagreement can create command-specific behavior unless
  the overlap support policy is shared and documented.

- GC correction, blacklist filtering, scaling weights, and tile fetch geometry
  probably still need to use aligned reference spans unless explicitly redesigned.
  That creates a split where midpoint placement is adjusted but some filters and
  weights remain aligned-coordinate based.

## Suggested first implementation boundary

Do not start by changing default midpoint semantics.

First implementation should probably add opt-in support with:

- independent clip and indel modes

- adjusted midpoint projection only when requested

- final-coordinate validation inside the aligned fragment span

- skip and count fragments whose adjusted midpoint projects outside the aligned
  fragment span

- deterministic insertion fallback

- counters for fallback and skipped projection cases

- tests for left/right soft clipping, insertion before and after midpoint,
  deletion before and after midpoint, insertion-at-midpoint fallback, deletion
  crossing the old aligned midpoint, and overlapping mates with agreeing and
  disagreeing CIGAR events

Open decision before implementation:

```text
Should insertion-at-midpoint count the insertion anchor, the nearest flanking
aligned base, or skip the fragment?

Should deletion-boundary ambiguity count the left flank, the right flank, the
flank closest to the original aligned midpoint, a deleted reference coordinate,
or skip the fragment?
```

Those choices should be made deliberately because they control the ambiguous
cases in reference-projected midpoint adjustment.
