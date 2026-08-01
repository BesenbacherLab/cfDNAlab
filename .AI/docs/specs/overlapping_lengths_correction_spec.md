# Average Overlapping Fragment Length Normalization Model Spec

`cfdna overlap-length-model` fits a sample-wide model of the relationship between positional
coverage and the average fragment length among fragments overlapping each eligible base. Its
purpose is to homogenize this relationship before comparing coverage across samples. The method
ports LIONHEART's two skewed Student-t mixture fits and is considered validated primarily for
100-220 bp fragments. It writes `<prefix>.overlap_length_model.zarr` for application by
`fcoverage --overlap-length-file`.

## Counting invariants

- A paired fragment length is `forward.pos` to `reverse.reference_end`. An unpaired fragment is
  `read.pos` to `read.reference_end`.
- Raw integer fragment depth drives the initial mixture's `1 / sqrt(depth)` spread. The refit and
  target use raw depth after bin-wise noise and skew division, `coverage > 0.5` filtering, and
  ties-to-even rounding, matching LIONHEART. GC weighting and genomic scaling affect only the
  observed coverage signal.
- A covered segment receives the fragment's full length from `forward.pos` to
  `reverse.reference_end` in the length-sum prefix, including when deletions, skipped regions, or
  the inter-mate gap are omitted from coverage.
- Blacklisted positions do not contribute to model sufficient statistics or an inferred fragment
  weight.
- Tiles aggregate only their cores. Per-bin signal sums/counts, the raw-depth histogram, and joint
  length-bin/depth counts with actual average-length sums are sufficient for both fits, so training
  performs one BAM sweep.

## Model and application

- The package stores all fitted curves, division factors, combined multiplicative weights, fit
  parameters, and compatibility settings. It stores no input paths or file fingerprints.
- `fcoverage` calculates a positional lookup weight from raw overlap depth and length sum, then
  assigns a fragment the base-pair-weighted mean across its original counted segments.
- The model scalar multiplies existing length-normalization and GC fragment weights.
  Genomic scaling remains positional after coverage accumulation.
- An `ignore_gap` mismatch is an error. Fragment-filter, pairing, trimming, and length-normalization
  differences warn and continue. GC/scaling/blacklist identity is documented rather than checked.
- Model-enabled fcoverage tiles use a two-maximum-fragment-length halo. Inference tracks the
  largest returned fragment start and finalizes positions before that start minus the maximum
  fragment length. It retains bounded 64 KiB prefix chunks until all queued fragments that
  reference them have been weighted, and returns fragments in exactly the inner iterator's order.
