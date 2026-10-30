**Title:** End-motif commands ignore `--min-length` and `--max-length`

On `main` at `3ca99fd3`, `end-motifs` and `interval-end-motifs` accept length bounds, but [`region_end_motifs()`](https://github.com/epifluidlab/FinaleToolkit/blob/3ca99fd3fdc63e22b69919ddf63e0d79426af089/src/finaletoolkit/frag/_end_motifs.py#L51-L180) counts directly from `AlignmentWrapper.fetch()` without applying them. The bounds reach the worker as `fraction_low`/`fraction_high`, but neither the worker nor the reader checks fragment lengths against them.

Could the inclusive length checks be added before the strand-specific counting branches?

Source inspection only, without runtime reproduction.
