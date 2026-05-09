# `cfdna coverage-weights`

Build genomic scaling factors that normalize large-scale coverage variation. The command measures average fragment coverage in stride bins, smooths it across larger bins, and writes multiplicative factors for downstream coverage-like features.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    setup["Define smoothing grid<br/>stride and large-bin size"]
    source["Run internal fcoverage<br/>fixed stride windows"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    filters["Keep usable fragments<br/>MAPQ, length, pairing, blacklist, GC rules"]
    coverage["Measure stride coverage<br/>average coverage per interval"]
    smooth["Smooth across large bins<br/>triangular overlap kernel"]
    normalize["Normalize genome-wide<br/>mean supported coverage becomes 1"]
    invert["Invert to multipliers<br/>low-coverage regions get higher factors"]
    output["Scaling-factor TSV<br/>coverage.scaling_factors.tsv"]
    consumer["Used by feature commands<br/>coverage-smoothing weights"]

    region_inputs["Optional region inputs<br/>chromosome selection and blacklists"]
    correction_inputs["Optional correction inputs<br/>GC package, GC tag, and 2bit reference"]
    run_inputs["Optional run controls<br/>read mode, MAPQ, fragment lengths, ignore gap"]

    input --> setup --> source --> fragments --> filters --> coverage --> smooth --> normalize --> invert --> output --> consumer

    region_inputs -.-> setup
    region_inputs -.-> filters
    correction_inputs -.-> filters
    correction_inputs -.-> coverage
    run_inputs -.-> fragments
    run_inputs -.-> filters
    run_inputs -.-> coverage

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,setup,source,fragments,filters,coverage,smooth,normalize,invert core;
    class region_inputs,correction_inputs,run_inputs optional;
    class output,consumer outputClass;
```

## Coverage Model

`coverage-weights` uses the same fragment creation and filtering behavior as `fcoverage`. In paired-end mode, coverage is based on the fragment span from the forward read position to the reverse read reference end. In `--reads-are-fragments` mode, each accepted read is counted as its own fragment.

By default, paired fragments cover the full inter-mate span. With `--ignore-gap`, uncovered sequence between non-overlapping mates is excluded, matching downstream `fcoverage --ignore-gap` analyses.

## Smoothing Model

The command first measures average coverage in fixed stride bins. It then applies a triangular overlap kernel derived from `--bin-size` and `--stride`, so each stride row reflects neighboring large-bin support.

After smoothing, supported rows are normalized to a genome-wide mean of 1 and inverted. Downstream commands multiply coverage-like features by these factors to reduce large-scale coverage variation.

## Output

The output is `<prefix>.coverage.scaling_factors.tsv`, or `coverage.scaling_factors.tsv` when no prefix is set. The TSV includes stride coordinates, raw stride average coverage, smoothed coverage, the multiplicative scaling factor, and metadata describing whether GC correction and inter-mate gap exclusion were used while building the weights.
