# `cfdna fragment-count-weights`

Build genomic scaling factors that normalize large-scale fragment-count variation. The command counts unit fragment mass in stride bins, smooths that mass across larger bins, and writes multiplicative factors for downstream count-like features.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    setup["Define smoothing grid<br/>stride and large-bin size"]
    source["Run internal fcoverage<br/>fixed stride windows"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    filters["Keep usable fragments<br/>MAPQ, length, pairing, blacklist, GC rules"]
    unit["Assign unit mass<br/>one fragment contributes about 1 total"]
    stride["Count stride bins<br/>fragment mass per genomic interval"]
    smooth["Smooth across large bins<br/>triangular overlap kernel"]
    normalize["Normalize genome-wide<br/>mean supported mass becomes 1"]
    invert["Invert to multipliers<br/>low-mass regions get higher factors"]
    output["Scaling-factor TSV<br/>fragment_counts.scaling_factors.tsv"]
    consumer["Used by feature commands<br/>fragment-count smoothing weights"]

    region_inputs["Optional region inputs<br/>chromosome selection and blacklists"]
    correction_inputs["Optional correction inputs<br/>GC package, GC tag, and 2bit reference"]
    run_inputs["Optional run controls<br/>read mode, MAPQ, fragment lengths"]

    input --> setup --> source --> fragments --> filters --> unit --> stride --> smooth --> normalize --> invert --> output --> consumer

    region_inputs -.-> setup
    region_inputs -.-> filters
    correction_inputs -.-> filters
    correction_inputs -.-> unit
    run_inputs -.-> fragments
    run_inputs -.-> filters

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,setup,source,fragments,filters,unit,stride,smooth,normalize,invert core;
    class region_inputs,correction_inputs,run_inputs optional;
    class output,consumer outputClass;
```

## Fragment-Mass Model

`fragment-count-weights` uses the same fragment creation and filtering behavior as `fcoverage`, but it counts in unit-mass mode. A fragment contributes approximately one total unit of mass, split across the stride bins covered by its span.

This differs from `coverage-weights`, where long fragments naturally contribute more total coverage because they cover more bases.

## Smoothing Model

The command first counts fragment mass in fixed stride bins. It then applies a triangular overlap kernel derived from `--bin-size` and `--stride`, so each stride row represents smoothed support from neighboring large bins.

After smoothing, supported rows are normalized to a genome-wide mean of 1 and inverted. Downstream commands multiply fragment counts by these factors to reduce large-scale count variation.

## Output

The output is `<prefix>.fragment_counts.scaling_factors.tsv`, or `fragment_counts.scaling_factors.tsv` when no prefix is set. The TSV includes stride coordinates, raw stride fragment mass, smoothed fragment mass, the multiplicative scaling factor, and metadata describing whether GC correction was used while building the weights.
