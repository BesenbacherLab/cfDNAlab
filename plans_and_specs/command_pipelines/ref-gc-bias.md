# `cfdna ref-gc-bias`

Build a reusable reference GC package for a genome assembly. The command samples reference start positions, counts GC content for each configured fragment length, and writes the expected GC-by-length matrix used by `cfdna gc-bias`.

## Pipeline

```mermaid
flowchart TD
    input["Reference genome<br/>2bit file"]
    setup["Define reference model<br/>chromosomes, length range, end offset"]
    regions["Select usable genome<br/>all bases or merged BED regions"]
    sample["Sample start positions<br/>same starts reused across lengths"]
    tiles["Scan reference in tiles<br/>with max-length sequence halos"]
    mask["Mask excluded bases<br/>blacklists and ambiguous sequence"]
    count["Count reference fragments<br/>GC by fragment length"]
    merge["Merge tile counts<br/>one genome-wide reference table"]
    smooth["Smooth raw GC counts<br/>optional length-row kernel"]
    percent["Convert to GC percentage<br/>with bin-width correction"]
    support["Build support masks<br/>impossible bins and sparse bins"]
    interpolate["Interpolate sparse bins<br/>optional per-length filling"]
    package["Reference GC package<br/>ref_gc_package.npz"]
    consumer["Used by gc-bias<br/>to fit sample-specific correction"]

    region_inputs["Optional region inputs<br/>BED include regions and blacklists"]
    run_inputs["Optional run controls<br/>sampling target, seed, threads"]

    input --> setup --> regions --> sample --> tiles --> mask --> count --> merge --> smooth --> percent --> support --> interpolate --> package --> consumer

    region_inputs -.-> regions
    region_inputs -.-> mask
    run_inputs -.-> sample
    run_inputs -.-> tiles

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,setup,regions,sample,tiles,mask,count,merge,smooth,percent,support,interpolate core;
    class region_inputs,run_inputs optional;
    class package,consumer outputClass;
```

## Reference Model

The command samples candidate start positions from the selected chromosomes and evaluates every configured fragment length from those starts. With `--by-bed`, overlapping and touching BED intervals are merged first, and a reference fragment contributes only when its full span fits inside the selected region.

`--end-offset` trims bases from both fragment ends before GC is counted. This lets the reference model focus on the fragment interior when end-proximal sequence composition should not drive GC correction.

## Support Model

Raw GC counts are converted into integer GC-percentage columns from 0 to 100. Because different raw GC counts can round into the same percentage bin, the command corrects for bin width before writing the package.

The package also stores support masks. One mask marks GC-percentage bins that are theoretically impossible for a length. The other marks bins with too little empirical reference support. By default, sparse bins are filled by interpolation so downstream correction can avoid avoidable zero-count gaps.

## Output

The output is `<prefix>.ref_gc_package.npz`, or `ref_gc_package.npz` when no prefix is set. This package is the reference-side input to `cfdna gc-bias` and can be reused for samples aligned to the same reference setup.
