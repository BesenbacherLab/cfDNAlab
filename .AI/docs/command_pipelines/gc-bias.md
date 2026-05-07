# `cfdna gc-bias`

Fit a sample-specific GC correction package from a BAM file. The command measures the sample's observed fragment GC distribution, compares it with a reusable reference GC package, and writes correction factors for downstream feature commands.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    reference["Reference inputs<br/>2bit genome and ref GC package"]
    setup["Load correction frame<br/>length range, GC bins, smoothing settings"]
    windows["Choose contribution windows<br/>global, fixed bins, or BED windows"]
    tiles["Scan BAM in tiles<br/>with fragment-length halos"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    gc["Measure fragment GC<br/>from the reference sequence"]
    assign["Assign fragments to windows<br/>overlap, containment, midpoint, or proportion"]
    normalize["Normalize each window<br/>mean-scale by usable ACGT"]
    merge["Merge window evidence<br/>including cross-tile windows"]
    observed["Observed sample table<br/>fragment length by GC percentage"]
    compare["Compare with reference<br/>sample bias relative to expected GC"]
    stabilize["Stabilize correction factors<br/>binning, interpolation, outlier handling"]
    package["GC correction package<br/>gc_bias_correction.npz"]
    consumer["Used by feature commands<br/>GC-aware fragment weights"]

    region_inputs["Optional region inputs<br/>BED windows and blacklists"]
    run_inputs["Optional run controls<br/>read mode, MAPQ, binning, outliers"]
    qc["Optional QC outputs<br/>intermediate arrays and plots"]

    input --> setup
    reference --> setup --> windows --> tiles --> fragments --> gc --> assign --> normalize --> merge --> observed --> compare --> stabilize --> package --> consumer

    region_inputs -.-> windows
    region_inputs -.-> gc
    run_inputs -.-> fragments
    run_inputs -.-> stabilize
    stabilize -.-> qc

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,reference,setup,windows,tiles,fragments,gc,assign,normalize,merge,observed,compare,stabilize core;
    class region_inputs,run_inputs,qc optional;
    class package,consumer outputClass;
```

## Correction Model

The reference GC package defines the expected fragment distribution for a genome setup: fragment lengths by GC percentage, plus the support masks and smoothing choices created by `cfdna ref-gc-bias`.

`gc-bias` measures the same length-by-GC table from the sample BAM. Fragments are counted in genomic windows so local coverage variation can be normalized before the command estimates the global sample GC profile.

## Window Model

With `--global`, all accepted fragments contribute to one sample-wide table. With fixed-size or BED windows, each window is normalized separately and then merged into the final observed table. Window assignment can count overlap proportionally, require any or all overlap, use the fragment midpoint, or require a minimum fragment-overlap proportion.

Blacklists remove problematic reference bases before fragment GC and window-normalization support are computed.

## Output

The main output is `<prefix>.gc_bias_correction.npz`, or `gc_bias_correction.npz` when no prefix is set. The package contains the correction matrix and metadata needed by downstream commands that apply GC-aware fragment weights. When requested, the command also writes intermediate arrays and QC plots for inspecting the fitted correction surface.
