# cfDNAlab | Rust Loaders <img src="https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png" align="right" height="155" />

Rust library helpers for loading files already written by `cfdna`.

These APIs open command outputs from disk and return typed Rust metadata, count containers, and selector builders. Use them when downstream Rust code needs to inspect or reuse cfDNAlab outputs without parsing TSV or Zarr files manually.

The loaders live under `cfdnalab::output_loaders` and are compiled with the matching command feature.

| Cargo feature   | Loader                    | Output                                                   |
| --------------- | ------------------------- | -------------------------------------------------------- |
| `cmd_midpoints` | `load_midpoints_output()` | `<prefix>.midpoint_profiles.zarr`                        |
| `cmd_ends`      | `load_ends_output()`      | `<prefix>.end_motifs.zarr`                               |
| `cmd_lengths`   | `load_lengths_output()`   | `<prefix>.length_counts.tsv`, optionally `.gz` or `.zst` |
| `cmd_fcoverage` | `load_fcoverage_output()` | non-positional aggregate `fcoverage` TSV outputs         |
| `cmd_ref_kmers` | `load_ref_kmers_output()` | `<prefix>.ref_kmer_counts.zarr`                          |

<br>

## Common Pattern

Each loader returns an output object with metadata and values. Start by checking what the file contains, then either read all values or pick the rows, groups, motifs, length bins, or positions you need.

Selectors keep the order you request and reject duplicate selectors on the same axis. Range selectors use half-open intervals with `Interval::new(start, end)?`.

```rust
use cfdnalab::{
    interval::Interval,
    output_loaders::load_lengths_output,
};

fn main() -> anyhow::Result<()> {
    // Load output file and check the available metadata
    let lengths = load_lengths_output("sample.length_counts.tsv.zst")?;
    println!("{}", lengths.output_metadata());

    // Select fragment lengths from 100-220bp
    let selected = lengths
        .select()
        .length_range(Interval::new(100, 221)?)
        .read()?;

    // Print the selected counts by row and length bin
    for (row_index, row_counts) in selected.counts().rows().enumerate() {
        for (length_bin, count) in selected.length_bins().iter().zip(row_counts) {
            println!(
                "{row_index}\t{}-{}\t{count}",
                length_bin.start(),
                length_bin.end()
            );
        }
    }

    Ok(())
}
```

<br>

## Midpoint Profiles

Midpoint profile stores contain a 3D count array with axes:

```text
group x length_bin x position
```

Opening the store reads metadata. Counts stay on disk until you call `read_all_counts()` or `select().read()`.

```rust
use cfdnalab::{
    interval::Interval,
    output_loaders::load_midpoints_output,
};

fn main() -> anyhow::Result<()> {
    // Load output file and check the available metadata
    let midpoints = load_midpoints_output("sample.midpoint_profiles.zarr")?;
    println!("{}", midpoints.output_metadata());

    // Check the axis sizes
    let groups = midpoints.group_metadata();
    let length_bins = midpoints.length_bins();
    let positions = midpoints.position_bins();

    println!(
        "{} groups, {} length bins, {} position bins",
        groups.len(),
        length_bins.len(),
        positions.len()
    );

    // Select groups LYL1 and GATA1 and fragment lengths from 100-220bp
    let selected = midpoints
        .select()
        .groups_by_name(&["LYL1", "GATA1"])
        .length_range(Interval::new(100, 221)?)
        .read()?;

    // Print a profile total for each selected group and length bin
    for (selected_group_index, group) in selected.groups().iter().enumerate() {
        for (selected_length_index, length_bin) in selected.length_bins().iter().enumerate() {
            let profile = selected
                .profile(selected_group_index, selected_length_index)
                .expect("selected profile indices should be in bounds");
            let profile_total = profile.iter().copied().sum::<f32>();

            println!("{}\t{:?}\t{profile_total}", group.name, length_bin.as_tuple());
        }
    }

    Ok(())
}
```

After `midpoints.select()`, use `groups()` for selecting group-axis indices, `groups_by_name()` for group labels, `length_bins()` for specific length-bin indices, `length_range()` for all length bins overlapping a half-open fragment length range, `positions()` for position-bin indices, and `position_range()` for interval-relative position ranges. On the loaded output, use `group_metadata()`, `length_bins()`, and `position_bins()` to inspect the full axes before selecting.

<br>

## End-Motif Counts

End-motif stores can be dense or sparse. Check `storage_mode()` before choosing how to access counts.

```rust
use cfdnalab::output_loaders::{
    EndMotifStorageMode,
    load_ends_output,
};

fn main() -> anyhow::Result<()> {
    // Load output file and check the available metadata
    let ends = load_ends_output("sample.end_motifs.zarr")?;
    println!("{}", ends.output_metadata());

    // Check row and motif metadata
    println!("{:?}", ends.row_mode());
    println!("{:?}", ends.motif_axis_kind());
    println!("{:?}", ends.motif_labels());

    // Select motif _AA
    let selected = ends
        .select()
        .motifs_by_label(&["_AA"])
        .read()?;

    // Print the selected counts
    match selected.storage_mode() {
        EndMotifStorageMode::Dense => {
            let counts = selected.dense_counts()?;
            for row in counts.rows() {
                println!("{:?}", row);
            }
        }
        EndMotifStorageMode::SparseCoo => {
            for entry in selected.sparse_counts()?.entries() {
                println!(
                    "{}\t{}\t{}",
                    entry.row_index, entry.motif_index, entry.count
                );
            }
        }
    }

    Ok(())
}
```

Windowed outputs provide `window_metadata()`, and you can select window rows with `select().windows(...)`. Grouped outputs provide `group_metadata()` and `group_index()`, and you can select grouped rows with `select().groups(...)` or `select().groups_by_name(...)`.

Sparse stores keep missing in-bounds cells as implicit zero counts. Use `sparse_counts()?.to_lookup_index()` for repeated random access or `to_dense_matrix()` only when the selected matrix is small enough to hold in memory.

Reference correction is available when the `cmd_ends` and `cmd_ref_kmers` features are both enabled. Load a reference k-mer store for the same k-mer size, windowing or grouping, motif settings, and reference genome, then pass it to `select_corrected_counts()`.

```rust
use cfdnalab::output_loaders::{
    load_ends_output,
    load_ref_kmers_output,
    UnsupportedReferencePolicy,
};

fn main() -> anyhow::Result<()> {
    let ends = load_ends_output("sample.end_motifs.zarr")?;
    let ref_kmers = load_ref_kmers_output("hg38.ref_kmers.zarr")?;

    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .motifs_by_label(&["_AA", "_GG"])
        .unsupported_reference_policy(UnsupportedReferencePolicy::KeepNaN)
        .read()?;

    let corrected_counts = corrected.to_dense_matrix()?;
    for row in corrected_counts.rows() {
        println!("{row:?}");
    }

    Ok(())
}
```

Reference correction divides each observed end-motif count by a reference-based correction factor for the matched row. This factor is computed from the motif frequencies in the reference k-mer output and normalized so a uniform reference composition leaves counts unchanged. Motifs that are common in the reference row are scaled down. Motifs that are rare in the reference row are scaled up. Only motifs with a positive reference frequency contribute to the row's correction support.

When motif labels contain both outside and inside bases, such as `AC_TG`, call `.two_sided_correction(...)` and choose how the two sides should be handled:

- `TwoSidedCorrectionMode::Joint` keeps full labels such as `AC_TG` and corrects each count with the matching reference k-mer, `ACTG`.
- `TwoSidedCorrectionMode::Split` keeps full labels such as `AC_TG`, but calculates the correction factor from the two sides separately. For `AC_TG`, separate correction factors are calculated for outside label `AC` and inside label `TG`. Those two correction factors are multiplied and applied to the observed `AC_TG` count. Use this when full two-sided motif labels should remain in the result, but the exact full reference k-mers are too sparse or the correction should treat outside and inside sequence composition separately.
- `TwoSidedCorrectionMode::Outside` returns outside labels such as `AC_`. For each outside label, all full motif counts with that outside label are summed first. For example, `AC_AA` and `AC_TG` both contribute to the `AC_` count. That summed count is corrected using the outside label `AC`.
- `TwoSidedCorrectionMode::Inside` returns inside labels such as `_TG`. For each inside label, all full motif counts with that inside label are summed first. For example, `AA_TG` and `AC_TG` both contribute to the `_TG` count. That summed count is corrected using the inside label `TG`.

For `Outside` and `Inside`, repeated side labels are deduplicated in their first loaded-motif occurrence order. The returned `EndMotifCountSelection::motif_labels()` and `motif_indices()` describe this corrected side axis, so use them to interpret matrix columns.

One-sided outputs do not accept an explicit mode.

Motif labels are matched to reference k-mers by removing `_`, for example `AT_CG` -> `ATCG`. Motif-group outputs are matched by group label. Both commands write forward-oriented motif labels, including right-end motifs from `cfdna ends`.

For `Split`, `Outside`, and `Inside`, side-specific reference frequencies are calculated from the loaded full-length reference k-mers. For example, the outside frequency for `AC` is the sum of frequencies for loaded k-mers with prefix `AC`, such as `ACTG` and `ACAA`. The inside frequency for `TG` is the corresponding sum over loaded k-mers with suffix `TG`. Separate shorter reference k-mer runs are not required.

A motifs file used for the reference output restricts these sums to the k-mers in that file. Without a motifs file, all k-mers in the reference output can contribute, including k-mers absent from the sample end-motif output.

By default, end-motif and reference k-mer rows must match exactly. A global reference k-mer store can be applied to every windowed or grouped end-motif row only when `.use_global_bias(true)` is set. That option requires a global reference store and is unnecessary when both outputs are global. Sample-observed motifs can be absent from the reference genome or have zero reference frequency in a row. Positive end-motif counts for those motifs are errors by default. Use `UnsupportedReferencePolicy::KeepNaN` to keep the selected shape and mark those cells as `NaN`.

<br>

## Length Counts

Length-count TSV outputs are loaded into row metadata, length-bin metadata, and a dense `DenseMatrix<f64>`. Rows can be `global`, windows, or groups.

```rust
use cfdnalab::{
    interval::Interval,
    output_loaders::{LengthOutputMode, load_lengths_output},
};

fn main() -> anyhow::Result<()> {
    // Load output file and check the available metadata
    let lengths = load_lengths_output("sample.length_counts.tsv.zst")?;
    println!("{}", lengths.output_metadata());

    // Check whether rows are global, windows, or groups
    match lengths.row_mode() {
        LengthOutputMode::Global => println!("global length counts"),
        LengthOutputMode::Windows => {
            for window in lengths.window_metadata()? {
                println!("{}\t{:?}", window.chrom, window.interval.as_tuple());
            }
        }
        LengthOutputMode::Groups => {
            for group in lengths.group_metadata()? {
                println!("{}\t{}", group.name, group.eligible_windows);
            }
        }
    }

    // Select fragment lengths from 100-220bp
    let selected = lengths
        .select()
        .length_range(Interval::new(100, 221)?)
        .read()?;

    // Print a total for each selected row
    for row_counts in selected.counts().rows() {
        let selected_total = row_counts.iter().copied().sum::<f64>();
        println!("{selected_total}");
    }

    Ok(())
}
```

Use `length_bin_for_length()` when you have a fragment length in bp and want the matching length-bin index. Use `length_bins_overlapping_range()` or `select().length_range()` when you want all whole length bins overlapping a half-open bp range.

Windowed outputs support `select().windows(&[...])`. Grouped outputs support `select().groups(&[...])` and `select().groups_by_name(&[...])`.

<br>

## fcoverage Aggregates

The `fcoverage` loader supports non-positional aggregate TSV outputs from `average`, `total`, and `summary_stats` modes.

**Note**: Positional bedGraph and per-window positional TSV outputs are intentionally out of scope.

```rust
use cfdnalab::output_loaders::{
    load_fcoverage_output,
    load_fcoverage_output_with_group_index,
};

fn main() -> anyhow::Result<()> {
    // Load windowed output and check the available metadata
    let windowed = load_fcoverage_output("sample.fcoverage.average.tsv.zst")?;
    println!("{}", windowed.output_metadata());

    // Print a value for each window
    for (window, value) in windowed
        .window_metadata()?
        .iter()
        .zip(windowed.values()?.iter().copied())
    {
        println!(
            "{}:{}-{}\t{value}",
            window.chrom,
            window.interval.start(),
            window.interval.end()
        );
    }

    // Load grouped output with the matching group index
    let grouped = load_fcoverage_output_with_group_index(
        "sample.fcoverage.total_on_unique_bases.tsv.zst",
        "sample.group_index.tsv",
    )?;

    // Select promoter and enhancer groups
    let selected = grouped
        .select()
        .groups_by_name(&["promoters", "enhancers"])
        .read()?;

    // Print the selected values by group
    for (group, value) in selected
        .group_metadata()?
        .iter()
        .zip(selected.values()?.iter().copied())
    {
        let group_name = group.name.as_deref().unwrap_or("<unnamed>");
        println!("{group_name}\t{value}");
    }

    Ok(())
}
```

The selected fcoverage rows are returned as `FCoverageSelection`, so handle scalar-value and summary-stat outputs separately.

```rust
use cfdnalab::output_loaders::{
    FCoverageSelection,
    load_fcoverage_output_with_group_index,
};

fn main() -> anyhow::Result<()> {
    // Load grouped summary stats output with the matching group index
    let fcoverage = load_fcoverage_output_with_group_index(
        "sample.fcoverage.summary_stats.tsv.zst",
        "sample.group_index.tsv",
    )?;

    // Select promoter group
    let selected = fcoverage
        .select()
        .groups_by_name(&["promoters"])
        .read()?;

    // Print scalar values or summary averages
    match selected {
        FCoverageSelection::Values(values) => {
            for value in values.values() {
                println!("{value}");
            }
        }
        FCoverageSelection::SummaryStats(stats) => {
            for stats_row in stats.stats() {
                println!("{}", stats_row.average);
            }
        }
    }

    Ok(())
}
```

Grouped fcoverage TSV files store numeric `group_idx` values. Use `load_fcoverage_output_with_group_index()` with the matching `group_index.tsv` when you want group names and `groups_by_name()` selection.

<br>

## Reference K-mer Frequencies

Reference k-mer stores can be dense or sparse. They store row-wise frequencies plus a row scaling factor that reconstructs counts.

For `ref-kmers` outputs written with `--motifs-file`, frequencies are normalized over the selected motifs or motif groups from that file. Unlisted k-mers are not part of the denominator, and the row scaling factor reconstructs selected k-mer or group counts.

With `--all-motifs`, the motif axis also keeps targets whose stored frequency is zero. Without a motifs file, those targets are all A/C/G/T k-mers for the configured `k`. With a motifs file, they are the motifs or motif groups listed in that file.

```rust
use cfdnalab::output_loaders::{
    load_ref_kmers_output,
    RefKmerStorageMode,
};

fn main() -> anyhow::Result<()> {
    // Load output file and check the available metadata
    let ref_kmers = load_ref_kmers_output("hg38.ref_kmer_counts.zarr")?;
    ref_kmers.ensure_reference_2bit_matches("hg38.2bit")?;
    println!("{}", ref_kmers.output_metadata());

    // Reconstruct a count from one frequency and its row scaling factor
    let motif_index = ref_kmers.motif_index("ACGT")?;
    let count = ref_kmers
        .count(0, motif_index)
        .expect("row and motif indices should be in bounds");
    println!("{count}");

    // Select grouped rows and motifs, then work with reconstructed counts
    let selected = ref_kmers
        .select()
        .groups_by_name(&["promoters", "enhancers"])
        .motifs_by_label(&["ACGT", "TGCA"])
        .read()?;
    let selected_counts = selected.to_dense_count_matrix()?;
    for (group, row_counts) in selected.group_metadata()?.iter().zip(selected_counts.rows()) {
        let count_total = row_counts.iter().copied().sum::<f64>();
        println!("{}\t{count_total}", group.name);
    }

    // Read the stored frequency data in its native mode
    match ref_kmers.storage_mode() {
        RefKmerStorageMode::Dense => {
            let frequencies = ref_kmers.dense_frequencies()?;
            println!("{:?}", frequencies.shape());
        }
        RefKmerStorageMode::SparseCoo => {
            let sparse_frequencies = ref_kmers.sparse_frequencies()?;
            for entry in sparse_frequencies.entries() {
                println!(
                    "{}\t{}\t{}",
                    entry.row_index, entry.motif_index, entry.frequency
                );
            }
            for entry in ref_kmers.sparse_count_entries()? {
                println!("{}\t{}\t{}", entry.row_index, entry.motif_index, entry.count);
            }
        }
    }

    Ok(())
}
```

Use `select().windows(...)`, `select().groups(...)`, or `select().groups_by_name(...)` for row subsets. Use `select().motifs(...)` or `select().motifs_by_label(...)` for motif subsets. Use `frequency()` or `frequency_for_motif()` when downstream code wants frequencies. Use `count()`, `count_for_motif()`, `sparse_count_entries()`, or `to_dense_count_matrix()` when it wants reconstructed counts.

<br>

## Error Handling

Loader methods return `OutputLoaderResult<T>`. Errors include path, line, array, or selector context where possible, so most applications can use `?` and add their own higher-level context at the call site.

```rust
use anyhow::Context;
use cfdnalab::output_loaders::load_lengths_output;

fn main() -> anyhow::Result<()> {
    // Add context to the load error
    let lengths = load_lengths_output("sample.length_counts.tsv.zst")
        .context("load length-count output for downstream analysis")?;

    // Check the count matrix shape
    println!("{:?}", lengths.counts().shape());
    Ok(())
}
```
