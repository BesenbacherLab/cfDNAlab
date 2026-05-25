#[cfg(not(feature = "cmd_fragment_kmers"))]
use crate::shared::kmers::kmer_codec::KmerSpec;
#[cfg(feature = "cmd_fragment_kmers")]
use crate::shared::{
    base::make_canonical,
    kmers::kmer_codec::{Kmer, KmerSpec},
};
use fxhash::FxHashMap;
#[cfg(feature = "cmd_fragment_kmers")]
use fxhash::FxHashSet;

/// Per-k map of “reference” counts
#[derive(Debug, Clone, PartialEq)]
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) struct DecodedCounts {
    pub(crate) counts: FxHashMap<u8, FxHashMap<String, f64>>, // k -> motif -> count
}

/// Prepare decoded counts for all kmer sizes in all windows.
///
/// Extracts motifs per kmer spec to allow future padding.
/// For kmers of size 1..6, this includes all possible motifs.
/// For larger kmer sizes, only the seen motifs is included as the number otherwise explodes.
///
/// * `windows`        – slice of per-window raw counts
/// * `canonical`      – canonical reverse complements when true
/// * `kmer_specs`     – validated specs for every k we want to keep
///
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn prepare_decoded_counts(
    windows: &[DecodedCounts],
    canonical: bool,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> (Vec<DecodedCounts>, FxHashMap<u8, Vec<String>>) {
    let n_windows = windows.len();

    // Initialise one empty DecodedCounts per window
    let mut out = vec![
        DecodedCounts {
            counts: FxHashMap::default()
        };
        n_windows
    ];

    let mut motifs_by_k: FxHashMap<u8, Vec<String>> = FxHashMap::default();

    // Loop over every k we validated
    for &k in kmer_specs.keys() {
        // Reference (match) bins for this k
        let (count_bins, motifs) =
            prepare_kmer_category(windows, kmer_specs, k as usize, canonical, k <= 6);

        // Insert into the corresponding window
        for i in 0..n_windows {
            out[i].counts.insert(k, count_bins[i].clone());
        }
        motifs_by_k.insert(k, motifs);
    }

    (out, motifs_by_k)
}

#[cfg(feature = "cmd_fragment_kmers")]
fn prepare_kmer_category(
    windows: &[DecodedCounts],
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    k: usize,
    canonical: bool,
    ensure_all: bool,
) -> (Vec<FxHashMap<String, f64>>, Vec<String>) {
    // Extract the raw maps
    let raw_bins = extract_bins(windows, k, canonical);

    // Build the (canonical) motif list once, if requested.
    let base_motifs: Vec<String> = if ensure_all {
        all_motifs(k, kmer_specs)
    } else {
        Vec::new()
    };

    // Build the (canonical) motif list *once* so we know what to pad with
    let mut motifs = collect_motifs(&raw_bins, base_motifs, canonical, ensure_all);
    motifs.sort_unstable();

    (raw_bins, motifs)
}

/// Collect per-window bins for the requested motif type and (optionally)
/// canonical them into strand-agnostic form.
///
/// * `windows` – slice of `DecodedCounts` (“one window” each).
/// * `k` – kmer-size to pull out of every `DecodedCounts`.
/// * `canonical` – if `true`, run the appropriate collapse_*_map helper.
///
/// Returns a fresh `Vec<FxHashMap<String, f64>>` – one map per window.
#[cfg(feature = "cmd_fragment_kmers")]
fn extract_bins(
    windows: &[DecodedCounts],
    k: usize, // pattern only; field values are ignored
    canonical: bool,
) -> Vec<FxHashMap<String, f64>> {
    windows
        .iter()
        .map(|dc| {
            // Pick the raw map for this window
            let raw: FxHashMap<String, f64> =
                dc.counts.get(&(k as u8)).cloned().unwrap_or_default();

            // Collapse if requested, otherwise return the raw map
            if canonical { collapse_map(raw) } else { raw }
        })
        .collect()
}

/// Collect motifs for a category, optionally ensuring the full universe and filtering 'N'
#[cfg(feature = "cmd_fragment_kmers")]
fn collect_motifs(
    windows: &[FxHashMap<String, f64>],
    base_motifs: Vec<String>,
    canonical: bool,
    ensure_all: bool,
) -> Vec<String> {
    // Universe of motifs to keep
    let set: FxHashSet<String> = if ensure_all {
        base_motifs.into_iter().collect()
    } else {
        windows.iter().flat_map(|m| m.keys().cloned()).collect()
    };

    // Strand-collapse if requested
    let collapsed_set = if canonical { collapse_set(&set) } else { set };

    // Convert to sorted Vec
    let mut v: Vec<String> = collapsed_set.into_iter().collect();
    v.sort_unstable();
    v
}

/// Collapse a map of kmer counts into canonical keys, summing counts
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn collapse_map(map: FxHashMap<String, f64>) -> FxHashMap<String, f64> {
    let mut out: FxHashMap<String, f64> = FxHashMap::default();
    out.reserve(map.len());

    for (kmer, count) in map {
        let canon = make_canonical(kmer, true, false);
        *out.entry(canon).or_insert(0.0) += count;
    }

    out
}

/// Collapse a set of reference kmers into a set of canonical keys
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn collapse_set(set: &FxHashSet<String>) -> FxHashSet<String> {
    set.iter()
        .map(|m| make_canonical(m.to_string(), true, false))
        .collect()
}

/// Return all possible reference motifs (4ᵏ) for a given k.
///
/// No motifs with 'N' are returned.
pub(crate) fn all_motifs(k: usize, specs: &FxHashMap<u8, KmerSpec>) -> Vec<String> {
    let spec = &specs[&(k as u8)];
    let max_code = 5u64.pow(k as u32) - 1; // no-N space
    (0..=max_code)
        .map(|c| spec.decode_kmer(c))
        .filter(|m| !m.contains('N'))
        .collect()
}

/// Split an aggregated `counts` map into per-k buckets.
///
/// * The `kmer_specs` dict tells us which k-values are valid and how to decode.
/// * Motifs that contain 'n' are discarded.
///
/// Returns one map for reference windows (“matches”) and one for mismatches.
#[cfg(feature = "cmd_fragment_kmers")]
pub(crate) fn split_and_decode_counts(
    counts: &FxHashMap<Kmer, f64>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> DecodedCounts {
    let mut count_bins: FxHashMap<u8, FxHashMap<String, f64>> = FxHashMap::default();

    for (&kmer, &cnt) in counts {
        // User-readable motif, e.g. "ACG"
        let motif = kmer.to_string(kmer_specs);

        // Drop N's
        if motif.contains('N') {
            continue;
        }

        let bucket = count_bins.entry(kmer.k).or_default();
        *bucket.entry(motif).or_insert(0.0) += cnt;
    }

    DecodedCounts { counts: count_bins }
}
