//! Motif extraction and reference-lookup helpers for the `ends` command.
//!
//! This module keeps the tile-local reference context, end-level motif
//! encoding, and per-window counting logic out of the top-level runner so
//! `ends.rs` can focus on orchestration.

use crate::{
    commands::ends::{
        config::EndsConfig,
        config_structs::{ClipStrategy, KmerSource, WindowMotifAssigner},
        counting::{EncodedEndMotifKey, EndCountsByWindow},
    },
    shared::{
        blacklist::{apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE},
        fragment::ends_fragment::{FragmentWithEnds, ResolvedFragmentEnd},
        interval::Interval,
        kmers::kmer_codec::{
            KmerCodes, KmerSpec, build_kmer_specs, build_left_aligned_codes_per_k,
        },
        reference::read_seq_in_range,
        tiled_run::Tile,
    },
};
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use std::sync::Arc;

/// Reference-backed motif resources for one tile.
///
/// This groups the per-tile state needed to validate and encode end motifs:
/// optional masked reference bases, optional radix-5 lookup tables for the
/// inside and outside halves, and the metadata needed to translate absolute
/// genomic motif starts into the preloaded tile-local reference slice.
pub(crate) struct TileMotifContext<'a> {
    /// Absolute genomic start of `reference_bases`
    reference_start: u64,
    /// Tile reference bases, already blacklist-masked when needed
    reference_bases: Option<Vec<u8>>,
    /// Spec for the inside half, if `k_inside > 0`
    inside_spec: Option<KmerSpec>,
    /// Spec for the outside half, if `k_outside > 0`
    outside_spec: Option<KmerSpec>,
    /// Precomputed masked-reference codes for inside lookups
    inside_codes: Option<Arc<KmerCodes>>,
    /// Precomputed masked-reference codes for outside lookups
    outside_codes: Option<Arc<KmerCodes>>,
    /// Blacklist intervals used to decide whether read-backed inside validation is needed
    blacklist_intervals: &'a [Interval<u64>],
    /// Chromosome length used for `sentinel_none` checks
    chrom_len: u64,
}

/// Identify which fragment end is being processed.
///
/// The left and right ends use different genomic lookup starts for both the
/// inside and outside halves, so this enum keeps that branching explicit.
#[derive(Clone, Copy)]
pub(crate) enum EndSide {
    Left,
    Right,
}

/// Track which fragment ends produced at least one motif count.
///
/// The same end can be counted into multiple windows, but the statistics for
/// `ends` should still count that end only once per fragment. This struct keeps
/// that per-fragment bookkeeping separate from the per-window accumulation map.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct CountedEndFlags {
    pub left_counted: bool,
    pub right_counted: bool,
}

impl CountedEndFlags {
    /// Merge count flags from another window of the same fragment.
    ///
    /// Parameters
    /// ----------
    /// - `other`:
    ///   Count flags collected from another candidate window for the same fragment
    pub(crate) fn merge(&mut self, other: Self) {
        self.left_counted |= other.left_counted;
        self.right_counted |= other.right_counted;
    }

    /// Return whether at least one end motif was counted for the fragment.
    ///
    /// Returns
    /// -------
    /// - `bool`:
    ///   `true` when either fragment end contributed at least one motif count
    pub(crate) fn any_counted(self) -> bool {
        self.left_counted || self.right_counted
    }

    /// Return how many distinct end motifs were counted for the fragment.
    ///
    /// Returns
    /// -------
    /// - `u64`:
    ///   Number of fragment ends that contributed at least one motif count
    pub(crate) fn counted_motif_total(self) -> u64 {
        self.left_counted as u64 + self.right_counted as u64
    }
}

/// Build an optional k-mer spec, treating `k=0` as an empty motif half.
///
/// Parameters
/// ----------
/// - `k`:
///   Requested k-mer size for one motif half
/// - `label`:
///   Human-readable side name used in error messages
///
/// Returns
/// -------
/// - `Result<Option<KmerSpec>>`:
///   `None` when `k=0`, otherwise the radix-5 codec spec for that half
pub(crate) fn build_optional_kmer_spec(k: usize, label: &str) -> Result<Option<KmerSpec>> {
    if k == 0 {
        return Ok(None);
    }

    let k_u8: u8 = k
        .try_into()
        .with_context(|| format!("{label} k-mer size {k} does not fit in u8"))?;
    let mut specs = build_kmer_specs(&[k_u8])?;
    Ok(specs.remove(&k_u8))
}

/// Compute the preloaded motif-reference span for one tile.
///
/// BAM fetch narrowing is independent of motif reference preload. For motif
/// lookups, `ends` needs the full tile fetch band plus any extra sequence that
/// raw clipping and outside-of-fragment k-mers may inspect beyond that band.
pub(crate) fn motif_reference_span_for_tile(
    tile: &Tile,
    chrom_len: u64,
    clip_strategy: ClipStrategy,
    max_soft_clips: u16,
    k_outside: usize,
) -> Result<Interval<u64>> {
    let raw_extra = if clip_strategy.uses_shifted_boundary() {
        u64::from(max_soft_clips)
    } else {
        0
    };
    let pad = raw_extra.saturating_add(k_outside as u64);
    let reference_start = (tile.fetch_start() as u64).saturating_sub(pad);
    let reference_end = (tile.fetch_end() as u64).saturating_add(pad).min(chrom_len);
    Interval::new(reference_start, reference_end).map_err(Into::into)
}

/// Prepare tile-local masked reference resources for motif encoding.
///
/// This loads and masks the tile reference slice only when the current run
/// actually needs it:
/// - always for outside motifs
/// - for reference-backed inside motifs
/// - for read-backed inside motifs when blacklist-driven motif validation is active
///
/// Parameters
/// ----------
/// - `opt`:
///   Full `ends` command configuration
/// - `tile`:
///   Tile currently being processed
/// - `reference_span`:
///   Tile-local reference span guaranteed to cover every valid motif lookup for this tile
/// - `chrom_len`:
///   Full chromosome length
/// - `blacklist_intervals`:
///   Merged blacklist intervals for this chromosome
/// - `inside_spec`:
///   Shared codec spec for the inside half, or `None` when `k_inside = 0`
/// - `outside_spec`:
///   Shared codec spec for the outside half, or `None` when `k_outside = 0`
///
/// Returns
/// -------
/// - `Result<TileMotifContext<'a>>`:
///   Tile-local reference resources ready for motif validation and encoding
pub(crate) fn build_tile_motif_context<'a>(
    opt: &EndsConfig,
    tile: &Tile,
    reference_span: Interval<u64>,
    chrom_len: u64,
    blacklist_intervals: &'a [Interval<u64>],
    inside_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
) -> Result<TileMotifContext<'a>> {
    let inside_spec = inside_spec.cloned();
    let outside_spec = outside_spec.cloned();

    if matches!(opt.clip.clip_strategy, ClipStrategy::RawAlignedBoundary)
        && matches!(opt.source_inside, KmerSource::Reference)
    {
        anyhow::bail!(
            "`--clip-strategy raw-aligned-boundary` cannot be combined with `--source-inside reference`"
        );
    }

    let needs_reference_bases = outside_spec.is_some()
        || (inside_spec.is_some()
            && (matches!(opt.source_inside, KmerSource::Reference)
                || !blacklist_intervals.is_empty()));
    let (reference_start, reference_end) = reference_span.as_tuple();

    if !needs_reference_bases {
        return Ok(TileMotifContext {
            reference_start,
            reference_bases: None,
            inside_spec,
            outside_spec,
            inside_codes: None,
            outside_codes: None,
            blacklist_intervals,
            chrom_len,
        });
    }

    let ref_2bit = opt
        .ref_2bit
        .as_ref()
        .context("Reference-backed motif extraction requires --ref-2bit")?;
    let mut reference_bases = read_seq_in_range(
        ref_2bit,
        &tile.chr,
        (reference_start as usize)..(reference_end as usize),
    )?;
    if !blacklist_intervals.is_empty() {
        apply_blacklist_mask_to_seq(&mut reference_bases, blacklist_intervals, reference_start);
    }

    let (inside_codes, outside_codes) = match (inside_spec.as_ref(), outside_spec.as_ref()) {
        (Some(inside_spec), Some(outside_spec)) if inside_spec.k == outside_spec.k => {
            let shared_codes =
                build_precomputed_reference_codes(Some(inside_spec), &reference_bases);
            (shared_codes.clone(), shared_codes)
        }
        _ => (
            build_precomputed_reference_codes(inside_spec.as_ref(), &reference_bases),
            build_precomputed_reference_codes(outside_spec.as_ref(), &reference_bases),
        ),
    };

    Ok(TileMotifContext {
        reference_start,
        reference_bases: Some(reference_bases),
        inside_spec,
        outside_spec,
        inside_codes,
        outside_codes,
        blacklist_intervals,
        chrom_len,
    })
}

/// Precompute one tile-local radix-5 code vector for a single motif half.
///
/// Parameters
/// ----------
/// - `spec`:
///   Codec spec for the requested half, or `None` when that half is empty
/// - `reference_bases`:
///   Tile-local reference slice, already blacklist-masked when needed
///
/// Returns
/// -------
/// - `Option<Arc<KmerCodes>>`:
///   Shared per-position codes for that `k`, or `None` when the half is empty
fn build_precomputed_reference_codes(
    spec: Option<&KmerSpec>,
    reference_bases: &[u8],
) -> Option<Arc<KmerCodes>> {
    let spec = spec?;
    let k: u8 = spec
        .k
        .try_into()
        .expect("validated k-mer size should fit into u8");
    let mut spec_map = FxHashMap::default();
    spec_map.insert(k, spec.clone());
    let mut codes_by_k = build_left_aligned_codes_per_k(reference_bases, &spec_map);
    Some(Arc::new(codes_by_k.remove(&k).expect(
        "missing precomputed k-mer codes after precomputation",
    )))
}

/// Count the relevant end motifs from one fragment into one output window.
///
/// Endpoint mode filters left and right ends independently against the current
/// window. All other assignment modes count every kept end in every candidate
/// window selected by the outer fragment-level overlap logic.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Sparse output counts updated in place
/// - `original_idx`:
///   Global output-row index for the current window
/// - `window_interval`:
///   Genomic coordinates of the current output window
/// - `fragment`:
///   Fragment with already-resolved ends
/// - `weight`:
///   Combined overlap, scaling, and GC weight for this count
/// - `motif_context`:
///   Tile-local reference resources
/// - `source_inside`:
///   Whether inside bases come from the read or the reference
/// - `assign_by`:
///   Window-assignment rule for deciding whether each end counts here
///
/// Returns
/// -------
/// - `Result<CountedEndFlags>`:
///   Which fragment ends contributed at least one motif count in this window
pub(crate) fn count_fragment_in_window(
    counts_by_window: &mut EndCountsByWindow,
    original_idx: u64,
    window_interval: Interval<u64>,
    fragment: &FragmentWithEnds,
    weight: f64,
    motif_context: &TileMotifContext<'_>,
    source_inside: KmerSource,
    assign_by: WindowMotifAssigner,
) -> Result<CountedEndFlags> {
    let mut counted_end_flags = CountedEndFlags::default();

    if let Some(left_end) = fragment.left_end.as_ref() {
        let count_left_end = match assign_by {
            WindowMotifAssigner::Endpoint => {
                window_interval.contains_point(left_end.boundary_pos as u64)
            }
            _ => true,
        };
        if count_left_end {
            if let Some(key) =
                maybe_encode_end_motif_key(left_end, EndSide::Left, motif_context, source_inside)?
            {
                counts_by_window
                    .entry(original_idx)
                    .or_default()
                    .incr_weighted(key, weight);
                counted_end_flags.left_counted = true;
            }
        }
    }

    if let Some(right_end) = fragment.right_end.as_ref() {
        let right_endpoint_pos = right_end
            .boundary_pos
            .checked_sub(1)
            .expect("right boundary must be > 0 for a valid half-open interval");
        let count_right_end = match assign_by {
            WindowMotifAssigner::Endpoint => {
                window_interval.contains_point(right_endpoint_pos as u64)
            }
            _ => true,
        };
        if count_right_end {
            if let Some(key) =
                maybe_encode_end_motif_key(right_end, EndSide::Right, motif_context, source_inside)?
            {
                counts_by_window
                    .entry(original_idx)
                    .or_default()
                    .incr_weighted(key, weight);
                counted_end_flags.right_counted = true;
            }
        }
    }

    Ok(counted_end_flags)
}

/// Encode one end motif if both halves are valid.
///
/// Invalid means either:
/// - no full reference k-mer exists at the requested coordinate
/// - the masked reference span contains a blacklisted base
///
/// In both cases the end is skipped by returning `None`.
///
/// Parameters
/// ----------
/// - `end`:
///   Resolved fragment end to encode
/// - `end_side`:
///   Whether this is the left or right fragment end
/// - `motif_context`:
///   Tile-local reference resources
/// - `source_inside`:
///   Whether inside bases come from the read or the reference
///
/// Returns
/// -------
/// - `Result<Option<EncodedEndMotifKey>>`:
///   The encoded motif key for a kept end, or `None` when the end should be skipped
fn maybe_encode_end_motif_key(
    end: &ResolvedFragmentEnd,
    end_side: EndSide,
    motif_context: &TileMotifContext<'_>,
    source_inside: KmerSource,
) -> Result<Option<EncodedEndMotifKey>> {
    let inside_code = encode_inside_code(end, end_side, motif_context, source_inside)?;
    let outside_code = encode_outside_code(end.boundary_pos as u64, end_side, motif_context)?;
    if motif_code_is_invalid(inside_code, motif_context.inside_spec.as_ref())
        || motif_code_is_invalid(outside_code, motif_context.outside_spec.as_ref())
    {
        return Ok(None);
    }

    Ok(Some(EncodedEndMotifKey {
        inside_code,
        outside_code,
        reverse_on_decode: matches!(end_side, EndSide::Right),
    }))
}

/// Check whether a code represents an invalid motif half.
///
/// Parameters
/// ----------
/// - `code`:
///   Encoded motif-half value
/// - `spec`:
///   Matching codec spec for that half, or `None` when the half is empty
///
/// Returns
/// -------
/// - `bool`:
///   `true` when the code is a sentinel and the end should be skipped
#[inline]
fn motif_code_is_invalid(code: u64, spec: Option<&KmerSpec>) -> bool {
    let Some(spec) = spec else {
        return false;
    };
    code == spec.sentinel_none() || code == spec.sentinel_n()
}

/// Encode the inside-fragment half for one end.
///
/// Read-backed mode validates blacklist overlap from the masked reference first,
/// then encodes the actual read bases. Reference-backed mode encodes directly
/// from the masked reference lookup.
///
/// Parameters
/// ----------
/// - `end`:
///   Resolved fragment end to encode
/// - `end_side`:
///   Whether this is the left or right fragment end
/// - `motif_context`:
///   Tile-local reference resources
/// - `source_inside`:
///   Whether inside bases come from the read or the reference
///
/// Returns
/// -------
/// - `Result<u64>`:
///   Encoded inside-half code or an invalid sentinel
fn encode_inside_code(
    end: &ResolvedFragmentEnd,
    end_side: EndSide,
    motif_context: &TileMotifContext<'_>,
    source_inside: KmerSource,
) -> Result<u64> {
    let Some(spec) = motif_context.inside_spec.as_ref() else {
        return Ok(0);
    };

    match source_inside {
        KmerSource::Read => {
            if let Some(masked_reference_code) =
                validate_blacklist_for_read_inside_code(end, end_side, spec, motif_context)?
            {
                if motif_code_is_invalid(masked_reference_code, Some(spec)) {
                    return Ok(masked_reference_code);
                }
            }
            Ok(spec.encode_kmer_bytes(&end.inside_bases))
        }
        KmerSource::Reference => {
            let start_pos = match end_side {
                EndSide::Left => end.boundary_pos as u64,
                EndSide::Right => {
                    let k = spec.k as u64;
                    if (end.boundary_pos as u64) < k {
                        return Ok(spec.sentinel_none());
                    }
                    end.boundary_pos as u64 - k
                }
            };
            get_reference_code(
                start_pos,
                spec,
                motif_context.inside_codes.as_deref(),
                motif_context,
            )
        }
    }
}

/// Validate the reference-addressable part of a read-backed inside motif against the
/// blacklist-masked reference.
///
/// This is only used for `source_inside=read`, where the emitted inside code still comes from the
/// read but blacklist filtering must remain genomic. `inside_reference_validation_bp` tells us how
/// much of `inside_bases` still maps to concrete reference positions; in
/// `raw-aligned-boundary`, clipped-only inside bases are intentionally ignored here because they
/// lie outside the aligned reference span.
///
/// Parameters
/// ----------
/// - `end`:
///   Resolved fragment end to validate
/// - `end_side`:
///   Whether this is the left or right fragment end
/// - `spec`:
///   Inside-half k-mer spec
/// - `motif_context`:
///   Tile-local reference resources
///
/// Returns
/// -------
/// - `Result<Option<u64>>`:
///   `None` when validation passed and the caller should keep the read-backed code, otherwise an
///   invalid sentinel that should be returned instead
fn validate_blacklist_for_read_inside_code(
    end: &ResolvedFragmentEnd,
    end_side: EndSide,
    spec: &KmerSpec,
    motif_context: &TileMotifContext<'_>,
) -> Result<Option<u64>> {
    if motif_context.blacklist_intervals.is_empty() {
        return Ok(None);
    }

    let validation_bp = end.inside_reference_validation_bp;
    if validation_bp == 0 {
        return Ok(None);
    }

    let start_pos = match end_side {
        EndSide::Left => end.boundary_pos as u64,
        EndSide::Right => {
            let validation_bp = validation_bp as u64;
            if (end.boundary_pos as u64) < validation_bp {
                return Ok(Some(spec.sentinel_none()));
            }
            end.boundary_pos as u64 - validation_bp
        }
    };

    let is_unmasked =
        masked_reference_span_is_valid(start_pos, validation_bp, spec, motif_context)?;
    if is_unmasked {
        Ok(None)
    } else {
        Ok(Some(spec.sentinel_n()))
    }
}

/// Check whether a reference-addressable genomic span stays inside the preloaded tile reference and
/// avoids blacklist masking.
fn masked_reference_span_is_valid(
    start_pos: u64,
    span_bp: usize,
    spec: &KmerSpec,
    motif_context: &TileMotifContext<'_>,
) -> Result<bool> {
    if start_pos + span_bp as u64 > motif_context.chrom_len {
        return Ok(false);
    }

    let reference_bases = motif_context
        .reference_bases
        .as_ref()
        .context("missing preloaded reference bases for blacklist validation")?;
    let local_start =
        try_reference_start_index(start_pos, span_bp, motif_context).with_context(|| {
            let loaded_end = motif_context.reference_start + reference_bases.len() as u64;
            format!(
                "motif reference lookup escaped preloaded tile span: start={start_pos}, k={}, loaded_reference_span=[{}, {})",
                spec.k, motif_context.reference_start, loaded_end
            )
        })?;

    Ok(!reference_bases[local_start..local_start + span_bp].contains(&BLACKLIST_BYTE))
}

/// Encode the outside-fragment half for one end from reference coordinates.
///
/// Parameters
/// ----------
/// - `boundary_pos`:
///   Assignment boundary for the end
/// - `end_side`:
///   Whether this is the left or right fragment end
/// - `motif_context`:
///   Tile-local reference resources
///
/// Returns
/// -------
/// - `Result<u64>`:
///   Encoded outside-half code or an invalid sentinel
fn encode_outside_code(
    boundary_pos: u64,
    end_side: EndSide,
    motif_context: &TileMotifContext<'_>,
) -> Result<u64> {
    let Some(spec) = motif_context.outside_spec.as_ref() else {
        return Ok(0);
    };

    let start_pos = match end_side {
        EndSide::Left => {
            let k = spec.k as u64;
            if boundary_pos < k {
                return Ok(spec.sentinel_none());
            }
            boundary_pos - k
        }
        EndSide::Right => boundary_pos,
    };

    get_reference_code(
        start_pos,
        spec,
        motif_context.outside_codes.as_deref(),
        motif_context,
    )
}

/// Look up one reference-backed motif half from the precomputed tile codes.
///
/// Parameters
/// ----------
/// - `start_pos`:
///   Genomic start position of the requested motif half
/// - `spec`:
///   Codec spec for the requested half
/// - `precomputed_codes`:
///   Tile-local code vector for this `k`
/// - `motif_context`:
///   Tile-local reference resources
///
/// Returns
/// -------
/// - `Result<u64>`:
///   Encoded reference-backed motif code or an invalid sentinel
fn get_reference_code(
    start_pos: u64,
    spec: &KmerSpec,
    precomputed_codes: Option<&KmerCodes>,
    motif_context: &TileMotifContext<'_>,
) -> Result<u64> {
    if start_pos + spec.k as u64 > motif_context.chrom_len {
        return Ok(spec.sentinel_none());
    }

    let codes = precomputed_codes
        .context("missing precomputed reference codes for a reference-backed motif lookup")?;
    let local_start = try_reference_start_index(start_pos, spec.k, motif_context).with_context(
        || {
            let loaded_end = motif_context.reference_start
                + motif_context.reference_bases.as_ref().map_or(0, Vec::len) as u64;
            format!(
                "motif reference lookup escaped preloaded tile span: start={start_pos}, k={}, loaded_reference_span=[{}, {})",
                spec.k, motif_context.reference_start, loaded_end
            )
        },
    )?;
    Ok(codes.get(local_start))
}

/// Translate a genomic motif start into a tile-local slice index when possible.
///
/// Parameters
/// ----------
/// - `start_pos`:
///   Genomic start position of the requested motif half
/// - `k`:
///   Requested motif length
/// - `motif_context`:
///   Tile-local reference resources
///
/// Returns
/// -------
/// - `Option<usize>`:
///   Tile-local start index, or `None` when the motif crosses the loaded tile slice
fn try_reference_start_index(
    start_pos: u64,
    k: usize,
    motif_context: &TileMotifContext<'_>,
) -> Option<usize> {
    if start_pos < motif_context.reference_start {
        return None;
    }

    let local_start = (start_pos - motif_context.reference_start) as usize;
    if local_start + k > motif_context.reference_bases.as_ref().map_or(0, Vec::len) {
        return None;
    }
    Some(local_start)
}

#[cfg(test)]
mod tests {
    include!("motifs_tests.rs");
}
