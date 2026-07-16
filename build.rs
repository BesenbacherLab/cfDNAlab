//! Defines internal cfg names for repeated feature groups.
//!
//! Cargo features are the public switches users enable, such as `cmd_lengths`
//! or `cmd_fcoverage`. Several source files need the same private groups of
//! those features, for example "commands that read indexed BAM files" or
//! "commands that load grouped BED files". This script defines those private
//! cfg names for rustc so source files can use one readable predicate without
//! adding more Cargo features or another dependency.

/// Define the internal cfg names used by the crate.
fn main() {
    define_any_feature_cfg(
        "output_loader_api",
        &[
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_ref_kmers",
        ],
    );
    define_any_feature_cfg(
        "has_cli_commands",
        &[
            "cmd_bam_to_bam",
            "cmd_bam_to_frag",
            "cmd_coverage_weights",
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_frag_to_bam",
            "cmd_fragment_count_weights",
            "cmd_fragment_kmers",
            "cmd_gc_bias",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_prepare_windows",
            "cmd_ref_kmers",
            "cmd_transitions",
            "cmd_visualize_positions",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "uses_kmers",
        &["cmd_ends", "cmd_fragment_kmers", "cmd_ref_kmers"],
    );
    define_any_feature_cfg(
        "reads_indexed_bam",
        &[
            "cmd_bam_to_bam",
            "cmd_bam_to_frag",
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_fragment_kmers",
            "cmd_gc_bias",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg("writes_bam_output", &["cmd_bam_to_bam", "cmd_frag_to_bam"]);
    define_any_feature_cfg(
        "loads_grouped_bed",
        &[
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_ref_kmers",
        ],
    );
    define_any_feature_cfg(
        "writes_text_outputs",
        &[
            "cmd_ends",
            "cmd_fragment_kmers",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_prepare_windows",
        ],
    );
    define_any_feature_cfg(
        "uses_progress_reporting",
        &[
            "cmd_bam_to_bam",
            "cmd_bam_to_frag",
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_fragment_kmers",
            "cmd_gc_bias",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_prepare_windows",
            "cmd_ref_kmers",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "uses_temp_chrom_names",
        &[
            "cmd_bam_to_frag",
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_frag_to_bam",
            "cmd_fragment_kmers",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_prepare_windows",
            "cmd_ref_kmers",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "uses_temp_chrom_name_map",
        &[
            "cmd_bam_to_frag",
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_frag_to_bam",
            "cmd_fragment_kmers",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_ref_kmers",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "uses_temp_chrom_paths",
        &["cmd_bam_to_frag", "cmd_frag_to_bam"],
    );
    define_any_feature_cfg(
        "uses_tile_window_helpers",
        &[
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_fragment_kmers",
            "cmd_gc_bias",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_ref_kmers",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "uses_tile_window_iter",
        &[
            "cmd_fcoverage",
            "cmd_fragment_kmers",
            "cmd_gc_bias",
            "cmd_midpoints",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "uses_bed_window_tier_helpers",
        &["cmd_ends", "cmd_lengths", "cmd_ref_kmers"],
    );
    define_any_feature_cfg(
        "uses_tile_bed_overlap_context",
        &["cmd_ends", "cmd_lengths"],
    );
    define_any_feature_cfg(
        "uses_temp_dirs",
        &[
            "cmd_bam_to_bam",
            "cmd_bam_to_frag",
            "cmd_coverage_weights",
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_frag_to_bam",
            "cmd_fragment_kmers",
            "cmd_gc_bias",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_prepare_windows",
            "cmd_ref_kmers",
            "cmd_transitions",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "checks_tile_bam_tid",
        &[
            "cmd_ends",
            "cmd_fcoverage",
            "cmd_lengths",
            "cmd_midpoints",
            "cmd_wps",
            "cmd_wps_peaks",
        ],
    );
    define_any_feature_cfg(
        "flattens_bed_windows",
        &["cmd_fcoverage", "cmd_gc_bias", "cmd_wps", "cmd_wps_peaks"],
    );
}

/// Define `alias` when at least one Cargo feature in `features` is enabled.
///
/// `cargo:rustc-check-cfg` tells rustc that the cfg name is intentional even
/// when this build does not enable it. `cargo:rustc-cfg` is printed only when
/// one of the listed features is active, matching an `any(feature = "...")`
/// gate without repeating that long predicate at every use site.
fn define_any_feature_cfg(alias: &str, features: &[&str]) {
    println!("cargo:rustc-check-cfg=cfg({alias})");
    if features.iter().any(|feature| feature_enabled(feature)) {
        println!("cargo:rustc-cfg={alias}");
    }
}

/// Return whether Cargo enabled a feature for this build.
fn feature_enabled(feature: &str) -> bool {
    let env_name = format!(
        "CARGO_FEATURE_{}",
        feature.to_ascii_uppercase().replace('-', "_")
    );
    std::env::var_os(env_name).is_some()
}
