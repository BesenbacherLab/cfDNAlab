#[cfg(feature = "cli")]
pub(crate) mod cli_app;
pub(crate) mod commands;
pub mod error;
pub(crate) mod shared;

pub use error::{Error, Result};

#[cfg(feature = "cli")]
pub fn run_cli() {
    cli_app::run_cli();
}

#[cfg(feature = "cli")]
pub fn build_docs_command() -> clap::Command {
    cli_app::build_docs_command()
}

pub mod reusable {
    pub mod interval {
        pub use crate::shared::interval::{
            IndexedInterval, Interval, Span, TouchingMergePolicy, merge_intervals,
            merge_sorted_intervals, push_merged_interval,
        };
    }

    pub mod blacklist {
        pub use crate::shared::blacklist::{
            BlacklistStrategy, apply_blacklist_mask_to_seq, compute_blacklist_overlap,
            is_blacklisted, load_blacklists,
        };
    }

    pub mod scale_genome {
        pub use crate::shared::scale_genome::{
            LoadedScalingFactors, ScalingBin, ScalingFactorsMetadata, ScalingGCMode, WindowScaling,
            apply_scaling_to_coverage_in_place, compute_per_window_scaling_over_fragment,
            compute_per_window_scaling_over_overlap, load_scaling_factors_tsv,
            scaling_gc_mode_for_run,
        };
    }

    pub mod gc_bias {
        pub use crate::commands::gc_bias::correct::{
            GCCorrector, GCLengthRange, LengthAgnosticGCCorrector,
            MarginalizeLengthsWeightingScheme, load_gc_corrector,
            load_length_agnostic_gc_corrector,
        };
        pub use crate::commands::gc_bias::counting::{
            GCPrefixes, build_gc_prefixes, get_gc_integer_percentage_for_window,
        };
        pub use crate::commands::gc_bias::package::GCCorrectionPackage;
    }

    pub mod bam {
        pub use crate::shared::bam::{Contigs, bam_contigs_info, bam_header_contigs};
    }

    pub mod reference {
        pub use crate::shared::reference::{
            ContigFootprintEntry, load_chrom_sizes, load_chrom_sizes_with_order, read_seq,
            read_seq_in_range, twobit_contig_footprint, twobit_contig_lengths, twobit_contig_names,
        };
    }

    pub mod overlaps {
        pub use crate::shared::overlaps::{
            OverlappingWindow, OverlappingWindows, create_overlapping_bins_by_size,
            find_overlapping_windows, fraction_overlap_of_a, half_open_intervals_overlap,
            overlap_len,
        };
    }

    pub mod positioning {
        pub use crate::shared::positioning::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
    }

    pub mod constants {
        pub use crate::shared::constants::{
            COVERAGE_WEIGHT_AUX_TAG, DEFAULT_MAX_SOFT_CLIPS, FRAGMENT_COUNT_WEIGHT_AUX_TAG,
            FRAGMENT_LENGTH_AUX_TAG, GC_WEIGHT_AUX_TAG, MAX_MAX_SOFT_CLIPS,
            MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION,
        };
    }

    pub mod parsing {
        pub use crate::commands::cli_common::{
            LengthBin, LengthBins, parse_length_bins, parse_output_prefix, parse_sam_aux_tag_name,
            resolve_length_bin_edges, validate_output_prefix, validate_sam_aux_tag_name,
        };
    }
}

pub mod run_like_cli {
    pub mod common {
        pub use crate::commands::cli_common::{
            ApplyGCArgFileOnly, ApplyGCArgs, AssignToWindowArgs, BaseSelectionArgs, ChromosomeArgs,
            ContigSource, ContigSourceKind, DistributionWindowSpec, DistributionWindowsArgs,
            FragmentLengthArgs, FragmentPositionSelectionArgs, GCWindowsArgs, IOCArgs, LogSpec,
            LoggingArgs, Ref2BitOptionalForGCArgs, Ref2BitRequiredArgs, ScaleGenomeArgs,
            UnpairedArgs, UnparsedPositionalSelectionSpec, WindowAssigner, WindowSpec, WindowsArgs,
        };
        pub use crate::shared::blacklist::BlacklistStrategy;
        pub use crate::shared::clip_mode::ClipMode;
        pub use crate::shared::indel_mode::{IndelMode, IndelMotifFilterPolicy};
        pub use crate::shared::positioning::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
    }

    #[cfg(feature = "cmd_bam_to_bam")]
    pub mod bam_to_bam {
        pub use crate::commands::bam_to_bam::config::BamToBamConfig;

        pub fn run_like_cli(config: &BamToBamConfig) -> anyhow::Result<()> {
            crate::commands::bam_to_bam::bam_to_bam::run(config)
        }
    }

    #[cfg(feature = "cmd_bam_to_frag")]
    pub mod bam_to_frag {
        pub use crate::commands::bam_to_frag::config::BamToFragConfig;

        pub fn run_like_cli(config: &BamToFragConfig) -> anyhow::Result<()> {
            crate::commands::bam_to_frag::bam_to_frag::run(config)
        }
    }

    #[cfg(feature = "cmd_coverage_weights")]
    pub mod coverage_weights {
        pub use crate::commands::coverage_weights::config::CoverageWeightsConfig;
        pub use crate::commands::coverage_weights::scaling_weights_config::ScalingWeightsArgs;

        pub fn run_like_cli(config: &CoverageWeightsConfig) -> anyhow::Result<()> {
            crate::commands::coverage_weights::coverage_weights::run(config)
        }
    }

    #[cfg(feature = "cmd_ends")]
    pub mod ends {
        pub use crate::commands::ends::config::EndsConfig;
        pub use crate::commands::ends::config_structs::{
            AssignMotifToWindowArgs, BaseQualityAggregation, BaseQualityComparisonOp,
            BaseQualityFilter, BaseQualityFilterScope, ClipStrategy, ClippingArgs, KmerSource,
            WindowMotifAssigner,
        };

        pub fn run_like_cli(config: &EndsConfig) -> anyhow::Result<()> {
            crate::commands::ends::ends::run(config)
        }
    }

    #[cfg(feature = "cmd_fcoverage")]
    pub mod fcoverage {
        pub use crate::commands::fcoverage::config::{FCoverageConfig, LengthNormalizationMode};
        pub use crate::commands::fcoverage::window_results::CoverageWindowAction;

        pub fn run_like_cli(config: &FCoverageConfig) -> anyhow::Result<()> {
            crate::commands::fcoverage::fcoverage::run(config)
        }
    }

    #[cfg(feature = "cmd_frag_to_bam")]
    pub mod frag_to_bam {
        pub use crate::commands::frag_to_bam::config::FragToBamConfig;

        pub fn run_like_cli(config: &FragToBamConfig) -> anyhow::Result<()> {
            crate::commands::frag_to_bam::frag_to_bam::run(config)
        }
    }

    #[cfg(feature = "cmd_fragment_count_weights")]
    pub mod fragment_count_weights {
        pub use crate::commands::fragment_count_weights::config::FragmentCountWeightsConfig;

        pub fn run_like_cli(config: &FragmentCountWeightsConfig) -> anyhow::Result<()> {
            crate::commands::fragment_count_weights::fragment_count_weights::run(config)
        }
    }

    #[cfg(feature = "cmd_fragment_kmers")]
    pub mod fragment_kmers {
        pub use crate::commands::fragment_kmers::config::{
            FragmentKmersConfig, FragmentKmersSharedArgs,
        };

        pub fn run_like_cli(config: &FragmentKmersConfig) -> anyhow::Result<()> {
            crate::commands::fragment_kmers::fragment_kmers::run(config)
        }
    }

    #[cfg(feature = "cmd_gc_bias")]
    pub mod gc_bias {
        pub use crate::commands::gc_bias::config::{GCConfig, OutlierMethodArg, OutlierScopeArg};

        pub fn run_like_cli(config: &GCConfig) -> anyhow::Result<()> {
            crate::commands::gc_bias::gc_bias::run(config)
        }
    }

    #[cfg(feature = "cmd_lengths")]
    pub mod lengths {
        pub use crate::commands::lengths::config::LengthsConfig;

        pub fn run_like_cli(config: &LengthsConfig) -> anyhow::Result<()> {
            crate::commands::lengths::lengths::run(config)
        }
    }

    #[cfg(feature = "cmd_midpoints")]
    pub mod midpoints {
        pub use crate::commands::midpoints::config::MidpointsConfig;
        pub use crate::commands::midpoints::smoothing::MidpointSmoothing;

        pub fn run_like_cli(config: &MidpointsConfig) -> anyhow::Result<()> {
            crate::commands::midpoints::midpoints::run(config)
        }
    }

    #[cfg(feature = "cmd_prepare_windows")]
    pub mod prepare_windows {
        pub use crate::commands::prepare_windows::config::{
            ComposeSpec, CoordinateSet, DedupKeep, DistSign, DistancePolicy, HeaderMode,
            MergeLabel, MergeScope, MissingScore, NearDirection, NearEdge, NearTiePolicy,
            OobPolicy, PrepareConfig,
        };
        pub use crate::commands::prepare_windows::near_file::NearDuplicatesPolicy;

        pub fn run_like_cli(config: &PrepareConfig) -> anyhow::Result<()> {
            crate::commands::prepare_windows::prepare_windows::run(config)
        }
    }

    #[cfg(feature = "cmd_ref_gc_bias")]
    pub mod ref_gc_bias {
        pub use crate::commands::ref_gc_bias::config::{RefGCBiasConfig, RefGCWindowsArgs};

        pub fn run_like_cli(config: &RefGCBiasConfig) -> anyhow::Result<()> {
            crate::commands::ref_gc_bias::ref_gc_bias::run(config)
        }
    }

    #[cfg(feature = "cmd_transitions")]
    pub mod transitions {
        pub use crate::commands::transitions::config::TransitionsConfig;

        pub fn run_like_cli(config: &TransitionsConfig) -> anyhow::Result<()> {
            crate::commands::transitions::transitions::run(config)
        }
    }

    #[cfg(feature = "cmd_visualize_positions")]
    pub mod visualize_positions {
        pub use crate::commands::visualize_positions::config::VisualizePositionsConfig;
        pub use crate::commands::visualize_positions::model::Style;

        pub fn run_like_cli(config: &VisualizePositionsConfig) -> anyhow::Result<()> {
            crate::commands::visualize_positions::visualize_positions::run(config)
        }
    }

    #[cfg(feature = "cmd_wps")]
    pub mod wps {
        pub use crate::commands::wps::config::{WPSConfig, WPSSharedConfig};

        pub fn run_like_cli(config: &WPSConfig) -> anyhow::Result<()> {
            crate::commands::wps::wps::run(config)
        }
    }

    #[cfg(feature = "cmd_wps_peaks")]
    pub mod wps_peaks {
        pub use crate::commands::wps_peaks::config::WPSPeaksConfig;
        pub use crate::commands::wps_peaks::window_peak_results::PeaksWindowAction;

        pub fn run_like_cli(config: &WPSPeaksConfig) -> anyhow::Result<()> {
            crate::commands::wps_peaks::wps_peaks::run(config)
        }
    }
}
