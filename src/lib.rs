#[cfg(feature = "cli")]
pub(crate) mod cli_app;
pub mod cli_command;
pub mod command_run;
pub(crate) mod commands;
pub mod error;
pub mod output_loaders;
pub(crate) mod shared;
#[cfg(any(feature = "testing", test))]
pub mod testing;

pub use cli_command::ToCliCommand;
pub use command_run::{CommandRunResult, RunOptions};
pub use error::{Error, Result};

#[cfg(feature = "cli")]
pub fn run_cli() {
    cli_app::run_cli();
}

#[cfg(feature = "cli")]
pub fn build_docs_command() -> clap::Command {
    cli_app::build_docs_command()
}

pub mod interval {
    pub use crate::shared::interval::{
        IndexedInterval, Interval, Span, TouchingMergePolicy, merge_intervals,
        merge_sorted_intervals, push_merged_interval,
    };
}

pub mod blacklist {
    pub use crate::shared::blacklist::{
        BlacklistStrategy, apply_blacklist_mask_to_seq, compute_blacklist_overlap, is_blacklisted,
        load_blacklists,
    };
}

pub mod scale_genome {
    pub use crate::shared::scale_genome::{
        LoadedScalingFactors, ScalingBin, ScalingFactorsMetadata, ScalingGCMode, WindowScaling,
        apply_scaling_to_coverage_in_place, compute_per_window_scaling_over_fragment,
        compute_per_window_scaling_over_overlap, load_scaling_factors_tsv, scaling_gc_mode_for_run,
    };
}

#[cfg(feature = "cmd_gc_bias")]
pub mod gc_bias {
    pub use crate::commands::gc_bias::correct::{
        GCCorrector, GCLengthRange, LengthAgnosticGCCorrector, MarginalizeLengthsWeightingScheme,
        load_gc_corrector, load_length_agnostic_gc_corrector,
    };
    pub use crate::commands::gc_bias::counting::{
        GCPrefixes, build_gc_prefixes, get_gc_integer_percentage_for_window,
    };
    pub use crate::commands::gc_bias::package::GCCorrectionPackage;
}

pub mod bam {
    pub use crate::shared::bam::{Contigs, bam_contigs_info, bam_header_contigs};
}

pub mod fragment {
    pub use crate::shared::fragment::minimal_fragment::{
        Fragment, MinimalReadInfo, collect_fragment, collect_fragment_from_single_read,
    };
    pub use crate::shared::read::{default_include_read_paired_end, default_include_read_unpaired};
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
        find_overlapping_windows, fraction_overlap_of_a, half_open_intervals_overlap, overlap_len,
    };
}

pub mod positioning {
    pub use crate::shared::positioning::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
}

pub mod indel_mode {
    pub use crate::shared::indel_mode::{IndelMode, IndelMotifFilterPolicy};
}

pub mod clip_mode {
    pub use crate::shared::clip_mode::ClipMode;
}

pub mod constants {
    pub use crate::shared::constants::{
        COVERAGE_WEIGHT_AUX_TAG, DEFAULT_MAX_SOFT_CLIPS, FRAGMENT_COUNT_WEIGHT_AUX_TAG,
        FRAGMENT_LENGTH_AUX_TAG, GC_CORRECTION_SCHEMA_VERSION, GC_WEIGHT_AUX_TAG,
        MAX_MAX_SOFT_CLIPS, MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION,
    };
}

pub mod parsing {
    pub use crate::commands::cli_common::{
        LengthBin, LengthBins, parse_length_bins, parse_output_prefix, parse_sam_aux_tag_name,
        resolve_length_bin_edges, validate_output_prefix, validate_sam_aux_tag_name,
    };
}

pub mod run_like_cli {
    pub use crate::command_run::{CommandRunResult, RunOptions};

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
        pub use crate::commands::bam_to_bam::bam_to_bam::{BamToBamRunResult, run_bam_to_bam};
        pub use crate::commands::bam_to_bam::config::BamToBamConfig;
    }

    #[cfg(feature = "cmd_bam_to_frag")]
    pub mod bam_to_frag {
        pub use crate::commands::bam_to_frag::bam_to_frag::{BamToFragRunResult, run_bam_to_frag};
        pub use crate::commands::bam_to_frag::config::BamToFragConfig;
    }

    #[cfg(feature = "cmd_coverage_weights")]
    pub mod coverage_weights {
        pub use crate::commands::coverage_weights::config::CoverageWeightsConfig;
        pub use crate::commands::coverage_weights::coverage_weights::{
            CoverageWeightsRunResult, run_coverage_weights,
        };
        pub use crate::commands::coverage_weights::scaling_weights_config::ScalingWeightsArgs;
    }

    #[cfg(feature = "cmd_ends")]
    pub mod ends {
        pub use crate::commands::ends::config::EndsConfig;
        pub use crate::commands::ends::config_structs::{
            AssignMotifToWindowArgs, BaseQualityAggregation, BaseQualityComparisonOp,
            BaseQualityFilter, BaseQualityFilterScope, ClipStrategy, ClippingArgs, KmerSource,
            WindowMotifAssigner,
        };
        pub use crate::commands::ends::ends::{EndsRunResult, run_ends};
    }

    #[cfg(feature = "cmd_fcoverage")]
    pub mod fcoverage {
        pub use crate::commands::fcoverage::config::{FCoverageConfig, LengthNormalizationMode};
        pub use crate::commands::fcoverage::fcoverage::{FCoverageRunResult, run_fcoverage};
        pub use crate::commands::fcoverage::window_results::CoverageWindowAction;
    }

    #[cfg(feature = "cmd_frag_to_bam")]
    pub mod frag_to_bam {
        pub use crate::commands::frag_to_bam::config::FragToBamConfig;
        pub use crate::commands::frag_to_bam::frag_to_bam::{FragToBamRunResult, run_frag_to_bam};
    }

    #[cfg(feature = "cmd_fragment_count_weights")]
    pub mod fragment_count_weights {
        pub use crate::commands::fragment_count_weights::config::FragmentCountWeightsConfig;
        pub use crate::commands::fragment_count_weights::fragment_count_weights::{
            FragmentCountWeightsRunResult, run_fragment_count_weights,
        };
    }

    #[cfg(feature = "cmd_fragment_kmers")]
    pub mod fragment_kmers {
        pub use crate::commands::fragment_kmers::config::{
            FragmentKmersConfig, FragmentKmersSharedArgs,
        };
        pub use crate::commands::fragment_kmers::fragment_kmers::{
            FragmentKmersRunResult, run_fragment_kmers,
        };
    }

    #[cfg(feature = "cmd_gc_bias")]
    pub mod gc_bias {
        pub use crate::commands::gc_bias::config::{GCConfig, OutlierMethodArg, OutlierScopeArg};
        pub use crate::commands::gc_bias::gc_bias::{GCBiasRunResult, run_gc_bias};
    }

    #[cfg(feature = "cmd_lengths")]
    pub mod lengths {
        pub use crate::commands::lengths::config::LengthsConfig;
        pub use crate::commands::lengths::lengths::{LengthsRunResult, run_lengths};
    }

    #[cfg(feature = "cmd_midpoints")]
    pub mod midpoints {
        pub use crate::commands::midpoints::config::MidpointsConfig;
        pub use crate::commands::midpoints::midpoints::{MidpointsRunResult, run_midpoints};
        pub use crate::commands::midpoints::smoothing::MidpointSmoothing;
    }

    #[cfg(feature = "cmd_prepare_windows")]
    pub mod prepare_windows {
        pub use crate::commands::prepare_windows::config::{
            ComposeSpec, CoordinateSet, DedupKeep, DistSign, DistancePolicy, HeaderMode,
            MergeLabel, MergeScope, MissingScore, NearDirection, NearEdge, NearTiePolicy,
            OobPolicy, PrepareConfig,
        };
        pub use crate::commands::prepare_windows::near_file::NearDuplicatesPolicy;
        pub use crate::commands::prepare_windows::prepare_windows::{
            PrepareWindowsRunResult, run_prepare_windows,
        };
    }

    #[cfg(feature = "cmd_ref_gc_bias")]
    pub mod ref_gc_bias {
        pub use crate::commands::ref_gc_bias::config::{RefGCBiasConfig, RefGCWindowsArgs};
        pub use crate::commands::ref_gc_bias::ref_gc_bias::{RefGCBiasRunResult, run_ref_gc_bias};
    }

    #[cfg(feature = "cmd_transitions")]
    pub mod transitions {
        pub use crate::commands::transitions::config::TransitionsConfig;
        pub use crate::commands::transitions::transitions::{
            TransitionsRunResult, run_transitions,
        };
    }

    #[cfg(feature = "cmd_visualize_positions")]
    pub mod visualize_positions {
        pub use crate::commands::visualize_positions::config::VisualizePositionsConfig;
        pub use crate::commands::visualize_positions::model::Style;
        pub use crate::commands::visualize_positions::visualize_positions::{
            VisualizePositionsRunResult, run_visualize_positions,
        };
    }

    #[cfg(feature = "cmd_wps")]
    pub mod wps {
        pub use crate::commands::wps::config::{WPSConfig, WPSSharedConfig};
        pub use crate::commands::wps::wps::{WPSRunResult, run_wps};
    }

    #[cfg(feature = "cmd_wps_peaks")]
    pub mod wps_peaks {
        pub use crate::commands::wps_peaks::config::WPSPeaksConfig;
        pub use crate::commands::wps_peaks::window_peak_results::PeaksWindowAction;
        pub use crate::commands::wps_peaks::wps_peaks::{WPSPeaksRunResult, run_wps_peaks};
    }
}
