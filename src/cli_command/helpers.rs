#![allow(
    dead_code,
    reason = "feature-limited builds compile this helper module but use only the enabled commands' helpers"
)]

use std::{
    ffi::OsString,
    path::{Path, PathBuf},
};

#[cfg(any(
    feature = "cmd_fragment_kmers",
    feature = "cmd_transitions",
    feature = "cmd_wps",
    feature = "cmd_wps_peaks"
))]
use crate::commands::cli_common::WindowsArgs;
use crate::commands::cli_common::{
    ApplyGCArgFileOnly, ApplyGCArgs, AssignToWindowArgs, ChromosomeArgs, DistributionWindowsArgs,
    FragmentLengthArgs, GCWindowsArgs, IOCArgs, LoggingArgs, Ref2BitRequiredArgs, ScaleGenomeArgs,
    UnpairedArgs, WindowAssigner,
};
#[cfg(any(
    feature = "cmd_fragment_kmers",
    feature = "cmd_transitions",
    feature = "cmd_visualize_positions"
))]
use crate::commands::cli_common::{BaseSelectionArgs, FragmentPositionSelectionArgs};
use crate::shared::{
    blacklist::BlacklistStrategy,
    clip_mode::ClipMode,
    indel_mode::{IndelMode, IndelMotifFilterPolicy},
    logging::LogSpec,
};

pub(crate) fn shell_quote(value: &OsString) -> String {
    let text = value.to_string_lossy();
    if text.is_empty() {
        return "''".to_string();
    }
    if text
        .chars()
        .all(|character| character.is_ascii_alphanumeric() || "_+-./:=,@%".contains(character))
    {
        text.into_owned()
    } else {
        format!("'{}'", text.replace('\'', "'\\''"))
    }
}

pub(crate) fn command_args(subcommand: &str) -> Vec<OsString> {
    vec![OsString::from("cfdna"), OsString::from(subcommand)]
}

pub(crate) fn push_flag(args: &mut Vec<OsString>, flag: &'static str) {
    args.push(OsString::from(flag));
}

pub(crate) fn push_bool(args: &mut Vec<OsString>, flag: &'static str, enabled: bool) {
    if enabled {
        push_flag(args, flag);
    }
}

pub(crate) fn push_value(args: &mut Vec<OsString>, flag: &'static str, value: impl ToString) {
    args.push(OsString::from(flag));
    args.push(OsString::from(value.to_string()));
}

pub(crate) fn push_path(args: &mut Vec<OsString>, flag: &'static str, path: &Path) {
    args.push(OsString::from(flag));
    args.push(path.as_os_str().to_os_string());
}

pub(crate) fn push_optional_path(
    args: &mut Vec<OsString>,
    flag: &'static str,
    path: Option<&Path>,
) {
    if let Some(path) = path {
        push_path(args, flag, path);
    }
}

pub(crate) fn push_path_values(
    args: &mut Vec<OsString>,
    flag: &'static str,
    paths: Option<&[PathBuf]>,
) {
    if let Some(paths) = paths.filter(|paths| !paths.is_empty()) {
        args.push(OsString::from(flag));
        args.extend(paths.iter().map(|path| path.as_os_str().to_os_string()));
    }
}

pub(crate) fn push_values<T: ToString>(args: &mut Vec<OsString>, flag: &'static str, values: &[T]) {
    if !values.is_empty() {
        args.push(OsString::from(flag));
        args.extend(values.iter().map(|value| OsString::from(value.to_string())));
    }
}

pub(crate) fn push_ioc(args: &mut Vec<OsString>, ioc: &IOCArgs) {
    push_path(args, "--bam", &ioc.bam);
    push_path(args, "--output-dir", &ioc.output_dir);
    push_value(args, "--n-threads", ioc.n_threads);
}

pub(crate) fn push_unpaired(args: &mut Vec<OsString>, unpaired: &UnpairedArgs) {
    push_bool(args, "--reads-are-fragments", unpaired.reads_are_fragments);
}

pub(crate) fn push_ref_2bit_required(args: &mut Vec<OsString>, ref_genome: &Ref2BitRequiredArgs) {
    push_path(args, "--ref-2bit", &ref_genome.ref_2bit);
}

pub(crate) fn push_output_prefix(args: &mut Vec<OsString>, output_prefix: &str) {
    push_value(args, "--output-prefix", output_prefix);
}

pub(crate) fn push_fragment_lengths(
    args: &mut Vec<OsString>,
    fragment_lengths: &FragmentLengthArgs,
) {
    push_value(
        args,
        "--min-fragment-length",
        fragment_lengths.min_fragment_length,
    );
    push_value(
        args,
        "--max-fragment-length",
        fragment_lengths.max_fragment_length,
    );
}

pub(crate) fn push_chromosomes(args: &mut Vec<OsString>, chromosomes: &ChromosomeArgs) {
    if let Some(chromosomes) = &chromosomes.chromosomes {
        push_values(args, "--chromosomes", chromosomes);
    }
    push_optional_path(
        args,
        "--chromosomes-file",
        chromosomes.chromosomes_file.as_deref(),
    );
}

pub(crate) fn push_scale_genome(args: &mut Vec<OsString>, scale_genome: &ScaleGenomeArgs) {
    push_optional_path(
        args,
        "--scaling-factors",
        scale_genome.scaling_factors.as_deref(),
    );
}

#[cfg(any(
    feature = "cmd_fragment_kmers",
    feature = "cmd_transitions",
    feature = "cmd_wps",
    feature = "cmd_wps_peaks"
))]
pub(crate) fn push_windows(args: &mut Vec<OsString>, windows: &WindowsArgs) {
    if let Some(by_size) = windows.by_size {
        push_value(args, "--by-size", by_size);
    }
    push_optional_path(args, "--by-bed", windows.by_bed.as_deref());
}

pub(crate) fn push_distribution_windows(
    args: &mut Vec<OsString>,
    windows: &DistributionWindowsArgs,
) {
    if let Some(by_size) = windows.by_size {
        push_value(args, "--by-size", by_size);
    }
    push_optional_path(args, "--by-bed", windows.by_bed.as_deref());
    push_optional_path(args, "--by-grouped-bed", windows.by_grouped_bed.as_deref());
}

pub(crate) fn push_gc_windows(args: &mut Vec<OsString>, windows: &GCWindowsArgs) {
    if let Some(by_size) = windows.by_size {
        push_value(args, "--by-size", by_size);
    } else if let Some(by_bed) = &windows.by_bed {
        push_path(args, "--by-bed", by_bed);
    } else if windows.global {
        push_flag(args, "--global");
    }
}

pub(crate) fn push_assign_to_window(args: &mut Vec<OsString>, assignment: &AssignToWindowArgs) {
    push_value(
        args,
        "--assign-by",
        window_assigner_value(&assignment.assign_by),
    );
}

pub(crate) fn push_apply_gc(args: &mut Vec<OsString>, gc: &ApplyGCArgs) {
    push_optional_path(args, "--gc-file", gc.gc_file.as_deref());
    if let Some(gc_tag) = &gc.gc_tag {
        push_value(args, "--gc-tag", gc_tag);
    }
    push_bool(args, "--neutralize-invalid-gc", gc.neutralize_invalid_gc);
}

pub(crate) fn push_apply_gc_file_only(args: &mut Vec<OsString>, gc: &ApplyGCArgFileOnly) {
    push_optional_path(args, "--gc-file", gc.gc_file.as_deref());
    push_bool(args, "--neutralize-invalid-gc", gc.neutralize_invalid_gc);
}

pub(crate) fn push_blacklist_common(
    args: &mut Vec<OsString>,
    blacklist: Option<&[PathBuf]>,
    blacklist_min_size: u64,
    blacklist_strategy: &BlacklistStrategy,
) {
    push_path_values(args, "--blacklist", blacklist);
    push_value(args, "--blacklist-min-size", blacklist_min_size);
    push_value(
        args,
        "--blacklist-strategy",
        blacklist_strategy_value(blacklist_strategy),
    );
}

pub(crate) fn push_logging(args: &mut Vec<OsString>, logging: &LoggingArgs) {
    let value = match &logging.log {
        LogSpec::Stdout => "stdout".to_string(),
        LogSpec::Quiet => "quiet".to_string(),
        LogSpec::File(None) => "file".to_string(),
        LogSpec::File(Some(path)) => format!("file={}", path.display()),
    };
    push_value(args, "--log", value);
}

#[cfg(any(
    feature = "cmd_fragment_kmers",
    feature = "cmd_transitions",
    feature = "cmd_visualize_positions"
))]
pub(crate) fn push_position_selection(
    args: &mut Vec<OsString>,
    selection: &FragmentPositionSelectionArgs,
) {
    push_values(
        args,
        "--frame",
        &selection
            .frame
            .iter()
            .map(|frame| frame.as_str())
            .collect::<Vec<_>>(),
    );
    push_values(args, "--positions", &selection.positions);
    push_values(args, "--step", &selection.step);
}

#[cfg(any(
    feature = "cmd_fragment_kmers",
    feature = "cmd_transitions",
    feature = "cmd_visualize_positions"
))]
pub(crate) fn push_base_selection(args: &mut Vec<OsString>, selection: &BaseSelectionArgs) {
    push_value(args, "--bases-from", selection.bases_from.as_str());
    push_value(
        args,
        "--mismatch-bases-from",
        selection.mismatch_bases_from.as_str(),
    );
}

pub(crate) fn indel_mode_value(mode: IndelMode) -> &'static str {
    match mode {
        IndelMode::Ignore => "ignore",
        IndelMode::Adjust => "adjust",
        IndelMode::Skip => "skip",
    }
}

pub(crate) fn clip_mode_value(mode: ClipMode) -> &'static str {
    match mode {
        ClipMode::Aligned => "aligned",
        ClipMode::Adjust => "adjust",
        ClipMode::Skip => "skip",
    }
}

pub(crate) fn indel_motif_filter_value(filter: IndelMotifFilterPolicy) -> &'static str {
    match filter {
        IndelMotifFilterPolicy::Auto => "auto",
        IndelMotifFilterPolicy::SkipAffectedEnd => "skip-affected-end",
        IndelMotifFilterPolicy::SkipAffectedFragment => "skip-affected-fragment",
    }
}

pub(crate) fn blacklist_strategy_value(strategy: &BlacklistStrategy) -> String {
    match strategy {
        BlacklistStrategy::Any => "any".to_string(),
        BlacklistStrategy::All => "all".to_string(),
        BlacklistStrategy::Midpoint => "midpoint".to_string(),
        BlacklistStrategy::Proportion(threshold) => format!("proportion={threshold}"),
    }
}

pub(crate) fn window_assigner_value(assigner: &WindowAssigner) -> String {
    match assigner {
        WindowAssigner::CountOverlap => "count-overlap".to_string(),
        WindowAssigner::Any => "any".to_string(),
        WindowAssigner::All => "all".to_string(),
        WindowAssigner::Midpoint => "midpoint".to_string(),
        WindowAssigner::Proportion(threshold) => format!("proportion={threshold}"),
    }
}
