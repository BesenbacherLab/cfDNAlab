use super::*;
use crate::shared::kmers::{
    kmer_codec::KmerSpec, process_counts::postprocess_selected_motif_counts,
};
use fxhash::FxHashMap;
use std::io::Write;

type SelectedCountsByWindow = FxHashMap<u64, FxHashMap<u32, f64>>;

fn store_positive_weight(weight: f64) -> anyhow::Result<bool> {
    anyhow::ensure!(weight >= 0.0, "negative weight {weight}");
    Ok(weight > 0.0)
}

fn allow_dense_output_size(_total_windows: usize, _total_motifs: usize) -> anyhow::Result<()> {
    Ok(())
}

fn target_labels(lookup: &SelectedMotifLookup) -> Vec<&str> {
    lookup.labels.iter().map(String::as_str).collect()
}

fn write_temp_motifs_file(contents: &str) -> anyhow::Result<tempfile::NamedTempFile> {
    let mut file = tempfile::NamedTempFile::new()?;
    write!(file, "{contents}")?;
    Ok(file)
}

fn spec_for(k: usize, label: &str) -> KmerSpec {
    build_optional_kmer_spec(k, label)
        .expect("valid k-mer size")
        .expect("non-zero k-mer size")
}

#[cfg(feature = "cmd_ref_kmers")]
fn ref_spec_for(kmer_size: usize) -> KmerSpec {
    spec_for(kmer_size, "kmer")
}

#[cfg(feature = "cmd_ref_kmers")]
fn ref_lookup_key(code: u64) -> EncodedMotifKey {
    EncodedMotifKey {
        inside_code: code,
        outside_code: 0,
        reverse_on_decode: false,
    }
}

#[cfg(feature = "cmd_ref_kmers")]
fn ref_reverse_lookup_key(code: u64) -> EncodedMotifKey {
    EncodedMotifKey {
        inside_code: code,
        outside_code: 0,
        reverse_on_decode: true,
    }
}

#[test]
fn parses_ungrouped_combined_motifs_into_ordered_targets_and_encoded_lookup() -> anyhow::Result<()> {
    // Arrange: two 2+2 motifs define two output motif targets in file order.
    let contents = "AC_GT\nTT_AA\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;

    // Assert: targets use the public combined motif labels directly.
    assert_eq!(lookup.column_kind, SelectedMotifColumnKind::Motif);
    assert_eq!(target_labels(&lookup), vec!["AC_GT", "TT_AA"]);

    let inside_spec = spec_for(2, "inside");
    let outside_spec = spec_for(2, "outside");

    // Assert: left-end storage keeps outside || inside without reverse-on-decode.
    let left_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"GT"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: false,
    };
    assert_eq!(lookup.target_for(left_key), Some(0));

    // Assert: right-end storage uses revcomp(inside) || revcomp(outside) and keeps
    // the reverse-on-decode state in the lookup key.
    let right_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"AC"),
        outside_code: outside_spec.encode_kmer_bytes(b"GT"),
        reverse_on_decode: true,
    };
    assert_eq!(lookup.target_for(right_key), Some(0));

    Ok(())
}

#[test]
fn parses_grouped_motifs_into_alphabetic_group_targets() -> anyhow::Result<()> {
    // Arrange: the first and third motifs share a group, while the middle motif creates
    // a second target. The dot and hyphen in group names are intentionally allowed.
    // Group labels are sorted alphabetically, so `group-one` comes before `group.one`.
    let contents = "AC_GT\tgroup.one\nTT_AA\tgroup-one\nGC_TA\tgroup.one\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;

    // Assert
    assert_eq!(lookup.column_kind, SelectedMotifColumnKind::MotifGroup);
    assert_eq!(target_labels(&lookup), vec!["group-one", "group.one"]);

    let inside_spec = spec_for(2, "inside");
    let outside_spec = spec_for(2, "outside");
    let first_group_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"GT"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: false,
    };
    let reused_group_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"TA"),
        outside_code: outside_spec.encode_kmer_bytes(b"GC"),
        reverse_on_decode: false,
    };
    let second_group_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"AA"),
        outside_code: outside_spec.encode_kmer_bytes(b"TT"),
        reverse_on_decode: false,
    };

    assert_eq!(lookup.target_for(first_group_key), Some(1));
    assert_eq!(lookup.target_for(reused_group_key), Some(1));
    assert_eq!(lookup.target_for(second_group_key), Some(0));

    Ok(())
}

#[test]
fn keeps_reverse_on_decode_state_distinct_when_mapping_groups() -> anyhow::Result<()> {
    // Arrange: these two final motifs intentionally share the same inside/outside codes
    // in different end states:
    // - right-end AC_GT stores inside=AC, outside=GT, reverse_on_decode=true
    // - left-end GT_AC stores inside=AC, outside=GT, reverse_on_decode=false
    let contents = "AC_GT\tgroup_a\nGT_AC\tgroup_b\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["group_a", "group_b"]);

    let inside_spec = spec_for(2, "inside");
    let outside_spec = spec_for(2, "outside");
    let shared_inside_code = inside_spec.encode_kmer_bytes(b"AC");
    let shared_outside_code = outside_spec.encode_kmer_bytes(b"GT");

    let right_state_for_first_row = EncodedMotifKey {
        inside_code: shared_inside_code,
        outside_code: shared_outside_code,
        reverse_on_decode: true,
    };
    let left_state_for_second_row = EncodedMotifKey {
        inside_code: shared_inside_code,
        outside_code: shared_outside_code,
        reverse_on_decode: false,
    };

    assert_eq!(lookup.target_for(right_state_for_first_row), Some(0));
    assert_eq!(lookup.target_for(left_state_for_second_row), Some(1));

    Ok(())
}

#[test]
fn parses_inside_only_motifs_with_omitted_or_leading_separator() -> anyhow::Result<()> {
    // Arrange: with no outside bases, users may write either `AC` or `_GT`.
    let contents = "AC\n_GT\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 0,
        },
    )?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["_AC", "_GT"]);

    let inside_spec = spec_for(2, "inside");
    let first_left_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"AC"),
        outside_code: 0,
        reverse_on_decode: false,
    };
    let first_right_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"GT"),
        outside_code: 0,
        reverse_on_decode: true,
    };
    assert_eq!(lookup.target_for(first_left_key), Some(0));
    assert_eq!(lookup.target_for(first_right_key), Some(0));

    Ok(())
}

#[test]
fn parses_outside_only_motifs_with_omitted_or_trailing_separator() -> anyhow::Result<()> {
    // Arrange: with no inside bases, users may write either `AC` or `GT_`.
    let contents = "AC\nGT_\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 0,
            k_outside: 2,
        },
    )?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["AC_", "GT_"]);

    let outside_spec = spec_for(2, "outside");
    let first_left_key = EncodedMotifKey {
        inside_code: 0,
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: false,
    };
    let first_right_key = EncodedMotifKey {
        inside_code: 0,
        outside_code: outside_spec.encode_kmer_bytes(b"GT"),
        reverse_on_decode: true,
    };
    assert_eq!(lookup.target_for(first_left_key), Some(0));
    assert_eq!(lookup.target_for(first_right_key), Some(0));

    Ok(())
}

#[test]
fn normalizes_lowercase_motif_bases_to_uppercase_output_labels() -> anyhow::Result<()> {
    // Arrange
    let contents = "ac_gt\tgroup_1\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["group_1"]);

    let inside_spec = spec_for(2, "inside");
    let outside_spec = spec_for(2, "outside");
    let left_key = EncodedMotifKey {
        inside_code: inside_spec.encode_kmer_bytes(b"GT"),
        outside_code: outside_spec.encode_kmer_bytes(b"AC"),
        reverse_on_decode: false,
    };
    assert_eq!(lookup.target_for(left_key), Some(0));

    Ok(())
}

#[test]
fn reads_motifs_file_from_path() -> anyhow::Result<()> {
    // Arrange
    let file = write_temp_motifs_file("AC_GT\n")?;

    // Act
    let lookup = parse_selected_end_motifs_file(file.path(), 2, 2)?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["AC_GT"]);

    Ok(())
}

#[test]
fn parses_crlf_terminated_motifs_file_rows() -> anyhow::Result<()> {
    // Arrange: Windows-authored TSV files commonly keep `\r` before `\n`. The parser should treat
    // that as the line ending, not as part of the motif or group label.
    let contents = "AC_GT\tgroup_one\r\nTT_AA\tgroup_two\r\n";

    // Act
    let lookup = parse_selected_motifs(
        contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;

    // Assert
    assert_eq!(lookup.column_kind, SelectedMotifColumnKind::MotifGroup);
    assert_eq!(target_labels(&lookup), vec!["group_one", "group_two"]);

    Ok(())
}

#[test]
fn shares_one_selected_subspace_when_inside_and_outside_large_k_match() -> anyhow::Result<()> {
    // Arrange: both motif halves use k=30, so they can share one selected half-code subset. The
    // full inside/outside motif pair is still filtered by the encoded lookup after each half has
    // been encoded.
    let contents = format!("{}_{}\n", "C".repeat(30), "A".repeat(30));

    // Act
    let lookup = parse_selected_motifs(
        &contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 30,
            k_outside: 30,
        },
    )?;

    // Assert
    match (&lookup.inside_spec, &lookup.outside_spec) {
        (
            Some(SelectedMotifHalfSpec::Subspace(inside_spec)),
            Some(SelectedMotifHalfSpec::Subspace(outside_spec)),
        ) => assert!(
            std::sync::Arc::ptr_eq(inside_spec, outside_spec),
            "inside and outside should share the same selected subspace"
        ),
        other => panic!("expected shared subspace specs, got {other:?}"),
    }

    Ok(())
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn parses_ref_kmers_with_forward_and_reverse_observation_lookup_states() -> anyhow::Result<()> {
    // Arrange: a single `_` can be used to reuse end-motif-style files. AACC is deliberately not
    // self reverse-complementary, so its reverse observation state is encoded as GGTT.
    let contents = "AA_CC\nTTAA\n";

    // Act
    let lookup =
        parse_selected_motifs(contents, SelectedMotifsFileKind::RefKmers { kmer_size: 4 })?;

    // Assert
    assert_eq!(lookup.column_kind, SelectedMotifColumnKind::Motif);
    assert_eq!(target_labels(&lookup), vec!["AACC", "TTAA"]);
    assert!(lookup.outside_spec.is_none());

    let spec = ref_spec_for(4);
    assert_eq!(
        lookup.target_for(ref_lookup_key(spec.encode_kmer_bytes(b"AACC"))),
        Some(0)
    );
    assert_eq!(
        lookup.target_for(ref_lookup_key(spec.encode_kmer_bytes(b"TTAA"))),
        Some(1)
    );

    assert_eq!(
        lookup.target_for(ref_reverse_lookup_key(spec.encode_kmer_bytes(b"GGTT"))),
        Some(0)
    );
    assert_eq!(
        lookup.target_for(ref_reverse_lookup_key(spec.encode_kmer_bytes(b"AACC"))),
        None
    );

    Ok(())
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn parses_grouped_ref_kmers_into_alphabetic_group_targets() -> anyhow::Result<()> {
    // Arrange: target order is alphabetic by group name. Later rows using an existing group
    // contribute to that already-created target.
    let contents = "AC_GT\tgroup.one\nTTAA\tgroup-one\nGG_CC\tgroup.one\n";

    // Act
    let lookup =
        parse_selected_motifs(contents, SelectedMotifsFileKind::RefKmers { kmer_size: 4 })?;

    // Assert
    assert_eq!(lookup.column_kind, SelectedMotifColumnKind::MotifGroup);
    assert_eq!(target_labels(&lookup), vec!["group-one", "group.one"]);

    let spec = ref_spec_for(4);
    assert_eq!(
        lookup.target_for(ref_lookup_key(spec.encode_kmer_bytes(b"ACGT"))),
        Some(1)
    );
    assert_eq!(
        lookup.target_for(ref_lookup_key(spec.encode_kmer_bytes(b"GGCC"))),
        Some(1)
    );
    assert_eq!(
        lookup.target_for(ref_lookup_key(spec.encode_kmer_bytes(b"TTAA"))),
        Some(0)
    );

    Ok(())
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn normalizes_lowercase_ref_kmers_to_uppercase_labels() -> anyhow::Result<()> {
    // Arrange
    let contents = "ac_gt\tgroup_1\n";

    // Act
    let lookup =
        parse_selected_motifs(contents, SelectedMotifsFileKind::RefKmers { kmer_size: 4 })?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["group_1"]);
    let spec = ref_spec_for(4);
    assert_eq!(
        lookup.target_for(ref_lookup_key(spec.encode_kmer_bytes(b"ACGT"))),
        Some(0)
    );

    Ok(())
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn reads_ref_kmers_motifs_file_from_path() -> anyhow::Result<()> {
    // Arrange
    let file = write_temp_motifs_file("AC_GT\n")?;

    // Act
    let lookup = parse_selected_ref_kmers_file(file.path(), 4)?;

    // Assert
    assert_eq!(target_labels(&lookup), vec!["ACGT"]);
    Ok(())
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn parses_crlf_terminated_ref_kmers_motifs_file_rows() -> anyhow::Result<()> {
    // Arrange: Windows-authored TSV files commonly keep `\r` before `\n`. The parser should treat
    // that as the line ending, not as part of the motif or group label.
    let contents = "ACGT\tgroup_one\r\nTTAA\tgroup_two\r\n";

    // Act
    let lookup =
        parse_selected_motifs(contents, SelectedMotifsFileKind::RefKmers { kmer_size: 4 })?;

    // Assert
    assert_eq!(lookup.column_kind, SelectedMotifColumnKind::MotifGroup);
    assert_eq!(target_labels(&lookup), vec!["group_one", "group_two"]);
    Ok(())
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn uses_selected_subspace_for_large_ref_kmers() -> anyhow::Result<()> {
    // Arrange: k = 30 cannot use the full radix-5 code space, so motifs-file counting must use a
    // byte-backed selected subspace. Only A^30 is listed, but T^30 must also be encoded because it
    // is the reverse complement observed in reference coordinates.
    let motif = "A".repeat(30);
    let reverse_motif = "T".repeat(30);
    let contents = format!("{motif}\n");

    // Act
    let lookup = parse_selected_motifs(
        &contents,
        SelectedMotifsFileKind::RefKmers { kmer_size: 30 },
    )?;

    // Assert
    let Some(SelectedMotifHalfSpec::Subspace(spec)) = lookup.inside_spec.as_ref() else {
        panic!("expected selected subspace for k = 30");
    };
    assert_eq!(spec.k, 30);
    let forward_code = spec.encode_kmer_bytes(motif.as_bytes());
    let reverse_code = spec.encode_kmer_bytes(reverse_motif.as_bytes());
    assert_eq!(lookup.target_for(ref_lookup_key(forward_code)), Some(0));
    assert_eq!(
        lookup.target_for(ref_reverse_lookup_key(reverse_code)),
        Some(0)
    );
    assert_eq!(
        spec.build_left_aligned_codes(reverse_motif.as_bytes()).get(0),
        reverse_code
    );
    Ok(())
}

#[test]
fn postprocesses_selected_counts_in_file_order_for_observed_targets() -> anyhow::Result<()> {
    // Arrange: target 2 is inserted before target 0 to make sure post-processing
    // compacts observed targets in file order, not hash-map insertion order.
    let lookup = parse_selected_motifs(
        "AC_GT\nTT_AA\nGC_TA\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;
    let mut counts_by_window = SelectedCountsByWindow::default();
    counts_by_window.entry(1).or_default().insert(2, 2.0);
    counts_by_window.entry(0).or_default().insert(0, 1.5);

    // Act
    let (all_bins, motif_order) = postprocess_selected_motif_counts(
        counts_by_window,
        3,
        &lookup.labels,
        false,
        store_positive_weight,
        allow_dense_output_size,
    )?;

    // Assert
    assert_eq!(motif_order, vec!["AC_GT", "GC_TA"]);
    assert_eq!(all_bins.len(), 3);
    assert_eq!(all_bins[0].get(&0), Some(&1.5));
    assert_eq!(all_bins[1].get(&1), Some(&2.0));
    assert!(all_bins[2].is_empty());

    Ok(())
}

#[test]
fn postprocesses_selected_all_motifs_with_unobserved_targets() -> anyhow::Result<()> {
    // Arrange
    let lookup = parse_selected_motifs(
        "AC_GT\nTT_AA\nGC_TA\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )?;
    let mut counts_by_window = SelectedCountsByWindow::default();
    counts_by_window.entry(0).or_default().insert(1, 2.0);

    // Act
    let (all_bins, motif_order) = postprocess_selected_motif_counts(
        counts_by_window,
        2,
        &lookup.labels,
        true,
        store_positive_weight,
        allow_dense_output_size,
    )?;

    // Assert
    assert_eq!(motif_order, vec!["AC_GT", "TT_AA", "GC_TA"]);
    assert_eq!(all_bins[0].get(&1), Some(&2.0));
    assert!(!all_bins[0].contains_key(&0));
    assert!(!all_bins[0].contains_key(&2));

    Ok(())
}

#[test]
fn postprocess_rejects_out_of_bounds_selected_window_idx() {
    // Arrange
    let lookup = parse_selected_motifs(
        "AC_GT\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect("valid lookup");
    let mut counts_by_window = SelectedCountsByWindow::default();
    counts_by_window.entry(2).or_default().insert(0, 1.0);

    // Act
    let error = postprocess_selected_motif_counts(
        counts_by_window,
        2,
        &lookup.labels,
        false,
        store_positive_weight,
        allow_dense_output_size,
    )
    .expect_err("out-of-bounds window index should fail");

    // Assert
    assert!(
        error.to_string().contains("out of bounds"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_empty_file() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("empty motifs file should fail");

    // Assert
    assert!(
        error.to_string().contains("at least one motif"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_blank_line_as_empty_motif() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_GT\n\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("blank motif line should fail");

    // Assert
    assert!(
        error.to_string().contains("line 2: motif label is empty"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_rows_with_more_than_two_columns() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_GT\tgroup\textra\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("three-column row should fail");

    // Assert
    assert!(
        error.to_string().contains("expected one or two"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_mixed_grouped_and_ungrouped_rows() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_GT\nTT_AA\tgroup\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("mixed column modes should fail");

    // Assert
    assert!(
        error.to_string().contains("one column for every row"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_empty_group_name() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_GT\t\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("empty group should fail");

    // Assert
    assert!(
        error.to_string().contains("group name is empty"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_group_names_with_whitespace() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_GT\tgroup one\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("group names with spaces should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid character"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_group_names_with_disallowed_punctuation() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_GT\tgroup/one\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("group names with slashes should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid character"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_motifs_containing_n() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AN_GT\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("N-containing motifs should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid base `N`"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_motifs_with_invalid_base_characters() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AX_GT\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("invalid base should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid base `X`"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_wrong_outside_length() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "A_GT\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("short outside side should fail");

    // Assert
    assert!(
        error.to_string().contains("outside length 1, expected 2"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_wrong_inside_length() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_G\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("short inside side should fail");

    // Assert
    assert!(
        error.to_string().contains("inside length 1, expected 2"),
        "unexpected error: {error}"
    );
}

#[test]
fn requires_separator_when_inside_and_outside_are_non_zero() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 2,
        },
    )
    .expect_err("combined motifs need separator");

    // Assert
    assert!(
        error.to_string().contains("must use `<outside>_<inside>`"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_extra_separator() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "A_C_G\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 1,
            k_outside: 1,
        },
    )
    .expect_err("extra separator should fail");

    // Assert
    assert!(
        error.to_string().contains("more than one `_`"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_trailing_separator_for_inside_only_motif() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "AC_\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 0,
        },
    )
    .expect_err("inside-only trailing separator should fail");

    // Assert
    assert!(
        error.to_string().contains("outside length 2, expected 0"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_leading_separator_for_outside_only_motif() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "_AC\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 0,
            k_outside: 2,
        },
    )
    .expect_err("outside-only leading separator should fail");

    // Assert
    assert!(
        error.to_string().contains("outside length 0, expected 2"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_duplicate_motifs_after_normalization() {
    // Arrange + Act: `AC` and `_AC` are the same inside-only public motif.
    let error = parse_selected_motifs(
        "AC\n_AC\n",
        SelectedMotifsFileKind::EndMotifs {
            k_inside: 2,
            k_outside: 0,
        },
    )
    .expect_err("normalized duplicate should fail");

    // Assert
    assert!(
        error.to_string().contains("duplicate motif `_AC`"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_empty_ref_kmers_motifs_file() {
    // Arrange + Act
    let error = parse_selected_motifs("", SelectedMotifsFileKind::RefKmers { kmer_size: 4 })
        .expect_err("empty motifs file should fail");

    // Assert
    assert!(
        error.to_string().contains("at least one k-mer"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_blank_ref_kmers_motifs_file_line_as_empty_motif() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\n\n",
        SelectedMotifsFileKind::RefKmers { kmer_size: 4 },
    )
    .expect_err("blank motif line should fail");

    // Assert
    assert!(
        error.to_string().contains("line 2: k-mer label is empty"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_ref_kmers_rows_with_more_than_two_columns() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\tgroup\textra\n",
        SelectedMotifsFileKind::RefKmers { kmer_size: 4 },
    )
    .expect_err("three-column row should fail");

    // Assert
    assert!(
        error.to_string().contains("expected one or two"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_mixed_grouped_and_ungrouped_ref_kmers_rows() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\nTTAA\tgroup\n",
        SelectedMotifsFileKind::RefKmers { kmer_size: 4 },
    )
    .expect_err("mixed column modes should fail");

    // Assert
    assert!(
        error.to_string().contains("one column for every row"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_empty_ref_kmers_group_name() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\t\n",
        SelectedMotifsFileKind::RefKmers { kmer_size: 4 },
    )
    .expect_err("empty group should fail");

    // Assert
    assert!(
        error.to_string().contains("group name is empty"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_ref_kmers_group_names_with_disallowed_punctuation() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\tgroup/one\n",
        SelectedMotifsFileKind::RefKmers { kmer_size: 4 },
    )
    .expect_err("group names with slashes should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid character"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_ref_kmers_containing_n_after_separator_collapse() {
    // Arrange + Act
    let error = parse_selected_motifs("AC_GN\n", SelectedMotifsFileKind::RefKmers { kmer_size: 4 })
        .expect_err("N-containing k-mer should fail");

    // Assert
    assert!(
        error.to_string().contains("invalid base `N`"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_wrong_ref_kmer_length_after_separator_collapse() {
    // Arrange + Act
    let error = parse_selected_motifs("AC_G\n", SelectedMotifsFileKind::RefKmers { kmer_size: 4 })
        .expect_err("short collapsed k-mer should fail");

    // Assert
    assert!(
        error.to_string().contains("length 3, expected 4"),
        "unexpected error: {error}"
    );
}

#[cfg(feature = "cmd_ref_kmers")]
#[test]
fn rejects_duplicate_ref_kmers_after_case_and_separator_normalization() {
    // Arrange + Act
    let error = parse_selected_motifs(
        "ACGT\nac_gt\n",
        SelectedMotifsFileKind::RefKmers { kmer_size: 4 },
    )
    .expect_err("normalized duplicate should fail");

    // Assert
    assert!(
        error.to_string().contains("duplicate k-mer `ACGT`"),
        "unexpected error: {error}"
    );
}

#[test]
fn rejects_encoded_key_conflicts_between_targets() {
    // Arrange: direct helper coverage for the defensive conflict check. Normal parsing
    // rejects duplicate motif rows before this path can usually be reached.
    let mut lookup = FxHashMap::default();
    let key = EncodedMotifKey {
        inside_code: 10,
        outside_code: 20,
        reverse_on_decode: true,
    };
    insert_lookup_key(&mut lookup, key, 0, "AC_GT", 1).expect("initial key insert should succeed");

    // Act
    let error = insert_lookup_key(&mut lookup, key, 1, "GT_AC", 2)
        .expect_err("conflicting target should fail");

    // Assert
    assert!(
        error.to_string().contains("already assigned"),
        "unexpected error: {error}"
    );
}
