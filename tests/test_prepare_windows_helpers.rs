#![cfg(feature = "cmd_prepare_windows")]

mod tests_prepare_windows_helpers {
    use anyhow::Result;
    use cfdnalab::commands::prepare_windows::chunk::apply_near_annotations;
    use cfdnalab::commands::prepare_windows::config::ComposeSpec;
    use cfdnalab::commands::prepare_windows::config::MergeScope;
    use cfdnalab::commands::prepare_windows::config::{CoordinateSet, PrepareConfig};
    use cfdnalab::commands::prepare_windows::filters::{
        MinPerKeyRuleState, MinPerWindowFilterData, collect_min_per_window_filter_data,
        normalize_min_per_rules, parse_exclude_rules, parse_min_per_rules, validate_available_keys,
        validate_compositions_available,
    };
    use cfdnalab::commands::prepare_windows::labels::{
        AtomicLabelPart, LabelKey, LabelSchema, LabelTuple, NO_NEAR_BIN_LABEL, NO_NEAR_LABEL,
        build_tuple_compositions,
    };
    use cfdnalab::commands::prepare_windows::near_file::{
        NearDuplicatesPolicy, NearIndex, Strand, load_near_index,
    };
    use cfdnalab::commands::prepare_windows::parsers::parse_distance_bins;
    use cfdnalab::commands::prepare_windows::postprocess::partition_safe_and_tail;
    use cfdnalab::commands::prepare_windows::prepare_windows::Window;
    use fxhash::FxHashSet;
    use std::io::Write;
    use std::sync::Arc;
    use tempfile::NamedTempFile;

    fn build_schema(specs: &[&str]) -> Result<LabelSchema> {
        let mut compose_specs: Vec<ComposeSpec> = Vec::with_capacity(specs.len());
        for spec in specs {
            let compose_spec = spec.parse::<ComposeSpec>().map_err(anyhow::Error::msg)?;
            compose_specs.push(compose_spec);
        }
        Ok(LabelSchema::new(&compose_specs)?)
    }

    fn all_available_parts() -> FxHashSet<AtomicLabelPart> {
        let mut parts: FxHashSet<AtomicLabelPart> = FxHashSet::default();
        parts.insert(AtomicLabelPart::Input);
        parts.insert(AtomicLabelPart::NearWindowSide);
        parts.insert(AtomicLabelPart::NearName);
        parts.insert(AtomicLabelPart::Bin);
        parts.insert(AtomicLabelPart::Cluster);
        parts
    }

    fn input_only_parts() -> FxHashSet<AtomicLabelPart> {
        let mut parts: FxHashSet<AtomicLabelPart> = FxHashSet::default();
        parts.insert(AtomicLabelPart::Input);
        parts
    }

    fn build_tuple(input: &str, bin: Option<&str>) -> LabelTuple {
        let mut tuple = LabelTuple::new(input.to_string());
        tuple.bin = bin.map(|value| value.to_string());
        tuple
    }

    fn build_window(chrom: &str, start: u32, end: u32, group_key: &str) -> Window {
        Window::from_bounds(
            Arc::from(chrom),
            start,
            end,
            start,
            end,
            Vec::new(),
            group_key.to_string(),
            None,
        )
        .expect("test window should be valid")
    }

    fn build_window_with_tuple(chrom: &str, start: u32, end: u32) -> Window {
        Window::from_bounds(
            Arc::from(chrom),
            start,
            end,
            start,
            end,
            vec![LabelTuple::new("A".to_string())],
            "A".to_string(),
            None,
        )
        .expect("test window should be valid")
    }

    fn assert_min_per_window_data(
        actual: &MinPerWindowFilterData,
        expected_before: &[Vec<String>],
        expected_after: &[Vec<String>],
        expected_kept: &[LabelTuple],
    ) {
        assert_eq!(actual.values_before_filter, expected_before);
        assert_eq!(actual.values_after_filter, expected_after);
        assert_eq!(actual.kept_tuples, expected_kept);
    }

    #[test]
    fn should_parse_compose_specs_and_resolve_keys() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin", "report=core,win-direction"])?;

        // Act
        let core_key = schema.resolve_key("core")?;
        let report_key = schema.resolve_key("report")?;

        // Assert
        assert!(matches!(core_key, LabelKey::Composition(0)));
        assert!(matches!(report_key, LabelKey::Composition(1)));
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_input_keyword() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["input=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_none_keyword() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["none=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_near_side_keyword() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["win-direction=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_near_name_keyword() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["near-name=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_bin_keyword() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["bin=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_cluster_keyword() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["cluster=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_parse_exclude_rules_for_atomic_and_composition_keys() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin"])?;
        let available_parts = all_available_parts();
        let specs = vec!["input=A".to_string(), "core=A.B".to_string()];

        // Act
        let rules = parse_exclude_rules(&specs, &schema, &available_parts)?;

        // Assert
        assert_eq!(rules.len(), 2);
        assert!(matches!(
            rules[0].key,
            LabelKey::Atomic(AtomicLabelPart::Input)
        ));
        assert!(matches!(rules[1].key, LabelKey::Composition(0)));
        assert_eq!(rules[0].term, "A");
        assert_eq!(rules[1].term, "A.B");
        Ok(())
    }

    #[test]
    fn should_reject_unknown_exclude_keys() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin"])?;
        let available_parts = all_available_parts();
        let specs = vec!["unknown=Z".to_string()];

        // Act
        let result = parse_exclude_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_exclude_rules_with_unavailable_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[])?;
        let available_parts = input_only_parts();
        let specs = vec!["win-direction=+".to_string()];

        // Act
        let result = parse_exclude_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_exclude_rules_with_unavailable_composition_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,win-direction"])?;
        let available_parts = input_only_parts();
        let specs = vec!["core=A.B".to_string()];

        // Act
        let result = parse_exclude_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_compose_with_unavailable_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,win-direction"])?;
        let available_parts = input_only_parts();

        // Act
        let result = validate_compositions_available(&schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_collect_min_per_values_when_rejections_and_missing_values_present() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let input_state = MinPerKeyRuleState::new(LabelKey::Atomic(AtomicLabelPart::Input), 1);
        let mut bin_state = MinPerKeyRuleState::new(LabelKey::Atomic(AtomicLabelPart::Bin), 1);
        bin_state.add_rejected_value("drop");
        let states = vec![input_state, bin_state];

        let tuple_one = build_tuple("A", Some("keep"));
        let tuple_two = build_tuple("B", Some("drop"));
        let tuple_three = build_tuple("A", None);
        let tuple_four = build_tuple("A", Some("keep"));
        let tuple_five = build_tuple("C", Some("keep"));
        let label_tuples = vec![
            tuple_one.clone(),
            tuple_two.clone(),
            tuple_three.clone(),
            tuple_four.clone(),
            tuple_five.clone(),
        ];

        let tuple_compositions: Vec<Vec<String>> = Vec::new();

        // Act
        let result =
            collect_min_per_window_filter_data(&label_tuples, &tuple_compositions, &states);

        // Assert
        let expected_before = vec![
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            vec!["drop".to_string(), "keep".to_string()],
        ];
        let expected_after = vec![
            vec!["A".to_string(), "C".to_string()],
            vec!["keep".to_string()],
        ];
        let expected_kept = vec![tuple_one, tuple_four, tuple_five];
        assert_min_per_window_data(&result, &expected_before, &expected_after, &expected_kept);
        Ok(())
    }

    #[test]
    fn should_collect_composition_values_when_composition_rule_present() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin"])?;
        let mut composition_state = MinPerKeyRuleState::new(LabelKey::Composition(0), 1);
        composition_state.add_rejected_value("A.y");
        let states = vec![composition_state];

        // core=input,bin yields values like "A.x" by joining input and bin with "."
        let tuple_one = build_tuple("A", Some("x"));
        let tuple_two = build_tuple("A", Some("y"));
        let tuple_three = build_tuple("B", Some("x"));
        let label_tuples = vec![tuple_one.clone(), tuple_two.clone(), tuple_three.clone()];
        let tuple_compositions = build_tuple_compositions(&label_tuples, &schema);

        // Act
        let result =
            collect_min_per_window_filter_data(&label_tuples, &tuple_compositions, &states);

        // Assert
        // Values are in core composition order, so input "A" with bin "y" becomes "A.y"
        let expected_before = vec![vec![
            "A.x".to_string(),
            "A.y".to_string(),
            "B.x".to_string(),
        ]];
        let expected_after = vec![vec!["A.x".to_string(), "B.x".to_string()]];
        let expected_kept = vec![tuple_one, tuple_three];
        assert_min_per_window_data(&result, &expected_before, &expected_after, &expected_kept);
        Ok(())
    }

    #[test]
    fn should_return_empty_values_when_window_has_no_tuples() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let input_state = MinPerKeyRuleState::new(LabelKey::Atomic(AtomicLabelPart::Input), 1);
        let bin_state = MinPerKeyRuleState::new(LabelKey::Atomic(AtomicLabelPart::Bin), 1);
        let states = vec![input_state, bin_state];
        let label_tuples: Vec<LabelTuple> = Vec::new();
        let tuple_compositions: Vec<Vec<String>> = Vec::new();

        // Act
        let result =
            collect_min_per_window_filter_data(&label_tuples, &tuple_compositions, &states);

        // Assert
        let expected_before = vec![Vec::new(), Vec::new()];
        let expected_after = vec![Vec::new(), Vec::new()];
        let expected_kept: Vec<LabelTuple> = Vec::new();
        assert_min_per_window_data(&result, &expected_before, &expected_after, &expected_kept);
        Ok(())
    }

    #[test]
    fn should_reject_exclude_rules_with_wrong_composition_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin"])?;
        let available_parts = all_available_parts();
        let specs = vec!["core=B".to_string()];

        // Act
        let result = parse_exclude_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_exclude_rules_with_wrong_nested_composition_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin", "report=core,win-direction"])?;
        let available_parts = all_available_parts();
        let specs = vec!["report=A.B".to_string()];

        // Act
        let result = parse_exclude_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_accept_exclude_rules_with_correct_nested_composition_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin", "report=core,win-direction"])?;
        let available_parts = all_available_parts();
        let specs = vec!["report=A.B.C".to_string()];

        // Act
        let rules = parse_exclude_rules(&specs, &schema, &available_parts)?;

        // Assert
        assert_eq!(rules.len(), 1);
        assert!(matches!(rules[0].key, LabelKey::Composition(1)));
        assert_eq!(rules[0].term, "A.B.C");
        Ok(())
    }

    #[test]
    fn should_parse_min_per_rules_for_atomic_and_composition_keys() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin"])?;
        let available_parts = all_available_parts();
        let specs = vec!["input=10".to_string(), "core=5".to_string()];

        // Act
        let rules = parse_min_per_rules(&specs, &schema, &available_parts)?;

        // Assert
        assert_eq!(rules.len(), 2);
        assert!(matches!(
            rules[0].key,
            LabelKey::Atomic(AtomicLabelPart::Input)
        ));
        assert!(matches!(rules[1].key, LabelKey::Composition(0)));
        assert_eq!(rules[0].min_count, 10);
        assert_eq!(rules[1].min_count, 5);
        Ok(())
    }

    #[test]
    fn should_reject_non_numeric_min_per_counts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[])?;
        let available_parts = all_available_parts();
        let specs = vec!["input=bad".to_string()];

        // Act
        let result = parse_min_per_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_min_per_rules_with_unavailable_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[])?;
        let available_parts = input_only_parts();
        let specs = vec!["near-name=10".to_string()];

        // Act
        let result = parse_min_per_rules(&specs, &schema, &available_parts);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_out_labels_with_unavailable_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[])?;
        let available_parts = input_only_parts();
        let out_labels = schema.resolve_keys(&["near-name".to_string()])?;

        // Act
        let result = validate_available_keys(&out_labels, &schema, &available_parts, "out-labels");

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_reject_merge_key_with_unavailable_parts() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[])?;
        let available_parts = input_only_parts();
        let merge_key = schema.resolve_key("win-direction")?;

        // Act
        let result = validate_available_keys(
            std::slice::from_ref(&merge_key),
            &schema,
            &available_parts,
            "merge-key",
        );

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_normalize_min_per_rules_by_key_and_zero_minimum() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[])?;
        let available_parts = all_available_parts();
        let specs = vec![
            "input=0".to_string(),
            "input=5".to_string(),
            "input=10".to_string(),
        ];
        let rules = parse_min_per_rules(&specs, &schema, &available_parts)?;

        // Act
        let normalized = normalize_min_per_rules(&rules, &schema);

        // Assert
        assert_eq!(normalized.len(), 1);
        assert!(matches!(
            normalized[0].key,
            LabelKey::Atomic(AtomicLabelPart::Input)
        ));
        assert_eq!(normalized[0].min_count, 10);
        Ok(())
    }

    #[test]
    fn should_normalize_min_per_rules_by_membership_ignoring_order() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["core=input,bin", "swap=bin,input"])?;
        let available_parts = all_available_parts();
        let specs = vec!["core=100".to_string(), "swap=200".to_string()];
        let rules = parse_min_per_rules(&specs, &schema, &available_parts)?;

        // Act
        let normalized = normalize_min_per_rules(&rules, &schema);

        // Assert
        assert_eq!(normalized.len(), 1);
        assert!(matches!(normalized[0].key, LabelKey::Composition(0)));
        assert_eq!(normalized[0].min_count, 200);
        Ok(())
    }

    #[test]
    fn should_normalize_min_per_rules_by_membership_across_atomic_and_composition() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&["inputonly=input"])?;
        let available_parts = all_available_parts();
        let specs = vec!["input=50".to_string(), "inputonly=75".to_string()];
        let rules = parse_min_per_rules(&specs, &schema, &available_parts)?;

        // Act
        let normalized = normalize_min_per_rules(&rules, &schema);

        // Assert
        assert_eq!(normalized.len(), 1);
        assert!(matches!(
            normalized[0].key,
            LabelKey::Atomic(AtomicLabelPart::Input)
        ));
        assert_eq!(normalized[0].min_count, 75);
        Ok(())
    }

    #[test]
    fn should_reject_compose_name_with_underscore() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let specs = ["input_only=input"];

        // Act
        let result = build_schema(&specs);

        // Assert
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn should_normalize_min_per_rules_by_membership_for_nested_compositions() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        let schema = build_schema(&[
            "core=input,bin",
            "report=core,win-direction",
            "swap=win-direction,bin,input",
        ])?;
        let available_parts = all_available_parts();
        let specs = vec!["report=120".to_string(), "swap=250".to_string()];
        let rules = parse_min_per_rules(&specs, &schema, &available_parts)?;

        // Act
        let normalized = normalize_min_per_rules(&rules, &schema);

        // Assert
        assert_eq!(normalized.len(), 1);
        assert!(matches!(normalized[0].key, LabelKey::Composition(1)));
        assert_eq!(normalized[0].min_count, 250);
        Ok(())
    }

    #[test]
    fn should_keep_tail_for_merge_gap_zero() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // Merge gap zero still allows overlaps across chunk boundaries
        // Even with only two windows, the last window must stay in the tail
        // The next chunk could add a window that overlaps or touches this one
        // Keeping it in the tail preserves correctness for cross-chunk merges
        let windows = vec![
            build_window("chr1", 10, 20, "A"),
            build_window("chr1", 30, 40, "A"),
        ];

        // Act
        let (safe_prefix, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Within,
            Some(0),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );

        // Assert
        assert_eq!(safe_prefix.len(), 1);
        assert_eq!(tail.len(), 1);
        assert_eq!(tail[0].start_for(CoordinateSet::Resized), 30);
        Ok(())
    }

    #[test]
    fn should_keep_tail_for_min_distance_zero() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // Minimum distance zero still requires checking the next chunk for overlaps
        // Even with only two windows, the last window must stay in the tail
        // The next chunk could add a window that overlaps or touches this one
        // Keeping it in the tail preserves correctness for cross-chunk spacing
        let windows = vec![
            build_window("chr1", 10, 20, "A"),
            build_window("chr1", 30, 40, "A"),
        ];

        // Act
        let (safe_prefix, tail) = partition_safe_and_tail(
            windows,
            Some(0),
            MergeScope::None,
            None,
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );

        // Assert
        assert_eq!(safe_prefix.len(), 1);
        assert_eq!(tail.len(), 1);
        assert_eq!(tail[0].start_for(CoordinateSet::Resized), 30);
        Ok(())
    }

    #[test]
    fn should_keep_tail_for_cluster_overlap() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // Clustering depends on overlap depth, so the last window must carry forward
        // Even with only two windows, the last window must stay in the tail
        // The next chunk could add a window that overlaps and changes cluster depth
        // Keeping it in the tail preserves correctness for cross-chunk clustering
        let windows = vec![
            build_window("chr1", 10, 20, "A"),
            build_window("chr1", 30, 40, "B"),
        ];

        // Act
        let (safe_prefix, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::None,
            None,
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            Some(2),
            CoordinateSet::Resized,
            None,
        );

        // Assert
        assert_eq!(safe_prefix.len(), 1);
        assert_eq!(tail.len(), 1);
        assert_eq!(tail[0].start_for(CoordinateSet::Resized), 30);
        Ok(())
    }

    #[test]
    fn should_keep_tail_for_cross_chunk_merge_gap_zero() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // Gap zero still allows a merge at the chunk boundary
        // Only groups whose final span reaches the boundary must stay in the tail
        // A ends at 30 while the boundary is 40 from B, so only B is kept
        let windows = vec![
            build_window("chr1", 10, 20, "A"),
            build_window("chr1", 20, 30, "A"),
            build_window("chr1", 40, 50, "B"),
        ];
        let merge_group_keys = vec!["A".to_string(), "A".to_string(), "B".to_string()];

        // Act
        let (safe_prefix, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Within,
            Some(0),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            Some(&merge_group_keys),
        );

        // Assert
        assert_eq!(safe_prefix.len(), 2);
        assert_eq!(tail.len(), 1);
        assert_eq!(tail[0].start_for(CoordinateSet::Resized), 40);
        Ok(())
    }

    #[test]
    fn should_keep_merge_chain_in_tail_when_gap_zero() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // Gap zero keeps the whole touching chain that reaches the boundary
        // Group A ends well before the boundary, so it is safe to write
        // Group B forms a chain (20-25-30-35) that can merge across chunks, so keep all of B
        let windows = vec![
            build_window("chr1", 0, 5, "A"),
            build_window("chr1", 20, 25, "B"),
            build_window("chr1", 25, 30, "B"),
            build_window("chr1", 30, 35, "B"),
        ];

        // Act
        let (safe_prefix, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Within,
            Some(0),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );

        // Assert
        assert_eq!(safe_prefix.len(), 1);
        assert_eq!(tail.len(), 3);
        assert_eq!(tail[0].start_for(CoordinateSet::Resized), 20);
        assert_eq!(tail[2].end_for(CoordinateSet::Resized), 35);
        Ok(())
    }

    #[test]
    fn should_load_near_without_strand_column() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // Missing strand column should default to '+' for every near interval
        let mut file = NamedTempFile::new()?;
        writeln!(file, "chr1\t10\t20\tGeneA")?;
        writeln!(file, "chr1\t30\t40\tGeneB")?;

        // Act
        let index = load_near_index(
            file.path(),
            '\t',
            false,
            None,
            Some(&[3]),
            false,
            NearDuplicatesPolicy::Error,
        )?;

        // Assert
        let chr1 = index.per_chrom.get("chr1").expect("chr1 near intervals");
        assert_eq!(chr1.intervals.len(), 2);
        assert_eq!(chr1.intervals[0].strand, Strand::Plus);
        assert_eq!(chr1.intervals[1].strand, Strand::Plus);
        assert_eq!(
            index.group_id_to_name,
            vec!["GeneA".to_string(), "GeneB".to_string()]
        );
        Ok(())
    }

    #[test]
    fn should_drop_windows_without_near_when_distance_max_set() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // No near intervals means no hit, so --distance-max drops windows
        let windows = vec![build_window_with_tuple("chr1", 10, 20)];
        let mut near_index = Some(NearIndex::default());
        let mut cfg = PrepareConfig::default();
        cfg.distance_max = Some(100);

        // Act
        let result =
            apply_near_annotations(windows, &mut near_index, &cfg, None, CoordinateSet::Resized);

        // Assert
        assert!(result.is_empty());
        Ok(())
    }

    #[test]
    fn should_label_no_near_when_bins_set_and_distance_max_unset() -> Result<()> {
        // Human verification status: unverified
        // Arrange
        // No near intervals should emit [NO-NEAR] bin and [NONE] labels
        let windows = vec![build_window_with_tuple("chr1", 10, 20)];
        let mut near_index = Some(NearIndex::default());
        let mut cfg = PrepareConfig::default();
        cfg.near_group_cols = vec!["3".to_string()];
        let bins = parse_distance_bins(&vec!["prox:<100".to_string()])?;

        // Act
        let result = apply_near_annotations(
            windows,
            &mut near_index,
            &cfg,
            Some(&bins),
            CoordinateSet::Resized,
        );

        // Assert
        assert_eq!(result.len(), 1);
        let tuple = &result[0].label_tuples[0];
        assert_eq!(tuple.near_side.as_deref(), Some(NO_NEAR_LABEL));
        assert_eq!(tuple.near_name.as_deref(), Some(NO_NEAR_LABEL));
        assert_eq!(tuple.bin.as_deref(), Some(NO_NEAR_BIN_LABEL));
        Ok(())
    }
}
