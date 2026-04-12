use super::*;

#[test]
fn parse_base_quality_filter_accepts_minimum_end_threshold() {
    // Arrange: `min in end >= 30` should parse into the smallest supported reduction on the
    // per-end scope with an inclusive threshold.
    //
    // Mental derivation:
    // - `min` -> `BaseQualityAggregation::Min`
    // - `end` -> `BaseQualityFilterScope::End`
    // - `>=` -> `BaseQualityComparisonOp::Ge`
    // - `30` -> floating threshold `30.0`
    let filter = "min in end >= 30"
        .parse::<BaseQualityFilter>()
        .expect("valid per-end filter should parse");

    // Assert
    assert_eq!(filter.aggregation, BaseQualityAggregation::Min);
    assert_eq!(filter.scope, BaseQualityFilterScope::End);
    assert_eq!(filter.op, BaseQualityComparisonOp::Ge);
    assert_eq!(filter.threshold, 30.0);
}

#[test]
fn parse_base_quality_filter_accepts_fragment_threshold_without_space_inside_comparison() {
    // Arrange: users may write the comparison as one token (`<25.5`) while still keeping the rest
    // of the expression space-separated. That should remain valid.
    //
    // Mental derivation:
    // - `mean` -> fragment-level mean score
    // - `<25.5` -> strict lower-than comparison with threshold 25.5
    let filter = "mean in fragment <25.5"
        .parse::<BaseQualityFilter>()
        .expect("comparison token without inner space should parse");

    // Assert
    assert_eq!(filter.aggregation, BaseQualityAggregation::Mean);
    assert_eq!(filter.scope, BaseQualityFilterScope::Fragment);
    assert_eq!(filter.op, BaseQualityComparisonOp::Lt);
    assert_eq!(filter.threshold, 25.5);
}

#[test]
fn parse_base_quality_filter_accepts_maximum_fragment_threshold() {
    // Arrange: `max in fragment < 20` is the strongest case for supporting `max`, because it
    // asks whether no end score reaches the threshold.
    //
    // Mental derivation:
    // - `max` -> `BaseQualityAggregation::Max`
    // - `fragment` -> `BaseQualityFilterScope::Fragment`
    // - `<` -> `BaseQualityComparisonOp::Lt`
    // - `20` -> floating threshold `20.0`
    let filter = "max in fragment < 20"
        .parse::<BaseQualityFilter>()
        .expect("valid max-based fragment filter should parse");

    // Assert
    assert_eq!(filter.aggregation, BaseQualityAggregation::Max);
    assert_eq!(filter.scope, BaseQualityFilterScope::Fragment);
    assert_eq!(filter.op, BaseQualityComparisonOp::Lt);
    assert_eq!(filter.threshold, 20.0);
}

#[test]
fn parse_base_quality_filter_trims_outer_whitespace_and_joins_comparison_tokens() {
    // Arrange: extra outer whitespace and splitting the comparison into two tokens should not
    // change the parsed meaning.
    //
    // Mental derivation:
    // - after trimming and token joining, the comparison becomes `<=40`
    // - the rest of the expression stays `mean in end`
    let filter = "  mean in end <= 40  "
        .parse::<BaseQualityFilter>()
        .expect("valid filter with extra whitespace should parse");

    // Assert
    assert_eq!(filter.aggregation, BaseQualityAggregation::Mean);
    assert_eq!(filter.scope, BaseQualityFilterScope::End);
    assert_eq!(filter.op, BaseQualityComparisonOp::Le);
    assert_eq!(filter.threshold, 40.0);
}

#[test]
fn parse_base_quality_filter_accepts_mixed_case_keywords() {
    // Arrange: the user-facing keywords are easier to work with if they are ASCII
    // case-insensitive. `MiN in FrAgMeNt >= 30` should therefore normalize cleanly.
    //
    // Mental derivation:
    // - `MiN` lowercases to `min`
    // - `FrAgMeNt` lowercases to `fragment`
    // - `>= 30` stays unchanged
    let filter = "MiN in FrAgMeNt >= 30"
        .parse::<BaseQualityFilter>()
        .expect("mixed-case keywords should parse");

    // Assert
    assert_eq!(filter.aggregation, BaseQualityAggregation::Min);
    assert_eq!(filter.scope, BaseQualityFilterScope::Fragment);
    assert_eq!(filter.op, BaseQualityComparisonOp::Ge);
    assert_eq!(filter.threshold, 30.0);
}

#[test]
fn parse_base_quality_filter_accepts_tabs_and_newlines_as_whitespace() {
    // Arrange: `split_whitespace()` means tabs and newlines are token separators just like spaces.
    //
    // Mental derivation:
    // - `mean` remains the aggregation
    // - `in`
    // - `end`
    // - `<=`
    // - `25`
    //   become the same token sequence regardless of whether the separators are spaces,
    //   tabs, or newlines.
    let filter = "mean\tin\nend\t<=\n25"
        .parse::<BaseQualityFilter>()
        .expect("tabs and newlines should parse like spaces");

    // Assert
    assert_eq!(filter.aggregation, BaseQualityAggregation::Mean);
    assert_eq!(filter.scope, BaseQualityFilterScope::End);
    assert_eq!(filter.op, BaseQualityComparisonOp::Le);
    assert_eq!(filter.threshold, 25.0);
}

#[test]
fn base_quality_filter_as_cli_expr_round_trips_to_the_same_filter() {
    // Arrange: formatting should preserve the exact semantics needed for a parse -> format ->
    // parse round-trip.
    let original = BaseQualityFilter {
        aggregation: BaseQualityAggregation::Mean,
        scope: BaseQualityFilterScope::Fragment,
        op: BaseQualityComparisonOp::Gt,
        threshold: 41.5,
    };

    // Act
    let reparsed = original
        .as_cli_expr()
        .parse::<BaseQualityFilter>()
        .expect("formatted filter should parse back");

    // Assert
    assert_eq!(reparsed, original);
}

#[test]
fn parse_base_quality_filter_errors_on_invalid_aggregation() {
    // Arrange: only `min`, `mean`, and `max` are part of the minimal grammar. `median` should
    // fail loudly.
    let err = "median in end >= 30"
        .parse::<BaseQualityFilter>()
        .expect_err("unsupported aggregation should fail");

    // Assert
    assert!(err.contains("Invalid base-quality aggregation"));
}

#[test]
fn parse_base_quality_filter_errors_when_in_keyword_is_missing() {
    // Arrange: the grammar requires the literal `in` between the aggregation and scope.
    let err = "min end >= 30"
        .parse::<BaseQualityFilter>()
        .expect_err("missing `in` should fail");

    // Assert
    assert!(err.contains("Invalid base-quality filter"));
}

#[test]
fn parse_base_quality_filter_errors_on_invalid_scope() {
    // Arrange: `motif` is intentionally not accepted because the base-quality filter grammar only
    // distinguishes `end` and `fragment`.
    let err = "mean in motif >= 30"
        .parse::<BaseQualityFilter>()
        .expect_err("unsupported scope should fail");

    // Assert
    assert!(err.contains("Invalid base-quality filter scope"));
}

#[test]
fn parse_base_quality_filter_errors_on_invalid_operator() {
    // Arrange: `==` is intentionally excluded from the minimal grammar.
    let err = "mean in end == 30"
        .parse::<BaseQualityFilter>()
        .expect_err("unsupported operator should fail");

    // Assert
    assert!(err.contains("Invalid base-quality filter"));
}

#[test]
fn parse_base_quality_filter_errors_on_negative_threshold() {
    // Arrange: base qualities are non-negative, so a negative threshold should be rejected during
    // parsing instead of being carried through as dead configuration.
    let err = "min in fragment >= -1"
        .parse::<BaseQualityFilter>()
        .expect_err("negative threshold should fail");

    // Assert
    assert!(err.contains("Base-quality threshold must be >= 0"));
}

#[test]
fn parse_base_quality_filter_errors_when_trailing_tokens_remain_after_threshold() {
    // Arrange: the threshold must be the last part of the expression. Extra tokens would make the
    // configuration too easy to misread.
    let err = "min in end >= 30 extra"
        .parse::<BaseQualityFilter>()
        .expect_err("trailing tokens should fail");

    // Assert
    assert!(err.contains("Invalid base-quality filter"));
}

#[test]
fn parse_base_quality_filter_errors_when_trailing_numeric_token_would_merge_into_threshold() {
    // Arrange: `30 5` must not silently merge into threshold `305`.
    let err = "min in end >= 30 5"
        .parse::<BaseQualityFilter>()
        .expect_err("trailing numeric token should fail");

    // Assert
    assert!(err.contains("Invalid base-quality filter"));
}
