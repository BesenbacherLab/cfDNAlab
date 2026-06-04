mod tests_parsers {

    use crate::commands::prepare_windows::parsers::{
        parse_cols_indices, parse_distance_bins, parse_record_line, parse_score_filter,
        resolve_column_indices,
    };
    use anyhow::Result;

    #[test]
    fn parse_distance_bins_and_match_labels() -> Result<()> {
        let bins = parse_distance_bins(&[
            "prox:<10".to_string(),
            "mid:10-20".to_string(),
            "dist:>20".to_string(),
        ])?;
        assert_eq!(bins.match_label(5), Some("prox"));
        assert_eq!(bins.match_label(15), Some("mid"));
        assert_eq!(bins.match_label(50), Some("dist"));
        assert_eq!(bins.match_label(-5), Some("prox"));
        Ok(())
    }

    #[test]
    fn parse_distance_bins_errors_on_invalid_expr() {
        let err = parse_distance_bins(&["bad".to_string()]).unwrap_err();
        assert!(format!("{err}").contains("Invalid distance bin spec"));
    }

    #[test]
    fn parse_distance_bins_prefers_first_matching_rule() -> Result<()> {
        let bins = parse_distance_bins(&["first:<=10".to_string(), "second:<=20".to_string()])?;
        assert_eq!(bins.match_label(5), Some("first"));
        assert_eq!(bins.match_label(15), Some("second"));
        Ok(())
    }

    #[test]
    fn parse_score_filter_evaluates_condition() -> Result<()> {
        let filter = parse_score_filter(">=1.5")?;
        assert!(filter.eval(2.0));
        assert!(!filter.eval(1.0));
        Ok(())
    }

    #[test]
    fn parse_score_filter_errors_on_invalid_operator() {
        let err = parse_score_filter("~=1.0").unwrap_err();
        assert!(format!("{err}").contains("Invalid score filter"));
    }

    #[test]
    fn resolve_indices_and_parse_record_line() -> Result<()> {
        let cols = resolve_column_indices("chrom=0,start=1,end=2", &["3".to_string()], Some("4"))?;
        let (chrom, start, end, group, score) =
            parse_record_line("chr1\t5\t10\tG\t3.5", '\t', &cols)?;
        assert_eq!(chrom, "chr1");
        assert_eq!((start, end), (5, 10));
        assert_eq!(group, "G");
        assert_eq!(score, Some(3.5));
        Ok(())
    }

    #[test]
    fn parse_record_line_handles_missing_group_columns() -> Result<()> {
        let cols = resolve_column_indices("chrom=0,start=1,end=2", &[], None)?;
        let (chrom, start, end, group, score) = parse_record_line("chr1\t0\t5", '\t', &cols)?;
        assert_eq!(chrom, "chr1");
        assert_eq!((start, end), (0, 5));
        assert!(group.is_empty());
        assert!(score.is_none());
        Ok(())
    }

    #[test]
    fn parse_record_line_errors_on_invalid_interval() {
        let cols = resolve_column_indices("chrom=0,start=1,end=2", &[], None).unwrap();
        let err = parse_record_line("chr1\t10\t5", '\t', &cols).unwrap_err();
        assert!(format!("{err}").contains("End must be greater than start"));
    }

    #[test]
    fn parse_cols_indices_requires_all_fields() {
        let err = parse_cols_indices("chrom=0,start=1").unwrap_err();
        assert!(format!("{err}").contains("cols: missing end="));
    }
}
