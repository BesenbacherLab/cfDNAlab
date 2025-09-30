mod tests {
    #[test]
    fn parse_tile_index_basic() {
        use cfdnalab::shared::tiled_run::parse_tile_index;
        assert_eq!(parse_tile_index("coverage.pos.chr1.12.tsv"), Some(12));
        assert_eq!(
            parse_tile_index("coverage.pos.chr10.000123.bedgraph.zst"),
            Some(123)
        );
        assert_eq!(
            parse_tile_index("coverage.part.chrX.7.part.tsv.zst"),
            Some(7)
        );
        assert_eq!(parse_tile_index("weird.noindex.zst"), None);
    }
}
