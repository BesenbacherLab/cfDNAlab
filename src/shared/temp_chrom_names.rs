use anyhow::{Result, bail};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub(crate) fn temp_chrom_token(index: usize) -> String {
    format!("chrom-{index:06}")
}

/// Bidirectional mapping between raw contig names and filesystem-safe temp-file tokens.
///
/// Raw contig names are biological identifiers, not path components. This mapper keeps those
/// identifiers out of intermediate filenames while preserving a reversible in-memory mapping for
/// reducers that need to associate temp files back to their source contig.
#[derive(Debug, Clone)]
pub(crate) struct TempChromNameMap {
    raw_to_token: HashMap<String, String>,
    token_to_raw: HashMap<String, String>,
}

impl TempChromNameMap {
    pub(crate) fn from_contigs(contigs: &[String]) -> Result<Self> {
        let mut raw_to_token = HashMap::with_capacity(contigs.len());
        let mut token_to_raw = HashMap::with_capacity(contigs.len());

        for (index, contig) in contigs.iter().enumerate() {
            if raw_to_token.contains_key(contig) {
                bail!(
                    "duplicate contig name '{}' cannot be mapped to a temp filename",
                    contig
                );
            }
            let token = temp_chrom_token(index);
            raw_to_token.insert(contig.clone(), token.clone());
            token_to_raw.insert(token, contig.clone());
        }

        Ok(Self {
            raw_to_token,
            token_to_raw,
        })
    }

    pub(crate) fn token_for(&self, contig: &str) -> Result<&str> {
        self.raw_to_token
            .get(contig)
            .map(String::as_str)
            .ok_or_else(|| anyhow::anyhow!("missing temp filename token for contig '{}'", contig))
    }
    pub(crate) fn raw_for(&self, token: &str) -> Result<&str> {
        self.token_to_raw
            .get(token)
            .map(String::as_str)
            .ok_or_else(|| anyhow::anyhow!("missing raw contig name for temp token '{}'", token))
    }

    pub(crate) fn path_with_suffix(
        &self,
        temp_dir: &Path,
        contig: &str,
        suffix: &str,
    ) -> Result<PathBuf> {
        Ok(temp_dir.join(format!("{}.{}", self.token_for(contig)?, suffix)))
    }
}

#[cfg(test)]
mod tests {
    use super::TempChromNameMap;

    #[test]
    fn maps_path_like_contigs_to_reversible_opaque_tokens() {
        // Manual expectations:
        // - Raw names must not appear in the generated token, even when they contain path syntax.
        // - The token stays reversible through the same in-memory map.
        // - A token-based suffix path is a single filename under the temp directory.
        let contigs = vec![
            "chr/with/slash".to_string(),
            "../chr2".to_string(),
            "chrom-000000".to_string(),
        ];
        let map = TempChromNameMap::from_contigs(&contigs).expect("valid contig map");

        assert_eq!(map.token_for("chr/with/slash").unwrap(), "chrom-000000");
        assert_eq!(map.token_for("../chr2").unwrap(), "chrom-000001");
        assert_eq!(map.raw_for("chrom-000002").unwrap(), "chrom-000000");

        let path = map
            .path_with_suffix(
                std::path::Path::new("/tmp/run"),
                "chr/with/slash",
                "frag.tmp",
            )
            .unwrap();
        assert_eq!(path, std::path::Path::new("/tmp/run/chrom-000000.frag.tmp"));
    }

    #[test]
    fn rejects_duplicate_raw_contigs() {
        let contigs = vec!["chr1".to_string(), "chr1".to_string()];
        let err = TempChromNameMap::from_contigs(&contigs)
            .expect_err("duplicate raw contigs should fail");
        assert!(
            err.to_string().contains("duplicate contig name 'chr1'"),
            "unexpected error: {err}"
        );
    }
}
