// TODO: Fix this. Just generated but chatty doesn't know frag files (finaledb) so it invents stuff

// ========================= frag file -> Fragment iterator (anyhow) =========================
//
// This module provides a streaming iterator over a "frag" file that yields `Fragment`.
// It supports these row formats:
//   - `tid  start  end`                (numeric TID)
//   - `chrom  start  end`              (BED3-style, with a name->tid mapper)
//   - Auto mode: try numeric TID first, else map `chrom` -> `tid`
//
// Lines starting with '#' or 'track ' or 'browser ' are ignored.
// Extra columns beyond the first three are ignored.
//
// Usage sketch:
//
// use std::sync::Arc;
//
// // Your mapper from chromosome name to tid (header-dependent):
// let mapper = Arc::new(move |chrom: &str| header_tid_by_name(chrom)); // -> Option<i32>
//
// let iter = FragFileIter::open("frags.txt", FragFormat::Auto(mapper))?;
// for frag in iter {
//     let frag = frag?;
//     // Use Fragment here
// }

use anyhow::{Context, Result, anyhow, ensure};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

use crate::utils::fragment::minimal_fragment::Fragment;

/// Function type for mapping chromosome names to `tid`
pub type NameToTidFn = Arc<dyn Fn(&str) -> Option<i32> + Send + Sync>;

/// Supported input formats
pub enum FragFormat {
    /// Columns: `tid  start  end`
    Tid3,
    /// Columns: `chrom  start  end` with provided mapper
    Bed3(NameToTidFn),
    /// Try `tid start end`, else `chrom start end` using mapper
    Auto(NameToTidFn),
}

/// Streaming iterator over a fragment file
pub struct FragFileIter<R: BufRead> {
    reader: R,
    buf: String,
    line_no: usize,
    fmt: FragFormat,
}

impl FragFileIter<BufReader<File>> {
    /// Open a path with a chosen format
    pub fn open(path: impl AsRef<Path>, fmt: FragFormat) -> Result<Self> {
        let f = File::open(&path).with_context(|| {
            format!(
                "failed to open fragment file: {}",
                path.as_ref().to_string_lossy()
            )
        })?;
        Ok(Self::new(BufReader::new(f), fmt))
    }
}

impl<R: BufRead> FragFileIter<R> {
    /// Construct from any `BufRead`
    pub fn new(reader: R, fmt: FragFormat) -> Self {
        Self {
            reader,
            buf: String::new(),
            line_no: 0,
            fmt,
        }
    }

    #[inline]
    fn parse_line(fmt: &FragFormat, line: &str, line_no: usize) -> Result<Fragment> {
        // Split on any whitespace; ignore trailing columns beyond the first three
        let mut it = line.split_whitespace();

        let c1 = it
            .next()
            .ok_or_else(|| anyhow!("line {}: missing column 1", line_no))?;
        let c2 = it
            .next()
            .ok_or_else(|| anyhow!("line {}: missing start", line_no))?;
        let c3 = it
            .next()
            .ok_or_else(|| anyhow!("line {}: missing end", line_no))?;

        // Helper to validate and build Fragment
        let mk = |tid: i32, start_s: &str, end_s: &str| -> Result<Fragment> {
            let start: u32 = start_s
                .parse()
                .with_context(|| format!("line {}: invalid start '{}'", line_no, start_s))?;
            let end: u32 = end_s
                .parse()
                .with_context(|| format!("line {}: invalid end '{}'", line_no, end_s))?;
            ensure!(
                end > start,
                "line {}: end <= start ({} <= {})",
                line_no,
                end,
                start
            );
            Ok(Fragment { tid, start, end })
        };

        match fmt {
            FragFormat::Tid3 => {
                let tid: i32 = c1
                    .parse()
                    .with_context(|| format!("line {}: invalid tid '{}'", line_no, c1))?;
                mk(tid, c2, c3)
            }
            FragFormat::Bed3(mapper) => {
                let tid = mapper(c1)
                    .ok_or_else(|| anyhow!("line {}: unknown contig '{}'", line_no, c1))?;
                mk(tid, c2, c3)
            }
            FragFormat::Auto(mapper) => {
                if let Ok(tid) = c1.parse::<i32>() {
                    mk(tid, c2, c3)
                } else {
                    let tid = mapper(c1)
                        .ok_or_else(|| anyhow!("line {}: unknown contig '{}'", line_no, c1))?;
                    mk(tid, c2, c3)
                }
            }
        }
    }
}

impl<R: BufRead> Iterator for FragFileIter<R> {
    type Item = Result<Fragment>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buf.clear();
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => return None, // EOF
                Ok(_) => {}
                Err(e) => {
                    return Some(Err(anyhow!(
                        "I/O error while reading fragment file at line {}: {}",
                        self.line_no + 1,
                        e
                    )));
                }
            }
            self.line_no += 1;

            let line = self.buf.trim_end_matches(&['\n', '\r'][..]);

            // Skip empties, comments, and common UCSC headers
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track ")
                || line.starts_with("browser ")
            {
                continue;
            }

            // Parse this line
            return Some(Self::parse_line(&self.fmt, line, self.line_no));
        }
    }
}
