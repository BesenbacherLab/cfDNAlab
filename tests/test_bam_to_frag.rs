mod tests_bam_to_frag {
    // tests/bam_to_frag_integration.rs

    use anyhow::{Context, Result};
    use flate2::read::GzDecoder;
    use rust_htslib::bam::index;
    use rust_htslib::bam::{
        self, Format, HeaderView, Writer,
        header::Header,
        record::{Cigar, Record},
    };
    use std::{
        fs::{self, File},
        io::Read as _,
        path::Path,
    };
    use tempfile::tempdir;

    // Bring your crate items into scope.
    use cfdnalab::commands::bam_to_frag::{bam_to_frag::run_inner, config::BamToFragConfig};
    use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs};

    #[test]
    fn bam_to_frag_smoke_two_chroms() -> Result<()> {
        // Temp working dir
        let work = tempdir().context("tempdir")?;
        let work_path = work.path();

        // Paths
        let bam_path = work_path.join("test.bam");
        let bai_path = work_path.join("test.bam.bai");
        let out_dir = work_path.join("out");
        fs::create_dir_all(&out_dir)?;

        // Create a tiny coordinate-sorted BAM with two chromosomes and three pairs:
        // chr1: pair A (R1 forward, R2 reverse, MAPQ=60) -> expect strand '+'
        // chr1: pair B (R1 forward, R2 reverse, MAPQ=0)  -> expect strand '+'
        // chr2: pair C (R1 reverse, R2 forward, min MAPQ=30) -> expect strand '-'
        write_test_bam(&bam_path)?;

        // Build BAI index
        index::build(
            bam_path.to_str().unwrap(),
            None,
            index::Type::BAI,
            1, // n_threads for indexing
        )
        .context("build BAI")?;
        assert!(bai_path.exists(), "BAI was not created");

        // Construct CLI config
        let ioc = IOCArgs {
            bam: bam_path.clone(),
            output_dir: out_dir.clone(),
            n_threads: 2,
        };

        // NOTE: if ChromosomeArgs has a different constructor in your crate, adjust this.
        // The intent here is "process all chromosomes".
        let chromosomes = default_chromosome_args();

        let mut cfg = BamToFragConfig::new(ioc, chromosomes);
        cfg.set_output_prefix("sample");
        cfg.set_min_mapq(0); // keep MAPQ=0 to test min mapq behavior
        cfg.set_require_proper_pair(false); // we never require proper pair
        cfg.set_blacklist(None);

        // Run the command
        let counters = run_inner(&cfg).context("run_inner failed")?;
        assert!(
            counters.base.counted_fragments >= 3,
            "Expected at least 3 fragments counted"
        );

        // Load merged frag file
        let frag_path = out_dir.join("sample.frag.tsv.gz");
        let recs = read_frag_gz(&frag_path)?;

        // Expect exactly 3 lines
        assert_eq!(
            recs.len(),
            3,
            "Unexpected number of fragment lines: {:?}",
            recs
        );

        // Parse into tuples for easier assertions
        // columns: chrom, start, end, mapq, strand
        let parsed: Vec<(String, u64, u64, u8, char)> = recs
            .iter()
            .map(|line| {
                let parts: Vec<&str> = line.split('\t').collect();
                assert_eq!(parts.len(), 5, "Bad line: {}", line);
                let chrom = parts[0].to_string();
                let start: u64 = parts[1].parse().unwrap();
                let end: u64 = parts[2].parse().unwrap();
                let mapq: u8 = parts[3].parse().unwrap();
                let strand: char = parts[4].chars().next().unwrap();
                (chrom, start, end, mapq, strand)
            })
            .collect();

        // Sort for stable assertions
        let mut parsed = parsed;
        parsed.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // chr1 pair A: R1 forward at 10002, R2 reverse at 10090 with 50M -> end=10140, min MQ = 60, strand '+'
        assert_eq!(parsed[0], ("chr1".into(), 10002, 10140, 60, '+'));

        // chr1 pair B: R1 forward at 10003, R2 reverse at 10087 with 50M -> end=10137, min MQ = 0, strand '+'
        assert_eq!(parsed[1], ("chr1".into(), 10003, 10137, 0, '+'));

        // chr2 pair C: R1 reverse at 20090, R2 forward at 20000 with 50M -> start=20000, end=20140, min MQ = 30, strand '-'
        assert_eq!(parsed[2], ("chr2".into(), 20000, 20140, 30, '-'));

        Ok(())
    }

    fn read_frag_gz(path: &Path) -> Result<Vec<String>> {
        let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mut gz = GzDecoder::new(f);
        let mut s = String::new();
        gz.read_to_string(&mut s)?;
        let lines = s.lines().map(|l| l.to_string()).collect();
        Ok(lines)
    }

    fn default_chromosome_args() -> ChromosomeArgs {
        // Adjust if your ChromosomeArgs differs; the goal is to include all contigs.
        // Many codebases implement Default for this.
        ChromosomeArgs::default()
    }

    /// Build a tiny BAM with two contigs and three inward-directed pairs, coordinate-sorted.
    /// All reads are primary, mapped, and paired. CIGAR is 50M.
    fn write_test_bam(path: &Path) -> Result<()> {
        // Header with two contigs
        let mut hdr = Header::new();
        hdr.push_record(
            bam::header::HeaderRecord::new(b"HD")
                .push_tag(b"VN", &"1.6")
                .push_tag(b"SO", &"coordinate"),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &100_000),
        );
        hdr.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr2")
                .push_tag(b"LN", &100_000),
        );

        // Create BAM writer
        let mut writer = Writer::from_path(path, &hdr, Format::BAM).context("create BAM writer")?;
        let header_view = HeaderView::from_header(&hdr);

        // Convenience closures
        let tid_chr1 = header_view.tid(b"chr1").expect("chr1 present") as i32;
        let tid_chr2 = header_view.tid(b"chr2").expect("chr2 present") as i32;

        // QNAMEs
        let q1 = b"pairA";
        let q2 = b"pairB";
        let q3 = b"pairC";

        // CIGAR 50M sequence and qual
        let cigar = vec![Cigar::Match(50)];
        let seq = b"ACGTN".repeat(10); // 50bp
        let qual = vec![30u8; 50];

        // ---- chr1 pair A: R1 forward @10002 (MAPQ 60), R2 reverse @10090 (MAPQ 60) ----
        let r1_a = make_rec(
            q1, tid_chr1, 10002, false, 60, &cigar, &seq, &qual, true, tid_chr1, 10090,
        );
        let r2_a = make_rec(
            q1, tid_chr1, 10090, true, 60, &cigar, &seq, &qual, false, tid_chr1, 10002,
        );

        // ---- chr1 pair B: R1 forward @10003 (MAPQ 0), R2 reverse @10087 (MAPQ 0) ----
        let r1_b = make_rec(
            q2, tid_chr1, 10003, false, 0, &cigar, &seq, &qual, true, tid_chr1, 10087,
        );
        let r2_b = make_rec(
            q2, tid_chr1, 10087, true, 0, &cigar, &seq, &qual, false, tid_chr1, 10003,
        );

        // ---- chr2 pair C: R1 reverse @20090 (MAPQ 30), R2 forward @20000 (MAPQ 40) ----
        let r1_c = make_rec(
            q3, tid_chr2, 20090, true, 30, &cigar, &seq, &qual, true, tid_chr2, 20000,
        );
        let r2_c = make_rec(
            q3, tid_chr2, 20000, false, 40, &cigar, &seq, &qual, false, tid_chr2, 20090,
        );

        // Write in coordinate order
        writer.write(&r1_a)?;
        writer.write(&r1_b)?;
        writer.write(&r2_b)?;
        writer.write(&r2_a)?;
        writer.write(&r2_c)?;
        writer.write(&r1_c)?;
        writer.finish()?;

        Ok(())
    }

    /// Construct a paired-end record.
    ///
    /// - `qname`: read name shared by the pair
    /// - `tid`, `pos`: target and 0-based start
    /// - `is_rev`: strand
    /// - `mapq`: mapping quality
    /// - `cigar`, `seq`, `qual`: alignment, bases, quals
    /// - `is_first_in_template`: true for R1, false for R2
    /// - `mtid`, `mpos`: mate reference and 0-based pos
    fn make_rec(
        qname: &[u8],
        tid: i32,
        pos: i64,
        is_rev: bool,
        mapq: u8,
        cigar: &[Cigar],
        seq: &[u8],
        qual: &[u8],
        is_first_in_template: bool,
        mtid: i32,
        mpos: i64,
    ) -> Record {
        let mut rec = Record::new();
        rec.set(qname, cigar, seq, qual);

        let mut flags: u16 = 0;
        flags |= 0x1; // paired
        if is_first_in_template {
            flags |= 0x40;
        } else {
            flags |= 0x80;
        }
        if is_rev {
            flags |= 0x10;
        }
        // No secondary/supplementary/duplicate/fail flags set.

        rec.set_tid(tid);
        rec.set_pos(pos);
        rec.set_flags(flags);
        rec.set_mapq(mapq);
        rec.set_mtid(mtid);
        rec.set_mpos(mpos);

        // Optional: template length (insert size). Sign convention is aligner-specific; not required by the code under test.
        // We skip setting isize.

        rec
    }
}
