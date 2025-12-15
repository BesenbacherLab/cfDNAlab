mod tests_gc_bias_window_logic {
    use anyhow::Result;
    use ndarray::{Array1, Array2};
    use tempfile::tempdir;

    use cfdnalab::commands::{
        cli_common::{ChromosomeArgs, IOCArgs},
        gc_bias::{
            config::GCConfig,
            counting::GCCounts,
            gc_bias::{process_window, stream_crossing_files},
        },
    };

    fn make_config(tmp: &tempfile::TempDir) -> GCConfig {
        let ioc = IOCArgs {
            bam: tmp.path().join("dummy.bam"),
            output_dir: tmp.path().join("out"),
            n_threads: 1,
        };
        let mut cfg = GCConfig::new(
            ioc,
            tmp.path().join("ref.2bit"),
            tmp.path().join("ref_gc"),
            ChromosomeArgs::default(),
        );
        cfg.set_min_window_acgt_pct(0);
        cfg
    }

    #[test]
    fn scales_window_by_mean_and_acgt_coverage() -> Result<()> {
        // Arrange: One length row (effective length 10 -> 11 GC bins). Only two bins set (2 and 4),
        // so mean = (2+4) / 11 = 0.54545...
        // Scale factor = (1/mean) * (num_acgt/avg_span) = (1/0.54545) * (40/100) = 0.73333...
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);

        let mut counts = GCCounts::new(10, 10, 0, (40, 50))?;
        counts.set(10, 0, 2.0);
        counts.set(10, 1, 4.0);

        // Act
        let scaled = process_window(counts, &cfg, Some(100.0))?.expect("window should be retained");

        // Assert
        let c0 = scaled.get(10, 0).unwrap();
        let c1 = scaled.get(10, 1).unwrap();
        assert!((c0 - 1.4666667).abs() < 1e-6);
        assert!((c1 - 2.9333334).abs() < 1e-6);
        Ok(())
    }

    // TODO: Validate this
    #[test]
    fn merges_crossing_files_and_scales_once_per_window() -> Result<()> {
        // Arrange: two crossing chunks for the same window idx=3, counts 2 and 3, acgt 20 and 30.
        // Merged counts=5, num_acgt=50 -> mean=5/11=0.45454..., scale=(1/0.45454)*(50/20)=5.5, final count=27.5.
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);
        let template = GCCounts::new(10, 10, 0, (0, 0))?;
        let counts_len = template.borrow_raw_counts().len();

        let file1 = tmp.path().join("cross.1.npz");
        let mut npz1 = ndarray_npy::NpzWriter::new(std::fs::File::create(&file1)?);
        npz1.add_array("idx", &Array1::from(vec![3u64]))?;
        npz1.add_array("acgt0", &Array1::from(vec![20u64]))?;
        npz1.add_array("acgt1", &Array1::from(vec![20u64]))?;
        let mut counts_arr1 = Array2::zeros((1, counts_len));
        counts_arr1[[0, 0]] = 2.0;
        npz1.add_array("counts", &counts_arr1)?;
        npz1.finish()?;

        let file2 = tmp.path().join("cross.2.npz");
        let mut npz2 = ndarray_npy::NpzWriter::new(std::fs::File::create(&file2)?);
        npz2.add_array("idx", &Array1::from(vec![3u64]))?;
        npz2.add_array("acgt0", &Array1::from(vec![30u64]))?;
        npz2.add_array("acgt1", &Array1::from(vec![30u64]))?;
        let mut counts_arr2 = Array2::zeros((1, counts_len));
        counts_arr2[[0, 0]] = 3.0;
        npz2.add_array("counts", &counts_arr2)?;
        npz2.finish()?;

        // Act
        let (merged, weight) = stream_crossing_files(vec![file1, file2], &template, &cfg, 20.0)?;

        // Assert
        assert_eq!(weight, 1, "one window should contribute once");
        let v = merged.get(10, 0).unwrap();
        assert!((v - 27.5).abs() < 1e-6);
        Ok(())
    }
}
