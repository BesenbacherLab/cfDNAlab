use crate::commands::gc_bias::correction_builder::{Correction, CorrectionMode};
use anyhow::Result;
use ndarray::{Array2, arr0};
use ndarray_npy::NpzWriter;
use std::fs::File;
use std::path::Path;

pub fn save_correction_npz(path: impl AsRef<Path>, corr: &Correction) -> Result<()> {
    let mut npz = NpzWriter::new(File::create(path)?);

    npz.add_array("weights", &corr.weights)?;
    npz.add_array("gc_min", &arr0(corr.gc_min as i64))?;
    npz.add_array("gc_max", &arr0(corr.gc_max as i64))?;
    npz.add_array("len_min", &arr0(corr.len_min as i64))?;
    npz.add_array("len_max", &arr0(corr.len_max as i64))?;
    npz.add_array("mode", &arr0(corr.mode as u8))?;

    match corr.mode {
        CorrectionMode::BySize => {
            let bs = corr.bin_size.expect("bin_size must be set for BySize");
            npz.add_array("bin_size", &arr0(bs as i64))?;
        }
        CorrectionMode::ByBed => {
            let wins = corr
                .windows
                .as_ref()
                .expect("windows must be set for ByBed");
            // Flatten Vec<(u64,u64)> -> Array2<u64> (N,2)
            let mut flat = Vec::with_capacity(wins.len() * 2);
            for (s, e) in wins {
                flat.push(*s);
                flat.push(*e);
            }
            let arr = Array2::from_shape_vec((wins.len(), 2), flat)?;
            npz.add_array("windows", &arr)?;
        }
        CorrectionMode::Global => {}
    }

    npz.finish()?;
    Ok(())
}

// Multiple objects in one .npz
pub fn save_corrections_npz(path: impl AsRef<Path>, list: &[Correction]) -> Result<()> {
    let mut npz = NpzWriter::new(File::create(path)?);
    for (i, corr) in list.iter().enumerate() {
        let p = |k: &str| format!("obj{}_{}", i, k);
        npz.add_array(&p("weights"), &corr.weights)?;
        npz.add_array(&p("gc_min"), &arr0(corr.gc_min as i64))?;
        npz.add_array(&p("gc_max"), &arr0(corr.gc_max as i64))?;
        npz.add_array(&p("len_min"), &arr0(corr.len_min as i64))?;
        npz.add_array(&p("len_max"), &arr0(corr.len_max as i64))?;
        npz.add_array(&p("mode"), &arr0(corr.mode as u8))?;
        match corr.mode {
            CorrectionMode::BySize => {
                let bs = corr.bin_size.expect("bin_size must be set for BySize");
                npz.add_array(&p("bin_size"), &arr0(bs as i64))?;
            }
            CorrectionMode::ByBed => {
                let wins = corr
                    .windows
                    .as_ref()
                    .expect("windows must be set for ByBed");
                let mut flat = Vec::with_capacity(wins.len() * 2);
                for (s, e) in wins {
                    flat.push(*s);
                    flat.push(*e);
                }
                let arr = Array2::from_shape_vec((wins.len(), 2), flat)?;
                npz.add_array(&p("windows"), &arr)?;
            }
            CorrectionMode::Global => {}
        }
    }
    npz.finish()?;
    Ok(())
}
