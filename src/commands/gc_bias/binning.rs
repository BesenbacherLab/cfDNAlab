use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use ndarray::{Array2, ArrayBase, ArrayView2, Axis, Data, Ix2};

#[derive(Debug, Clone)]
pub struct BinnedAxis {
    pub index_to_bin: FxHashMap<usize, usize>,
    pub bin_to_indices: FxHashMap<usize, Vec<usize>>,
    pub num_bins: usize,
}

pub enum CollapseAggregation {
    Sum,
    Mean,
}

pub fn bin_greedily_by_mass<S>(
    counts: &ArrayBase<S, Ix2>,
    axis: usize,
    min_mass_pct: f64,
) -> Result<BinnedAxis>
where
    S: Data<Elem = f64>,
{
    ensure!(axis < 2, "axis must be 0 or 1");
    ensure!(
        (0.0..=100.0).contains(&min_mass_pct),
        "min_mass_pct must be within 0..=100"
    );

    // Sum along the other axis to get per-index mass
    let masses = match axis {
        0 => counts.sum_axis(Axis(1)),
        _ => counts.sum_axis(Axis(0)),
    };

    let total_mass: f64 = masses.iter().sum();
    if total_mass == 0.0 {
        return Ok(BinnedAxis {
            index_to_bin: FxHashMap::default(),
            bin_to_indices: FxHashMap::default(),
            num_bins: 0,
        });
    }

    let min_mass = total_mass * (min_mass_pct / 100.0);
    let mut bins: Vec<Vec<usize>> = Vec::new();
    let mut running_mass = 0.0;
    let mut current_bin_indices: Vec<usize> = Vec::new();

    for (idx, &mass) in masses.iter().enumerate() {
        running_mass += mass;
        current_bin_indices.push(idx);

        if running_mass >= min_mass {
            bins.push(current_bin_indices.clone());
            current_bin_indices.clear();
            running_mass = 0.0;
        }
    }

    if !current_bin_indices.is_empty() {
        if bins.is_empty() {
            bins.push(current_bin_indices);
        } else {
            bins.last_mut().unwrap().extend(current_bin_indices);
        }
    }

    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();

    for (bin_idx, indices) in bins.iter().enumerate() {
        bin_to_indices.insert(bin_idx, indices.clone());
        for &idx in indices {
            index_to_bin.insert(idx, bin_idx);
        }
    }

    let num_bins = bin_to_indices.len();

    Ok(BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins,
    })
}

pub fn collapse_counts_by_bins<S>(
    counts: &ArrayBase<S, Ix2>,
    axis: usize,
    bins: &BinnedAxis,
    agg: CollapseAggregation,
    mass_counts: Option<ArrayView2<'_, f64>>,
) -> Result<Array2<f64>>
where
    S: Data<Elem = f64>,
{
    ensure!(axis < 2, "axis must be 0 or 1");
    if let Some(mass) = mass_counts.as_ref() {
        ensure!(
            mass.dim() == counts.dim(),
            "mass_counts must have same shape as counts"
        );
        if matches!(agg, CollapseAggregation::Sum) {
            bail!("mass_counts provided for Sum aggregation; weighted sums are unsupported");
        }
    }

    let (n_rows, n_cols) = counts.dim();
    match axis {
        0 => {
            let weights = mass_counts.as_ref().map(|m| m.sum_axis(Axis(1)));
            let mut out = Array2::<f64>::zeros((bins.num_bins, n_cols));
            for bin_idx in 0..bins.num_bins {
                if let Some(indices) = bins.bin_to_indices.get(&bin_idx) {
                    let mut denom = 0.0;
                    let mut count = 0usize;
                    for &row_idx in indices {
                        let source = counts.row(row_idx);
                        let mut dest = out.row_mut(bin_idx);
                        match agg {
                            CollapseAggregation::Sum => {
                                dest += &source;
                            }
                            CollapseAggregation::Mean => {
                                if let Some(ref weights_vec) = weights {
                                    let weight = weights_vec[row_idx];
                                    denom += weight;
                                    dest.scaled_add(weight, &source);
                                } else {
                                    dest += &source;
                                    count += 1;
                                }
                            }
                        }
                    }
                    if matches!(agg, CollapseAggregation::Mean) {
                        let mut dest = out.row_mut(bin_idx);
                        if weights.is_some() {
                            if denom > 0.0 {
                                dest /= denom;
                            } else if !indices.is_empty() {
                                dest /= indices.len() as f64;
                            }
                        } else if count > 0 {
                            dest /= count as f64;
                        }
                    }
                }
            }
            Ok(out)
        }
        _ => {
            let weights = mass_counts.as_ref().map(|m| m.sum_axis(Axis(0)));
            let mut out = Array2::<f64>::zeros((n_rows, bins.num_bins));
            for bin_idx in 0..bins.num_bins {
                if let Some(indices) = bins.bin_to_indices.get(&bin_idx) {
                    let mut denom = 0.0;
                    let mut count = 0usize;
                    for &col_idx in indices {
                        let source = counts.column(col_idx);
                        let mut dest = out.column_mut(bin_idx);
                        match agg {
                            CollapseAggregation::Sum => {
                                dest += &source;
                            }
                            CollapseAggregation::Mean => {
                                if let Some(ref weights_vec) = weights {
                                    let weight = weights_vec[col_idx];
                                    denom += weight;
                                    dest.scaled_add(weight, &source);
                                } else {
                                    dest += &source;
                                    count += 1;
                                }
                            }
                        }
                    }
                    if matches!(agg, CollapseAggregation::Mean) {
                        let mut dest = out.column_mut(bin_idx);
                        if let Some(_) = weights {
                            if denom > 0.0 {
                                dest /= denom;
                            } else if !indices.is_empty() {
                                dest /= indices.len() as f64;
                            }
                        } else if count > 0 {
                            dest /= count as f64;
                        }
                    }
                }
            }
            Ok(out)
        }
    }
}

pub fn compute_bin_edges(bins: &BinnedAxis, start_value: u32, max_value: u32) -> Result<Vec<u32>> {
    ensure!(
        bins.num_bins > 0,
        "Bin definition must contain at least one bin"
    );
    let mut edges = Vec::with_capacity(bins.num_bins + 1);
    for bin_idx in 0..bins.num_bins {
        let indices = bins
            .bin_to_indices
            .get(&bin_idx)
            .context("Missing indices for bin")?;
        let min_idx = indices
            .iter()
            .min()
            .copied()
            .context("Bin indices cannot be empty")?;
        edges.push(start_value + min_idx as u32);
    }
    edges.push(max_value);
    Ok(edges)
}

pub fn bins_from_edges(edges: &[u32]) -> Result<BinnedAxis> {
    ensure!(
        edges.len() >= 2,
        "Bin edges must contain at least a start and end entry"
    );
    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();
    let base = edges[0];

    for (bin_idx, window) in edges.windows(2).enumerate() {
        let edge_start = window[0];
        let edge_end = window[1];
        ensure!(
            edge_start >= base,
            "Edge values must be >= the first edge ({}). Found {}",
            base,
            edge_start
        );
        ensure!(
            edge_end >= edge_start,
            "Bin edges must be non-decreasing. Found {} then {}",
            edge_start,
            edge_end
        );
        let start_idx = (edge_start - base) as usize;
        let mut end_idx = (edge_end - base) as usize;
        let is_last_bin = bin_idx == edges.len() - 2;
        if is_last_bin {
            // Last edge is inclusive, so add one to get the exclusive bound.
            end_idx += 1;
        }
        ensure!(
            end_idx > start_idx,
            "Edge interval [{}, {}{} should contain at least one value",
            edge_start,
            edge_end,
            if is_last_bin { "]" } else { ")" }
        );
        let indices: Vec<usize> = (start_idx..end_idx).collect();
        bin_to_indices.insert(bin_idx, indices.clone());
        for idx in indices {
            index_to_bin.insert(idx, bin_idx);
        }
    }

    Ok(BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins: edges.len() - 1,
    })
}
