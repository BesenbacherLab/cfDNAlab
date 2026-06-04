use anyhow::{Context, Result, anyhow};

pub(crate) fn parse_lengths(list: Option<&str>, range: Option<&str>) -> Result<Vec<u32>> {
    if let Some(list) = list {
        let mut lengths = Vec::new();
        for item in list.split(',') {
            let trimmed = item.trim();
            if trimmed.is_empty() {
                continue;
            }
            let value = trimmed
                .parse::<u32>()
                .with_context(|| format!("\"{}\" is not a positive integer", trimmed))?;
            if value == 0 {
                return Err(anyhow!(
                    "fragment lengths must be positive (example: 120,150)"
                ));
            }
            lengths.push(value);
        }
        return Ok(lengths);
    }

    if let Some(range) = range {
        return parse_length_range(range);
    }

    Ok(Vec::new())
}

fn parse_length_range(spec: &str) -> Result<Vec<u32>> {
    let parts: Vec<&str> = spec.split(':').collect();
    if parts.len() < 2 || parts.len() > 3 {
        return Err(anyhow!(
            "length range must follow MIN:MAX[:STEP] (example: 100:200:25)"
        ));
    }
    let min = parts[0]
        .trim()
        .parse::<u32>()
        .with_context(|| format!("\"{}\" is not a positive integer", parts[0].trim()))?;
    let max = parts[1]
        .trim()
        .parse::<u32>()
        .with_context(|| format!("\"{}\" is not a positive integer", parts[1].trim()))?;
    if min == 0 || max == 0 {
        return Err(anyhow!(
            "length range values must be positive (example: 100:200:25)"
        ));
    }
    if min > max {
        return Err(anyhow!(
            "length range requires MIN <= MAX (example: 120:200:10)"
        ));
    }
    let step = if parts.len() == 3 {
        parts[2]
            .trim()
            .parse::<u32>()
            .with_context(|| format!("\"{}\" is not a positive integer", parts[2].trim()))?
    } else {
        20
    };
    if step == 0 {
        return Err(anyhow!(
            "length range step must be positive (example: 80:200:20)"
        ));
    }
    let mut result = Vec::new();
    let mut current = min;
    while current <= max {
        result.push(current);
        match current.checked_add(step) {
            Some(next) if next > current => current = next,
            _ => break,
        }
    }
    if *result.last().unwrap_or(&0) != max && (max - result.last().unwrap_or(&min)) % step != 0 {
        // Ensure the upper bound is included if it falls on the lattice.
        let mut last = *result.last().unwrap();
        while let Some(next) = last.checked_add(step) {
            if next > max {
                break;
            }
            result.push(next);
            last = next;
        }
        if *result.last().unwrap() != max {
            result.push(max);
        }
    }
    Ok(result)
}
