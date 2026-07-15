pub(crate) fn benchmark_frag_input_loading(
    frag_path: &Path,
    read_in_background: bool,
) -> Result<u64> {
    let column_layout = resolve_frag_column_layout(frag_path, None, false, false)?;
    let reader = if read_in_background {
        open_text_reader_in_background(frag_path)
    } else {
        open_text_reader(frag_path)
    }?;
    let mut non_empty_lines_seen = 0_u64;
    let mut parsed_fragments = 0_u64;

    for (line_index, line_result) in reader.lines().enumerate() {
        let line_number = line_index as u64 + 1;
        let line = line_result.with_context(|| format!("Reading line {line_number}"))?;
        if line.trim().is_empty() {
            continue;
        }
        non_empty_lines_seen += 1;
        if column_layout.skip_first_non_empty_line && non_empty_lines_seen == 1 {
            continue;
        }
        std::hint::black_box(parse_frag_line(
            &line,
            line_number,
            &column_layout.indices,
        )?);
        parsed_fragments += 1;
    }

    Ok(parsed_fragments)
}
