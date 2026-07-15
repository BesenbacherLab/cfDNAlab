use super::*;
use std::{fs, io::Read};
use tempfile::TempDir;

fn sorted_entry_names(path: &Path) -> Result<Vec<String>> {
    let mut names = fs::read_dir(path)?
        .map(|entry_result| {
            entry_result.map(|entry| entry.file_name().to_string_lossy().into_owned())
        })
        .collect::<std::io::Result<Vec<_>>>()?;
    names.sort();
    Ok(names)
}

#[test]
fn background_text_reader_preserves_the_complete_input_stream() -> Result<()> {
    // Arrange
    let root = TempDir::new()?;
    let input_path = root.path().join("windows.bed");
    let expected = "chr1\t10\t20\tgroup_a\nchr2\t30\t40\tgroup_b\n";
    fs::write(&input_path, expected)?;

    // Act
    let mut reader = open_text_reader_in_background(&input_path)?;
    let mut observed = String::new();
    reader.read_to_string(&mut observed)?;

    // Assert
    assert_eq!(observed, expected);
    Ok(())
}

#[test]
fn final_output_files_moves_all_recorded_files_after_temp_writes() -> Result<()> {
    let root = TempDir::new()?;
    let output_dir = root.path().join("output");
    fs::create_dir(&output_dir)?;
    let run_temp_dir = output_dir.join("run_tmp");
    fs::create_dir(&run_temp_dir)?;

    let mut final_outputs = FinalOutputFiles::new(&run_temp_dir)?;
    let final_output_temp_dir = run_temp_dir.join("final_outputs");

    assert_eq!(sorted_entry_names(&output_dir)?, ["run_tmp"]);
    assert_eq!(sorted_entry_names(&run_temp_dir)?, ["final_outputs"]);
    assert!(sorted_entry_names(&final_output_temp_dir)?.is_empty());

    let files = [
        ("counts.npy", "counts\n"),
        ("group_index.tsv", "0\tgroup_a\n"),
        ("settings.json", "{\"minimum_length\":30}\n"),
    ];

    let mut expected_final_paths = Vec::new();
    let mut expected_temp_paths = Vec::new();
    for (file_name, contents) in files {
        let final_path = output_dir.join(file_name);
        let temp_path = final_outputs.temp_path_for(&final_path)?;

        assert_eq!(temp_path, final_output_temp_dir.join(file_name));
        assert!(!final_path.exists());

        fs::write(&temp_path, contents)?;

        expected_final_paths.push((final_path, contents));
        expected_temp_paths.push(temp_path);
    }
    final_outputs.record_temp_files_with_same_names_in(
        expected_temp_paths.iter().cloned(),
        &output_dir,
    )?;

    assert_eq!(sorted_entry_names(&output_dir)?, ["run_tmp"]);
    assert_eq!(
        sorted_entry_names(&final_output_temp_dir)?,
        ["counts.npy", "group_index.tsv", "settings.json"]
    );
    for (final_path, _) in &expected_final_paths {
        assert!(!final_path.exists());
    }
    for temp_path in &expected_temp_paths {
        assert!(temp_path.exists());
    }

    final_outputs.move_into_place()?;

    assert_eq!(
        sorted_entry_names(&output_dir)?,
        ["counts.npy", "group_index.tsv", "run_tmp", "settings.json"]
    );
    assert!(sorted_entry_names(&final_output_temp_dir)?.is_empty());
    for (final_path, contents) in expected_final_paths {
        assert_eq!(fs::read_to_string(final_path)?, contents);
    }
    for temp_path in expected_temp_paths {
        assert!(!temp_path.exists());
    }

    Ok(())
}

#[test]
fn move_output_file_replaces_existing_directory_output() -> Result<()> {
    let root = TempDir::new()?;
    let temp_path = root.path().join("temp.zarr");
    let final_path = root.path().join("sample.zarr");
    fs::create_dir(&temp_path)?;
    fs::create_dir(&final_path)?;
    fs::write(temp_path.join("zarr.json"), "new\n")?;
    fs::write(final_path.join("zarr.json"), "old\n")?;

    move_output_file_into_place(&temp_path, &final_path)?;

    assert!(!temp_path.exists());
    assert_eq!(fs::read_to_string(final_path.join("zarr.json"))?, "new\n");
    Ok(())
}

#[test]
fn move_output_file_rejects_directory_replacement_without_filename() -> Result<()> {
    let root = TempDir::new()?;
    let temp_path = root.path().join("temp.zarr");
    fs::create_dir(&temp_path)?;

    let err = move_output_file_into_place(&temp_path, Path::new("."))
        .expect_err("directory replacement without a filename should fail");

    assert!(
        err.to_string().contains("no final path component"),
        "unexpected error: {err:#}"
    );
    assert!(temp_path.exists());
    Ok(())
}

#[test]
fn move_output_file_rejects_directory_replacement_without_allowed_extension() -> Result<()> {
    let root = TempDir::new()?;
    let temp_path = root.path().join("temp.zarr");
    let final_path = root.path().join("sample.output_dir");
    fs::create_dir(&temp_path)?;
    fs::create_dir(&final_path)?;

    let err = move_output_file_into_place(&temp_path, &final_path)
        .expect_err("directory replacement without an allowed extension should fail");

    assert!(
        err.to_string().contains("extension is not one of"),
        "unexpected error: {err:#}"
    );
    assert!(temp_path.exists());
    assert!(final_path.exists());
    Ok(())
}
