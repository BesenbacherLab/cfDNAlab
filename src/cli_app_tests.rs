//! Tests for shared terminal help formatting.

#[test]
fn subcommand_help_identifies_the_full_command_path() {
    let mut command = super::build_terminal_command();
    let ref_kmers = command
        .find_subcommand_mut("ref-kmers")
        .expect("ref-kmers should be available when cmd_ref_kmers is enabled");

    let help = ref_kmers.render_long_help().to_string();

    assert!(help.lines().any(|line| line == "Command: cfdna ref-kmers"));
}

#[test]
fn top_level_help_identifies_the_program() {
    let help = super::build_terminal_command()
        .render_long_help()
        .to_string();

    assert!(help.lines().any(|line| line == "Command: cfdna"));
}
