#![cfg(output_loader_api)]
//! Public API tests for the Rust output-loader error type.
//!
//! These tests exercise the error surface the way a downstream Rust app would:
//! importing the public type, returning `OutputLoaderResult`, converting
//! contextual `anyhow` errors with `?`, and inspecting display text and sources.

use anyhow::Context as _;
use cfdnalab::output_loaders::{OutputLoaderError, OutputLoaderResult};
use std::{error::Error, io};

/// Verify the public error type has the standard app-facing error traits.
#[test]
fn output_loader_error_can_be_stored_as_thread_safe_standard_error() {
    // Arrange, Act, Assert:
    // Downstream apps should be able to store loader errors behind ordinary
    // error trait objects and move them across threads if their app does that.
    let loader_error = OutputLoaderError::from(anyhow::anyhow!("thread-safe loader failure"));
    let boxed_error: Box<dyn Error + Send + Sync> = Box::new(loader_error);

    assert_eq!(boxed_error.to_string(), "thread-safe loader failure");
}

/// Verify direct conversion from anyhow preserves the visible error message.
#[test]
fn output_loader_error_can_be_built_from_anyhow_error() {
    // Arrange:
    // Public loader internals and downstream wrappers convert contextual
    // anyhow errors into the stable loader error type.
    let source_error = anyhow::anyhow!("malformed output table");

    // Act
    let loader_error = OutputLoaderError::from(source_error);

    // Assert
    assert_eq!(loader_error.to_string(), "malformed output table");
    assert_eq!(
        loader_error.as_anyhow().to_string(),
        "malformed output table"
    );
    assert!(
        Error::source(&loader_error).is_none(),
        "plain anyhow messages have no lower-level source"
    );
}

/// Verify downstream apps can return OutputLoaderResult and use `?` on anyhow errors.
#[test]
fn output_loader_result_accepts_contextual_anyhow_errors_from_app_code() {
    // Arrange:
    // This models app code that wraps lower-level IO or parsing failures with
    // app-specific context before returning the loader error type.
    fn read_user_output() -> anyhow::Result<()> {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing required count column",
        ))
        .context("parse user-selected output")?;
        Ok(())
    }

    fn app_entrypoint() -> OutputLoaderResult<()> {
        read_user_output().context("load output in downstream app")?;
        Ok(())
    }

    // Act
    let loader_error = app_entrypoint().expect_err("app output loading should fail");

    // Assert
    assert_eq!(loader_error.to_string(), "load output in downstream app");
    assert_eq!(
        loader_error.as_anyhow().root_cause().to_string(),
        "missing required count column"
    );

    let chain = loader_error
        .as_anyhow()
        .chain()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    assert_eq!(
        chain,
        vec![
            "load output in downstream app",
            "parse user-selected output",
            "missing required count column",
        ]
    );

    let standard_error: &dyn Error = &loader_error;
    assert_eq!(standard_error.to_string(), "load output in downstream app");
    assert_eq!(
        standard_error
            .source()
            .expect("contextual error should expose a source")
            .to_string(),
        "parse user-selected output"
    );
}

/// Verify downstream apps can propagate loader errors into their own anyhow result.
#[test]
fn output_loader_result_can_be_propagated_into_anyhow_result_by_app_code() {
    // Arrange:
    // This models the common external-app shape: call a cfDNAlab loader from
    // code that returns anyhow::Result and add app-level context at the callsite.
    fn library_loader_boundary() -> OutputLoaderResult<()> {
        Err(anyhow::Error::new(io::Error::new(
            io::ErrorKind::InvalidData,
            "corrupt row offsets",
        ))
        .context("read cfDNAlab output")
        .into())
    }

    fn app_entrypoint() -> anyhow::Result<()> {
        library_loader_boundary().context("run downstream report")?;
        Ok(())
    }

    // Act
    let app_error = app_entrypoint().expect_err("downstream app should see loader failure");

    // Assert
    assert_eq!(app_error.to_string(), "run downstream report");
    assert_eq!(app_error.root_cause().to_string(), "corrupt row offsets");

    let chain = app_error
        .chain()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    assert_eq!(
        chain,
        vec![
            "run downstream report",
            "read cfDNAlab output",
            "corrupt row offsets",
        ]
    );
}
