use std::fmt::{Display, Formatter};

/// Error returned by Rust output loaders.
///
/// Loader errors keep the detailed file, line, path, and schema context from
/// the parser that failed. The concrete error type is stable for downstream
/// code, while the message remains the main way to inspect malformed output
/// files.
#[derive(Debug)]
pub struct OutputLoaderError {
    source: anyhow::Error,
}

impl OutputLoaderError {
    /// Build a loader error from a plain message.
    pub(crate) fn message(message: impl Into<String>) -> Self {
        Self {
            source: anyhow::anyhow!(message.into()),
        }
    }

    /// Return the underlying contextual error.
    pub fn as_anyhow(&self) -> &anyhow::Error {
        &self.source
    }
}

impl Display for OutputLoaderError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.source, formatter)
    }
}

impl std::error::Error for OutputLoaderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source.source()
    }
}

impl From<anyhow::Error> for OutputLoaderError {
    fn from(source: anyhow::Error) -> Self {
        Self { source }
    }
}

/// Result alias returned by Rust output loaders.
pub type OutputLoaderResult<T> = std::result::Result<T, OutputLoaderError>;
