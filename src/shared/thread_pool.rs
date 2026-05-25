use anyhow::{Result, anyhow};
use std::sync::OnceLock;

static INIT: OnceLock<Result<(), String>> = OnceLock::new();

/// Initialize the global Rayon thread pool once, ignoring subsequent requests with
/// different sizes (tests may invoke multiple commands sequentially).
pub(crate) fn init_global_pool(num_threads: usize) -> Result<()> {
    let res = INIT.get_or_init(|| {
        match rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
        {
            Ok(_) => Ok(()),
            Err(err) => {
                let msg = err.to_string();
                if msg.contains("already been initialized") {
                    Ok(())
                } else {
                    Err(msg)
                }
            }
        }
    });

    match res {
        Ok(()) => Ok(()),
        Err(msg) => Err(anyhow!(msg.clone())),
    }
}
