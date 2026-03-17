# Contributing

We welcome contributions! We have a few guidelines for this project.

Most importantly in the age of AI: While we totally support including AI agents in the process of making your contribution, this project requires human validation of all code to ensure scientific integrity. That means, we expect you to have manually verified your contributions and that we will be doing so as well. Since that takes a long time for large PRs, please keep them to a readable length. If it takes a month to read, it will not be accepted. This should not keep you from being ambitious though, if the idea and code is great, we want to see it.

Further, please skim through the code of existing commands. Many things are sort of standardized across commands. E.g., use similar config files, fragment iterators (you can build custom *Fragment types but use the same setup for iterating over them), etc. If you need help converting your idea to our setup, just reach out.

On copyright: Please be sure not to contribute code that others have the right to. If you copy from other work (that has a permissable license allowing this), be transparent about in the PR message, give credit, etc.

## Testing

Test all important logic using **mentally derived expectations**. Do not compute expectations (unless necessary, in which case this must be clearly marked) or see passing tests as a goal in themselves. Testing is for detecting bugs and bad assumptions.

Use `tests/` for behavior users can rely on:

- Public crate API
- CLI behavior and output files
- End-to-end workflows
- Regressions in released behavior

Use module-local `*_tests.rs` files for important internal logic that should stay private or `pub(crate)`.
Include the tests in the footer of the module with:

```rust
#[cfg(test)]
mod tests {
    include!("interpolation_tests.rs");
}
```

Do not make items `pub` just to make a test compile.

## API boundaries

This repo is CLI-first. Keep the public Rust API small and intentional.

Choose visibility based on support intent, not test convenience:

- private for local implementation details
- `pub(crate)` for crate-internal shared helpers
- `pub` only for supported external API

Before making something `pub`, ask:

- Is this a real reusable workflow, config type, or result type?
- Would we be comfortable supporting this name and behavior later?

If the answer is unclear, keep it private or `pub(crate)`.
