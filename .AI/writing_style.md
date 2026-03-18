# Writing Style

## Code Style

- **Write for humans.** Prefer clear, direct language over buzzwords or slang. Comments should help a typical programmer understand intent and edge cases.
- **Comment generously.** Explain *why* as well as *what*. If you remove outdated comments, **add** updated ones so net documentation never decreases.
- **Keep code simple.** Choose the simplest readable approach unless a more complex solution brings substantial gains (speed/memory). If complexity seems warranted, briefly note the trade-off or ask before proceeding.
- **No inline tests.** Public behavior tests live under `tests/`. Internal logic that should stay private or `pub(crate)` goes in sibling `*_tests.rs` files included from the owning module, not as inline test bodies mixed into production code.
- **Descriptive variable names.** Variables should have descriptive names so the code is readable. Never use single-letter variables. E.g. a "window start" position can be called "window_start", "win_start", but never "ws" or "s".

## CLI help

Help strings are defined via docstrings in the config files. This needs to be useful for any newcomer or experienced user.

## Docstrings

Docstrings should read like a short tutorial, then details, then structured sections. You may also add examples when they are relevant.

Bullet points in CLI-facing documentation (config files) should have a newline between them, otherwise CLI collapses the sentences.

Reduce the number of semi-colons in docstrings and comments. Use comma or dot instead.

In-line comments start with title-cased first word and does not have a terminal dot *in the end*. E.g. `// A comment`

Never use "…" or similar non-ascii symbols. They don't work in the terminal. Use "...", "->", ">=", etc.

Adapt to my language. If I use "center", don't use "centre". Don't use words/phrases that humans rarely use when unnecessary, like "emit", "bubbles up"/"bubbling".

**Order**

1. **Summary (pedagogical):** What this does and when to use it. (The pedagogical part is implicit, don't add "friendly summary" etc.)
2. **Technical details:** Key behavior, assumptions, edge cases, complexity notes.
3. **Parameters**
4. **Returns**

**Template**

```python
def fn(...):
    """
    Short, friendly summary that teaches the idea in plain language.

    Technical details that note important behavior, invariants, and caveats.
    Mention performance characteristics if relevant.

    Parameters
    ----------
    - `arg1`:
        What it is and how it is used.
    - `arg2`:
        Constraints, defaults, and special cases.

    Returns
    -------
    - `out`:
        What is returned and how to interpret it.
    """
```

### Clap

Don't use `long = "name"`, just `long`. It is automatically filled in.

Clap already specifies default values, so don't add "default is xx".

## Additionals

Don't spend time reordering imports manually. Just let autoformatting do that for us.

DO NOT make conclusions about the code, if you haven't re-read it! If I ask for another review, expect the code to have changed! You cannot just read half a page and then make claims about the whole code base. Otherwise you will not catch new errors.

If a file is uploaded, ALWAYS look at it before answering.

When making larger refactors, do NOT remove comments unless they are no longer true. If it's there, it's because I find it relevant. Keep them.

DO NOT USE SEMICOLONS ";" IN DOCSTRINGS!
DO NOT BE LAZY WITH COMMENTS!
