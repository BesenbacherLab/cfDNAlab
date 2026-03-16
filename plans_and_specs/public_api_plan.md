# Public API plan

This plan lists the initial set of public elements for rust devs. 

---

Before marking something `pub`, ask:

- Would we be comfortable supporting this name and behavior in a future release?
- Does this represent a real workflow, config type, or result type that another crate should use?
- Would exposing this make later cleanup harder?

---

## Public structs and methods

Note: Many of these are feature-based, so they should be pub even if users need to add the feature to get access to them.

 - All config structs to enable running the commands from rust
 - The shared cli_common structs needed to create/change config structs
 - All command runners
 - Main fragment + iteration classes for making it easy to iterate fragments in downstream software - although there are many different kinds of fragments (analysis-dependent), so this could get messy
 - Code for applying GC correction in an analysis
 - Code for applying genomic smoothing in an analysis
 - Code for applying blacklisting in an analysis

Candidates:
 - Main tile-run helpers? I would want to tile downstream software as well
 - Plotters? These are used in the commands as well, so perhaps they need to be public?
