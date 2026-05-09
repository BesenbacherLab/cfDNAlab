#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"

cargo run --bin gen_cli_docs \
  --features cli,docs_gen,cmd_bam_to_bam,cmd_bam_to_frag,cmd_frag_to_bam,cmd_coverage_weights,cmd_fragment_count_weights,cmd_ends,cmd_fcoverage,cmd_gc_bias,cmd_lengths,cmd_midpoints,cmd_ref_gc_bias \
  -- \
  --out-dir "${repo_root}/website/docs/generated/cli" \
  --scope release

"${repo_root}/website/scripts/generate_release_notes.sh"
