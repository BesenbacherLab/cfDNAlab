#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
changelog_file="${repo_root}/CHANGELOG.md"
release_notes_file="${repo_root}/website/docs/generated/release-notes.md"

if [[ ! -f "${changelog_file}" ]]; then
  echo "Expected changelog file at ${changelog_file}" >&2
  exit 1
fi

mkdir -p "$(dirname "${release_notes_file}")"

{
  echo "<!-- AUTO-GENERATED FILE - DO NOT EDIT -->"
  echo "<!-- Source: CHANGELOG.md -->"
  echo "<!-- Generated path: website/docs/generated/release-notes.md -->"
  echo
  echo "# Release notes"
  echo
  awk '
    NR == 1 && $0 ~ /^# / { next }
    $0 ~ /^<br ?\/?>$/ { next }
    $0 == "" {
      if (!blank_line) {
        print
      }
      blank_line = 1
      next
    }
    {
      blank_line = 0
      print
    }
  ' "${changelog_file}"
} > "${release_notes_file}"

echo "Generated ${release_notes_file}"
