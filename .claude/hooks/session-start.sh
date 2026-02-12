#!/bin/bash
set -euo pipefail

# Only run in remote (Claude Code on the web) environments
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# Install Perl::Critic for linting Perl code
if ! perl -e 'use Perl::Critic' &>/dev/null; then
  apt-get update -qq || true
  apt-get install -y -qq libperl-critic-perl
fi
