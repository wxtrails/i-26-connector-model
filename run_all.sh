#!/usr/bin/env bash
cd "$(dirname "$0")" && exec uv run python run_all.py "$@"
