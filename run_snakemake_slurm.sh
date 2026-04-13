#!/usr/bin/env bash
set -euo pipefail
snakemake   --use-conda   --jobs 20   --latency-wait 60   --rerun-incomplete  --default-resources threads=5 mem_mb=8000 time_min=120  --cluster "sbatch -c {threads} --mem={resources.mem_mb}M -t {resources.time_min}"   -p "$@"
