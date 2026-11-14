#!/usr/bin/env bash
# Reproduce the complete AlphaGenome 2 GHSR analysis.
#
# Requires:
#   * Python >= 3.11 with requirements.txt installed
#   * ALPHA_GENOME_API_KEY exported in the environment
#   * outbound HTTPS to gdmscience.googleapis.com, rest.ensembl.org and
#     api.ncbi.nlm.nih.gov
#
# Usage:
#   bash scripts/run_all.sh            # full run, including the ~25 min ISM step
#   FAST=1 bash scripts/run_all.sh     # skip the ISM and sequence-length steps

set -euo pipefail

cd "$(dirname "$0")/.."
export PYTHONPATH="${PYTHONPATH:-}:src"

if [[ -z "${ALPHA_GENOME_API_KEY:-}" ]]; then
  echo "ERROR: ALPHA_GENOME_API_KEY is not set." >&2
  echo "Export the AlphaGenome API key before running; never store it in a file." >&2
  exit 1
fi

step() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

step "Phase I - variant annotation and coordinate verification"
python -m ghsr_ag2.phase1_annotation

step "Phase II - condition selection and data-quality audit"
python -m ghsr_ag2.phase2_conditions

step "Phase III - AlphaGenome 2 variant-effect prediction"
python -m ghsr_ag2.phase3_scoring

step "Phase III-b - REF vs ALT base-resolution tracks"
python -m ghsr_ag2.phase3b_tracks

step "Phase IV - robustness and quality control"
if [[ "${FAST:-0}" == "1" ]]; then
  echo "FAST=1: skipping the ISM reference distribution and the"
  echo "sequence-context sensitivity analysis."
  python -m ghsr_ag2.phase4_qc
else
  python -m ghsr_ag2.phase4_qc --ism --seqlen
fi

step "Phase IV-b - cross-sequence-context stability"
if [[ "${FAST:-0}" == "1" ]]; then
  echo "FAST=1: skipped (requires the sequence-context runs from Phase IV)."
else
  python -m ghsr_ag2.phase4b_stability
fi

step "Phase V - integration with the existing experimental findings"
python -m ghsr_ag2.phase5_integration

step "Phase VI - figures"
python -m ghsr_ag2.phase6_figures

step "CSV exports of the raw model output"
python scripts/export_csv.py

step "Done"
echo "Tables:  results/tables/"
echo "Figures: results/figures/"
echo "Raw model output: data/raw/phase3/phase3_raw_scores_main.parquet"
echo "Provenance:      data/raw/provenance_*.json"
