#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG="$SCRIPT_DIR/config"
PIPELINE="$SCRIPT_DIR/pipeline"

tail -n +2 "$CONFIG/datasets.tsv" | while IFS=$'\t' read -r ds h5ad ct_col gene_col species n_cells tier notes; do
  [[ "$tier" != "1" ]] && continue
  tail -n +2 "$CONFIG/gwas_scenarios.tsv" | while IFS=$'\t' read -r scen gs_src top_n desc; do
    OUT_BASE="$SCRIPT_DIR/results/$ds/$scen"
    echo "=== $ds | $scen ==="
    python "$PIPELINE/00_pseudobulk.py" --h5ad "$h5ad" --ct-col "$ct_col" --gene-col "$gene_col" --out "$OUT_BASE/pseudobulk"
  done
done
echo "Done."
