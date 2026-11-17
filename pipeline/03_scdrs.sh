#!/usr/bin/env bash
# 03_scdrs.sh — scDRS scoring + per-cell-type statistics
set -euo pipefail
H5AD=""; GS=""; CT_COL="cell_type"; OUTDIR=""; H5AD_SPECIES="mouse"; GS_SPECIES="human"
N_CTRL=1000; PYTHON="/home/rihabaltaf/envs/scdrs/bin/python3"
while [[ $# -gt 0 ]]; do case "$1" in
    --h5ad)         H5AD="$2";         shift 2 ;;
    --gs)           GS="$2";           shift 2 ;;
    --ct-col)       CT_COL="$2";       shift 2 ;;
    --outdir)       OUTDIR="$2";       shift 2 ;;
    --h5ad-species) H5AD_SPECIES="$2"; shift 2 ;;
    --gs-species)   GS_SPECIES="$2";   shift 2 ;;
    --n-ctrl)       N_CTRL="$2";       shift 2 ;;
    --python)       PYTHON="$2";       shift 2 ;;
    *) echo "Unknown: $1"; exit 1 ;;
esac; done
[[ -z "$H5AD" || -z "$GS" || -z "$OUTDIR" ]] && { echo "ERROR: --h5ad --gs --outdir required"; exit 1; }
mkdir -p "$OUTDIR" "$OUTDIR/mc_stats"
echo "[03_scdrs] h5ad=$H5AD  species=$H5AD_SPECIES  gs-species=$GS_SPECIES"
echo "[03_scdrs] Scoring cells ..."
/home/rihabaltaf/envs/scdrs/bin/scdrs compute_score \
    --h5ad-file     "$H5AD" \
    --h5ad-species  "$H5AD_SPECIES" \
    --gs-file       "$GS" \
    --gs-species    "$GS_SPECIES" \
    --flag-filter-data  True \
    --flag-raw-count    True \
    --n-ctrl        "$N_CTRL" \
    --flag-return-ctrl-raw-score  True \
    --flag-return-ctrl-norm-score True \
    --out-folder    "$OUTDIR"
echo "[03_scdrs] score_cell complete"
echo "[03_scdrs] Computing per-cell-type stats ..."
for SCORE_FILE in "$OUTDIR"/*.score.gz; do
    [[ -f "$SCORE_FILE" ]] || continue
    TRAIT=$(basename "$SCORE_FILE" .score.gz)
    echo "  mc_stats: $TRAIT"
    /home/rihabaltaf/envs/scdrs/bin/scdrs perform_downstream \
        --h5ad-file    "$H5AD" \
        --score-file   "$SCORE_FILE" \
        --out-folder   "$OUTDIR/mc_stats" \
        --group-analysis "$CT_COL" \
        --flag-filter-data True \
        || echo "  WARNING: perform_downstream failed for $TRAIT"
done
echo "[03_scdrs] Done ✓"
