#!/usr/bin/env bash
# run_tier1.sh — Tier 1 benchmark: TMS brain x all methods, in parallel
# Usage: bash run_tier1.sh [--dry-run]
set -euo pipefail

BENCH=/home/rihabaltaf/benchmarking
PIPE=/home/rihabaltaf/cell-type-gwas-benchmark/pipeline
LOGS=$BENCH/logs/tier1
PYTHON=/home/rihabaltaf/envs/scdrs/bin/python3
RSCRIPT=/home/rihabaltaf/envs/r_ewce/bin/Rscript
H5AD=$BENCH/data/tms_facs_brain.h5ad
GS=$BENCH/data/benchmark_traits_scdrs.gs
CT_COL=cell_ontology_class
SPECIES=mouse
EXPORT_DIR=$BENCH/tier1_data
MAGMA_GENES_DIR=$BENCH/magma/genes_out
OUT_EWCE=$BENCH/results/ewce
OUT_CEPO=$BENCH/results/cepo_magma
OUT_SCDRS=$BENCH/results/scdrs
OUT_FGSEA=$BENCH/results/fgsea
OUT_CCT=$BENCH/results/cct
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true
mkdir -p "$LOGS"
for d in "$OUT_EWCE" "$OUT_CEPO" "$OUT_SCDRS" "$OUT_FGSEA" "$OUT_CCT"; do mkdir -p "$d"; done

echo "=== Step 1: Export data ==="
if $DRY_RUN; then
    echo "[DRY] 00_pseudobulk.py --h5ad $H5AD --ct-col $CT_COL --gs $GS --out $EXPORT_DIR"
else
    "$PYTHON" "$PIPE/00_pseudobulk.py" --h5ad "$H5AD" --ct-col "$CT_COL" --gs "$GS" --out "$EXPORT_DIR" \
        2>&1 | tee "$LOGS/00_export.log"
fi

echo ""; echo "=== Step 2: Methods (parallel) ==="
launch(){
    local name="$1"; shift; local log="$LOGS/${name}.log"
    if $DRY_RUN; then echo "[DRY] $name: $*"; echo "0"; return; fi
    echo "[launch] $name -> $log"
    nohup "$@" > "$log" 2>&1 & echo $!
}
PID_EWCE=$(launch "01_ewce" "$RSCRIPT" "$PIPE/01_ewce.R" --datadir "$EXPORT_DIR" --outdir "$OUT_EWCE" --reps 1000 --species "$SPECIES")
PID_CEPO=$(launch "02_cepo_magma" "$RSCRIPT" "$PIPE/02_cepo_magma.R" --datadir "$EXPORT_DIR" --genes-dir "$MAGMA_GENES_DIR" --outdir "$OUT_CEPO" --species "$SPECIES")
PID_SCDRS=$(launch "03_scdrs" bash "$PIPE/03_scdrs.sh" --h5ad "$H5AD" --gs "$GS" --ct-col "$CT_COL" --outdir "$OUT_SCDRS" --h5ad-species "$SPECIES" --gs-species human --python "$PYTHON")
PID_FGSEA=$(launch "04_fgsea" bash -c "
    for i in \$(seq 1 60); do [[ -f '$OUT_EWCE/ctd_specificity.tsv' ]] && break; sleep 60; done
    $RSCRIPT $PIPE/04_fgsea.R --specificity '$OUT_EWCE/ctd_specificity.tsv' --gs '$GS' --outdir '$OUT_FGSEA'")
echo "PIDs: EWCE=$PID_EWCE  CEPO=$PID_CEPO  scDRS=$PID_SCDRS  fGSEA=$PID_FGSEA"
echo "Monitor: tail -f $LOGS/01_ewce.log"
if ! $DRY_RUN; then
    for pid in "$PID_EWCE" "$PID_CEPO" "$PID_SCDRS" "$PID_FGSEA"; do
        wait "$pid" && echo "PID $pid OK" || echo "PID $pid FAILED"
    done
fi

echo ""; echo "=== Step 3: CCT ==="
CCT_ARGS="--outdir $OUT_CCT"
[[ -f "$OUT_EWCE/ewce_results.csv"        ]] && CCT_ARGS="$CCT_ARGS --ewce  $OUT_EWCE/ewce_results.csv"
[[ -f "$OUT_FGSEA/fgsea_results.csv"      ]] && CCT_ARGS="$CCT_ARGS --fgsea $OUT_FGSEA/fgsea_results.csv"
[[ -f "$OUT_CEPO/cepo_magma_results.csv"  ]] && CCT_ARGS="$CCT_ARGS --cepo  $OUT_CEPO/cepo_magma_results.csv"
[[ -d "$OUT_SCDRS/mc_stats"               ]] && CCT_ARGS="$CCT_ARGS --scdrs $OUT_SCDRS/mc_stats"
if $DRY_RUN; then echo "[DRY] 05_cct.R $CCT_ARGS"
else "$RSCRIPT" "$PIPE/05_cct.R" $CCT_ARGS 2>&1 | tee "$LOGS/05_cct.log"; fi
echo "=== Done ==="
