#!/bin/bash
# 03_scdrs.sh — scDRS scoring + downstream
# Usage: bash 03_scdrs.sh <h5ad> <gs_file> <ct_col> <out_dir>
H5AD=$1; GS=$2; CT_COL=$3; OUT=$4
mkdir -p "$OUT/scores" "$OUT/downstream"
echo "[scdrs] h5ad=$H5AD gs=$GS ct_col=$CT_COL out=$OUT"
