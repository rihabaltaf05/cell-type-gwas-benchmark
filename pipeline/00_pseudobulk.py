#!/usr/bin/env python3
"""
00_pseudobulk.py — Export h5ad to MTX format for R-based pipeline methods.
"""
import argparse, os, sys
import scipy.io, scipy.sparse as sp
import pandas as pd

CHUNK = 5000

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad",    required=True)
    p.add_argument("--ct-col",  required=True)
    p.add_argument("--gs",      required=True)
    p.add_argument("--out",     required=True)
    p.add_argument("--use-raw", action="store_true")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    import scanpy as sc
    print(f"[00] Loading {args.h5ad} (backed='r') ...", flush=True)
    adata = sc.read_h5ad(args.h5ad, backed='r').to_memory()
    print(f"     {adata.shape[0]} cells x {adata.shape[1]} genes", flush=True)

    ct_col = args.ct_col
    if ct_col not in adata.obs.columns:
        print(f"ERROR: '{ct_col}' not in obs. Have: {list(adata.obs.columns)}")
        sys.exit(1)

    if args.use_raw:
        if adata.raw is None: print("ERROR: --use-raw but adata.raw is None"); sys.exit(1)
        print("     Using adata.raw.X", flush=True)
        X = adata.raw.X; genes = list(adata.raw.var_names)
    else:
        print("     Using adata.X", flush=True)
        X = adata.X; genes = list(adata.var_names)

    X_T = sp.csr_matrix(X.T)
    print(f"     Writing counts.mtx ({X_T.shape[0]} genes x {X_T.shape[1]} cells) ...", flush=True)
    scipy.io.mmwrite(os.path.join(args.out, "counts.mtx"), X_T)
    pd.Series(genes).to_csv(os.path.join(args.out, "genes.tsv"), index=False, header=False)
    pd.Series(list(adata.obs_names)).to_csv(os.path.join(args.out, "barcodes.tsv"), index=False, header=False)

    meta = adata.obs[[ct_col]].copy(); meta.columns = ["cell_type"]
    meta.to_csv(os.path.join(args.out, "cell_metadata.csv"))
    ct_counts = meta["cell_type"].value_counts()
    print(f"     {meta['cell_type'].nunique()} cell types:", flush=True)
    for ct, n in ct_counts.items(): print(f"       {ct}: {n}", flush=True)

    print(f"[00] Parsing gene lists from {args.gs} ...", flush=True)
    gene_lists = {}
    with open(args.gs) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split('\t')
            if len(parts) < 2: continue
            gene_lists[parts[0]] = [g.split(':')[0] for g in parts[1].replace(',', ' ').split() if ':' in g]
    with open(os.path.join(args.out, "gene_lists.tsv"), 'w') as f:
        for trait, g in gene_lists.items():
            f.write(f"{trait}\t{','.join(g)}\n")
    print(f"     {len(gene_lists)} traits written", flush=True)

    bg = sorted(set(genes))
    pd.Series(bg).to_csv(os.path.join(args.out, "background_genes.tsv"), index=False, header=False)
    print(f"     Background: {len(bg)} genes", flush=True)
    print("[00] Export done ✓", flush=True)

if __name__ == "__main__":
    main()
