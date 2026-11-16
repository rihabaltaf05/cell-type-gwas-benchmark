#!/usr/bin/env python3
"""Compute mean expression per cell type. Usage: python 00_pseudobulk.py --h5ad PATH --ct-col COL --out DIR"""
import argparse, os
if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--ct-col", required=True)
    p.add_argument("--gene-col", default="var_names")
    p.add_argument("--out", required=True)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    print(f"[pseudobulk] h5ad={args.h5ad} ct_col={args.ct_col} out={args.out}")
