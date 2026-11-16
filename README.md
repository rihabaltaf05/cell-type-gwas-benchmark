# Cell-Type × GWAS Benchmarking Pipeline

Benchmarking EWCE, Cepo-MAGMA, scDRS, fgsea, and Cauchy CCT for trait-cell-type
association mapping. Replicates and extends Li et al. 2025 medRxiv.

Novel angle: method robustness under degraded GWAS power.

## Datasets (Tier 1)
| ID | Description | Cells | Species |
|---|---|---|---|
| tms_facs_brain | Tabula Muris Senis FACS brain | 20,666 | Mouse |

## Methods
| Step | Method | Script |
|---|---|---|
| 00 | Pseudobulk (shared) | pipeline/00_pseudobulk.py |
| 01 | EWCE | pipeline/01_ewce.R |
| 02 | Cepo + MAGMA | pipeline/02_cepo_magma.R |
| 03 | scDRS | pipeline/03_scdrs.sh |
| 04 | fgsea | pipeline/04_fgsea.R |
| 05 | Cauchy CCT | pipeline/05_cct.R |

## GWAS Power Scenarios
| ID | Genes per trait |
|---|---|
| full_power | ~1000 |
| top500 | 500 |
| top200 | 200 |
| top100 | 100 |

## Usage
bash run_all.sh
