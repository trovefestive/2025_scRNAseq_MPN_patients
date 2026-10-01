#!/usr/bin/env python3
"""
Stage 3 (downstream-analysis-plan.md): per-sample Scrublet doublet detection.

For each SRR, reads Cell Ranger's filtered_feature_bc_matrix/, restricts to the
QC-passing barcodes from Stage 2 (seurat_qc/<SRR>_barcodes.txt), runs Scrublet,
and writes per-cell doublet_score + predicted_doublet to
seurat_qc/<SRR>_scrublet.csv. Doublets are NOT dropped here -- that happens in
apply_doublet_filter.R, pre-merge, per the plan.
"""
import os
import sys
import gzip
import scipy.io
import scipy.sparse
import numpy as np
import pandas as pd
import scrublet as scr

PROJ_DIR = os.environ.get("PROJ_DIR", os.getcwd())  # run from the project root
os.chdir(PROJ_DIR)

with open("runs.txt") as f:
    runs = [l.strip() for l in f if l.strip()]

summary_rows = []

for srr in runs:
    print(f"=== {srr} ===", flush=True)
    mtx_dir = os.path.join("cellranger_count", srr, "outs", "filtered_feature_bc_matrix")
    barcode_file = os.path.join("seurat_qc", f"{srr}_barcodes.txt")

    if not os.path.isdir(mtx_dir) or not os.path.exists(barcode_file):
        print(f"  missing input for {srr}, skipping", flush=True)
        continue

    # Cell Ranger mtx: genes x cells -- transpose to cells x genes for scrublet
    mat = scipy.io.mmread(os.path.join(mtx_dir, "matrix.mtx.gz")).tocsc().T.tocsr()
    with gzip.open(os.path.join(mtx_dir, "barcodes.tsv.gz"), "rt") as f:
        all_barcodes = [l.strip() for l in f]

    qc_barcodes = set(pd.read_csv(barcode_file, header=None)[0])
    keep_idx = [i for i, bc in enumerate(all_barcodes) if bc in qc_barcodes]
    kept_barcodes = [all_barcodes[i] for i in keep_idx]
    mat_qc = mat[keep_idx, :]

    n_cells = mat_qc.shape[0]
    print(f"  {n_cells} QC-passing cells (of {len(all_barcodes)} total)", flush=True)

    if n_cells < 50:
        print(f"  WARNING: only {n_cells} cells, scrublet may be unreliable; running anyway", flush=True)

    scrub = scr.Scrublet(mat_qc, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=min(30, n_cells - 1)
    )

    out_df = pd.DataFrame({
        "barcode": kept_barcodes,
        "doublet_score": doublet_scores,
        "predicted_doublet": predicted_doublets,
    })
    out_df.to_csv(os.path.join("seurat_qc", f"{srr}_scrublet.csv"), index=False)

    n_doublets = int(np.sum(predicted_doublets))
    threshold = scrub.threshold_ if hasattr(scrub, "threshold_") else np.nan
    summary_rows.append({
        "SRR": srr,
        "n_cells": n_cells,
        "n_predicted_doublets": n_doublets,
        "pct_doublets": round(100 * n_doublets / n_cells, 2),
        "scrublet_threshold": threshold,
    })
    print(f"  {n_doublets}/{n_cells} predicted doublets ({100*n_doublets/n_cells:.2f}%)", flush=True)

summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join("seurat_qc", "scrublet_summary.csv"), index=False)
print("\n=== Scrublet summary ===")
print(summary.to_string(index=False))
print("\nStage 3a (Scrublet doublet detection) complete.")
