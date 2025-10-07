#!/usr/bin/env python3
# Usage: python scripts/run_aucell_gmtmanual.py <dataset> <cpus>
# e.g.,   python scripts/run_aucell_gmtmanual.py pbmc 32

import argparse, os, time
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix

# ---- CLI ----
p = argparse.ArgumentParser(description="Run AUCell (omicverse) classification with manual GMT parsing")
p.add_argument("dataset", choices=["pbmc", "neuro", "cancer"])
p.add_argument("cpus", type=int)
args = p.parse_args()

# ---- Threads ----
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"]:
    os.environ[v] = str(args.cpus)
os.environ.setdefault("UCX_TLS", "tcp,sm,self")

# ---- Paths ----
REPO = Path(__file__).resolve().parents[1]
ds_dir = REPO / "rebuttal_datasets" / "AUCell"
runtime_dir = REPO / "running_time"; runtime_dir.mkdir(parents=True, exist_ok=True)
pred_dir    = REPO / "predictions";   pred_dir.mkdir(parents=True, exist_ok=True)

test_candidates = [
    ds_dir / f"AUCell_test_{args.dataset}.h5ad",
    ds_dir / f"{args.dataset}.h5ad",
]
test_h5ad = next((p for p in test_candidates if p.exists()), None)
if test_h5ad is None:
    raise FileNotFoundError(f"Missing test AnnData. Tried: {', '.join(map(str, test_candidates))}")

train_gmt  = ds_dir / f"AUCell_train_{args.dataset}.gmt"
pred_out   = pred_dir / f"AUCell_{args.dataset}.csv"
summary_fn = runtime_dir / f"AUCell_{args.dataset}.txt"

if not train_gmt.exists():
    raise FileNotFoundError(f"Missing training GMT: {train_gmt}")

# ---- Original GMT parser (description field optional) ----
def read_gmt(path: Path) -> dict:
    gs = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            genes = parts[2:] if len(parts) > 2 else parts[1:]
            genes = [g for g in genes if g]
            if genes:
                gs[name] = genes
    if not gs:
        raise ValueError(f"No signatures parsed from {path}")
    return gs

gmt = read_gmt(train_gmt)

# ---- Load data ----
adata = sc.read_h5ad(test_h5ad)
if not issparse(adata.X):
    adata.X = csr_matrix(adata.X)

# ---- AUCell scoring (EXACT same logic as your snippet) ----
import omicverse as ov

start_time = time.time()
ov.single.pathway_aucell(
    adata,
    pathway_names=list(gmt.keys()),
    pathways_dict=gmt
)

desired_cols = [f"{k}_aucell" for k in list(gmt.keys())]
existing_cols = [c for c in desired_cols if c in adata.obs.columns]
if not existing_cols:
    raise RuntimeError("AUCell produced no score columns in adata.obs.")

scores_df = adata.obs[existing_cols].copy()
scores_df.rename(columns=dict(zip(existing_cols, [c[:-7] for c in existing_cols])), inplace=True)

sorted_scores_idx = np.argsort(-scores_df.values, axis=1)
top_score_idx = sorted_scores_idx[:, 0]
second_top_score_idx = sorted_scores_idx[:, 1]  # kept for parity

adata.obs['AUCell_classification'] = [scores_df.columns[i] for i in top_score_idx]

if 'Cell type' in adata.obs.columns and pd.api.types.is_categorical_dtype(adata.obs['Cell type']):
    categories_order = adata.obs['Cell type'].cat.categories
    adata.obs['AUCell_classification'] = pd.Categorical(
        adata.obs['AUCell_classification'],
        categories=categories_order,
        ordered=True
    )

execution_time = time.time() - start_time

# ---- Save predictions ----
out = pd.DataFrame({
    "cell_id": adata.obs_names,
    "prediction": adata.obs["AUCell_classification"].astype(str).values,
})
out.to_csv(pred_out, index=False)
print(f"[OK] Saved predictions -> {pred_out}")

# ---- Append runtime ----
need_header = (not summary_fn.exists()) or (summary_fn.stat().st_size == 0)
with open(summary_fn, "a") as f:
    if need_header:
        f.write("TIME(s)\tCPUs\n")
    f.write(f"{execution_time}\t{args.cpus}\n")
print(f"[OK] Appended runtime -> {summary_fn}  (elapsed: {execution_time:.3f}s, cpus={args.cpus})")

