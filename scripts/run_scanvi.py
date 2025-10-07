#!/usr/bin/env python3
# Usage: python scripts/run_scanvi.py <dataset> <cpus>
# Example: python scripts/run_scanvi.py neuro 64

import argparse, os, time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import anndata as ad
import h5py
from scipy.sparse import issparse, csr_matrix


# ---------------- CLI (same minimal interface) ----------------
p = argparse.ArgumentParser(description="Run scANVI semi-supervised classification")
p.add_argument("dataset", choices=["pbmc", "neuro", "cancer"])
p.add_argument("cpus", type=int)
args = p.parse_args()

# ---------------- Threads (match -c) ----------------
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"]:
    os.environ[v] = str(args.cpus)
os.environ.setdefault("UCX_TLS", "tcp,sm,self")  # quiet IB noise on non-IB nodes
torch.set_num_threads(args.cpus)

# ---------------- Paths (script lives in repo/scripts) ----------------
REPO = Path(__file__).resolve().parents[1]
ds_dir = REPO / "rebuttal_datasets" / "scANVI"
runtime_dir = REPO / "running_time"; runtime_dir.mkdir(parents=True, exist_ok=True)
pred_dir    = REPO / "predictions";    pred_dir.mkdir(parents=True, exist_ok=True)

train_h5ad = ds_dir / f"scANVI_train_{args.dataset}.h5ad"
test_h5ad  = ds_dir / f"scANVI_test_{args.dataset}.h5ad"
pred_out   = pred_dir / f"scANVI_{args.dataset}.csv"
summary_fn = runtime_dir / f"scANVI_{args.dataset}.txt"

if not train_h5ad.exists():
    raise FileNotFoundError(f"Missing training AnnData: {train_h5ad}")
if not test_h5ad.exists():
    raise FileNotFoundError(f"Missing test AnnData: {test_h5ad}")

# ---------------- H5AD sanitization (fix reader hiccups) ----------------
def sanitize_h5ad(path: Path):
    """Remove known problematic entries so AnnData can read cleanly."""
    with h5py.File(path, "r+") as f:
        if "uns" in f and "log1p" in f["uns"]:
            del f["uns"]["log1p"]
            print(f"[fix] Removed /uns/log1p from {path}")

sanitize_h5ad(train_h5ad)
sanitize_h5ad(test_h5ad)

# ---------------- Load ----------------
train = sc.read_h5ad(train_h5ad)
test  = sc.read_h5ad(test_h5ad)

# Keep matrices sparse for memory friendliness
if not issparse(train.X): train.X = csr_matrix(train.X)
if not issparse(test.X):  test.X  = csr_matrix(test.X)

# ---------------- Label column ----------------
# Prefer 'CELLTYPE'; fall back if needed.
label_candidates = ["CELLTYPE", "Cell type", "celltype", "labels", "label"]
label_col = next((c for c in label_candidates if c in train.obs.columns), None)
if label_col is None:
    raise KeyError(f"No label column found in training obs. Tried: {', '.join(label_candidates)}")

# Make sure the "Unknown" category exists for test labels
train.obs["SCANVI_LABELS"] = train.obs[label_col].astype("category")
test.obs["SCANVI_LABELS"]  = pd.Categorical(["Unknown"] * test.n_obs, categories=list(train.obs["SCANVI_LABELS"].cat.categories) + ["Unknown"])

# Add a batch/split column before concatenation
train.obs["batch"] = "train"
test.obs["batch"]  = "test"

# Align genes by outer-joining features (missing genes -> 0)
adata = ad.concat([train, test], join="outer", label=None, fill_value=0)

# ---------------- scANVI ----------------
import scvi
scvi.settings.seed = 1

# Ensure SCANVI_LABELS exist and set test cells to "Unknown"
truth_candidates = ["CELLTYPE","Cell type","celltype","labels","label"]
if "SCANVI_LABELS" not in adata.obs.columns:
    src = next((c for c in truth_candidates if c in adata.obs.columns), None)
    if src is not None:
        adata.obs["SCANVI_LABELS"] = adata.obs[src].astype(str)
    else:
        # fallback: create empty labels; SCANVI will treat all as unlabeled
        adata.obs["SCANVI_LABELS"] = "Unknown"

# Your original split logic: batch == "test"
is_test = (adata.obs["batch"].astype(str).values == "test")
# mark test cells as unlabeled for SCANVI
adata.obs.loc[is_test, "SCANVI_LABELS"] = "Unknown"

# Pretrain SCVI on all cells (CPU)
scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
scvi_model = scvi.model.SCVI(adata, n_latent=30)

t0 = time.perf_counter()
scvi_model.train(
    max_epochs=20,
    accelerator="cpu",
    devices=1,
    batch_size=1024,
    plan_kwargs={"weight_decay": 0.0},
    enable_progress_bar=False,
)

# Train SCANVI with labels (train) + unlabeled (test)
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key="SCANVI_LABELS",
    unlabeled_category="Unknown",
)
scanvi_model.train(
    max_epochs=20,
    accelerator="cpu",
    devices=1,
    batch_size=1024,
    enable_progress_bar=False,
)
elapsed = time.perf_counter() - t0

# ---------------- Predict on test split ONLY ----------------
adata_test = adata[is_test].copy()
pred_test = scanvi_model.predict(adata_test)  # ndarray of labels for test only

out = pd.DataFrame({
    "cell_id": adata_test.obs_names,
    "prediction": np.asarray(pred_test).astype(str)
})

# Add truth if present in the (subset) AnnData
truth_col = next((c for c in truth_candidates if c in adata_test.obs.columns), None)
if truth_col is not None:
    out["truth"] = adata_test.obs[truth_col].astype(str).values

# Preserve the test cells' original order
out = out.set_index("cell_id").loc[adata_test.obs_names].reset_index()

out.to_csv(pred_out, index=False)
print(f"[OK] Saved predictions -> {pred_out}")

# ---------------- Append runtime (single per-dataset file) ----------------
need_header = (not summary_fn.exists()) or (summary_fn.stat().st_size == 0)
with open(summary_fn, "a") as f:
    if need_header:
        f.write("TIME(s)\tCPUs\n")
    f.write(f"{elapsed}\t{args.cpus}\n")
print(f"[OK] Appended runtime -> {summary_fn}  (elapsed: {elapsed:.3f}s, cpus={args.cpus})")

