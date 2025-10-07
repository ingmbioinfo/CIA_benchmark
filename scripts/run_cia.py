#!/usr/bin/env python3
# CIA classify runner — same interface as Celltypist:
# Usage: python scripts/run_cia.py <dataset> <cpus>
# Example: python scripts/run_cia.py pbmc 32

import argparse, os, time, io
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix

# ---------------- CLI (only dataset, cpus) ----------------
p = argparse.ArgumentParser(description="Run CIA classification")
p.add_argument("dataset", choices=["pbmc", "cancer", "neuro"])
p.add_argument("cpus", type=int)
args = p.parse_args()

# ---------------- Threads to match -c ----------------
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"]:
    os.environ[v] = str(args.cpus)
# Optional: avoid UCX/IB noise on nodes without IB
os.environ.setdefault("UCX_TLS", "tcp,sm,self")

# ---------------- Paths (relative; script is in repo/scripts) ----------------
REPO = Path(__file__).resolve().parents[1]
ds_dir = REPO / "rebuttal_datasets" / "CIA"
runtime_dir = REPO / "running_time"; runtime_dir.mkdir(parents=True, exist_ok=True)
pred_dir = REPO / "predictions";    pred_dir.mkdir(parents=True, exist_ok=True)

test_h5ad  = ds_dir / f"CIA_test_{args.dataset}.h5ad"
train_gmt  = ds_dir / f"CIA_train_{args.dataset}.gmt"
pred_out   = pred_dir / f"CIA_{args.dataset}.csv"
summary_fn = runtime_dir / f"CIA_{args.dataset}.txt"

# ---------------- Import CIA ----------------
from cia import investigate  # will raise if env not set

# ---------------- Helpers ----------------
def ensure_sparse_and_raw(adata):
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    if adata.raw is None:
        adata.raw = adata  # snapshot current X/var as raw
    return adata

def assert_not_html(path: Path):
    with open(path, "rb") as f:
        head = f.read(1024).lower()
    if b"<!doctype html" in head or b"<html" in head:
        raise ValueError(
            f"{path} looks like an HTML page (probably a GitHub 'blob'). "
            "Download the **Raw** file and try again."
        )

def pad_gmt_for_pandas_c_engine(src: Path, dst: Path):
    """
    Read a GMT (tab-separated, variable-length rows) and write a padded TSV
    where all rows have the same number of fields. This avoids pandas C-engine
    'Expected N fields, saw M' errors inside CIA.
    """
    lines = []
    with open(src, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.strip() == "":
                continue
            parts = line.split("\t")
            # require at least 3 fields: name, desc, >=1 gene
            if len(parts) < 3:
                continue
            lines.append(parts)
    if not lines:
        raise ValueError(f"No valid signature lines found in {src}")
    max_len = max(len(r) for r in lines)
    with open(dst, "w", encoding="utf-8") as out:
        for r in lines:
            if len(r) < max_len:
                r = r + [""] * (max_len - len(r))
            out.write("\t".join(r) + "\n")
    return dst

# ---------------- Load data ----------------
if not test_h5ad.exists():
    raise FileNotFoundError(f"Missing test AnnData: {test_h5ad}")
if not train_gmt.exists():
    raise FileNotFoundError(f"Missing training GMT: {train_gmt}")

assert_not_html(train_gmt)
padded_gmt = runtime_dir / f"{train_gmt.stem}__padded.gmt"
pad_gmt_for_pandas_c_engine(train_gmt, padded_gmt)

adata = sc.read_h5ad(test_h5ad)
ensure_sparse_and_raw(adata)

# ---------------- Classify ----------------
t0 = time.perf_counter()
investigate.CIA_classify(
    data=adata,
    signatures_input=str(padded_gmt),  # pass padded file to keep pandas C-engine happy
    n_cpus=args.cpus,
)
elapsed = time.perf_counter() - t0

# ---------------- Save predictions ----------------
label_col = "CIA prediction"  # CIA default
if label_col not in adata.obs.columns:
    # fallback: try to find something that looks like CIA's output
    cand = [c for c in adata.obs.columns if "cia" in c.lower() and "pred" in c.lower()]
    if not cand:
        raise KeyError(f"Prediction column '{label_col}' not found in adata.obs.")
    label_col = cand[0]

df = pd.DataFrame({
    "cell_id": adata.obs_names,
    "prediction": adata.obs[label_col].astype(str).values,
}, index=None)
df.to_csv(pred_out, index=False)
print(f"[OK] Saved predictions -> {pred_out}")

# ---------------- Append runtime (single per-dataset file) ----------------
need_header = (not summary_fn.exists()) or (summary_fn.stat().st_size == 0)
with open(summary_fn, "a") as f:
    if need_header:
        f.write("TIME(s)\tCPUs\n")
    f.write(f"{elapsed}\t{args.cpus}\n")
print(f"[OK] Appended runtime -> {summary_fn}  (elapsed: {elapsed:.3f}s, cpus={args.cpus})")

