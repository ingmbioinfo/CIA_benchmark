#!/usr/bin/env python3

import argparse, os, time, pathlib, sys
import numpy as np
import scanpy as sc
import scBalance as sb

# ------------------------ CLI ------------------------
ap = argparse.ArgumentParser(description="Train + predict with scBalance on selected dataset")
ap.add_argument("--dataset", required=True, choices=["neuro","pbmc","cancer"],
                help="Which dataset split to use")
ap.add_argument("--cpus", type=int, required=True, help="Number of CPUs/threads to use")
args = ap.parse_args()

# Threads for BLAS/NumExpr/etc.
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = str(args.cpus)

# ------------------------ Resolve relative paths ------------------------
script_dir = pathlib.Path(__file__).resolve().parent
repo_dir   = script_dir.parent  # one level up from scripts/

data_dir = repo_dir / "rebuttal_datasets" / "scBalance"
train_path = data_dir / f"scBalance_train_{args.dataset}.h5ad"
test_path  = data_dir / f"scBalance_test_{args.dataset}.h5ad"

if not train_path.exists() or not test_path.exists():
    print(f"[ERROR] Expected files not found:\n  {train_path}\n  {test_path}", file=sys.stderr)
    sys.exit(2)

runtime_dir = repo_dir / "running_time"
pred_dir    = repo_dir / "predictions"
runtime_dir.mkdir(parents=True, exist_ok=True)
pred_dir.mkdir(parents=True, exist_ok=True)

runtime_file = runtime_dir / f"scBalance_{args.dataset}_{args.cpus}.txt"
pred_file    = pred_dir    / f"scBalance_{args.dataset}.csv"
# ------------------------ Helpers ------------------------
def maybe_raw_to_adata(ad):
    if getattr(ad, "raw", None) is not None and ad.raw is not None:
        try:
            return ad.raw.to_adata()
        except Exception:
            return ad
    return ad

# ------------------------ Load data ------------------------
#np.random.seed(3)
train = sc.read(str(train_path))
test  = sc.read(str(test_path))

common_genes = train.var_names.intersection(test.var_names)
X_train = train.to_df()[common_genes]
X_test = test.to_df()[common_genes]

# Labels with numeric codes
y_train = train.obs[["Cell type"]].rename(columns={"Cell type": "Label"})
y_test = test.obs[["Cell type"]].rename(columns={"Cell type": "Label"})

label_mapping = dict(
    enumerate(train.obs["Cell type"].astype("category").cat.categories)
)

y_train["Label"] = y_train["Label"].astype("category").cat.codes
y_test["Label"] = y_test["Label"].astype("category").cat.codes

# cast data to float32
X_test = X_test.astype('float32')
X_train = X_train.astype('float32')

# --- Run scBalance and measure time ---
start = time.time()
pred_label = sb.scBalance(
    X_test, X_train, y_train, processing_unit="cpu"
)
end = time.time()
runtime = end - start

# ------------------------ Running time output (single file) ------------------------
# Summary file per dataset, tab-separated with header "TIME(s)\tCPUs"
summary_file = runtime_dir / f"scBalance_{args.dataset}.txt"

# Create header if file doesn't exist or is empty
need_header = (not summary_file.exists()) or (summary_file.stat().st_size == 0)
with open(summary_file, "a") as f:
    if need_header:
        f.write("TIME(s)\tCPUs\n")
    f.write(f"{runtime}\t{args.cpus}\n")

print(f"[OK] Appended runtime to: {summary_file}")

# Save CSV: index,scBalance
df_results = y_test.copy()
df_results["pred_label"] = pred_label
df_results["pred_label"] = df_results["pred_label"].map(label_mapping)
df_results.drop(columns=['Label'], inplace=True)

df_results.to_csv(pred_file)

print(f"[OK] Saved predictions to: {pred_file}")

