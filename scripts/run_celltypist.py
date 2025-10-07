import argparse, os, time, pathlib, sys
import numpy as np
import scanpy as sc
import celltypist

# ------------------------ CLI ------------------------
ap = argparse.ArgumentParser(description="Train + predict with Celltypist on selected dataset")
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

data_dir = repo_dir / "rebuttal_datasets" / "Celltypist"
train_path = data_dir / f"Celltypist_train_{args.dataset}.h5ad"
test_path  = data_dir / f"Celltypist_test_{args.dataset}.h5ad"

if not train_path.exists() or not test_path.exists():
    print(f"[ERROR] Expected files not found:\n  {train_path}\n  {test_path}", file=sys.stderr)
    sys.exit(2)

runtime_dir = repo_dir / "running_time"
pred_dir    = repo_dir / "predictions"
runtime_dir.mkdir(parents=True, exist_ok=True)
pred_dir.mkdir(parents=True, exist_ok=True)

runtime_file = runtime_dir / f"Celltypist_{args.dataset}_{args.cpus}.txt"
pred_file    = pred_dir    / f"Celltypist_{args.dataset}.csv"
# ------------------------ Helpers ------------------------
def maybe_raw_to_adata(ad):
    if getattr(ad, "raw", None) is not None and ad.raw is not None:
        try:
            return ad.raw.to_adata()
        except Exception:
            return ad
    return ad

# ------------------------ Load data ------------------------
np.random.seed(3)
train = sc.read(str(train_path))
test  = sc.read(str(test_path))

train = maybe_raw_to_adata(train)
test  = maybe_raw_to_adata(test)

# Fixed label key
label_key = "Cell type"
if label_key not in train.obs.columns:
    print(f"[ERROR] '{label_key}' not found in train.obs columns: {list(train.obs.columns)}", file=sys.stderr)
    sys.exit(3)


# ------------------------ Train ------------------------
t0 = time.perf_counter()
model = celltypist.train(
    train,
    use_SGD=True,
    labels=label_key,
    n_jobs=args.cpus,
    feature_selection=True
)
elapsed = time.perf_counter() - t0

# ------------------------ Running time output (single file) ------------------------
# Summary file per dataset, tab-separated with header "TIME(s)\tCPUs"
summary_file = runtime_dir / f"Celltypist_{args.dataset}.txt"

# Create header if file doesn't exist or is empty
need_header = (not summary_file.exists()) or (summary_file.stat().st_size == 0)
with open(summary_file, "a") as f:
    if need_header:
        f.write("TIME(s)\tCPUs\n")
    f.write(f"{elapsed}\t{args.cpus}\n")

print(f"[OK] Appended runtime to: {summary_file}")

# ------------------------ Predict ------------------------
pred = celltypist.annotate(test, model=model)
yhat = getattr(pred, "predicted_labels", None)
if yhat is None:
    yhat = getattr(pred, "pred", None)
if yhat is None:
    raise RuntimeError("celltypist.annotate returned no labels (neither .predicted_labels nor .pred).")

test.obs["Prediction celltypist"] = yhat

# Save CSV: index,celltypist
with open(pred_file, "w") as f:
    f.write("index,celltypist\n")
    for idx, lbl in zip(test.obs_names, test.obs["Prediction celltypist"]):
        f.write(f"{idx},{lbl}\n")

print(f"[OK] Saved predictions to: {pred_file}")
