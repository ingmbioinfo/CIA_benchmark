#!/usr/bin/env python3
import os
import sys
import glob
import pandas as pd

# --- deps ---------------------------------------------------------------------
try:
    import scanpy as sc
except Exception:
    sys.stderr.write("ERROR: scanpy is required (conda/pip install scanpy)\n")
    raise

# Prefer cia.report; fallback to cia.reports if older version
try:
    from cia.report import compute_classification_metrics
except Exception:
    from cia.reports import compute_classification_metrics  # type: ignore

# --- helpers ------------------------------------------------------------------
def classifier_label_from_filename(path: str) -> str:
    """
    Map filename stem to a readable classifier label.
    - 'CIA_R_*.csv' -> 'CIA_R'
    - 'CIA_*.csv'   -> 'CIA_Python'
    - otherwise: token before first underscore (e.g., 'AUCell', 'Celltypist', 'SingleR', 'scANVI', 'scBalance')
    """
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]  # e.g. "CIA_R_neuro" or "CIA_neuro"
    if stem.startswith("CIA_R"):
        return "CIA_R"
    if stem.startswith("CIA"):
        return "CIA_Python"
    return stem.split("_", 1)[0]

def read_prediction_firstcol(path: str) -> pd.Series:
    """Read CSV with first column as index; keep only the first data column (drop the rest)."""
    df = pd.read_csv(path, index_col=0)
    if df.shape[1] == 0:
        raise ValueError(f"{os.path.basename(path)} has no data columns.")
    if df.shape[1] > 1:
        # keep only the first non-index column
        df = df.iloc[:, [0]]
    s = df.iloc[:, 0].copy()
    s.index = s.index.astype(str)
    return s

def find_ref_obs(obs: pd.DataFrame) -> str:
    for col in ("Cell type", "CELLTYPE"):
        if col in obs.columns:
            return col
    raise KeyError("Neither 'Cell type' nor 'CELLTYPE' found in adata.obs")

# --- main ---------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        sys.stderr.write(
            "Usage: python3 perfomances.py <dataset> [unassigned_label]\n"
            "  <dataset> in {neuro, pbmc, cancer}\n"
            "  [unassigned_label] optional, e.g. 'Unknown' (default: '')\n"
        )
        sys.exit(2)

    dataset = sys.argv[1].strip().lower()
    unassigned_label = sys.argv[2] if len(sys.argv) >= 3 else ""

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.normpath(os.path.join(script_dir, ".."))

    adata_path = os.path.join(repo_root, "rebuttal_datasets", "CIA", f"CIA_test_{dataset}.h5ad")
    pred_dir  = os.path.join(repo_root, "predictions")
    out_dir   = os.path.join(repo_root, "performances")
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.isfile(adata_path):
        sys.stderr.write(f"ERROR: missing test set: {adata_path}\n")
        sys.exit(1)

    print(f"[INFO] Reading: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    adata.obs.index = adata.obs.index.astype(str)

    ref_col = find_ref_obs(adata.obs)
    print(f"[INFO] Using ref_obs = '{ref_col}'")

    # Detect if ALL adata indices end with "-1"
    adata_all_dash1 = adata.obs.index.str.endswith("-1").all()

    files = sorted(glob.glob(os.path.join(pred_dir, f"*_{dataset}.csv")))
    files = [f for f in files if not os.path.basename(f).startswith(("merged_", "performances_", "metrics_"))]
    if not files:
        sys.stderr.write(f"ERROR: no prediction files '*_{dataset}.csv' in {pred_dir}\n")
        sys.exit(1)

    print(f"[INFO] Found {len(files)} prediction file(s):")
    pred_cols = []
    for f in files:
        lbl = classifier_label_from_filename(f)
        print(f"  - {os.path.basename(f)}  ->  {lbl}")

        s = read_prediction_firstcol(f)

        # If adata indices all end with "-1" but prediction indices don't, append "-1"
        if adata_all_dash1 and not s.index.str.endswith("-1").all():
            s.index = s.index.astype(str) + "-1"

        s = s.reindex(adata.obs.index)
        adata.obs[lbl] = s
        pred_cols.append(lbl)

    print(f"[INFO] Running CIA metrics on columns: {pred_cols}")
    report_df = compute_classification_metrics(
        data=adata,
        classification_obs=pred_cols,
        ref_obs=ref_col,
        unassigned_label=unassigned_label
    )

    out_path = os.path.join(out_dir, f"performances_{dataset}.csv")
    report_df=report_df.round(2)
    report_df.to_csv(out_path, index=True)
    print(f"[OK] Wrote metrics to: {out_path}")
    print(f"[OK] Shape: {report_df.shape[0]} x {report_df.shape[1]}")

if __name__ == "__main__":
    main()
