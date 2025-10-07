#!/usr/bin/env python3
import sys, os, glob
from math import floor, ceil

# ---------- percentiles (linear interpolation, NumPy-like) ----------
def pct(vals, p):
    if not vals:
        return float("nan")
    x = sorted(vals)
    n = len(x)
    if n == 1:
        return x[0]
    pos = (p / 100.0) * (n - 1)
    lo = int(floor(pos))
    hi = int(ceil(pos))
    if lo == hi:
        return x[lo]
    frac = pos - lo
    return x[lo] + frac * (x[hi] - x[lo])

def float_str(x):
    return f"{x:.15g}"

def to_label(filename):
    """
    Map file basename to classifier label.

    Expects filenames like:
      <LABEL>_<dataset>.txt
    where <LABEL> may contain underscores (e.g., CIA_R, CIA_python).

    We remove the final "_<dataset>" chunk by splitting from the RIGHT.
    """
    base = os.path.basename(filename)
    name, _ext = os.path.splitext(base)          # e.g., CIA_R_pbmc
    if "_" in name:
        label = name.rsplit("_", 1)[0]           # e.g., CIA_R
    else:
        label = name

    # Optional normalization (only for exact matches)
    low = label.lower()
    if low == "cia":           # plain "CIA" means Python implementation in your setup
        return "CIA_Python"
    if low == "cia_python":
        return "CIA_Python"
    if low == "cia_r":
        return "CIA_R"

    return label

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 make_summary.py <dataset>\n  <dataset> in {neuro, pbmc, cancer}", file=sysstderr)
        sys.exit(2)

    dataset = sys.argv[1].strip().lower()

    # Resolve ../running_time relative to THIS script file (not CWD)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_dir = os.path.normpath(os.path.join(script_dir, "..", "running_time"))

    if not os.path.isdir(run_dir):
        raise SystemExit(f"running_time directory not found at: {run_dir}")

    pattern = os.path.join(run_dir, f"*_{dataset}.txt")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No files matching '*_{dataset}.txt' in {run_dir}")

    # (label, cpu) -> list of times
    buckets = {}

    for f in files:
        label = to_label(f)
        with open(f) as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                # Skip header lines like "TIME(s) CPUs"
                head = s.lower().replace("\t", " ")
                if head.startswith("time"):
                    continue
                parts = s.split()
                if len(parts) < 2:
                    continue
                try:
                    t = float(parts[0])
                    cpu = int(parts[1])
                except ValueError:
                    continue
                buckets.setdefault((label, cpu), []).append(t)

    # Preferred label order; unseen labels fall to the end alphabetically
    pref_order = ["AUCell", "Celltypist", "CIA_Python", "CIA_R", "Garnett", "SingleR", "scANVI"]
    labels_present = sorted(
        {lab for (lab, _) in buckets.keys()},
        key=lambda L: (pref_order.index(L) if L in pref_order else 999, L)
    )

    out_path = os.path.join(run_dir, f"{dataset}_summary.txt")
    with open(out_path, "w") as out:
        out.write("CPUs Median IQR IQR_lower IQR_upper classifier\n")
        for lab in labels_present:
            cpus = sorted({cpu for (L, cpu) in buckets.keys() if L == lab})
            for cpu in cpus:
                times = buckets[(lab, cpu)]
                q1 = pct(times, 25)
                med = pct(times, 50)
                q3 = pct(times, 75)
                iqr = q3 - q1
                out.write(
                    f"{cpu} {float_str(med)} {float_str(iqr)} "
                    f"{float_str(q1)} {float_str(q3)} {lab}\n"
                )

    print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()

