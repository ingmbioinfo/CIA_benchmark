#!/usr/bin/env bash
# Run from the repo's scripts/ directory.
# Writes ONLY into existing ../rebuttal_datasets/{classifier} folders.
# Does NOT create directories and does NOT overwrite existing files.

set -euo pipefail
trap 'ret=$?; echo "[FAIL] line $LINENO: $BASH_COMMAND (exit $ret)"; exit $ret' ERR

# If run via Slurm, ensure paths are relative to the submit directory
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Base target: one level up from scripts/
BASE_DIR="$(cd .. && pwd)/rebuttal_datasets"
[[ -d "$BASE_DIR" ]] || { echo "ERROR: missing base dir: $BASE_DIR"; exit 2; }

have() { command -v "$1" >/dev/null 2>&1; }

fetch() { # fetch <url> <dest>
  local url="$1" dest="$2"
  if [[ -f "$dest" ]]; then
    echo "[SKIP] already downloaded: $dest"
    return 0
  fi
  if have curl; then
    curl -fL --retry 5 --retry-connrefused --continue-at - -o "$dest" "$url"
  elif have wget; then
    wget -c --tries=10 --timeout=30 -O "$dest" "$url"
  else
    echo "ERROR: need curl or wget." >&2; return 1
  fi
}

decompress_xz_keep() { # keep .xz; skip if output exists
  local path="$1"
  local out="${path%.xz}"
  if [[ -f "$out" ]]; then
    echo "[SKIP] already extracted: $out"
    return 0
  fi
  if have xz; then
    xz -dk "$path"
  elif have unxz; then
    unxz -k "$path"
  else
    echo "ERROR: need xz or unxz to decompress $path" >&2; return 1
  fi
  echo "[OK] extracted: $(basename "$out")"
}

# ---------------- URLs ----------------
readarray -t URLS <<'EOF'
https://zenodo.org/records/17241856/files/Celltypist_test_cancer.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/Celltypist_test_neuro.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/Celltypist_test_pbmc.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/Celltypist_train_cancer.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/Celltypist_train_neuro.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/Celltypist_train_pbmc.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/CIA_R_test_cancer.rds.xz?download=1
https://zenodo.org/records/17241856/files/CIA_R_test_neuro.rds.xz?download=1
https://zenodo.org/records/17241856/files/CIA_R_test_pbmc.rds.xz?download=1
https://zenodo.org/records/17241856/files/garnett_test_cancer.rds.xz?download=1
https://zenodo.org/records/17241856/files/garnett_test_neuro.rds.xz?download=1
https://zenodo.org/records/17241856/files/garnett_test_pbmc.rds.xz?download=1
https://zenodo.org/records/17241856/files/garnett_train_cancer.rds.xz?download=1
https://zenodo.org/records/17241856/files/garnett_train_neuro.rds.xz?download=1
https://zenodo.org/records/17241856/files/garnett_train_pbmc.rds.xz?download=1
https://zenodo.org/records/17241856/files/scANVI_test_cancer.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scANVI_test_neuro.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scANVI_test_pbmc.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scANVI_train_cancer.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scANVI_train_neuro.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scANVI_train_pbmc.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scBalance_test_cancer.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scBalance_test_neuro.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scBalance_test_pbmc.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scBalance_train_cancer.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scBalance_train_neuro.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/scBalance_train_pbmc.h5ad.xz?download=1
https://zenodo.org/records/17241856/files/SingleR_test_cancer.rds.xz?download=1
https://zenodo.org/records/17241856/files/SingleR_test_neuro.rds.xz?download=1
https://zenodo.org/records/17241856/files/SingleR_test_pbmc.rds.xz?download=1
https://zenodo.org/records/17241856/files/SingleR_train_cancer.rds.xz?download=1
https://zenodo.org/records/17241856/files/SingleR_train_neuro.rds.xz?download=1
https://zenodo.org/records/17241856/files/SingleR_train_pbmc.rds.xz?download=1
EOF

# ---------------- main ----------------
download_and_extract() {
  local url="$1"
  local fname="${url##*/}"; fname="${fname%%\?*}"   # strip ?download=1
  local tool="${fname%_test_*}"; tool="${tool%_train_*}"
  [[ "$tool" == "$fname" ]] && tool="${fname%%_*}"

  local dest_dir="$BASE_DIR/$tool"
  [[ -d "$dest_dir" ]] || { echo "ERROR: expected dir not found: $dest_dir"; exit 3; }

  local compressed="$dest_dir/$fname"
  local extracted="${compressed%.xz}"

  echo "[INFO] → $tool : $fname"
  fetch "$url" "$compressed"
  decompress_xz_keep "$compressed" && rm -f "$compressed"


  # Special rule: copy Celltypist test .h5ad to AUCell & CIA
  if [[ "$fname" =~ ^Celltypist_test_(cancer|neuro|pbmc)\.h5ad\.xz$ ]]; then
    local ds="${BASH_REMATCH[1]}"
    local aucell_dir="$BASE_DIR/AUCell"
    local cia_dir="$BASE_DIR/CIA"
    [[ -d "$aucell_dir" ]] || { echo "ERROR: expected dir not found: $aucell_dir"; exit 4; }
    [[ -d "$cia_dir"   ]] || { echo "ERROR: expected dir not found: $cia_dir"; exit 4; }

    local aucell_target="$aucell_dir/AUCell_test_${ds}.h5ad"
    local cia_target="$cia_dir/CIA_test_${ds}.h5ad"

    [[ -f "$aucell_target" ]] || { cp -n "$extracted" "$aucell_target" && echo "[OK] copied → $aucell_target"; }
    [[ -f "$cia_target"   ]] || { cp -n "$extracted" "$cia_target"   && echo "[OK] copied → $cia_target"; }
  fi
}

for url in "${URLS[@]}"; do
  [[ -z "$url" ]] && continue
  download_and_extract "$url"
done

echo "[DONE] Wrote into existing dirs under: $BASE_DIR"

