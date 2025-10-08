#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-$(pwd)/rebuttal_datasets}"
mkdir -p "$BASE_DIR"

# ---- all URLs ----
read -r -d '' URLS <<'EOF'
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

echo "[INFO] Base directory: $BASE_DIR"
echo "[INFO] Will keep both compressed (.xz) and extracted files."

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found." >&2; exit 1; }; }
need_cmd curl
need_cmd xz

download_and_extract() {
  local url="$1"
  local filename
  filename="$(basename "${url%%\?*}")"   # strip ?download=1
  local tool="${filename%%_*}"           # prefix before first underscore
  local dest_dir="$BASE_DIR/$tool"
  mkdir -p "$dest_dir"

  local compressed_path="$dest_dir/$filename"
  local extracted_path="${compressed_path%.xz}"

  echo "[INFO] -> $tool : $filename"
  if [[ ! -f "$compressed_path" ]]; then
    curl -fL --retry 5 --retry-connrefused --continue-at - -o "$compressed_path" "$url"
  else
    echo "[SKIP] already downloaded: $compressed_path"
  fi

  if [[ ! -f "$extracted_path" ]]; then
    # -d decompress, -k keep .xz, -T0 use all cores
    xz -T0 -dk "$compressed_path"
    echo "[OK] extracted: $(basename "$extracted_path")"
  else
    echo "[SKIP] already extracted: $extracted_path"
  fi

  # Special rule: Celltypist_test_{dataset}.h5ad -> also AUCell & CIA
  if [[ "$filename" =~ ^Celltypist_test_(cancer|neuro|pbmc)\.h5ad\.xz$ ]]; then
    local ds="${BASH_REMATCH[1]}"
    mkdir -p "$BASE_DIR/AUCell" "$BASE_DIR/CIA"

    local aucell_target="$BASE_DIR/AUCell/AUCell_test_${ds}.h5ad"
    local cia_target="$BASE_DIR/CIA/CIA_test_${ds}.h5ad"

    if [[ ! -f "$aucell_target" ]]; then
      cp -n "$extracted_path" "$aucell_target"
      echo "[OK] copied to: AUCell/AUCell_test_${ds}.h5ad"
    else
      echo "[SKIP] AUCell copy exists: $aucell_target"
    fi

    if [[ ! -f "$cia_target" ]]; then
      cp -n "$extracted_path" "$cia_target"
      echo "[OK] copied to: CIA/CIA_test_${ds}.h5ad"
    else
      echo "[SKIP] CIA copy exists: $cia_target"
    fi
  fi
}

# Iterate URLs
while IFS= read -r url; do
  [[ -z "$url" ]] && continue
  download_and_extract "$url"
done <<< "$URLS"

echo "[DONE] All files processed under: $BASE_DIR"

