# CIA Benchmark Repository

A reproducible benchmark comparing single‑cell classifiers on **PBMC**, **neuro**, and **cancer** datasets, with scripts, conda environments, and Slurm wrappers to run locally or on HPC.

> **Execution convention:** all commands below assume your current working directory is `CIA_benchmark_repository/scripts/` and paths are resolved relative to that.

---

## Table of contents

- [Repository layout](#repository-layout)
- [Top‑level folders](#top-level-folders)
- [Run convention](#run-convention)
- [Datasets (`rebuttal_datasets/`)](#datasets-rebuttal_datasets)
  - [Naming scheme](#naming-scheme)
- [Scripts (`scripts/`)](#scripts-scripts)
  - [Runner example: `run_cia.py`](#runner-example-run_cia.py)
  - [Slurm wrappers: `run_{classifier}.sbatcher`](#slurm-wrappers-run_classifiersbatcher)
  - [Aggregators](#aggregators)
- [Conda environments (`envs/`)](#conda-environments-envs)
- [Predictions (`predictions/`)](#predictions-predictions)
- [Running times (`running_time/`)](#running-times-running_time)
- [Performances (`performances/`)](#performances-performances)
- [Reproduce the benchmark](#reproduce-the-benchmark)
- [Notes on reproducibility](#notes-on-reproducibility)

---

## Repository layout

```
CIA_benchmark_repository/
├── rebuttal_datasets/
│   ├── AUCell/
│   ├── Celltypist/
│   ├── CIA/
│   ├── CIA_R/
│   ├── garnett/
│   ├── scANVI/
│   ├── scBalance/
│   └── SingleR/
├── scripts/
├── envs/
├── predictions/
├── running_time/
└── performances/
```

## Top‑level folders

- **`rebuttal_datasets/`** — benchmark inputs per classifier; train/reference and test sets for **pbmc**, **neuro**, **cancer**.
- **`scripts/`** — runnable bash/Python/R pipelines. **Assume `PWD = scripts/`.**
- **`envs/`** — conda environment YAMLs for reproducible runs.
- **`predictions/`** — classifier outputs (predicted labels).
- **`running_time/`** — per‑run timing logs and per‑dataset summaries.
- **`performances/`** — aggregated performance metrics tables.

## Run convention

All examples assume:

```bash
cd CIA_benchmark_repository/scripts
bash <script>.sh   # or: python <script>.py / Rscript <script>.R
```

---

## Datasets (`rebuttal_datasets/`)

After downloading (see **Reproduce the benchmark**), you should have:

```
rebuttal_datasets/
├── AUCell/
├── Celltypist/
├── CIA/
├── CIA_R/
├── garnett/
├── scANVI/
├── scBalance/
└── SingleR/
```

Each classifier folder contains both train (or reference, e.g. `.gmt` where applicable) and test datasets for `pbmc`, `neuro`, and `cancer` benchmarks.

### Naming scheme

For every `{classifier} ∈ {AUCell, Celltypist, CIA, CIA_R, garnett, scANVI, scBalance, SingleR}` and
every `{dataset} ∈ {pbmc, neuro, cancer}`, files follow:

```
{classifier}_train_{dataset}.{ext}
{classifier}_test_{dataset}.{ext}
```

`{ext}` depends on the classifier (`.h5ad`, `.rds`, `.gmt`, …).

---

## Scripts (`scripts/`)

- **`download_data.sh`** — downloads the ready‑to‑use datasets into the proper classifier directories under `../rebuttal_datasets/{Classifier}/`, following the naming above.
- **`run_{classifier}.py|R`** — classifier runners (AUCell, Celltypist, CIA, CIA_R, Garnett, scANVI, scBalance, SingleR).
- **`run_{classifier}.sbatcher`** — Slurm submission wrappers that mirror the runners.

### Runner example: `run_cia.py`

**Usage**
```bash
python run_cia.py <dataset> <cpus>
# <dataset>: pbmc | neuro | cancer
```

**Behaviour (example: `python run_cia.py pbmc 32`)**

- **Inputs**: `../rebuttal_datasets/CIA/CIA_train_pbmc.gmt` (reference) and `../rebuttal_datasets/CIA/CIA_test_pbmc.h5ad` (test).
- **Predictions**: `../predictions/CIA_pbmc.csv`.
- **Runtime log**: appends to `../running_time/CIA_pbmc.txt` with header `TIME(s)  CPUs` (here `CPUs=32`).

Other runners follow the same contract:

- **Predictions**: `../predictions/{CLF}_{dataset}.csv`
- **Runtime log**: `../running_time/{CLF}_{dataset}.txt` (tab‑separated; header: `TIME(s)  CPUs`)
- **Optional per‑dataset roll‑up**: `../running_time/{dataset}_summary.txt` created by aggregators.

### Slurm wrappers: `run_{classifier}.sbatcher`

**Usage**
```bash
sbatch -c <cpus> run_cia.sbatcher <dataset> <cpus>
# example:
sbatch -c 32 run_cia.sbatcher pbmc 32
```

**Behaviour**

- **Args**: `<dataset>` ∈ `{pbmc, neuro, cancer}`; `<cpus>` is the thread count. If omitted, defaults to `$SLURM_CPUS_PER_TASK` (or `1`).
- **Thread binding**: exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` to your requested value so libraries honour `-c`.
- **Working dir**: `cd "$SLURM_SUBMIT_DIR"` to run next to the runner.
- **Launch**: `python|Rscript run_{classifier}.* <dataset> <cpus>`.
- **Logging**: prints dataset/CPUs, job id, and node to STDOUT.

### Aggregators

#### `make_summary.py` (runtime roll‑ups)

**Usage**
```bash
python3 make_summary.py <dataset>   # pbmc | neuro | cancer
```

**What it does**

- Scans `../running_time/*_<dataset>.txt` and groups rows by `(classifier_label, CPUs)`.
- Normalises labels (e.g. `CIA_R_* → CIA_R`, `CIA_* → CIA_Python`).
- Expects per‑file lines like `TIME(s)  CPUs` header then `<seconds> <cpus>` rows.
- Computes **Q1**, **median**, **Q3**, **IQR** (NumPy‑like linear interpolation).
- Writes `../running_time/<dataset>_summary.txt` with columns:
  `CPUs  Median  IQR  IQR_lower  IQR_upper  classifier`.

#### `performances.py` (classification metrics)

**Usage**
```bash
python3 performances.py <dataset> [unassigned_label]
# <dataset>: pbmc | neuro | cancer
# optional unassigned_label defaults to ""
```

**What it does**

- Loads test AnnData from `../rebuttal_datasets/CIA/CIA_test_<dataset>.h5ad`.
- Detects the ground‑truth column in `adata.obs` (`"Cell type"` or `"CELLTYPE"`).
- Reads all `../predictions/*_<dataset>.csv` (skips `merged_`, `performances_`, `metrics_` prefixes), keeping only the **first** data column from each file.
- Normalises classifier labels from filenames: `CIA_R_* → CIA_R`; `CIA_* → CIA_Python`; others use the stem before the first underscore.
- **Index harmonisation**: if all `adata.obs.index` end with `-1` but a prediction index does not, it appends `-1` to the prediction indices before reindexing.
- Adds one column per classifier to `adata.obs`, then calls `cia.report.compute_classification_metrics` with `classification_obs=<those columns>`, `ref_obs=<detected>`, and optional `unassigned_label`.
- Writes the aggregated table to `../performances/performances_<dataset>.csv`.

---

## Conda environments (`envs/`)

Create each environment once, then activate before running the corresponding scripts.

```bash
conda env create -f envs/<ENV>.yml
conda activate <env-name>   # the name is defined inside the YAML
```

**Available environments and what they run**

- `AUCell.yml` — AUCell (Python)
  - Local: `python run_aucell.py <dataset> <cpus>`
  - Slurm: `sbatch -c <cpus> run_aucell.sbatcher <dataset> <cpus>`

- `CIA_and_Celltypist.yml` — shared env for CIA (Python), Celltypist, and the aggregators
  - CIA (Python): `python run_cia.py <dataset> <cpus>` / `sbatch -c <cpus> run_cia.sbatcher <dataset> <cpus>`
  - Celltypist: `python run_celltypist.py <dataset> <cpus>` / `sbatch -c <cpus> run_celltypist.sbatcher <dataset> <cpus>`
  - Aggregators: `python make_summary.py <dataset>` and `python performances.py <dataset>`

- `CIA_R.yml` — CIA (R implementation)
  - Local: `Rscript run_cia_R.R <dataset> <cpus>`
  - Slurm: `sbatch -c <cpus> run_cia_R.sbatcher <dataset> <cpus>`

- `garnett.yml` — Garnett (R)
  - Local: `Rscript run_garnett.R <dataset> <reduction> <cpus>`
  - Slurm: `sbatch -c <cpus> run_garnett.sbatcher <dataset> <reduction> <cpus>`
  - **Note**: `reduction=UMAP` was used for **neuro** and **cancer** datasets, and `reduction=t-SNE` for **PBMC**.

- `scANVI.yml` — scANVI (scvi‑tools, Python)
  - Local: `python run_scanvi.py <dataset> <cpus>`
  - Slurm: `sbatch -c <cpus> run_scanvi.sbatcher <dataset> <cpus>`

- `scBalance.yml` — scBalance (Python)
  - Local: `python run_scbalance.py <dataset> <cpus>`
  - Slurm: `sbatch -c <cpus> run_scbalance.sbatcher <dataset> <cpus>`

- `SingleR.yml` — SingleR (R)
  - Local: `Rscript run_singler.R <dataset> <cpus>`
  - Slurm: `sbatch -c <cpus> run_singler.sbatcher <dataset> <cpus>`

---

## Predictions (`predictions/`)

Classifier outputs live here. One CSV per classifier × dataset.

**Naming**

```
{Classifier}_{dataset}.csv
# Classifier ∈ {AUCell, Celltypist, CIA, CIA_R, garnett, scANVI, scBalance, SingleR}
# dataset   ∈ {pbmc, neuro, cancer}
```

Examples: `AUCell_pbmc.csv`, `Celltypist_neuro.csv`, `CIA_cancer.csv`, `CIA_R_pbmc.csv`, `garnett_neuro.csv`, `scANVI_cancer.csv`, `scBalance_pbmc.csv`, `SingleR_neuro.csv`.

**CSV schema (strict, minimal)**

- **Index (column 0)**: cell barcodes (strings). Performances script reindexes to the test AnnData.
- **First data column**: predicted label (any column name). Additional columns, if present, are ignored by `performances.py`.

---

## Running times (`running_time/`)

Per‑run wall‑clock timings, one file per classifier × dataset, plus per‑dataset roll‑ups from `make_summary.py`.

**Per‑run logs**

```
{Classifier}_{dataset}.txt   # tab‑separated
```

Each file has a single header and one row per run:

```
TIME(s)  CPUs
<seconds>  <cpus>
<seconds>  <cpus>
...
```

**Example**

```
TIME(s)  CPUs
12.84    32
12.55    32
13.01    32
```

**Roll‑ups**

```
{dataset}_summary.txt
# Columns: CPUs  Median  IQR  IQR_lower  IQR_upper  classifier
```

Example row:
```
32   12.64  0.46  12.55  13.01  CIA_Python
```

---

## Performances (`performances/`)

Aggregated classification metrics produced by `scripts/performances.py`:

```
performances_cancer.csv  performances_neuro.csv  performances_pbmc.csv
```

Each CSV has one row per classifier present in `../predictions/*_<dataset>.csv` and columns:

- **SE — Sensitivity (Recall)**: `TP / (TP + FN)`
- **SP — Specificity**: `TN / (TN + FP)`
- **PR — Precision**: `TP / (TP + FP)`
- **ACC — Accuracy**: `(TP + TN) / (TP + TN + FP + FN)`
- **F1 — F1‑score**: `2 × (Precision × Recall) / (Precision + Recall)`

**Example (excerpt)**

```csv
,SE,SP,PR,ACC,F1
AUCell,0.5355,0.9821,0.5355,0.9656,0.5355
CIA_Python,0.7754,0.9914,0.7754,0.9834,0.7754
Celltypist,0.7857,0.9918,0.7857,0.9841,0.7857
SingleR,0.7298,0.9903,0.7435,0.9807,0.7366
garnett,0.0965,0.9957,0.4643,0.9624,0.1598
scANVI,0.7776,0.9914,0.7776,0.9835,0.7776
scBalance,0.7774,0.9914,0.7774,0.9835,0.7774
```

---

## Reproduce the benchmark

1) **Download data**

```bash
cd scripts
bash download_data.sh
```
This populates `../rebuttal_datasets/` with train/test sets for all classifiers and datasets.

2) **Install conda environments**

```bash
conda env create -f envs/<ENV>.yml
```

3) **Run classifiers**

For each classifier, activate its conda env. For each dataset (`pbmc`, `neuro`, `cancer`), launch **10 runs** with the corresponding `sbatcher` at each CPU setting:

```
ncpus ∈ {1, 2, 8, 16, 32, 64}
```

For **scANVI**, run **3 times** at `ncpus ∈ {16, 32, 64}`.

**Example (CIA, PBMC, 32 CPUs; repeat 10×):**
```bash
sbatch -c 32 run_cia.sbatcher pbmc 32
```

4) **Aggregate**

Activate the `CIA_and_Celltypist` env, then:

- **Roll up runtimes:**
  ```bash
  python make_summary.py pbmc
  python make_summary.py neuro
  python make_summary.py cancer
  ```

- **Compute performances:**
  ```bash
  python performances.py pbmc
  python performances.py neuro
  python performances.py cancer
  ```

Outputs will appear in `../running_time/*_summary.txt` and `../performances/performances_<dataset>.csv`.

---

## Notes on reproducibility

- **Conda + Slurm required for apples‑to‑apples timings.** Reported runtimes depend on conda environments and the Slurm scheduler on our HPC. Absolute times will vary across machines.
- **Direct launch is fine if Slurm is unavailable.** Run any classifier via its `run_{classifier}.*` script once the env is active. Predictions and metrics will match; runtime logs/summaries will differ.
- **Memory for cancer datasets.** Some cancer runs needed large memory (e.g. `--mem=400G`). Adjust your `sbatcher` headers accordingly.
