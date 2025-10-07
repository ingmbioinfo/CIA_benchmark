#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_singler.R <dataset> <cpus>
# Example: Rscript scripts/run_singler.R pbmc 32

suppressPackageStartupMessages({
  library(SingleR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript scripts/run_singler.R <dataset> <cpus>")
dataset <- args[[1]]
cpus    <- as.integer(args[[2]])

# Thread env to match -c
Sys.setenv(
  OMP_NUM_THREADS        = cpus,
  OPENBLAS_NUM_THREADS   = cpus,
  MKL_NUM_THREADS        = cpus,
  NUMEXPR_NUM_THREADS    = cpus,
  VECLIB_MAXIMUM_THREADS = cpus,
  UCX_TLS                = "tcp,sm,self"
)

# Resolve repo root (script lives in repo/scripts)
args_full   <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_full[grep("^--file=", args_full)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_dir    <- normalizePath(file.path(script_dir, ".."))

# Paths
ds_dir      <- file.path(repo_dir, "rebuttal_datasets", "SingleR")
runtime_dir <- file.path(repo_dir, "running_time"); dir.create(runtime_dir, showWarnings = FALSE, recursive = TRUE)
pred_dir    <- file.path(repo_dir, "predictions");  dir.create(pred_dir,  showWarnings = FALSE, recursive = TRUE)

train_rds   <- file.path(ds_dir, sprintf("SingleR_train_%s.rds", dataset))
test_rds    <- file.path(ds_dir, sprintf("SingleR_test_%s.rds",  dataset))
pred_out    <- file.path(pred_dir, sprintf("SingleR_%s.csv", dataset))
summary_fn  <- file.path(runtime_dir, sprintf("SingleR_%s.txt", dataset))

if (!file.exists(train_rds)) stop("Missing training RDS: ", train_rds)
if (!file.exists(test_rds))  stop("Missing test RDS: ", test_rds)

# Load data
ref  <- readRDS(train_rds)  # expect SingleCellExperiment / SummarizedExperiment
test <- readRDS(test_rds)

# ---- Robust detection of 'cell type' column in training ----------------------
canon <- function(x) gsub("[^[:alnum:]]", "", tolower(x))  # lower + remove spaces/_/punct
cn_ref <- colnames(colData(ref))
hit <- which(canon(cn_ref) == "celltype")
if (length(hit) == 0) {
  stop(
    "Training object must have a cell-type column in colData(ref). ",
    "Accepted variants include: 'celltype', 'CELLTYPE', 'Cell type'. ",
    "Found: ", paste(cn_ref, collapse = ", ")
  )
}
label_col <- cn_ref[hit[1]]  # use the first match

# Run SingleR
set.seed(1)
bp <- BiocParallel::MulticoreParam(workers = cpus)
t0 <- proc.time()[["elapsed"]]
pred <- SingleR(
  test        = test,
  ref         = ref,
  labels      = colData(ref)[[label_col]],
  de.method   = "wilcox",
  num.threads = cpus,
  BPPARAM     = bp
)
elapsed <- proc.time()[["elapsed"]] - t0

# Use pruned.labels when available; else labels
pred_vec <- if ("pruned.labels" %in% colnames(pred)) pred[, "pruned.labels"] else pred[, "labels"]

# Save predictions
out <- data.frame(
  cell_id    = colnames(test),
  prediction = as.character(pred_vec),
  stringsAsFactors = FALSE
)
write.table(out, pred_out, row.names = FALSE, quote = FALSE, sep = ',')
cat(sprintf("[OK] Saved predictions -> %s (label_col: %s)\n", pred_out, label_col))

# Append runtime (single per-dataset file)
needs_header <- !file.exists(summary_fn) || isTRUE(file.info(summary_fn)$size == 0)
con <- file(summary_fn, open = if (file.exists(summary_fn)) "a" else "w"); on.exit(close(con), add = TRUE)
if (needs_header) writeLines("TIME(s)\tCPUs", con)
writeLines(sprintf("%.6f\t%d", elapsed, cpus), con)
cat(sprintf("[OK] Appended runtime -> %s  (elapsed: %.3fs, cpus=%d)\n", summary_fn, elapsed, cpus))
