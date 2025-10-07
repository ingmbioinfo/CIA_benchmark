#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_cia_R.R <dataset> <cpus>
# Example: Rscript scripts/run_cia_R.R pbmc 32

suppressPackageStartupMessages({
  library(CIA)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/run_cia_R.R <dataset> <cpus>")
}
dataset <- args[[1]]
cpus    <- as.integer(args[[2]])

# Thread env (match -c)
Sys.setenv(
  OMP_NUM_THREADS       = cpus,
  OPENBLAS_NUM_THREADS  = cpus,
  MKL_NUM_THREADS       = cpus,
  NUMEXPR_NUM_THREADS   = cpus,
  VECLIB_MAXIMUM_THREADS= cpus,
  UCX_TLS               = "tcp,sm,self"
)

# Resolve repo root from this script's path
args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_full[grep("^--file=", args_full)])
script_dir  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
repo_dir    <- normalizePath(file.path(script_dir, ".."))

# Paths (relative to repo root)
ds_dir      <- file.path(repo_dir, "rebuttal_datasets", "CIA_R")
runtime_dir <- file.path(repo_dir, "running_time"); dir.create(runtime_dir, showWarnings = FALSE, recursive = TRUE)
pred_dir    <- file.path(repo_dir, "predictions");  dir.create(pred_dir,  showWarnings = FALSE, recursive = TRUE)

test_rds    <- file.path(ds_dir, sprintf("CIA_R_test_%s.rds",  dataset))
train_gmt   <- file.path(ds_dir, sprintf("CIA_R_train_%s.gmt", dataset))
pred_out    <- file.path(pred_dir, sprintf("CIA_R_%s.csv", dataset))
summary_fn  <- file.path(runtime_dir, sprintf("CIA_R_%s.txt", dataset))

if (!file.exists(test_rds))  stop("Missing test RDS: ", test_rds)
if (!file.exists(train_gmt)) stop("Missing training GMT: ", train_gmt)

# Guard against HTML "blob" instead of raw GMT
head_bytes <- tolower(paste(readLines(train_gmt, n = 3, warn = FALSE), collapse = " "))
if (grepl("<!doctype html|<html", head_bytes)) {
  stop(train_gmt, " looks like an HTML page. Download the RAW .gmt file.")
}

# Load data
so  <- readRDS(test_rds)
gmt <- CIA::load_signatures(train_gmt)

# Classify (only dataset, cpus — defaults for the rest)
t0 <- proc.time()[["elapsed"]]
so <- CIA::CIA_classify(
  data              = so,
  signatures_input  = gmt,
  n_cpus            = cpus,
  similarity_threshold = 0
)
elapsed <- proc.time()[["elapsed"]] - t0

# Save predictions
pred <- data.frame(
  cell_id    = colnames(so),
  prediction = as.character(so@meta.data[["CIA_prediction"]]),
  stringsAsFactors = FALSE
)
write.table(pred, pred_out, row.names = FALSE, quote=FALSE, sep=',' )
cat(sprintf("[OK] Saved predictions -> %s\n", pred_out))

# Append runtime (single per-dataset file)
needs_header <- !file.exists(summary_fn) || (file.info(summary_fn)$size %||% 0) == 0
con <- file(summary_fn, open = if (file.exists(summary_fn)) "a" else "w")
on.exit(close(con), add = TRUE)
if (needs_header) writeLines("TIME(s)\tCPUs", con)
writeLines(sprintf("%f\t%d", elapsed, cpus), con)
cat(sprintf("[OK] Appended runtime -> %s  (elapsed: %.3fs, cpus=%d)\n", summary_fn, elapsed, cpus))

`%||%` <- function(a, b) if (is.null(a)) b else a

