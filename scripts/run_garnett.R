#!/usr/bin/env Rscript

# ------------------------ Load libraries ------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(monocle3)
  library(Seurat)
  library(garnett)
})

# ------------------------ CLI ------------------------
option_list <- list(
  make_option(c("--dataset"), type="character",
              help="Which dataset split to use", metavar="DATASET"),
  make_option(c("--cpus"), type="integer",
              help="Number of CPUs/threads to use", metavar="N"),
  make_option(c("--reduction"), type="character", default="UMAP",
              help="Reduction method to use [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list,
                           description="Train + predict with Garnett on selected dataset")
opt <- parse_args(opt_parser)

# ------------------------ Validate arguments ------------------------
if (is.null(opt$dataset) || !(opt$dataset %in% c("neuro","pbmc","cancer"))) {
  print_help(opt_parser)
  stop("Error: --dataset must be one of: neuro, pbmc, cancer", call.=FALSE)
}
if (is.null(opt$cpus)) {
  print_help(opt_parser)
  stop("Error: --cpus is required", call.=FALSE)
}

# ------------------------ Resolve relative paths ------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
repo_dir <- normalizePath(file.path(script_dir, ".."))

data_dir   <- file.path(repo_dir, "rebuttal_datasets", "garnett")
train_path <- file.path(data_dir, sprintf("garnett_train_%s.rds", opt$dataset))
test_path  <- file.path(data_dir, sprintf("garnett_test_%s.rds",  opt$dataset))
signature_path <- file.path(data_dir, sprintf("garnett_train_%s.txt",  opt$dataset))

# ------------------------ Check required files ------------------------
if (!file.exists(train_path) || !file.exists(test_path) || !file.exists(signature_path)) {
  message(sprintf("[ERROR] Expected files not found:\n  %s\n  %s\n  %s",
                  train_path, test_path, signature_path))
  quit(status = 2)
}

# ------------------------ Create output dirs ------------------------
runtime_dir <- file.path(repo_dir, "running_time")
pred_dir    <- file.path(repo_dir, "predictions")
dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

runtime_file <- file.path(runtime_dir, sprintf("garnett_%s.txt", opt$dataset))
pred_file    <- file.path(pred_dir,    sprintf("garnett_%s.csv", opt$dataset))

# ------------------------ Load training data ------------------------
cds <- readRDS(train_path)
colData(cds)$garnett_cluster <- clusters(cds, reduction_method=opt$reduction)

# ------------------------ Train classifier ------------------------
start_time <- Sys.time()
my_classifier <- train_cell_classifier(
  cds = cds,
  marker_file = signature_path,
  db = org.Hs.eg.db::org.Hs.eg.db,
  cds_gene_id_type = "SYMBOL",
  marker_file_gene_id_type = "SYMBOL",
  min_observations = 50,
  max_training_samples = 800,
  cores = opt$cpus
)
end_time <- Sys.time()

training_time <- end_time - start_time
rm(cds)

# ------------------------ Load test data & classify ------------------------
cds <- readRDS(test_path)
start_time <- Sys.time()
cds <- classify_cells(
  cds,
  my_classifier,
  db = org.Hs.eg.db::org.Hs.eg.db,
  cluster_extend = TRUE,
  cds_gene_id_type = "SYMBOL"
)
end_time <- Sys.time()
testing_time <- end_time - start_time

# ------------------------ Compute times in seconds ------------------------
training_sec <- as.numeric(training_time, units = "secs")
testing_sec  <- as.numeric(testing_time,  units = "secs")
total_sec    <- training_sec + testing_sec

# ------------------------ Export runtime to file ------------------------
if (!file.exists(runtime_file)) {
  write("TIME(s)\tCPUs", file = runtime_file)
}

write(sprintf("%.3f\t%d", total_sec, opt$cpus),
      file = runtime_file,
      append = TRUE)

print(sprintf("[OK] Recorded total runtime: %.3f s (training %.3f s + inference %.3f s)",
              total_sec, training_sec, testing_sec))
print(paste("[OK] Appended runtime to:", runtime_file))

# ------------------------ Export classification results ------------------------
write.csv(as.data.frame(colData(cds)['cell_type']),
          file = pred_file,
          row.names = TRUE)

print(paste("[OK] Wrote predictions to:", pred_file))

