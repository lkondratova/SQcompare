#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(digest)
  library(readr)
  library(dplyr)
  library(tools)
})

# Argument Parser
option_list <- list(
  make_option(c("--collapse_ISM"), type = "logical", default = FALSE,
              help = "Collapse ISM isoforms into existant FSM (if available) [default: FALSE]"),
  make_option(c("--input_files"), type = "character", help = "TSV file with paths to SQANTI3 outputs")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate and Read Input
if (is.null(opt$input_files)) {
  stop("Please provide the --input_files argument (a TSV file)")
}

input_df <- read_tsv(opt$input_files, col_names = FALSE, comment = "#")

if (ncol(input_df) < 2) {
  stop("Input TSV must have at least two columns: classification.txt and junctions.txt")
}

colnames(input_df)[1:2] <- c("classification_path", "junctions_path")
if (ncol(input_df) >= 3) {
  colnames(input_df)[3] <- "expression_path"
}

# Initialize Storage
samples_info <- list()
n_samples <- nrow(input_df)

# Process Each Sample
for (i in seq_len(n_samples)) {
  class_file <- input_df$classification_path[i]
  junction_file <- input_df$junctions_path[i]
  expr_file <- if ("expression_path" %in% colnames(input_df)) input_df$expression_path[i] else NA
  
  # Check if files exist
  if (!file.exists(class_file)) stop(paste("Classification file not found:", class_file))
  if (!file.exists(junction_file)) stop(paste("Junctions file not found:", junction_file))
  if (!is.na(expr_file) && !file.exists(expr_file)) stop(paste("Expression file not found:", expr_file))
  
  # Get sample name from classification filename
  sample_name <- file_path_sans_ext(basename(class_file))
  sample_name <- sub("_classification$", "", sample_name)

  # Generate hash from combined paths
  hash_input <- paste(class_file, junction_file, expr_file, sep = "|")
  hash_id <- digest(hash_input, algo = "md5")

  # Read files
  classification_data <- read_tsv(class_file, show_col_types = FALSE)
  junction_data <- read_tsv(junction_file, show_col_types = FALSE)
  expression_data <- if (!is.na(expr_file)) read_tsv(expr_file, show_col_types = FALSE) else NULL
  
  # Store
  samples_info[[sample_name]] <- list(
    hash = hash_id,
    classification = classification_data,
    junctions = junction_data,
    expression = expression_data
  )
}

# Summary Output
cat(sprintf("Loaded %d samples.\n", n_samples))
cat("Sample names:\n")
cat(paste0(" - ", names(samples_info), collapse = "\n"), "\n")
cat("\n--collapse_ISM =", opt$collapse_ISM, "\n")

# Save Intermediate Object
saveRDS(samples_info, file = "sqanti3_samples_parsed.rds")