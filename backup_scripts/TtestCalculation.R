# Calculating t-test results
rm(list = ls())
source("Scripts/CalculationTools.R")

# Function to process a single file with error handling
process_single_file <- function(file, output_dir, dataset_name) {
  tryCatch({
    # Generate paths
    meta_file <- paste0("Data/MetaData/metadata_", 
                       gsub(basename(file), pattern=".txt$", replacement=""), 
                       ".txt")
    
    # Process data
    countsT <- LoadCountsT(file)
    meta <- LoadMeta(meta_file)
    
    # Process data
    countsT <- countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
    meta <- meta[intersect(colnames(countsT), rownames(meta)), , drop = F]
    
    rownames(meta) <- colnames(countsT)
    colnames(meta) <- "conditions"
    
    results <- calcTtest(countsT, meta)
    
    # Create output directory if it doesn't exist
    output_path <- file.path(output_dir, dataset_name, "ttest")
    if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
    
    # Write results
    write.table(results, 
                file.path(output_path,
                         paste0(gsub(basename(file), pattern=".txt$", replacement=""), "_t.txt")),
                sep = "\t", row.names = T)
    
    list(status = "success", file = basename(file))
  }, error = function(e) {
    # Log error and return error status
    message(sprintf("Error processing %s: %s", basename(file), e$message))
    list(status = "error", file = basename(file), error = e$message)
  })
}

# Function to process a single dataset
processDataset <- function(input_dir, output_dir, dataset_name) {
  all_files <- list.files(file.path(input_dir), pattern = "*.txt", full.names = TRUE)
  
  cat(sprintf("\nProcessing %s dataset...\n", dataset_name))
  
  # Process files
  results <- lapply(all_files, function(file) {
    result <- process_single_file(file, output_dir, dataset_name)
    cat(sprintf("  %-20s ... %s\n", basename(file), result$status))
    result
  })
  
  # Print summary
  successful <- sum(sapply(results, function(x) x$status == "success"))
  cat(sprintf("Completed %d/%d files\n", successful, length(all_files)))
  
  invisible(NULL)
}

# Define dataset configurations
datasets <- list(
  list(input_dir = "Data/CountsTable/RDPNorm", 
       output_dir = "Results/PkgResults",
       name = "RDP"),
  list(input_dir = "Data/CountsTable/dada2Norm",
       output_dir = "Results/PkgResults",
       name = "dada2"),
  list(input_dir = "Data/CountsTable/WGSNorm",
       output_dir = "Results/PkgResults",
       name = "WGS")
)

# Process all datasets
invisible(lapply(datasets, function(ds) {
  processDataset(ds$input_dir, ds$output_dir, ds$name)
}))


