#shuffle counts in every sample 100x and test all pval for resampled whole taxon of counts tables
rm(list = ls())
set.seed(9527)
source("Scripts/CalculationTools.R")

# Function to process a single file
process_single_file <- function(file, dataset_name) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  
  # Generate file paths
  raw_file <- paste0("Data/CountsTable/", dataset_name, "Raw/", name, ".txt")
  meta_file <- paste0("Data/MetaData/metadata_", name, ".txt")
  
  # Load data using CalculationTools functions
  RawcountsT <- LoadCountsT(raw_file)
  meta <- LoadMeta(meta_file)

  # Filter out low counts taxa
  if (all(which(rowMeans(RawcountsT)<2) == 0)) {
    RawcountsT <- RawcountsT
  } else {
    RawcountsT <- RawcountsT[-c(which(rowMeans(RawcountsT)<2)), , drop = F]
  }

  # Align data
  RawcountsT <- RawcountsT[, intersect(colnames(RawcountsT), rownames(meta)), drop = F]
  meta <- meta[intersect(colnames(RawcountsT), rownames(meta)), , drop = F]

  rownames(meta) <- colnames(RawcountsT)
  colnames(meta) <- "conditions"

  # Resample
  RawcountsT <- resampleWholeTaxonRNBINOM(RawcountsT, meta, 1)

  # Initialize p-value matrices
  ttestpvals <- c()
  wilcoxpvals <- c()
  deseq2pvals <- c()
  edgerpvals <- c()
  
  # Run 100 iterations
  for (i in 1:100) {
    shuffledCountsT <- apply(RawcountsT, 2, sample)
    NormCountsT <- normFun(shuffledCountsT)

    ttestresults <- calcTtest(NormCountsT, meta)
    ttestpvals <- cbind(ttestpvals, ttestresults$pval)
    wilcoxresults <- calcWilcox(NormCountsT, meta)
    wilcoxpvals <- cbind(wilcoxpvals, wilcoxresults$pval)
    deseq2results <- calcDESeq2(shuffledCountsT, meta)
    deseq2pvals <- cbind(deseq2pvals, deseq2results$pval)
    edgerresults <- calcEdgeR(shuffledCountsT, meta)
    edgerpvals <- cbind(edgerpvals, edgerresults$pval)
  }

  # Set row names
  rownames(ttestpvals) <- rownames(RawcountsT)
  rownames(wilcoxpvals) <- rownames(RawcountsT)
  rownames(deseq2pvals) <- rownames(RawcountsT)
  rownames(edgerpvals) <- rownames(RawcountsT)

  # Handle NA values
  ttestpvals[is.na(ttestpvals)] <- 1
  wilcoxpvals[is.na(wilcoxpvals)] <- 1

  # Save results
  save_dir <- paste0("Results/NBResampleWholeTaxonShuffleInSampleDump/", dataset_name, "/")
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  write.table(ttestpvals, paste0(save_dir, name, "_t.txt"), sep = "\t", row.names = T)
  write.table(wilcoxpvals, paste0(save_dir, name, "_wilcox.txt"), sep = "\t", row.names = T)
  write.table(deseq2pvals, paste0(save_dir, name, "_deseq2.txt"), sep = "\t", row.names = T)
  write.table(edgerpvals, paste0(save_dir, name, "_edgeR.txt"), sep = "\t", row.names = T)
  
  return(list(status = "success", file = name))
}

# Function to process a dataset
processDataset <- function(dataset_name) {
  all_files <- list.files(paste0("Data/CountsTable/", dataset_name, "Raw"), 
                         pattern = "*.txt", full.names = TRUE)
  
  cat(sprintf("\nProcessing %s dataset...\n", dataset_name))
  
  # Process files
  results <- lapply(all_files, function(file) {
    result <- process_single_file(file, dataset_name)
    cat(sprintf("  %-20s ... %s\n", basename(file), result$status))
    result
  })
  
  # Print summary
  successful <- sum(sapply(results, function(x) x$status == "success"))
  cat(sprintf("Completed %d/%d files\n", successful, length(all_files)))
  
  invisible(NULL)
}

# Dataset configuration
datasets <- list(
  list(name = "RDP"),
  list(name = "dada2"),
  list(name = "WGS"),
)

# Process all datasets
invisible(lapply(datasets, function(ds) processDataset(ds$name)))
