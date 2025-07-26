#resampled whole taxon shuffled significant pvalue visualization
rm(list = ls())
source("Scripts/PlotTools.R")
library(ggExtra)
library(patchwork)
library(ggplot2)

makePlot <- function(fractT, classifier, name, shuffleType) {
  df_long <- stack(fractT)
  
  p <- ggplot(df_long, aes(x = ind, y = values, fill = ind)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("red", "tan2", "purple", "cornflowerblue")) +
    theme_classic() +
    labs(title = paste0("(", classifier, "-", name, ")\n", "NB-Resample Whole Taxon-\n", shuffleType),
         x = "DAA methods",
         y = "Fraction of significant results") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylim(0, 0.8)
  
  return(p)
}

calcFraction <- function(pvalT) {
  fractions <- apply(pvalT, 2,
                     function(column){
                       sum(column < 0.05)/length(column)
                     })
  return(fractions)
}

# Function to process a single file for a specific shuffle type
process_single_file <- function(file, dataset_name, shuffle_type) {
  name <- gsub(basename(file), pattern=".txt$", replacement="")
  
  # Read p-value files
  ttestpvals <- read.table(paste0("Results/", shuffle_type$dir, "/", dataset_name, "/", name, "_t.txt"), 
                           header = T, row.names = 1, sep = "\t")
  wilcoxpvals <- read.table(paste0("Results/", shuffle_type$dir, "/", dataset_name, "/", name, "_wilcox.txt"), 
                            header = T, row.names = 1, sep = "\t")
  deseq2pvals <- read.table(paste0("Results/", shuffle_type$dir, "/", dataset_name, "/", name, "_deseq2.txt"), 
                            header = T, row.names = 1, sep = "\t")
  edgerpvals <- read.table(paste0("Results/", shuffle_type$dir, "/", dataset_name, "/", name, "_edgeR.txt"), 
                           header = T, row.names = 1, sep = "\t")
  
  # Calculate significant fractions
  sigT <- calcFraction(ttestpvals)
  sigW <- calcFraction(wilcoxpvals)
  sigD <- calcFraction(deseq2pvals)
  sigE <- calcFraction(edgerpvals)
  
  fractionT <- cbind(sigT, sigW, sigD, sigE)
  colnames(fractionT) <- c("t-test", "Wilcoxon", "DESeq2", "edgeR")
  fractionT <- data.frame(fractionT, check.names = F)
  
  # Create plot
  p <- makePlot(fractionT, dataset_name, name, shuffle_type$name)
  return(p)
}

# Function to process a dataset type
processDatasetType <- function(dataset_name) {
  # Define shuffle types and their directories
  shuffle_types <- list(
    list(name = "Shuffle sample tags", dir = "NBResampleWholeTaxonShuffleSampleTagDump"),
    list(name = "Shuffle counts in sample", dir = "NBResampleWholeTaxonShuffleInSampleDump"),
    list(name = "Shuffle counts in taxon", dir = "NBResampleWholeTaxonShuffleInTaxonDump"),
    list(name = "Shuffle counts table", dir = "NBResampleWholeTaxonShuffleCountsTDump")
  )
  
  cat(sprintf("\nProcessing %s dataset...\n", dataset_name))
  
  # Get all files for this dataset
  all_files <- list.files(paste0("Data/CountsTable/", dataset_name, "Raw"), 
                          pattern = "*.txt", full.names = TRUE)
  
  # Process all shuffle types for this dataset
  all_shuffle_plots <- list()
  for (shuffle_idx in seq_along(shuffle_types)) {
    shuffle_type <- shuffle_types[[shuffle_idx]]
    
    # Process all files for this shuffle type
    file_plots <- list()
    for (file_idx in seq_along(all_files)) {
      file <- all_files[[file_idx]]
      p <- process_single_file(file, dataset_name, shuffle_type)
      file_plots[[file_idx]] <- p
      cat(sprintf("  %-20s ... processed\n", basename(file)))
    }
    
    all_shuffle_plots[[shuffle_idx]] <- file_plots
  }
  
  # Create combined plot with shuffle types as rows and files as columns
  combined_plots <- NULL
  for (shuffle_idx in seq_along(all_shuffle_plots)) {
    shuffle_plots <- all_shuffle_plots[[shuffle_idx]]
    
    # Combine plots for this shuffle type (one row)
    row_plots <- NULL
    for (plot_idx in seq_along(shuffle_plots)) {
      if (is.null(row_plots)) {
        row_plots <- wrap_elements(shuffle_plots[[plot_idx]])
      } else {
        row_plots <- row_plots + wrap_elements(shuffle_plots[[plot_idx]])
      }
    }
    
    # Add this row to the combined plot
    if (is.null(combined_plots)) {
      combined_plots <- row_plots
    } else {
      combined_plots <- combined_plots / row_plots
    }
  }
  
  # Create output directory
  plot_dir <- paste0("Plots/", dataset_name, "/")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Save plot
  png(paste0(plot_dir, "NBResampleWholeTaxonShuffledFractionOfSignificantResults(", dataset_name, ").png"), 
      width = 1200 * length(all_files), height = 4800, res = 300)
  print(combined_plots)
  dev.off()
  
  cat(sprintf("Completed %s dataset\n", dataset_name))
  invisible(NULL)
}

# Dataset configuration
datasets <- list(
  list(name = "RDP"),
  list(name = "dada2"),
  list(name = "WGS")
)

# Process all datasets
invisible(lapply(datasets, function(ds) processDatasetType(ds$name)))
