# Plot log10 pvalue plots
rm(list = ls())
source("Scripts/PlotTools.R")
library(ggplot2)
library(dplyr)
library(patchwork)

# Define dataset configurations
datasets <- list(
  list(
    name = "RDP",
    input_dir = "Data/CountsTable/RDPRaw",
    results_dir = "Results/PkgResults/RDP",
    plot_dir = "Plots/RDP"
  ),
  list(
    name = "dada2",
    input_dir = "Data/CountsTable/dada2Raw",
    results_dir = "Results/PkgResults/dada2",
    plot_dir = "Plots/dada2"
  ),
  list(
    name = "WGS",
    input_dir = "Data/CountsTable/WGSRaw",
    results_dir = "Results/PkgResults/WGS",
    plot_dir = "Plots/WGS"
  ),
  list(
    name = "RNAseq",
  )
)

# Function to read results file safely
read_results_safely <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  read.table(file_path, header = T, row.names = 1)
}

# Function to create all comparison plots for a single dataset
create_comparison_plots <- function(name, results_dir, file_name) {
  # Read all results
  tRES <- read_results_safely(file.path(results_dir, "ttest", paste0(file_name, "_t.txt")))
  wRES <- read_results_safely(file.path(results_dir, "Wilcoxon", paste0(file_name, "_wilcox.txt")))
  dRES <- read_results_safely(file.path(results_dir, "DESeq2", paste0(file_name, "_deseq2.txt")))
  eRES <- read_results_safely(file.path(results_dir, "edgeR", paste0(file_name, "_edger.txt")))
  
  # Process results
  tRES <- processTtestRes(tRES)
  wRES <- processWilcoxonRes(wRES)
  dRES <- processDESeq2Res(dRES)
  eRES <- processEdgeRRes(eRES)
  
  # Create all pairwise comparison plots
  plots <- list(
    # tvw
    p1 = Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", 
                      paste0("(", name, "-", file_name, ") t-test vs. Wilcoxon")),
    # tvd
    p2 = Log10PvPPlot(tRES, dRES, "t-test", "DESeq2",
                      paste0("(", name, "-", file_name, ") t-test vs. DESeq2")),
    # tve
    p3 = Log10PvPPlot(tRES, eRES, "t-test", "edgeR",
                      paste0("(", name, "-", file_name, ") t-test vs. edgeR")),
    # wvd
    p4 = Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2",
                      paste0("(", name, "-", file_name, ") Wilcoxon vs. DESeq2")),
    # wve
    p5 = Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR",
                      paste0("(", name, "-", file_name, ") Wilcoxon vs. edgeR")),
    # dve
    p6 = Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR",
                      paste0("(", name, "-", file_name, ") DESeq2 vs. edgeR"))
  )
  
  # Combine all plots
  combined_plots <- plots$p1 + plots$p2 + plots$p3 + plots$p4 + plots$p5 + plots$p6
  combined_plots + plot_layout(ncol = 3)
}

# Process each dataset
for (ds in datasets) {
  # Create output directory if it doesn't exist
  if (!dir.exists(ds$plot_dir)) {
    dir.create(ds$plot_dir, recursive = TRUE)
  }
  
  # Get all input files
  all_names <- gsub(
    basename(list.files(ds$input_dir, pattern = "*.txt", full.names = TRUE)),
    pattern = ".txt$",
    replacement = ""
  )
  
  # Process each file in the dataset
  for (name in all_names) {
    # Create plot
    png(
      paste0(ds$plot_dir, "/Log10PvPPlots(", name, ").png"),
      width = 5000,
      height = 2400,
      res = 300
    )
    par(mar = c(5, 6, 4, 1) + .1)
    
    # Generate and print plots
    tryCatch({
      print(create_comparison_plots(ds$name, ds$results_dir, name))
    }, error = function(e) {
      message(paste("Error processing", name, "in", ds$name, "dataset:", e$message))
    })
    
    dev.off()
  }
}
