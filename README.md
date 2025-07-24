# Differential Analysis Project

This project contains R scripts for performing differential analysis on various types of biological data including:
- RDP (Ribosomal Database Project)
- DADA2 (Divisive Amplicon Denoising Algorithm)
- WGS (Whole Genome Sequencing)
- RNA-seq (RNA Sequencing)

## Project Structure

```
DAProject/
├── Data/
│   ├── CountsTable/      # Contains normalized and raw count data
│   └── MetaData/         # Contains metadata files
├── Plots/                # Output directory for generated plots
├── Results/              # Output directory for analysis results
└── Scripts/             # R scripts for analysis
    ├── CalculationTools.R    # Common functions and utilities
    ├── Log10PvPPlots.R      # Plotting functions
    ├── PlotTools.R          # Additional plotting utilities
    ├── TtestCalculation.R   # T-test analysis
    ├── WilcoxonCalculation.R # Wilcoxon test analysis
    ├── DESeq2Calculation.R   # DESeq2 differential expression analysis
    └── edgeRCalculation.R    # edgeR differential expression analysis
```

## Analysis Methods

The project implements multiple statistical methods for differential analysis:
1. T-test
2. Wilcoxon test
3. DESeq2
4. edgeR

Each analysis method is implemented in its own script, with common functions sourced from `CalculationTools.R`. 