# Differential Analysis Project

Comprehensive differential analysis pipeline for biological data including RDP, DADA2, WGS, and RNA-seq datasets.

## Quick Start

```bash
# Run the complete manuscript figure generation pipeline
./run_pipeline.sh
```

## Project Structure

```
DAProject/
├── Data/
│   ├── CountsTable/      # Raw and normalized count data
│   └── MetaData/         # Metadata files
├── Plots/                # Generated plots
├── Results/              # Analysis results
└── Scripts/             # Analysis scripts
    ├── CalculationTools.R    # Common functions
    ├── PlotTools.R          # Plotting utilities
    ├── *Calculation.R       # Statistical test scripts
    ├── Shuffle*Calculation.R # Shuffling analysis scripts
    └── NB*Calculation.R     # NB resampling scripts
```

## Analysis Methods

### Statistical Tests
- **T-test** - Parametric comparison
- **Wilcoxon** - Non-parametric comparison  
- **DESeq2** - Differential expression
- **edgeR** - Differential expression

### Shuffling Analysis
- **Sample Tags** - Random sample reassignment
- **Counts in Sample** - Within-sample shuffling
- **Counts in Taxon** - Within-taxon shuffling
- **Counts Table** - Complete table shuffling

### NB Resampling
- **Whole Taxon Resampling** - Negative binomial resampling
- **Combined with Shuffling** - Applied to resampled data

## Pipeline Steps

1. **Statistical Tests** - T-test, Wilcoxon, DESeq2, edgeR
2. **Log10 P-value Plots** - Method comparisons
3. **Shuffle Analysis** - Robustness assessment
4. **NB Resample Analysis** - Distribution simulation
5. **Verification** - Output validation

## Output Files

### Results
- `Results/PkgResults/[dataset]/[method]/[file]_[method].txt`

### Plots
- `Plots/[dataset]/Log10PvPPlots([file]).png`
- `Plots/[dataset]/FractionOfSignificantResults([dataset]).png`
- `Plots/[dataset]/NBResampleWholeTaxonShuffledFractionOfSignificantResults([dataset]).png`

## Dependencies

### R Packages
- `DESeq2`, `edgeR` - Differential expression
- `coin` - Non-parametric tests
- `ggplot2`, `patchwork`, `ggExtra` - Plotting

### System
- R (≥3.6)
- Unix-like environment
- Sufficient disk space

## Individual Scripts

```bash
# Statistical tests
Rscript -e "source('Scripts/TtestCalculation.R')"
Rscript -e "source('Scripts/WilcoxonCalculation.R')"
Rscript -e "source('Scripts/DESeq2Calculation.R')"
Rscript -e "source('Scripts/edgeRCalculation.R')"

# Plotting
Rscript -e "source('Scripts/Log10PvPPlots.R')"
Rscript -e "source('Scripts/ShuffleFourWaysPvalSigCountsBoxplots.R')"
Rscript -e "source('Scripts/NBResampleWholeTaxonShuffleFourWaysPvalSigCountsBoxPlots.R')"
```

## Features

- **Modular Design** - Common functions in `CalculationTools.R`
- **Error Handling** - Robust error recovery
- **Progress Reporting** - Clear execution feedback
- **Automated Pipeline** - Single command execution
- **Comprehensive Output** - Statistical results and visualizations 