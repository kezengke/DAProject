# Differential Analysis Project

Comprehensive differential analysis pipeline for biological data including RDP, DADA2, WGS, and RNA-seq datasets.

## Quick Start

```bash
# Run the complete manuscript figure generation pipeline (all datasets)
./run_pipeline.sh

# Run the pipeline without RNAseq datasets (faster for testing)
./run_pipeline_no_rnaseq.sh
```

### Pipeline Options

- **Full Pipeline** (`./run_pipeline.sh`): Processes all 4 types of datasets (RDP, dada2, WGS, RNAseq)
- **Fast Pipeline** (`./run_pipeline_no_rnaseq.sh`): Processes 3 types datasets (RDP, dada2, WGS) - much faster for testing

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

