#!/bin/bash

# Manuscript Figure Generation Pipeline
# This script runs all analysis and plotting scripts in the correct order
# to generate all figures for the manuscript

echo "=========================================="
echo "Starting Manuscript Figure Generation Pipeline"
echo "=========================================="

# Set error handling
set -e  # Exit on any error

# Function to print section headers
print_section() {
    echo ""
    echo "=========================================="
    echo "$1"
    echo "=========================================="
}

# Function to check if a script completed successfully
check_success() {
    if [ $? -eq 0 ]; then
        echo "SUCCESS: $1 completed successfully"
    else
        echo "ERROR: $1 failed"
        exit 1
    fi
}

# Step 1: Run statistical test calculation scripts
print_section "Step 1: Running Statistical Test Calculations"

echo "Running t-test calculations..."
Rscript -e "source('Scripts/TtestCalculation.R')"
check_success "T-test calculations"

echo "Running Wilcoxon calculations..."
Rscript -e "source('Scripts/WilcoxonCalculation.R')"
check_success "Wilcoxon calculations"

echo "Running DESeq2 calculations..."
Rscript -e "source('Scripts/DESeq2Calculation.R')"
check_success "DESeq2 calculations"

echo "Running edgeR calculations..."
Rscript -e "source('Scripts/edgeRCalculation.R')"
check_success "edgeR calculations"

# Step 2: Generate Log10 P-value plots
print_section "Step 2: Generating Log10 P-value Plots"

echo "Creating Log10 P-value plots..."
Rscript -e "source('Scripts/Log10PvPPlots.R')"
check_success "Log10 P-value plots"

# Step 3: Run shuffle calculation scripts
print_section "Step 3: Running Shuffle Calculations"

echo "Running shuffle sample tags calculations..."
Rscript -e "source('Scripts/ShuffleSampleTags100xPvalCalculation.R')"
check_success "Shuffle sample tags calculations"

echo "Running shuffle counts in sample calculations..."
Rscript -e "source('Scripts/ShuffleInSample100xPvalCalculation.R')"
check_success "Shuffle counts in sample calculations"

echo "Running shuffle counts in taxon calculations..."
Rscript -e "source('Scripts/ShuffleInTaxon100xPvalCalculation.R')"
check_success "Shuffle counts in taxon calculations"

echo "Running shuffle counts table calculations..."
Rscript -e "source('Scripts/ShuffleCountsT100xPvalCalculation.R')"
check_success "Shuffle counts table calculations"

# Step 4: Generate shuffle fraction of significant counts boxplots
print_section "Step 4: Generating Shuffle Fraction of Significant Counts Boxplots"

echo "Creating shuffle fraction of significant counts boxplots..."
Rscript -e "source('Scripts/ShuffleFourWaysPvalSigCountsBoxplots.R')"
check_success "Shuffle fraction boxplots"

# Step 5: Run NB resample calculation scripts
print_section "Step 5: Running NB Resample Calculations"

echo "Running NB resample shuffle sample tags calculations..."
Rscript -e "source('Scripts/NBResampleWholeTaxonShuffleSampleTagPvalCalculation.R')"
check_success "NB resample shuffle sample tags calculations"

echo "Running NB resample shuffle counts in sample calculations..."
Rscript -e "source('Scripts/NBResampleWholeTaxonShuffleInSamplePvalCalculation.R')"
check_success "NB resample shuffle counts in sample calculations"

echo "Running NB resample shuffle counts in taxon calculations..."
Rscript -e "source('Scripts/NBResampleWholeTaxonShuffleInTaxonPvalCalculation.R')"
check_success "NB resample shuffle counts in taxon calculations"

echo "Running NB resample shuffle counts table calculations..."
Rscript -e "source('Scripts/NBResampleWholeTaxonShuffleCountsTPvalCalculation.R')"
check_success "NB resample shuffle counts table calculations"

# Step 6: Generate NB resample fraction of significant counts boxplots
print_section "Step 6: Generating NB Resample Fraction of Significant Counts Boxplots"

echo "Creating NB resample fraction of significant counts boxplots..."
Rscript -e "source('Scripts/NBResampleWholeTaxonShuffleFourWaysPvalSigCountsBoxPlots.R')"
check_success "NB resample fraction boxplots"

# Step 7: Summary and verification
print_section "Step 7: Summary and Verification"

echo "Checking generated plots..."

# Check Log10 P-value plots
echo "Log10 P-value plots:"
for dataset in RDP dada2 WGS RNAseq; do
    if [ -d "Plots/$dataset" ]; then
        count=$(ls Plots/$dataset/Log10PvPPlots*.png 2>/dev/null | wc -l)
        echo "   $dataset: $count plots"
    else
        echo "   $dataset: No plots directory found"
    fi
done

# Check shuffle fraction boxplots
echo "Shuffle fraction boxplots:"
for dataset in RDP dada2 WGS RNAseq; do
    if [ -f "Plots/$dataset/FractionOfSignificantResults($dataset).png" ]; then
        echo "   $dataset: SUCCESS - Generated"
    else
        echo "   $dataset: ERROR - Missing"
    fi
done

# Check NB resample fraction boxplots
echo "NB resample fraction boxplots:"
for dataset in RDP dada2 WGS RNAseq; do
    if [ -f "Plots/$dataset/NBResampleWholeTaxonShuffledFractionOfSignificantResults($dataset).png" ]; then
        echo "   $dataset: SUCCESS - Generated"
    else
        echo "   $dataset: ERROR - Missing"
    fi
done

echo ""
echo "=========================================="
echo "Manuscript Figure Generation Complete!"
echo "=========================================="
echo ""
echo "Generated figures:"
echo "• Log10 P-value comparison plots (4 datasets × multiple files each)"
echo "• Shuffle fraction of significant counts boxplots (4 datasets)"
echo "• NB resample fraction of significant counts boxplots (4 datasets)"
echo ""
echo "All figures are saved in the Plots/ directory structure." 