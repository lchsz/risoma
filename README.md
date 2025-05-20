# risoma: Comprehensive Analysis of microRNA Isoforms (isomiRs)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

The `risoma` package provides a comprehensive solution for analyzing microRNA isoforms (isomiRs) with enhanced sensitivity and flexibility. It implements an innovative region-specific threshold strategy to improve detection of complex variants while integrating essential functionalities for isomiR research.

## Key Features

-   **Region-Specific Variant Detection**\
    Independent control of InDel/SNP thresholds across 5' end, seed region, and 3' end

-   **Integrated Analysis Pipeline**\
    miRBase integration → isomiR detection → expression quantification → differential analysis → visualization

-   **Enhanced Sensitivity**\
    Improved "seed-and-extension" algorithm detects 28.6% more unique variants compared to existing tools

-   **Flexible Data Structures**\
    S4 classes (`Isomir`/`IsomirDataSet`) for cross-sample data management

-   **Interactive Visualization**\
    Built-in functions for heatmaps, expression profiles, and sequence alignments

## Installation

### System Requirements

1.  **R Version**: Ensure R (≥ 4.4) is installed.\

2.  **Compilation Tools**:

    -   **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
    -   **macOS**: Install Xcode Command Line Tools (`xcode-select --install` in Terminal).
    -   **Linux**: Install development tools (e.g., `build-essential` for Ubuntu/Debian).

### R Dependencies

``` r
# CRAN Packages
install.packages(c("dplyr", "ggplot2", "pheatmap", "Rcpp", "stringr"))

# Install suggested packages (optional for vignettes/tests)  
install.packages(c("knitr", "rmarkdown", "testthat")) 

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install(c("apeglm", "Biostrings", "DESeq2", "mirbase.db", "msa"))
```

### Install `risoma` from Source

``` bash
# Clone the Git repository:
git clone https://github.com/lchsz/risoma

cd risoma
R CMD build .
R CMD INSTALL risoma_0.1.0.tar.gz
```

## Quick Start

### Basic Analysis Workflow

``` r
library(risoma)

# Load miRNA reference data
vvi_mirnas <- load_mirnas("vvi") # Vitis vinifera miRNAs

# Detect isomiRs from FASTQ files
sample_info <- system.file("extdata", "sample_info.csv", package = "risoma")
fq_dir <- system.file("extdata", "fastqs", package = "risoma")
isomir_data <- detect_isomirs(sample_info, fq_dir, vvi_mirnas, min_tpm = 5)

# Explore detected isoforms
get_isoform_num(isomir_data)  # Isoforms per reference miRNA
head(calc_group_expr_num(isomir_data)) # Expression statistics
```

### Differential Expression Analysis

``` r
# Prepare data for DESeq2
deg_data <- get_deg_data(isomir_data, 
                        treatment = "Flower_FB", 
                        control = "Inflorescence_WD")

# Run DESeq2 workflow
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = deg_data$expr,
  colData = deg_data$sample_info,
  design = ~ group
)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, alpha = 0.05)
```

### Tissue Specificity Analysis

``` r
tsi_df <- calc_tsi(isomir_data)

ggplot2::ggplot(tsi_df, aes(x = type, y = tsi)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1) +
  labs(title = "TSI Distribution by Isoform Type")
```

## Documentation

```         
browseVignettes("risoma")
```
