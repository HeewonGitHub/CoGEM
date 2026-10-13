# CoGEM

This repository contains the R code for CoGEM, a framework for identifying phenotype-associated gene modules using communicability-guided random forests.

---

## Overview

CoGEM uses protein-protein interaction (PPI) information to guide feature sampling in a random forest. Gene importance is calibrated against independently permuted shadow data, and the selected gene module is evaluated with a phenotype-permutation test.

The analysis is organized into six steps:

1. Set analysis parameters and load the input data.
2. Construct a PPI network and calculate gene communicability.
3. Fit a communicability-guided random forest and obtain out-of-bag predictions.
4. Calculate observed gene importance from out-of-bag permutation loss.
5. Compare observed importance with shadow-based empirical null values and select genes.
6. Test the selected module using joint module permutation and phenotype permutation.

---

## Requirements

The script uses R and the following packages:

```r
library(OmnipathR)
library(readr)
library(dplyr)
```

`OmnipathR::bioplex_all()` retrieves BioPlex interactions through OmniPathR, so an internet connection is required when the script runs.

---

## Input Data

Place the following files in the working directory before running `sampleR_CoGEM.R`:

- `EXPRESSION.csv`: numeric gene-expression matrix with samples in rows and genes in columns.
- `Phenotype.csv`: one binary phenotype value (`0` or `1`) per sample, in the same sample order as `EXPRESSION.csv`.

Gene column names in the expression matrix must correspond to gene symbols used in the BioPlex interaction data. The script uses `colnames(X)` to match expression features to PPI genes. Because `read.table()` is called without `header = TRUE`, check the imported column names and adjust the input file or import settings if needed.

---

## Step-by-Step Usage

Before sourcing the script, set `it` to an integer run index. It is used in random-number seeds.

```r
it <- 1
source("sampleR_CoGEM.R")
```

### Step 1. Parameters and Input Data

The script loads the required packages and defines the forest, shadow, permutation, and gene-selection settings. It then retrieves BioPlex interactions and reads the expression and phenotype files.

| Parameter | Default | Description |
|---|---:|---|
| `ntree` | 199 | Number of trees in each forest. |
| `B` | 199 | Number of shadow replicates used to calculate gene-level empirical p-values. |
| `Omega` | 199 | Number of phenotype permutations used to calculate the module p-value. |
| `mtry` | 10 | Number of candidate genes considered at each node. |
| `kC` | 50 | Maximum number of communicability-linked candidate genes per gene. |
| `max_depth` | 8 | Maximum tree depth. |
| `min_leaf` | 5 | Minimum number of samples allowed in each child node. |
| `alpha` | 0.01 | Cutoff applied to gene-level empirical p-values for selecting genes. |

Larger values of `B` and `Omega` provide finer empirical p-value resolution and increase computation time. The script selects genes using `p_gene < alpha`; it does not apply a multiple-testing correction at this step.

### Step 2. Network and Communicability

The PPI data are restricted to genes in the expression matrix and represented as an undirected adjacency matrix. The script calculates the degree-normalized adjacency matrix, obtains its matrix exponential, and sets self-communicability to zero. For each gene, it retains up to `kC` genes with the highest communicability and normalizes their values to form feature-sampling probabilities.

### Step 3. Communicability-Guided Random Forest

The script grows bootstrap trees using Gini impurity to select splits. Candidate genes at the root are sampled uniformly. At each non-root node, candidate genes are sampled according to their communicability-based probabilities relative to the parent split gene. The fitted forest is used to calculate out-of-bag (OOB) predictions, an OOB classification table, and the OOB Brier loss.

### Step 4. Observed Gene Importance

For each tree, the script identifies genes used for splitting and permutes each such gene across that tree's OOB samples. Gene importance is the increase in Brier loss after permutation, averaged over the forest. The script prints the 20 genes with the largest observed importance values.

### Step 5. Shadow-Calibrated Gene Identification

For each of the `B` shadow replicates, every gene column is independently permuted across samples. A forest is fitted to each shadow matrix, and OOB permutation importance is calculated. The observed importance is compared with the corresponding shadow importance values to obtain a gene-level empirical p-value:

\[
p_j = \frac{1 + \sum_{b=1}^{B} I(I_{j}^{\mathrm{null},(b)} \geq I_j^{\mathrm{obs}})}{B+1}.
\]

Genes with `p_gene < alpha` are retained for the module-level test in Step 6.

### Step 6. Module-Level Phenotype Association Test

The selected genes are held fixed as the candidate module. First, the script calculates the observed module statistic as the increase in OOB Brier loss after jointly permuting the module genes across samples. It then permutes the phenotype `Omega` times, refits a forest for each permuted phenotype, and recalculates the module statistic. The empirical module p-value is calculated as

\[
p_{\mathrm{module}} = \frac{1 + \sum_{r=1}^{\Omega} I(Q_r^{\mathrm{null}} \geq Q^{\mathrm{obs}})}{\Omega+1}.
\]

---

## Expected Output

The script prints the OOB classification table and Brier loss, the 20 highest-importance genes, and the selected module size, observed module statistic, and empirical module p-value. Intermediate results, including `gene_importance`, `p_gene`, `candidate_genes`, and `Q_null`, remain available in the R environment. The script does not write these results to separate files.

---

## Notes

- The phenotype must be binary and coded as `0` and `1` for the Gini split and Brier-loss calculations.
- `it` must be defined before execution.
- The default values of `B` and `Omega` are intended for a sample run; increase them for finer empirical p-value resolution, allowing for longer run times.
