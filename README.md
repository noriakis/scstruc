
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scstruc

The package for analysing the gene regulatory network based on Bayesian
network structure of single-cell transcriptomics data. The function
works primarily with `SingleCellExperiment`. Multiple algorithms
tailored for single-cell transcriptomics data are prepared. The inferred
networks are validated based on the causal relationships between genes.

<!-- badges: start -->

[![R-CMD-check](https://github.com/noriakis/scstruc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noriakis/scstruc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

Using `devtools`:

``` r
devtools::install_github("noriakis/scstruc")
```

## Examples

Based on `SingleCellExperiment`, the network is inferred and plotted.

``` r
library(scran)
library(scstruc)
library(ggraph)
library(bnlearn)
sce <- mockSCE()
sce <- logNormCounts(sce)
included_genes <- sample(row.names(sce), 100)
gs <- scstruc(sce, included_genes, changeSymbol=FALSE)
plotNet(gs$net, gs$data, showText=FALSE)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="3600" style="display: block; margin: auto;" />

Using bootstrapping, the averaged network is obtained. This time, `L1MB`
algorithm with the selection of neighbors based on BIC is used.

``` r
library(glmnet)
gs2 <- scstruc(sce, included_genes, algorithm="glmnet_BIC", boot=TRUE,
               changeSymbol=FALSE, R=20)
#> Bootstrapping specified
plotAVN(gs2$net, sizeRange=c(1,3))
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="2400" style="display: block; margin: auto;" />
