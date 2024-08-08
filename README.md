
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scstruc

The package for analysing the gene regulatory network based on the
Bayesian network structure of single-cell transcriptomics data. The
function works with `SingleCellExperiment` and `SpatialExperiment`.

## Installation

Using `devtools`:

``` r
devtools::install_github("noriakis/scstruc")
```

For using the CCDr algorithm in the latest R version, please install
`ccdrAlgorithm` package from the following repository, which removes the
`||` usage error.

``` r
devtools::install_github("noriakis/ccdrAlgorithm")
```

## Examples

``` r
library(scran)
library(scstruc)
library(bnlearn)
sce <- mockSCE()
sce <- logNormCounts(sce)
included_genes <- sample(row.names(sce), 100)
gs <- scstruc(sce, included_genes, changeSymbol=FALSE)
#> Using default bnlearn algorithm
plotNet(gs$net, gs$data, showText=FALSE)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="3600" style="display: block; margin: auto;" />

Using the bootstrapping, the averaged network is obtained.

``` r
library(glmnet)
gs2 <- scstruc(sce, included_genes, algorithm="glmnet_BIC.boot",
               changeSymbol=FALSE, algorithm.args=list("R"=20))
#> Bootstrapping specified
plotAVN(gs2$net, sizeRange=c(1,3))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="2400" style="display: block; margin: auto;" />

Celltype-specific bootstrapped networks.

``` r
sce <- mockSCE()
sce2 <- mockSCE()
sces <- list(sce, sce2)
sces <- lapply(sces, function(x) logNormCounts(x))
sces <- lapply(sces, function(x) {
  colLabels(x) <- sample(c("Celltype_1","Celltype_2"), ncol(x), replace=TRUE)
  return(x)})

included_genes <- sample(row.names(sce), 50)
booted <- lapply(sces, function(x) {
  scstruc(x, included_genes,
    changeSymbol=FALSE,
    labelName="Celltype_1",
    label="label",
    algorithm="glmnet_CV.boot",
    algorithm.args=list("R"=10))
})
#> Bootstrapping specified
#> Bootstrapping specified
ce <- coreBootEdges(booted)
ce
#> # A tbl_graph: 47 nodes and 73 edges
#> #
#> # A directed multigraph with 1 component
#> #
#> # Node Data: 47 × 1 (active)
#>    name     
#>    <chr>    
#>  1 Gene_0008
#>  2 Gene_1713
#>  3 Gene_0059
#>  4 Gene_1113
#>  5 Gene_0065
#>  6 Gene_0078
#>  7 Gene_0124
#>  8 Gene_0918
#>  9 Gene_1040
#> 10 Gene_0184
#> # ℹ 37 more rows
#> #
#> # Edge Data: 73 × 5
#>    from    to strength direction net  
#>   <int> <int>    <dbl>     <dbl> <chr>
#> 1     1     2      0.6     0.667 1    
#> 2     3     4      0.9     0.778 1    
#> 3     5     6      0.6     0.75  1    
#> # ℹ 70 more rows
```
