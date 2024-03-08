
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scstruc

The package for analysing the gene regulatory network based on the
structure of single-cell transcriptomics data, inferred by Bayesian
network. The function works with `SingleCellExperiment` and
`SpatialExperiment`.

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
gs <- globalStruc(sce, included_genes, return_bn=TRUE, return_data=TRUE, change_symbol=FALSE, reg=FALSE)
#> Using default MMHC
fitted <- bn.fit(gs[[1]], gs[[2]])
ggraph(bn_fit_to_igraph(fitted), layout="fr") + 
  geom_edge_diagonal(aes(color=coef),
      arrow=arrow(type="closed", length=unit(2,"mm")),
      start_cap=circle(1,"mm"), end_cap=circle(1,"mm"))+
  geom_node_point()+
  theme_graph()
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="2400" style="display: block; margin: auto;" />

Using the bootstrapping, the averaged network is obtained.

``` r
library(glmnet)
gs2 <- globalStruc(sce, included_genes, return_bn=TRUE, penalty="glmnet_BIC.boot",
                  change_symbol = FALSE, algorithm.args=list("R"=100))
#> Bootstrapping specified
plotAVN(gs2)
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
booted <- lapply(sces, function(x) perLabelStruc(x, label_name = "Celltype_1",
     algorithm.args=list("R"=10), barcode_column="row",
     candidate_genes=included_genes, nonzero=0.5,
     penalty="glmnet_CV.boot"))
#> Celltype_1 
#> Bootstrapping specified
#> Celltype_1 
#> Bootstrapping specified
ce <- coreEdges(booted, "Celltype_1")
ce
#> # A tbl_graph: 33 nodes and 36 edges
#> #
#> # A directed simple graph with 4 components
#> #
#> # Node Data: 33 × 1 (active)
#>    name     
#>    <chr>    
#>  1 Gene_0074
#>  2 Gene_0918
#>  3 Gene_0109
#>  4 Gene_0306
#>  5 Gene_1059
#>  6 Gene_0280
#>  7 Gene_1234
#>  8 Gene_1917
#>  9 Gene_0356
#> 10 Gene_1562
#> # ℹ 23 more rows
#> #
#> # Edge Data: 36 × 4
#>    from    to strength direction
#>   <int> <int>    <dbl>     <dbl>
#> 1     1     2      0.7     0.714
#> 2     3     4      0.7     0.786
#> 3     3     2      0.9     0.667
#> # ℹ 33 more rows
```
