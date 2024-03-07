
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
