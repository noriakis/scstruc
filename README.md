
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

## Examples

``` r
library(scran)
library(scstruc)
library(bnlearn)
sce <- mockSCE()
sce <- logNormCounts(sce)
included_genes <- sample(row.names(sce), 100)
gs <- globalStruc(sce, included_genes, return_bn=TRUE, return_data=TRUE)
fitted <- bn.fit(gs[[1]], gs[[2]])
ggraph(bn_fit_to_igraph(fitted), layout="fr") + 
  geom_edge_diagonal(aes(color=coef),
      arrow=arrow(type="closed", length=unit(2,"mm")),
      start_cap=circle(1,"mm"), end_cap=circle(1,"mm"))+
  geom_node_point()+
  theme_graph()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="2400" style="display: block; margin: auto;" />
