
# DGOplot

<!-- badges: start -->
<!-- badges: end -->

The goal of DGOplot is to perform enrichment analses on Disease Ontology and Gene Ontology and plot the result into a bar plot and a gene association network.

## Installation

You can install the latest version of DGOplot using:
``` r
require("devtools")
install_github("JoelleJee/DGOplot")
library("DGOplot")
```
## Overview

An overview of the package is illustrated below.
![](C:\Users\Joelle\Desktop\BCB410\Assignment\Jee_Y_A1.png)
![](C:\Users\Joelle\Desktop\BCB410\DGOnetplot.jpg)
![](C:\Users\Joelle\Desktop\BCB410\DGObarplot.png)

## Contributions

The author of the package is Yoonsun Jee. The functions available within this package include:

```r
library("DGOplot")
lsf.str("package:DGOplot")
```
- enrichDGO
- DGObarplot
- DGOnetplot

The enrichDGO function performs the DO and GO enrichment analyses using functions from 
DOSE and clusterProfiler R packages. The DOSE R package is used for DO enrichment analysis while the 
clusterProfiler R package is used for GO enrichment analysis.

The DGObarplot function takes in the result from enrichDGO function and plots it into a double bar graph
using ggplot2 R package.

The DGOnetPlot function uses the igraph R package and qgraph R package to graph a gene association network
from the enrichDGO results.

## Example

This a simple example using the dataset provided from DOSE R package.

``` r
library(DGOplot)
## basic example code
data(geneList)
gene <- names(geneList)[abs(geneList) > 2]
x <- enrichDGO(gene, universe = names(geneList))
DGObarPlot(x, showCategory=5)
DGOnetPlot(x, showCategory=5)

```

