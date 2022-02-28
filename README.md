# BASS (Bayesian Analytics for Spatial Segmentation)
## Overview
<p align="justify"> BASS performs transcriptomic analyses at two anatomic scales. At the single-cell scale, we perform cell type clustering and cluster cells into multiple cell types. At the tissue regional scale, we segment the tissue into spatial domains in a de novo fashion and characterize the cell type composition in each detected domain. We carry out the two analyses at different scales jointly in a coherent fashion based on a Bayesian hierarchical model, which allows us to seamlessly integrate gene expression information with spatial information to improve the effectiveness of both analyses. Importantly, our method allows for multi-sample integrative analysis of spatial transcriptomic data measured on multiple tissue sections in the same anatomic region. Integrative analysis of multiple spatial transcriptomic data allows us to borrow critical biological information across tissue sections to further enhance analytic performance. </p>

## Installation
Install BASS R package maintained in github through the "devtools" package.
```r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/BASS")
library(BASS)
```

## How to use `BASS`
Follow the [tutorial](https://zhengli09.github.io/BASS/)