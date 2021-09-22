# BASS (Bayesian Analytics for Spatial Segmentation)
## Overview
<p align="justify"> BASS performs transcriptomic analyses at two different scales. At the single cell scale, we perform cell type clustering and cluster cells into multiple cell types. At the tissue regional scale, we segment tissues into tissue structures in a de novo fashion and characterize the cell type composition in each detected tissue region separately. Importantly, the analyses at the two different scales are performed in a coherent way through a Bayesian hierarchical model, which allows us to seamlessly integrate gene expression information with spatial information to improve analyses at both scales. In addition, data of multiple tissue sections can be integrated together by BASS to further enhance the model performance. </p>

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