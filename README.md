# BASS: Bayesian Analytics for Spatial Segmentation
## Overview
![BASS Overview](https://github.com/zhengli09/BASS-Analysis/blob/master/docs/BASS_workflow.png)
<p align="justify"> Spatial transcriptomic studies perform gene expression 
profiling on tissues with spatial localization information. With technological 
advances, many spatial transcriptomic studies are reaching single-cell spatial 
resolution and are capable of collecting data from multiple tissue sections. 
Here, we develop a computational method, BASS, that enables multi-scale and 
multi-sample analysis for single-cell resolution spatial transcriptomics. 
Specifically, BASS performs multi-scale analyses in the form of cell type 
clustering at the single-cell scale and spatial domain detection at the tissue 
regional scale. The two analytic tasks are carried out in a coherent fashion 
through a Bayesian hierarchical modeling framework. For both analytic tasks, 
BASS properly accounts for the spatial correlation structure and seamlessly 
integrates gene expression information with spatial localization information to 
enhance analytic performance. In addition, BASS is capable of performing 
multi-sample analysis via joint modeling of multiple tissue sections/samples, 
facilitating the cross-sample integration of spatial transcriptomics. </p>

## Installation
BASS is implemented as an R (>= 4.0.3) package with underlying efficient C++ 
code interfaced through Rcpp and RcppArmadillo. BASS depends on a few other 
R packages that include GIGrvg, Matrix, harmony, label.switching, mclust, Rcpp, 
RcppArmadillo, RcppDist, SPARK, scran, and scater. Please refer to the package 
[DESCRIPTION](https://github.com/zhengli09/BASS/blob/master/DESCRIPTION) file 
for details. Dependent packages are supposed to be automatically installed while 
installing BASS.

Install the BASS R package maintained in github through the `devtools` package.
```r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/BASS")
library(BASS)
?BASS
```

## How to use `BASS`
Check our [vignettes](https://zhengli09.github.io/BASS-Analysis/).

## Citing the work
If you find the `BASS` package, any of the source code, or processed data 
in this repository and in the [BASS-analysis](https://github.com/zhengli09/BASS-Analysis) 
repository useful for your work, please cite:

> Zheng Li, Xiang Zhou# (2022). Multi-scale and multi-sample analysis enables 
> accurate cell type clustering and spatial domain detection in spatial 
> transcriptomic studies. Under revision.

Visit our [group website](http://www.xzlab.org) for more statistical tools on 
analyzing spatial transcriptomic data.

## Release Notes
### v1.1.0:
* Changed functional programming to OO programming
* Optimized the implementation for Swendsen-Wang algorithm
* Updated threshold and message for checking the convergence of the spatial 
parameter beta