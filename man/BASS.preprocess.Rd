% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BASS.R
\name{BASS.preprocess}
\alias{BASS.preprocess}
\title{Pre-process gene expression data}
\usage{
BASS.preprocess(
  BASS,
  doLogNormalize = TRUE,
  geneSelect = c("sparkx", "hvgs"),
  doPCA = TRUE,
  scaleFeature = TRUE,
  doBatchCorrect = TRUE,
  nHVG = 2000,
  nSE = 3000,
  nPC = 20
)
}
\arguments{
\item{BASS}{A BASS object created by \code{createBASSObject}.}

\item{doLogNormalize}{Logical. Whether to perform the library size 
normalization followed by a log2-transformation (after adding a 
pseudo-count of 1).}

\item{geneSelect}{Character. Perform feature selection by either selecting 
spatially variable genes with SPARK-X (Zhu et al., 2021) (\code{sparkx}) 
or highly variable genes with the scran package (\code{hvgs}).}

\item{doPCA}{Logical. Whether to perform PCA on the normalized expression 
matrix to extract \code{nPC} low-dimensional expression features.}

\item{scaleFeature}{Logical. Whether to center and standardize each gene 
to have a mean of 0 and standard deviation of 1 before performing PCA.}

\item{doBatchCorrect}{Logical. Whether to perform batch effect adjustment 
with Harmony (Korsunsky et al., 2019) to align expression data from 
different tissue sections.}

\item{nHVG}{Numeric. Number of highly variable genes to select if 
\code{geneSelect} is set to be \code{hvgs}.}

\item{nSE}{Numeric. Number of spatially variable genes to select if 
\code{geneSelect} is set to be \code{sparkx}.}

\item{nPC}{Numeric. Number of top PCs to extract if \code{doPCA} is set to 
be \code{TRUE}.}
}
\value{
A BASS object with the processed gene expression data stored in the 
  \code{X_run} slot. Refer to \linkS4class{BASS} for a complete list of slots 
  in the object.
}
\description{
Pre-process the gene expression data by first combining cells/spots across 
all L tissue sections, conduct library size normalization followed by a 
log2-transformation (after adding a pseudo-count of 1), select feature genes,  
perform dimension reduction with PCA on the normalized expression matrix to 
extract J low-dimensional expression features, and perform batch effect 
adjustment with Harmony (Korsunsky et al., 2019) to align expression data 
from different tissue sections.
}
