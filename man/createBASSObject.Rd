% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BASS.R
\name{createBASSObject}
\alias{createBASSObject}
\title{Create a BASS object}
\usage{
createBASSObject(
  X,
  xy,
  C,
  R,
  init_method = c("kmeans", "mclust"),
  psi0 = 1,
  n0 = 1,
  alpha0 = 1,
  k = 4,
  burnin = 2000,
  nsample = 5000,
  beta_method = c("SW", "fix"),
  beta = 1,
  step_size = 0.1,
  tol = 0.001,
  potts_burnin = 10,
  potts_nsample = 10
)
}
\arguments{
\item{X}{An L-list of gene expression count matrices for L tissue sections.
For the lth tissue section, the matrix has dimension PxNl representing P 
common genes shared across tissue sections and Nl cells/spots in the
lth tissue section.}

\item{xy}{An L-list of spatial coordinates matrices for L tissue sections.
For the lth tissue section, the matrix has dimension Nlx2 representing the 
same Nl cells/spots as in the expression count matrix and their x and y 
spatial coordinates.}

\item{C}{Numeric. Number of cell types.}

\item{R}{Numeric. Number of spatial domains.}

\item{init_method}{Character. Initialize the cell type clusters and spatial 
domains with either k-means clustering (\code{kmeans}) or Gaussian mixture 
model (\code{mclust}).}

\item{psi0}{Numeric. psi0*I is the scale matrix of the inv-Wishart prior on 
the variance-covariance matrix, where I is an identity matrix.}

\item{n0}{Numeric. Degrees of freedom of the inv-Wishart prior on the 
variance-covariance matrix.}

\item{alpha0}{Numeric. Concentration parameter of the Dirichlet prior
specified on the cell type composition vector in each spatial domain.}

\item{k}{Numeric. Minimum number of neighbors for each cell/spot based on 
the Euclidean distance.}

\item{burnin}{Numeric. Number of burn-in iterations in the MCMC.}

\item{nsample}{Numeric. Number of posterior samples in the MCMC.}

\item{beta_method}{Character. Fix the cell-cell interaction parameter 
to be \code{beta} (\code{fix}) or estimate the parameter based on the 
data at hand (\code{SW}) .}

\item{beta}{Numeric. Pre-specified cell-cell interaction parameter if 
\code{beta_method} is set to be \code{fix} or otherwise its initial value.}

\item{step_size}{Numeric. Step size of a uniform random walk.}

\item{tol}{Numeric. Threshold of convergence for the cell-cell interaction 
parameter.}

\item{potts_burnin}{Numeric. Number of burn-in iterations in Potts sampling.}

\item{potts_nsample}{Numeric. Number of Potts samples to approximate the 
partition ratio.}
}
\value{
A BASS object. Refer to \linkS4class{BASS} for a complete list of
  slots in the object.
}
\description{
Create a BASS object to store data, settings, hyper-parameters, MCMC 
controlling parameters, and results for the spatial transcriptomic analysis.
}
