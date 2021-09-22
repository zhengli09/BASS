#ifndef _PARTITION_RATIO_H_
#define _PARTITION_RATIO_H_
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]

// Use Monte Carlo method to approximate the ratio of the partition function in log scale
double compPartitionRatio(
  const double betaProposed,
  const double beta,
  const std::map<int, std::set<int>> &V,
  const int K, // number of distint clusters
  const int N, // number of pixels/spots
  const int M, // Monte Carlo iterations
  const int B = 100 // number of burn-in steps in Potts sampling
  );

// Compute number of homogeneous cliques in the Potts model
double compNHC(
  const std::map<int, std::set<int>> &V,
  const Rcpp::IntegerVector &z
  );

// Compute number of homogeneous cliques in the Potts model
double compNHC(
  const std::map<int, std::set<int>> &V,
  const Rcpp::IntegerVector &z
  );


// A helper function in runSW that merges two labels of patches
void mergePatches(
  Rcpp::IntegerVector &patches,
  const int patch1,
  const int patch2
  );

// Sampling from a simple Potts model using Swendsen-Wang algorithm
Rcpp::IntegerMatrix smpPottsSW(
  const double beta,
  const std::map<int, std::set<int>> &V,
  const int K, // number of distint clusters
  const int N, // number of pixels/spots
  const int M, // number of sampling steps
  const int B = 100 // number of burn-in steps
  );

// Working function of Swendsen-Wang algorithm
void runSW(
  const double beta,
  const std::map<int, std::set<int>> &V,
  const int K,
  Rcpp::IntegerVector &z,
  Rcpp::IntegerVector &patches,
  Rcpp::LogicalVector &isBond
  );

// Test Swendsen-Wang algorithm
Rcpp::List testPottsSampling(
  const Rcpp::NumericMatrix coord, // coordinates of pixels
  const double r, // distance between spots to define neighbors
  const double beta,
  const int K, // number of distint clusters
  const int M // number of sampling steps 
  );

// Test whether the interaction paramter beta in Potts model can be estimated
Rcpp::List testPostBeta(
  const Rcpp::IntegerVector x,
  const Rcpp::NumericMatrix coord,
  const double r,
  const int K,
  double beta, // initial value for beta
  const double epsilon,
  const int M,
  const int S
  );
#endif //_PARTITION_RATIO_H_