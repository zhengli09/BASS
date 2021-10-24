// Author: Zheng Li
// Date: 2020-08-13
// Functions used in Gibbs sampler
#ifndef _GIBBS_FUNCTIONS_H_
#define _GIBBS_FUNCTIONS_H_
#include <RcppDist.h>

// Initialize parameters
void initParams(
  arma::mat &mu, // J x C
  Rcpp::IntegerVector &z, // N x 1
  Rcpp::IntegerVector &c, // N x 1
  double &beta,
  arma::vec &lambda, // J x 1
  arma::vec &d, // J x 1
  arma::vec &R2, // J x 1
  Rcpp::String initMethod, // initialization method for c and z
  const arma::mat &X, // N x J 
  const int R
  );

// Draw samples from a Dirichlet distribution
Rcpp::NumericMatrix rdirichlet(
  const int numSmps, 
  const Rcpp::NumericVector alphas
  );

// An array of arrays that contain cells of each cell type
std::map<int, std::set<int>> updateT(
  const Rcpp::IntegerVector &c // N x 1 vector of cell types for N cells
  );

// An array of arrays that contain cells of each spatial cluster
std::map<int, std::set<int>> updateS(
 const Rcpp::IntegerVector &z // N x 1 vector of spatial clusters for N cell
  );

// An array of arrays that contain cells belonging to a unique cell type
// and spatial cluster combination
std::map<std::string, std::set<int>> updateU(
  const std::map<int, std::set<int>> &S, 
  const std::map<int, std::set<int>> &T
  );

// A graph represented by an Adjacency list, with neighbor defined as
// cells with distance less than the threshold value r 
std::map<int, std::set<int>> getNeighbors(
  const double r, // threshold distance for two cells to be neighbors
  const Rcpp::NumericMatrix &coord // N x 2 coordinate matrix
  );

// Construct a KNN graph violently
std::map<int, std::set<int>> constructKNN(
  const int K, // number of neighbors for each spot/cell
  const Rcpp::NumericMatrix &coord // N x 2 coordinate matrix
  );

// Get average number of neighbors per cell
double avgNeighbors(
  const double r, 
  const Rcpp::NumericMatrix &coord
  );

// Update mean parameter of gene expression
void updateMu(
  const Rcpp::NumericMatrix &X, // N x J gene expression matrix, 
                                // X(i, j) is expression value of gene j of cell i
  std::map<int, std::set<int>> &T,
  const Rcpp::NumericMatrix &sigma2, // J x C matrix
  const double kappa, // hyper-parameter in mean prior
  Rcpp::NumericMatrix &mu // J x C matrix
  );

// Update mean paramter of gene expression. A hierarchical prior is placed on the 
// prior mean and prior variance. A flat normal prior is placed on the prior mean
// and a shrinkage prior on the prior variance.
// Reference: Model-based clustering based on sparse finite Gaussian mixtures.
void updateMuNG(
  const arma::mat &X, // N x J gene expression matrix, 
                      // X(i, j) is expression value of gene j of cell i
  std::map<int, std::set<int>> &T,
  const Rcpp::NumericMatrix &sigma2, // J x C matrix
  const arma::vec &d, // J x 1 vector of prior mean of mu
  const arma::vec &lambda, // J x 1 vector of scaling factors of prior variance of mu
  const arma::vec &R2, // J x 1 vector of squared range of each feature
  arma::mat &mu // J x C matrix
  );

// Update mean parameter under EEE covariance structure
void updateMuEEE(
  const arma::mat &X, // N x J gene expression matrix, 
                      // X(i, j) is expression value of gene j of cell i
  std::map<int, std::set<int>> &T,
  const arma::mat &Sigmainv, // J x J matrix
  const arma::vec &d, // J x 1 vector of prior mean of mu
  const arma::vec &lambda, // J x 1 vector of scaling factors of prior variance of mu
  const arma::vec &R2, // J x 1 vector of squared range of each feature
  arma::mat &mu // J x C matrix
  );

// Update scaling factors of prior variance of mu
void updateLambda(
  double v1, // shape paramter of a gamma distribution
  double v2, // rate paramter of a gamma distribution
  const arma::mat &mu, // J x C matrix
  const arma::vec &d, // J x 1 vector of prior mean of mu
  const arma::vec &R2, // J x 1 vector of squared range of each feature
  arma::vec &lambda // J x 1 vector of scaling factors of prior variance of mu
  );

// Update prior mean of mu
void updateD(
  const arma::mat &mu, // J x C matrix
  const arma::vec &lambda, // J x 1 vector of scaling factors of prior variance of mu
  const arma::vec &R2, // J x 1 vector of squared range of each feature
  arma::vec &d // J x 1 vector of prior mean of mu
  );

// Update variance parameter of gene expression
void updateSigma2(
  const Rcpp::NumericMatrix &X, 
  std::map<int, std::set<int>> &T, 
  const Rcpp::NumericMatrix &mu, 
  const double a, 
  const double b, 
  Rcpp::NumericMatrix &sigma2
  );

// Update variance parameter under EII covariance structure
void updateSigma2EII(
  const arma::mat &X,
  std::map<int, std::set<int>> &T,
  const arma::mat &mu,
  const double a,
  const double b,
  Rcpp::NumericMatrix &sigma2
  );

// Update variance parameter under EEE covariance structure
void updateSigmaEEE(
  const arma::mat &X,
  const Rcpp::IntegerVector &c,
  const arma::mat &mu,
  const arma::mat &W0inv,
  const int n0,
  arma::mat &Sigmainv
  );

// Update probability of cell types given spatial clusters
void updatePi(
  std::map<std::string, std::set<int>> &U, 
  const double alpha0, // concentration parameter in Dirichlet distribution
  Rcpp::NumericMatrix &pi // C x K matrix and pi[c, k] represents the 
                          // probability of being cell type c given spatial
                          // cluster k 
  );

// Update cell types
void updateC(
  const arma::mat &X, // N x J
  const arma::mat &mu, // J x C
  const Rcpp::NumericMatrix &sigma2, // J x C
  const Rcpp::NumericMatrix &pi, // C x K
  const Rcpp::IntegerVector &z, // N x 1
  Rcpp::IntegerVector &c // N x 1
  );

// Update cell types with the extended covariance structure:
// Equal unconstrained covariance matrix for each cluster
void updateCEEE(
  const arma::mat &X, // N x J
  const arma::mat &mu, // J x C
  const arma::mat &Sigmainv, // J x J
  const Rcpp::NumericMatrix &pi, // C x K
  const Rcpp::IntegerVector &z, // N x 1
  Rcpp::IntegerVector &c // N x 1
  );

// Update tissue structure labels
void updateZ(
  std::map<int, std::set<int>> &V, 
  const Rcpp::NumericMatrix &pi, // C x K
  const Rcpp::IntegerVector &c, // N x 1
  const double beta, // interaction parameter in Potts
  Rcpp::IntegerVector &z // N x 1
  );

// Update tissue structure labels with Swendsen-Wang algorithm
void updateZSW(
  std::map<int, std::set<int>> &V,
  const Rcpp::NumericMatrix &pi, // C x K
  const Rcpp::IntegerVector &c, // N x 1
  const double beta, // interaction parameter in Potts
  Rcpp::IntegerVector &z // N x 1
  );

// Update Potts interaction parameter
bool updateBeta(
  const Rcpp::IntegerVector &z, 
  const std::map<int, std::set<int>> &V, 
  const int K, 
  const double epsilon, // uniform random walk parameter
  const int M, // Monte Carlo iterations to approximate partition ratio
  const int B, // number of burn-in steps in Potts sampling
  const double betaMax, // Upper bound of beta
  double &beta
  );

// Update Potts interaction parameter in a fast way
bool updateBetaFast(
  const Rcpp::IntegerVector &z, 
  const std::map<int, std::set<int>> &V, 
  const double epsilon, // uniform random walk parameter
  const Rcpp::NumericMatrix &NHC, // number of homogeneous cliques under different beta
  double &beta
  );

// Check convergence of beta
bool checkBetaConvergence(
  const Rcpp::NumericVector &burninBeta,
  int iter,
  double tol = 0.1
  );

// Compute Pseudo-likelihood of a Potts model in log scale
double evalLogSL(
  const double beta, 
  const Rcpp::IntegerVector &z, 
  std::map<int, std::set<int>> &V, 
  const int K
  );

// debugging functions
void printVs(std::vector<std::map<int, std::set<int>>> Vs);

// update all parameters (a wrapping function)
void updateAll(
  const arma::mat &X,
  const int L,
  const double alpha0,
  const Rcpp::String covStruc,
  const double a,
  const double b,
  const arma::mat &W0inv,
  const int n0,
  const double kappa,
  const double beta,
  const arma::vec &R2,
  const Rcpp::IntegerVector &smpIdx,
  std::vector<std::map<int, std::set<int>>> &Vs,
  Rcpp::IntegerVector &zl,
  Rcpp::IntegerVector &cl,
  Rcpp::NumericMatrix &pi,
  Rcpp::NumericMatrix &sigma2,
  arma::mat &Sigmainv,
  arma::mat &mu,
  Rcpp::IntegerVector &c,
  Rcpp::IntegerVector &z,
  arma::vec &lambda,
  arma::vec &d,
  std::map<int, std::set<int>> &T,
  std::map<int, std::set<int>> &S,
  std::map<std::string, std::set<int>> &U
  );
#endif