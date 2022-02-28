// Zheng Li
// 2020-08-13
// Gibbs sampler of a full bayesian spatial domain segmentation model
// Bayesian Analytics for Spatial Segmentation (BASS)
// Change log:
// 2021-01-02: Cannot use std::vector.pushback(Rcpp::NumericMatrix), this is pass-by-reference
//             and all the values will be the last one. Instead, use cube data structure in 
//             Armadillo.
// 2021-04-18: Extend BASS to multiple samples
// 2021-04-18: Wrapped parameter updates into one function "updateAll"

#include "compPartitionRatio.h"
#include "GibbsFuncs.h"
#include <RcppDist.h>
#include <string>
using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
Rcpp::List Gibbs(
  const Rcpp::NumericMatrix Xin, // N x J gene expression matrix
  Rcpp::NumericMatrix xy, // N x 2 coordinates
  const Rcpp::IntegerVector Ns, // L x 1 vector of sample sizes
  const int C, // number of cell types
  const int R, // number of spatial domains
  Rcpp::String initMethod, // initialization method for c and z
  Rcpp::String covStruc, // covariance structure for expression features
  const double kappa, // prior parameter of gene expression (not used)
  const double alpha0, // concentration parameter of Dirichlet distribution
  const double a, // shape parameter of inverse gamma prior
  const double b, // scale parameter of inverse gamma prior
  const Rcpp::NumericMatrix W0in, // scale matrix of Wishart prior
  const int n0, // degrees of freedom of Wishart prior
  const int k, // number of neighbors for kNN graph
  const int warmUp, // number of burn-in iterations
  const int numSamples, // number of posterior samples
  Rcpp::String betaEstApproach, // approach to estimate beta
  const double betaIn, // fixed value of beta if the betaEstApproach is set to be "FIXED"
  const double betaMax, // upper bound of beta
  const double epsilon, // uniform random walk step size
  const double betaTol, // tolerance for checking beta convergence
  const int M, // Monte Carlo iterations to approximate partition ratio
  const int B // number of burn-in steps in Potts sampling
  )
{
  arma::mat X = Rcpp::as<arma::mat>(Xin);
  arma::mat W0inv = arma::inv(Rcpp::as<arma::mat>(W0in));
  int N = X.n_rows; // total number of cells in L samples
  int J = X.n_cols; // number of genes shared across L samples
  int L = Ns.length(); // number of samples
  Rcpp::IntegerVector smpIdx = Rcpp::cumsum(Ns);
  smpIdx.push_front(0);

  std::map<int, std::set<int>> S; // spatial domain s -> indeces of cells belonging to the domain
  std::map<int, std::set<int>> T; // cell type c -> indeces of cells belonging to the type
  std::map<std::string, set<int>> U; // "sc" -> indeces of cells in spatial domain s and of cell type c
  std::vector<std::map<int, std::set<int>>> Vs; // L-vector of adjacency list for each sample

  arma::mat mu(J, C); // J x C mean parameter of expression feature
  Rcpp::NumericMatrix sigma2(J, C); // J x C variance of expression of feature j in cell type c
  arma::mat Sigmainv(J, J); // J x J covariance matrix (EEE model)
  Rcpp::NumericMatrix pi(C, R); // C x R probability of being cell type c in spatial domain r
  Rcpp::IntegerVector z(N); // N x 1 vector of spatial domain labels
  Rcpp::IntegerVector zl; // Ns x 1 vector of spatial domain labels (sub-vector of z)
  Rcpp::IntegerVector initZ(N); // N x 1 vector of initial spatial domains from k-means
  Rcpp::IntegerVector c(N); // N x 1 vector of cell types
  Rcpp::IntegerVector cl; // Ns x 1 vector of cell types (sub-vector of c)
  Rcpp::IntegerVector initC(N); // N x 1 vector of initial cell types from k-means
  arma::vec lambda(J); // J x 1 vector of scaling factors of prior variance of mu
  arma::vec d(J); // J x 1 vector of prior mean of mu
  arma::vec R2(J); // J x 1 vector of squared range of each gene
  double beta = betaIn; // scalar, interaction paramter of Potts model
  double acceptBetaP = 0.0; // scalar, acceptance probability of beta
  bool isConverged = false; // convergence check of beta

  // posterior samples
  arma::cube posMu(J, C, numSamples, arma::fill::zeros); // posterior samples of mu
  arma::cube posSigma2(J, C, numSamples, arma::fill::zeros); // posterior samples of sigma2
  arma::cube posSigma(J, J, numSamples, arma::fill::zeros); // posterior samples of Sigma
  arma::cube posPi(C, R, numSamples, arma::fill::zeros); // posterior samples of pi
  Rcpp::NumericMatrix posZ(N, numSamples); // posterior samples of spatial domains
  Rcpp::NumericMatrix posC(N, numSamples); // posterior samples of cell types 
  Rcpp::NumericVector burninBeta; // burn-in samples of beta for checking convergence
  arma::mat posLambda(J, numSamples); // posterior samples of lambda
  arma::mat posD(J, numSamples); // posterior samples of d

  // parameter initialization
  initParams(mu, z, c, lambda, d, R2, initMethod, X, R);
  // Note: assignment by "=" is copy by reference, not value 
  // initC = c;
  initC.assign(c.begin(), c.end());
  initZ.assign(z.begin(), z.end());
  S = updateS(z);
  T = updateT(c);
  U = updateU(S, T);
  // knn graph for each sample
  for(int l = 0; l < L; l++)
  {
    Rcpp::NumericMatrix xyl = xy(Range(smpIdx(l), smpIdx(l+1)-1), _);
    Vs.push_back(constructKNN(k, xyl));
  }

  // 1.handling beta
  if(betaEstApproach == "FIXED")
  {
    cout << "Fixing beta to be: " << beta << endl;
  }
  else if(betaEstApproach == "ACCUR_EST" || betaEstApproach == "FAST_EST")
  {
    if(L > 1)
    {
      cout << "Estimating beta in the first tissue section" << endl;
    }
    cout << "Estimating beta..." << endl;
    int iter = 0;
    while(!isConverged)
    {
      updateAll(X, L, alpha0, covStruc, a, b, W0inv, n0, kappa, beta, R2, smpIdx,
        Vs, zl, cl, pi, sigma2, Sigmainv, mu, c, z, lambda, d, T, S, U);

      if(betaEstApproach == "ACCUR_EST")
      {
        acceptBetaP += updateBeta(z[Range(0, smpIdx(1)-1)], Vs[0], R, epsilon, 
          M, B, betaMax, beta);
      }

      if(iter >= 200 && iter % 100 == 0)
      {
        isConverged = checkBetaConvergence(burninBeta, iter, betaTol);
      }

      burninBeta.push_back(beta);
      iter += 1;
    }

    beta = Rcpp::mean(burninBeta[Rcpp::seq(iter - 100, iter - 1)]);
    cout << "Parameter beta has converged: fix beta to be " << beta << endl;
  }


  // 2.warming up
  for(int iter = 0; iter < warmUp; iter++)
  {
    if((iter + 1) % 100 == 0)
    {
      cout << "\rWarming up iteration: " << iter + 1 << "/" << warmUp << flush;
    }
    // update all parameters
    updateAll(X, L, alpha0, covStruc, a, b, W0inv, n0, kappa, beta, R2, smpIdx,
      Vs, zl, cl, pi, sigma2, Sigmainv, mu, c, z, lambda, d, T, S, U);
  }
  cout << endl;


  // 3.sampling
  for(int iter = 0; iter < numSamples; iter++)
  {
    if((iter + 1) % 100 == 0)
    {
      cout << "\rSampling iteration: " << iter + 1 << "/" << numSamples << flush;
    }
    // update all parameters
    updateAll(X, L, alpha0, covStruc, a, b, W0inv, n0, kappa, beta, R2, smpIdx,
      Vs, zl, cl, pi, sigma2, Sigmainv, mu, c, z, lambda, d, T, S, U);

    // posteriors
    posMu.slice(iter) = mu;
    if(covStruc == "EII")
    {
      posSigma2.slice(iter) = Rcpp::as<arma::mat>(sigma2);
    }
    else if(covStruc == "EEE")
    {
      posSigma.slice(iter) = arma::inv(Sigmainv);
    }
    posPi.slice(iter) = Rcpp::as<arma::mat>(pi);
    posZ(_, iter) = z;
    posC(_, iter) = c;
    posLambda.col(iter) = lambda;
    posD.col(iter) = d;
  }
  cout << endl << "done" << endl;

  Rcpp::List res = Rcpp::List::create(
    _["mu"] = posMu, 
    _["sigma2"] = posSigma2, 
    _["Sigma"] = posSigma,
    _["pi"] = posPi, 
    _["z"] = posZ, 
    _["init_z"] = initZ,
    _["c"] = posC,
    _["init_c"] = initC,
    _["burninBeta"] = burninBeta,
    _["beta"] = beta,
    _["acceptBetaP"] = acceptBetaP / burninBeta.length(),
    _["lambda"] = posLambda,
    _["d"] = posD,
    _["R2"] = R2
    );
  return res;
}
