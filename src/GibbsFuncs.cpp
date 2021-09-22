// Author: Zheng Li
// Date: 2020-08-13
// Functions used in Gibbs sampler
#include <RcppArmadillo.h>
#include "compPartitionRatio.h"
#include "GibbsFuncs.h"
#include <string>
using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

void initParams(
  Rcpp::NumericMatrix &mu, // J x C
  Rcpp::IntegerVector &z, // N x 1
  Rcpp::IntegerVector &c, // N x 1
  double &beta,
  Rcpp::NumericVector &lambda, // J x 1
  Rcpp::NumericVector &d, // J x 1
  Rcpp::NumericVector &R2, // J x 1
  Rcpp::String initMethod,
  const Rcpp::NumericMatrix &X, // N x J 
  const int R
  )
{
  int J = mu.nrow();
  int C = mu.ncol();
  Rcpp::List kmResC;
  Rcpp::List kmResR;
  Rcpp::List mclustResC;
  Rcpp::List mclustResR;
  Rcpp::List mclustResParam;

  if(initMethod == "kmeans")
  {
    Rcpp::Function kmeans("kmeans");
    // 1. Initialize c
    kmResC = kmeans(_["x"] = X, _["centers"] = C, _["iter.max"] = 100, _["nstart"] = 1000);
    c = kmResC["cluster"];
    // 2. Initialize mu
    mu = Rcpp::transpose(Rcpp::as<NumericMatrix>(kmResC["centers"]));
    // 3. Initialize z
    kmResR = kmeans(_["x"] = X, _["centers"] = R, _["iter.max"] = 100, _["nstart"] = 1000);
    z = kmResR["cluster"];
  } 
  else if(initMethod == "mclust")
  {
    Rcpp::Function Mclust("Mclust");
    // 1. Initialize c
    mclustResC = Mclust(_["data"] = X, _["G"] = C, _["modelNames"] = "EEE", 
      _["verbose"] = false);
    c = mclustResC["classification"];
    // 2. Initialize mu
    mclustResParam = mclustResC["parameters"];
    mu = Rcpp::as<NumericMatrix>(mclustResParam["mean"]);
    // 3. Initialize z
    mclustResR = Mclust(_["data"] = X, _["G"] = R, _["modelNames"] = "EEE", 
      _["verbose"] = false);
    z = mclustResR["classification"];
  }
  c = c - 1;
  z = z - 1;

  // 4. Initialize beta
  beta = 1.0;

  // 5. Initialize lambda to 1 (no shrinkage)
  lambda.fill(1.0);

  // 6. Initialize d to 0
  d.fill(0.0);

  // 7. Compute R2: squared range of each feature
  for(int j = 0; j < J; j++)
  {
    Rcpp::NumericVector range = Rcpp::range(X(_, j));
    R2(j) = std::pow(range(1) - range(0), 2.0);
  }
}


std::map<int, std::set<int>> updateT(
  const Rcpp::IntegerVector &c
  )
{
  int N = c.length();
  std::map<int, std::set<int>> T;
  for(int n = 0; n < N; n++)
  {
    T[c(n)].insert(n);
  }
  return T;
}


std::map<int, std::set<int>> updateS(
 const Rcpp::IntegerVector &z
  )
{
  int N = z.length();
  std::map<int, set<int>> S;
  for(int n = 0; n < N; n++)
  {
    S[z(n)].insert(n);
  }
  return S;
}


std::map<std::string, std::set<int>> updateU(
  const std::map<int, std::set<int>> &S, 
  const std::map<int, std::set<int>> &T
  )
{
  string key;
  std::map<std::string, std::set<int>> U;
  // traverse S and T to construct U
  for(auto s: S)
  {
    for(auto t: T)
    {
      std::set<int> u; // store each combination
      key = to_string(s.first) + to_string(t.first);
      // intersect sets S[k] and T[c]
      std::set_intersection(
      s.second.begin(), s.second.end(),
      t.second.begin(), t.second.end(),
      std::inserter(u, u.begin())
      );
      U[key] = u;
    }
  }
  return U;
}


std::map<int, std::set<int>> getNeighbors(
  const double r,
  const Rcpp::NumericMatrix &coord
  )
{
  int N = coord.nrow();
  std::map<int, set<int>> V;
  double dist; // distance between cell n and m
  for(int n = 0; n < N; n++)
  {
    for(int m = 0; m < N; m++)
    {
      dist = std::sqrt(Rcpp::sum(Rcpp::pow(coord(n, _) - coord(m, _), 2.0)));
      if(m != n && dist <= r)
      {
        V[n].insert(m);
      }
    }
  }
  return V;
}


std::map<int, std::set<int>> constructKNN(
  const int K,
  const Rcpp::NumericMatrix &coord 
  )
{
  int N = coord.nrow();
  if(N < K)
  {
    cout << "ERROR: Number of cells fewer than number of neighbours K" << endl;
    exit(1);
  }
  std::map<int, set<int>> V;
  double dist;
  std::vector<std::pair<double, int>> distAll;

  for(int n = 0; n < N; n++)
  {
    for(int m = 0; m < N; m++)
    {
      dist = std::sqrt(Rcpp::sum(Rcpp::pow(coord(n, _) - coord(m, _), 2.0)));
      distAll.push_back(std::make_pair(dist, m));
    }
    std::sort(distAll.begin(), distAll.end()); // sorting based on first value of pair objects
    for(int k = 1; k < (K + 1); k++) // skip itself
    {
      V[n].insert(distAll[k].second);
    }
    distAll.clear();
  }
  return V;
}


// [[Rcpp::export]]
double avgNeighbors(
  const double r, 
  const Rcpp::NumericMatrix &coord
  )
{
  int N = coord.nrow();
  double avgNeighbors = 0;
  std::map<int, std::set<int>> V;
  
  V = getNeighbors(r, coord);
  for(auto v: V)
  {
    avgNeighbors += v.second.size();
  }
  avgNeighbors = avgNeighbors / N;
  return avgNeighbors;
}


void updateMu(
  const Rcpp::NumericMatrix &X,
  std::map<int, std::set<int>> &T,
  const Rcpp::NumericMatrix &sigma2,
  const double kappa,
  Rcpp::NumericMatrix &mu
  )
{
  int J = mu.nrow();
  int C = mu.ncol();
  double sumXij;

  for(int j = 0; j < J; j++)
  {
    for(int c = 0; c < C; c++)
    {
      sumXij = 0;
      for(auto i: T[c])
      {
        sumXij += X(i, j);
      }
      mu(j ,c) = R::rnorm(
        sumXij / (T[c].size() + kappa),
        std::sqrt(sigma2(j, c) / (T[c].size() + kappa))
        );
    }
  }
}


void updateMuNG(
  const Rcpp::NumericMatrix &X,
  std::map<int, std::set<int>> &T,
  const Rcpp::NumericMatrix &sigma2,
  const Rcpp::NumericVector &d,
  const Rcpp::NumericVector &lambda,
  const Rcpp::NumericVector &R2,
  Rcpp::NumericMatrix &mu
  )
{
  int J = mu.nrow();
  int C = mu.ncol();
  double sumXij;

  for(int j = 0; j < J; j++)
  {
    for(int c = 0; c < C; c++)
    {
      sumXij = 0;
      for(auto i: T[c])
      {
        sumXij += X(i, j);
      }
      mu(j, c) = R::rnorm(
        (lambda(j) * R2(j) * sumXij + d(j) * sigma2(j, c)) / 
          (T[c].size() * lambda(j) * R2(j) + sigma2(j, c)),
        std::sqrt((lambda(j) * R2(j) * sigma2(j, c)) / 
          (T[c].size() * lambda(j) * R2(j) + sigma2(j, c)))
        );
    }
  } 
}


void updateLambda(
  double v1,
  double v2,
  const Rcpp::NumericMatrix &mu,
  const Rcpp::NumericVector &d,
  const Rcpp::NumericVector &R2,
  Rcpp::NumericVector &lambda
  )
{
  Rcpp::Function rgig("rgig");
  int J = mu.nrow();
  int C = mu.ncol();
  double aj = 2.0 * v2;
  double bj;
  double pC = v1 - C / 2.0;

  for(int j = 0; j < J; j++)
  {
    bj = 0;
    for(int c = 0; c < C; c++)
    {
      bj += std::pow((mu(j, c) - d(j)), 2.0) / R2(j);
    }
    SEXP gig = rgig(_["n"] = 1, _["lambda"] = pC, _["chi"] = bj, _["psi"] = aj);
    lambda(j) = *REAL(gig);
  }
}


void updateD(
  const Rcpp::NumericMatrix &mu,
  const Rcpp::NumericVector &lambda,
  const Rcpp::NumericVector &R2,
  Rcpp::NumericVector &d
  )
{
  int J = mu.nrow();
  int C = mu.ncol();
  double sumMujc;

  for(int j = 0; j < J; j++)
  {
    sumMujc = 0;
    for(int c = 0; c < C; c++)
    {
      sumMujc += mu(j, c);
    }
    d(j) = R::rnorm(sumMujc / C, std::sqrt(lambda(j) * R2(j) / C));
  }
}


void updateSigma2(
  const Rcpp::NumericMatrix &X, 
  std::map<int, std::set<int>> &T, 
  const Rcpp::NumericMatrix &mu, 
  const double a, 
  const double b, 
  Rcpp::NumericMatrix &sigma2
  )
{
  int J = sigma2.nrow();
  int C = sigma2.ncol();
  double sumXMinusMu2;

  for(int j = 0; j < J; j++)
  {
    for(int c = 0; c < C; c++)
    {
      sumXMinusMu2 = 0;
      for(auto i: T[c])
      {
        sumXMinusMu2 += std::pow(X(i, j) - mu(j, c), 2.0);
      }
      sigma2(j, c) = R::rgamma(a + T[c].size() / 2.0, 1.0 / (b + sumXMinusMu2 / 2.0));
      sigma2(j, c) = 1.0 / sigma2(j, c);
    }
  }
}


void updateSigma2EII(
  const Rcpp::NumericMatrix &X,
  std::map<int, std::set<int>> &T,
  const Rcpp::NumericMatrix &mu,
  const double a,
  const double b,
  Rcpp::NumericMatrix &sigma2
  )
{
  int N = X.nrow();
  int J = sigma2.nrow();
  int C = sigma2.ncol();
  double sigma2EII;
  double sumXMinusMu2 = 0;

  for(int j = 0; j < J; j++)
  {
    for(int c = 0; c < C; c++)
    {
      for(auto i: T[c])
      {
        sumXMinusMu2 += std::pow(X(i, j) - mu(j, c), 2.0);
      }
    }
  }
  sigma2EII = R::rgamma(a + N * J / 2.0, 1.0 / (b + sumXMinusMu2 / 2.0));
  sigma2EII = 1.0 / sigma2EII;
  std::fill(sigma2.begin(), sigma2.end(), sigma2EII);
}


void updatePi(
  std::map<std::string, std::set<int>> &U, 
  const double alpha0,
  Rcpp::NumericMatrix &pi
  )
{
  int C = pi.nrow();
  int R = pi.ncol();
  std::string key;
  Rcpp::NumericVector alphasNew(C);

  for(int r = 0; r < R; r++)
  {
    for(int c = 0; c < C; c++)
    {
      key = to_string(r) + to_string(c);
      alphasNew(c) = alpha0 + U[key].size();
    }
    pi(_, r) = rdirichlet(1, alphasNew)(0, _);
  }
}


void updateC(
  const Rcpp::NumericMatrix &X,
  const Rcpp::NumericMatrix &mu,
  const Rcpp::NumericMatrix &sigma2,
  const Rcpp::NumericMatrix &pi,
  const Rcpp::IntegerVector &z,
  Rcpp::IntegerVector &c
  )
{
  int N = X.nrow();
  int C = pi.nrow();
  double sumXMinusMu2;
  double logSigma2;
  Rcpp::IntegerVector cLabels = seq(0, C - 1);
  Rcpp::NumericVector probCK(C); // posterior probability of c_i = c given z_i = k
  
  for(int n = 0; n < N; n++)
  {
    for(int c = 0; c < C; c++)
    {
      sumXMinusMu2 = -Rcpp::sum(Rcpp::pow(X(n, _) - mu(_, c), 2.0) / (2.0 * sigma2(_, c)));
      logSigma2 = -0.5 * Rcpp::sum(Rcpp::log(sigma2(_, c)));
      probCK(c) = std::log(pi(c, z(n))) + logSigma2 + sumXMinusMu2;
    }
    probCK = probCK - Rcpp::max(probCK);
    probCK = Rcpp::exp(probCK);
    c(n) = Rcpp::sample(cLabels, 1, true, probCK)[0];
  }
}


void updateZ(
  std::map<int, std::set<int>> &V, 
  const Rcpp::NumericMatrix &pi,
  const Rcpp::IntegerVector &c,
  const double beta,
  Rcpp::IntegerVector &z
  )
{
  int N = z.length();
  int R = pi.ncol();
  Rcpp::NumericVector Nzi(R); // number of neighbours with the same label
  Rcpp::IntegerVector kLabels = seq(0, R - 1);
  Rcpp::NumericVector probKC(R); // posterior probability of z_i = z given c_i = c

  for(int n = 0; n < N; n++)
  {
    Nzi.fill(0.0);
    for(auto i: V[n])
    {
      Nzi(z(i)) += 1;
    }
    probKC = pi(c(n), _) * Rcpp::exp(beta * Nzi);
    z(n) = Rcpp::sample(kLabels, 1, true, probKC)[0];
  }
}


void updateZSW(
  std::map<int, std::set<int>> &V,
  const Rcpp::NumericMatrix &pi,
  const Rcpp::IntegerVector &c,
  const double beta,
  Rcpp::IntegerVector &z
  )
{
  int N = z.length();
  int R = pi.ncol();
  Rcpp::IntegerVector patches(N, -1);
  Rcpp::LogicalVector isBond(N, false);
  std::map<int, std::set<int>> patchSet;
  Rcpp::IntegerVector labels = Rcpp::seq(0, R - 1);
  Rcpp::IntegerVector newLabels;
  Rcpp::NumericVector probR(R); // probability of being each tissue structure 
  int patch = 0;

  // sample bond variable and form patches
  for(int n = 0; n < N; n++)
  {
    if(!isBond(n))
    {
      patches(n) = patch;
      isBond(n) = true;
    }

    for(auto i: V[n])
    {
      if(z(i) == z(n) && R::runif(0, 1) <= (1.0 - std::exp(-0.5 * beta)))
      {
        if(!isBond(i))
        {
          patches(i) = patches(n);
          isBond(i) = true;
        }
        else if(patches(i) != patches(n))
        {
          mergePatches(patches, patches(i), patches(n));
        }
      }
    }

    if(patch == patches(n))
    {
      patch += 1;
    }
  }

  // create patch set for easy access
  for(int n = 0; n < N; n++)
  {
    patchSet[patches(n)].insert(n);
  }

  // update z
  for(int p = 0; p < patch; p++)
  {
    probR.fill(0.0);
    for(auto i: patchSet[p])
    {
      probR += Rcpp::log(pi(c(i), _));
    }
    probR = Rcpp::exp(probR - Rcpp::max(probR));
    newLabels.push_back(Rcpp::sample(labels, 1, true, probR)[0]);
  }
  for(int n = 0; n < N; n++)
  {
    z(n) = newLabels(patches(n));
  }
}


bool updateBeta(
  const Rcpp::IntegerVector &z, 
  const std::map<int, std::set<int>> &V, 
  const int R, 
  const double epsilon,
  const int M,
  const int B,
  const double betaMax,
  double &beta
  )
{
  int N = z.length();
  double betaProposed;
  double partitionRatio;
  double A; // Acceptance probability
  bool doAccept = false;
  
  betaProposed = std::abs(R::runif(-epsilon, epsilon) + beta);
  betaProposed = std::min(betaProposed, betaMax);
  partitionRatio = compPartitionRatio(betaProposed, beta, V, R, N, M, B);
  A = std::exp(partitionRatio + (betaProposed - beta) * compNHC(V, z));
  A = std::min(1.0, A);

  if(R::runif(0, 1) <= A)
  {
    beta = betaProposed;
    doAccept = true;
  }
  return doAccept;
}


bool updateBetaFast(
  const Rcpp::IntegerVector &z,
  const std::map<int, std::set<int>> &V, 
  const double epsilon,
  const Rcpp::NumericMatrix &NHC,
  double &beta
  )
{
  int M = NHC.ncol(); // number of Potts samples
  Rcpp::NumericVector exps(M); // exponents
  double betaProposed;
  double offset;
  double partitionRatio;
  double A; // Acceptance probability
  bool doAccept = false;

  betaProposed = std::min(4.0, std::abs(R::runif(-epsilon, epsilon) + beta));
  // round betaProposed to the first decimal point
  betaProposed = std::round(betaProposed * 10.0) / 10.0;

  exps = (beta - betaProposed) * NHC(int(betaProposed * 10), _);
  offset = Rcpp::max(exps);
  partitionRatio = -std::log(M) + offset + std::log(Rcpp::sum(Rcpp::exp(exps - offset)));

  A = std::exp(partitionRatio + (betaProposed - beta) * compNHC(V, z));
  A = std::min(1.0, A);

  if(R::runif(0, 1) <= A)
  {
    beta = betaProposed;
    doAccept = true;
  }
  return doAccept;
}


bool checkBetaConvergence(
  const Rcpp::NumericVector &burninBeta,
  int iter,
  double tol
  )
{
  Rcpp::NumericVector betas_prev;
  Rcpp::NumericVector betas_next;
  double beta_prev;
  double beta_next;
  bool isConverged = false;
  
  betas_prev = burninBeta[Rcpp::seq(iter - 200, iter - 101)];
  betas_next = burninBeta[Rcpp::seq(iter - 100, iter - 1)];
  beta_prev = Rcpp::mean(betas_prev);
  beta_next = Rcpp::mean(betas_next);
  
  if(std::abs(beta_prev - beta_next) < tol)
  {
    isConverged = true;
  }
  return isConverged;
}


double evalLogSL(
  const double beta, 
  const Rcpp::IntegerVector &z, 
  std::map<int, std::set<int>> &V, 
  const int K
  )
{
  int N = z.length();
  double logSL = 0.0;
  Rcpp::NumericVector Nzi(K);

  for(int n = 0; n < N; n++)
  {
    Nzi.fill(0.0);
    for(auto i: V[n])
    {
      Nzi(z(i)) += 1;
    }
    logSL += std::log(std::exp(beta * Nzi(z(n))) / Rcpp::sum(Rcpp::exp(beta * Nzi)));
  }
  return logSL;
}


Rcpp::NumericMatrix rdirichlet(
  const int numSmps, 
  const Rcpp::NumericVector alphas
  )
{
    int distSize = alphas.length();
    double normalizeTerm;
    // each row will be a draw from a Dirichlet
    Rcpp::NumericMatrix smps (numSmps, distSize);

    for(int i = 0; i < numSmps; ++i)
    {
        // loop through the smps and draw Gamma variables
        for(int j = 0; j < distSize; ++j)
        {
            smps(i, j) = R::rgamma(alphas[j], 1.0);
        }
        // normalize
        normalizeTerm = Rcpp::sum(smps(i, _));
        for (int j = 0; j < distSize; ++j)
        {
            smps(i, j) = smps(i, j) / normalizeTerm;
        }
    }
    return smps;
}


// [[Rcpp::export]]
Rcpp::List testFastBetaEst(
  const Rcpp::IntegerVector x,
  const Rcpp::NumericMatrix coord,
  const double r,
  double beta,
  const double epsilon,
  const int S,
  const Rcpp::NumericMatrix NHC 
  )
{
  int N = x.length();
  map<int, set<int>> V;
  double dist;
  double nAccept = 0;
  Rcpp::NumericVector posBeta (S);

  for(int n = 0; n < N; n++)
  {
    for(int m = 0; m < N; m++)
    {
      dist = std::sqrt(Rcpp::sum(Rcpp::pow(coord(n, _) - coord(m, _), 2)));
      if(m != n && dist <= r)
      {
        V[n].insert(m);
      }
    }
  }

  for(int s = 0; s < S; s++)
  {
    nAccept += updateBetaFast(x, V, epsilon, NHC, beta);
    posBeta(s) = beta;
  }
  
  Rcpp::List res = Rcpp::List::create(
    _["posBetas"] = posBeta,
    _["acceptProb"] = nAccept / S
    );
  return res;
}


void printVs(std::vector<std::map<int, std::set<int>>> Vs)
{
  int L = Vs.size();
  int Nl;
  for(int l = 0; l < L; l++)
  {
    cout << "***** Tissue slice " << l << " *****" << endl;
    Nl = Vs[l].size();
    for(int n = 0; n < Nl; n++)
    {
      cout << "  " << n << ": "; 
      for(auto i: Vs[l][n])
      {
        cout << i;
      }
      cout << endl;
    }
  }
}


void updateAll(
  const Rcpp::NumericMatrix &X,
  const int L,
  const double alpha0,
  const double a,
  const double b,
  const double kappa,
  const double beta,
  const Rcpp::NumericVector &R2,
  const Rcpp::IntegerVector &smpIdx,
  std::vector<std::map<int, std::set<int>>> &Vs,
  Rcpp::IntegerVector &zl,
  Rcpp::IntegerVector &cl,
  Rcpp::NumericMatrix &pi,
  Rcpp::NumericMatrix &sigma2,
  Rcpp::NumericMatrix &mu,
  Rcpp::IntegerVector &c,
  Rcpp::IntegerVector &z,
  Rcpp::NumericVector &lambda,
  Rcpp::NumericVector &d,
  std::map<int, std::set<int>> &T,
  std::map<int, std::set<int>> &S,
  std::map<std::string, std::set<int>> &U
  )
{
  updatePi(U, alpha0, pi);
  updateSigma2EII(X, T, mu, a, b, sigma2);
  updateMuNG(X, T, sigma2, d, lambda, R2, mu);
  updateC(X, mu, sigma2, pi, z, c);
  for(int l = 0; l < L; l++)
  {
    zl =  z[Range(smpIdx(l), smpIdx(l+1)-1)];
    cl = c[Range(smpIdx(l), smpIdx(l+1)-1)];
    updateZSW(Vs[l], pi, cl, beta, zl);
    z[Range(smpIdx(l), smpIdx(l+1)-1)] = zl;
  }
  updateLambda(0.5, 0.5, mu, d, R2, lambda);
  updateD(mu, lambda, R2, d);
  T = updateT(c);
  S = updateS(z);
  U = updateU(S, T);
}