#include "compPartitionRatio.h"
#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

double compPartitionRatio(
  const double betaProposed,
  const double beta,
  const std::map<int, std::set<int>> &V,
  const int R,
  const int N,
  const int M,
  const int B
  )
{
  Rcpp::IntegerMatrix smps(N, M); // Each column is a sample from Potts
  Rcpp::NumericVector exps(M); // exponents
  double offset;
  double partitionRatio = 0;

  smps = smpPottsSW(betaProposed, V, R, N, M, B);
  for(int m = 0; m < M; m++)
  {
    exps(m) = (beta - betaProposed) * compNHC(V, smps(_, m));
    // partitionRatio += std::exp((beta - betaProposed) * compNHC(V, smps(_, m))); // Numeric problem
  }
  // offset = Rcpp::min(exps);
  offset = Rcpp::max(exps);
  partitionRatio = -std::log(M) + offset + std::log(Rcpp::sum(Rcpp::exp(exps - offset))); 
  return partitionRatio;
}


Rcpp::IntegerMatrix smpPottsSW(
  const double beta,
  const std::map<int, std::set<int>> &V,
  const int R,
  const int N, 
  const int M,
  const int B
  )
{
  Rcpp::IntegerVector z(N);
  Rcpp::IntegerVector patches(N, -1);
  Rcpp::LogicalVector isBond(N, false);
  Rcpp::IntegerMatrix smps(N, M);
  Rcpp::IntegerVector labels = Rcpp::seq(0, R - 1);

  // randomly initialize z
  z = Rcpp::sample(labels, N, true);
  // burn-in
  for(int b = 0; b < B; b++)
  {
    runSW(beta, V, R, z, patches, isBond);
    patches.fill(-1);
    isBond.fill(false);
  }
  // sampling
  for(int m = 0; m < M; m++)
  {
    runSW(beta, V, R, z, patches, isBond);
    smps(_, m) = z;
    patches.fill(-1);
    isBond.fill(false);
  }
  return smps;
}


void runSW(
  const double beta,
  const std::map<int, std::set<int>> &V,
  const int R,
  Rcpp::IntegerVector &z,
  Rcpp::IntegerVector &patches,
  Rcpp::LogicalVector &isBond
  )
{
  int N = z.length();
  int patch = 1;
  Rcpp::IntegerVector labels = Rcpp::seq(0, R - 1);
  Rcpp::IntegerVector newLabels;

  for(int n = 0; n < N; n++)
  {
    if(!isBond(n))
    {
      patches(n) = patch;
      isBond(n) = true;
    }

    if(V.find(n) != V.end())
    {
      for(set<int>::iterator it = V.at(n).begin(); it != V.at(n).end(); it++)
      {
        // Because each edge is visted twice, we adjust the probability of 
        // making a bond by the following.
        if(z(*it) == z(n) && R::runif(0, 1) <= (1.0 - std::exp(-0.5 * beta)))
        {
          if(!isBond(*it))
          {
            patches(*it) = patches(n);
            isBond(*it) = true;
          }
          else if(patches(*it) != patches(n))
          {
            mergePatches(patches, patches(*it), patches(n));
          }
        }
      }
    }

    if(patch == patches(n))
    {
      patch += 1;
    }
  }

  newLabels = Rcpp::sample(labels, patch, true);
  for(int n = 0; n < N; n++)
  {
    z(n) = newLabels(patches(n) - 1);
  }
}


void mergePatches(
  Rcpp::IntegerVector &patches,
  const int patch1,
  const int patch2
  )
{
  int N = patches.length();
  int mergedPatch = std::min(patch1, patch2);
  int droppedPatch = std::max(patch1, patch2);
  for(int n = 0; n < N; n++)
  {
    if(patches(n) == droppedPatch)
    {
      patches(n) = mergedPatch;
    }
  }
}


double compNHC(
  const std::map<int, std::set<int>> &V,
  const Rcpp::IntegerVector &z
  )
{
  int N = z.length();
  double NHC = 0.0;

  for(int n = 0; n < N; n++)
  {
    if(V.find(n) != V.end())
    {
      for(set<int>::iterator it = V.at(n).begin(); it != V.at(n).end(); it++)
      { 
        if(z(*it) == z(n))
        {
          NHC += 1.0;
        }
      }
    }
  }
  NHC = NHC / 2.0;
 return NHC;
}


// [[Rcpp::export]]
Rcpp::List testPottsSampling(
  const Rcpp::NumericMatrix coord,
  const double r,
  const double beta,
  const int R,
  const int M
  )
{
  // construct a graph
  int N = coord.nrow();
  Rcpp::IntegerMatrix smps (N, M);
  Rcpp::NumericVector NHC (M);
  map<int, set<int>> V;
  double dist; // distance between pixel n and m

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
  // Show V
  // for(int n = 0; n < N; n++)
  // {
  //   cout << n << ":{";
  //   for(set<int>::iterator it = & V.at(n).begin(); it != & V.at(n).end(); it++)
  //   {
  //     cout << *it << " ";
  //   }
  //   cout << "}" << endl;
  // }
  
  smps = smpPottsSW(beta, V, R, N, M);
  for(int m = 0; m < M; m++)
  {
    NHC(m) = compNHC(V, smps(_, m));
  }

  Rcpp::List res = Rcpp::List::create(
    _["samples"] = smps,
    _["NHCs"] = NHC
    );
  return res;
}


// [[Rcpp::export]]
Rcpp::List testPostBeta(
  const Rcpp::IntegerVector x,
  const Rcpp::NumericMatrix coord,
  const int R,
  double beta,
  const double epsilon,
  const int M,
  const int S
  )
{
  int N = x.length();
  map<int, set<int>> V;
  double dist;
  std::vector<std::pair<double, int>> distAll;
  double betaProposed;
  double partitionRatio;
  double A;
  Rcpp::NumericVector As (S);
  double nAccept = 0;
  Rcpp::NumericVector posBeta (S);

  for(int n = 0; n < N; n++)
  {
    for(int m = 0; m < N; m++)
    {
      dist = std::sqrt(Rcpp::sum(Rcpp::pow(coord(n, _) - coord(m, _), 2.0)));
      distAll.push_back(std::make_pair(dist, m));
    }
    std::sort(distAll.begin(), distAll.end()); // sorting based on first value of pair objects
    for(int r = 1; r < (R + 1); r++) // skip itself
    {
      V[n].insert(distAll[r].second);
    }
    distAll.clear();
  }

  for(int s = 0; s < S; s++)
  {
    betaProposed = beta + R::runif(-epsilon, epsilon);
    partitionRatio = compPartitionRatio(betaProposed, beta, V, R, N, M);
    A = std::exp(partitionRatio + (betaProposed - beta) * compNHC(V, x));
    As(s) = A;
    A = std::min(1.0, A);
  
    if(R::runif(0, 1) <= A)
    {
      beta = betaProposed;
      nAccept += 1;
    }
    posBeta(s) = beta;
  }

  Rcpp::List res = Rcpp::List::create(
    _["posBetas"] = posBeta,
    _["acceptProb"] = nAccept / S,
    _["As"] = As,
    _["NHC"] = compNHC(V, x)
    );
  return res;
}