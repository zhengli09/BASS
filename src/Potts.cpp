// Author: Zheng Li
// Date: 2022-04-16
// Potts model

#include "Potts.h"
#include <RcppArmadillo.h>
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

void PottsModel::initialize_model(const sp_mat &V, int R)
{
  this->V = V;
  this->R = R;
  N = V.n_cols;
  // initialization
  z = floor(randu<vec>(N) * R);
  u = sp_mat(N, N);
  visited = vec(N);
  groups = vec(N);
}

void PottsModel::update_u()
{
  u.zeros();
  for(int n = 0; n < N; n++)
  {
    for(sp_mat::col_iterator it = V.begin_col(n);
      it != V.end_col(n); ++it)
    {
      if(it.row() > n && z(n) == z(it.row()))
      {
        u(it.row(), n) = (randu<double>() * exp(beta)) > 1;
      }
    }
  }
  u += u.t();
}

void PottsModel::depth_first_search(int n)
{
  visited(n) = 1;
  groups(n) = ngroup - 1;
  for(sp_mat::col_iterator it = V.begin_col(n);
    it != V.end_col(n); ++it)
  {
    if(visited(it.row()) == 0 && u(it.row(), n) == 1)
    {
      depth_first_search(it.row());
    }
  }
}

void PottsModel::update_groups()
{
  groups.zeros();
  visited.zeros();
  ngroup = 0;
  for(int n = 0; n < N; n++)
  {
    if(visited(n) == 0)
    {
      ngroup += 1;
      depth_first_search(n);
    }
  }
}

void PottsModel::compute_kernel()
{
  kernel = 0;
  for(int n = 0; n < N; n++)
  {
    for(sp_mat::col_iterator it = V.begin_col(n);
      it != V.end_col(n); ++it)
    {
      if(it.row() > n && z(it.row()) == z(n))
      {
        kernel += 1;
      }
    }
  }
}

vec PottsModel::run(double beta, int burnin, int nsample)
{
  this->beta = beta;
  int total_iter = burnin + nsample;
  vec kernels(nsample);
  for(int iter = 0; iter < total_iter; iter++)
  {
    update_u();
    update_groups();
    for(int k = 0; k < ngroup; k++)
    {
      uvec idx = find(groups == k);
      z(idx).fill(floor(randu<double>() * R));
    }
    if(iter >= burnin)
    {
      compute_kernel();
      kernels(iter-burnin) = kernel;
    }
  }
  return kernels;
}

