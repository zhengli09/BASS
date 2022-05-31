// Author: Zheng Li
// Date: 2022-04-16
// Potts model
#ifndef _POTTS_MODEL_H_
#define _POTTS_MODEL_H_
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

class PottsModel
{
  private:
    vec z; // labels
    double beta; // interaction parameter
    sp_mat V; // neighborhood graph
    int R; // number of states
    int N; // sample size
    // SW-algorithm relevant variables
    sp_mat u; // bond variables
    vec groups;
    vec visited;
    int ngroup;
    // profile
    double kernel; // \sum_{i}\sum_{i'~i}I(z_i'==z_i)
    
  public:
    void initialize_model(const sp_mat &V, int R);
    void update_u();
    void depth_first_search(int n);
    void update_groups();
    void compute_kernel();
    vec run(double beta, int burnin, int nsample);
};
#endif
