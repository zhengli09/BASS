// Author: Zheng Li
// Date: 2022-04-14
// Purpose:
// BASS implementation

#include "Potts.h"
#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::plugins(cpp11)]]

vec rdirichlet(const vec alpha)
{
  int dist_dim = alpha.n_elem;
  double normalize_const;
  vec smp(dist_dim);
  for(int j = 0; j < dist_dim; j++)
  {
    smp(j) = randg(distr_param(alpha(j), 1.0));
  }
  normalize_const = accu(smp);
  smp /= normalize_const;
  return smp;
}


class BASSModel
{
  private:
    struct STData
    {
      int N; // total number of cells/spots
      int J; // number of expression features
      int L; // number of tissue sections
      vec ns; // sample size in each tissue section
      mat X; // JxN data matrix
      mat xy; // Nx2 spatial coordinates
      vec R2; // squared range of each feature
      vec section_idx; // begin and end index of each tissue section
    } dat;
    
    struct HyperParas
    {
      int C; // number of cell type clusters
      int R; // number of spatial domains
      double alpha0; // concentration parameter of the Dirichlet prior
      mat Psi0; // scale matrix of the inv-Wishart prior
      int n0; // df of the inverse-Wishart prior
      double nu1; // prior parameter of the normal-gamma prior
      double nu2; // prior parameter of the normal-gamma prior
      double beta_min; // lower bound of beta
      double beta_max; // upper bound of beta
      int k; // minimum number of neighbors of each cell/spot
    } hyper_paras;
    
    struct BASSParas
    {
      uvec c; // cell type labels
      uvec z; // spatial domain labels
      mat pi; // CxR cell type proportions in each spatial domain
      mat mu; // JxC mean parameters of expression
      mat Sigma; // JxJ covariance matrix of expression
      mat Sigma_inv; // inverse of Sigma
      vec lambda; // Jx1 vector of local shrinkage parameters for each feature
      vec d; // Jx1 prior mean parameters of mu across cell type clusters
      double loglik;
    } paras;
    
    struct PottsParas
    {
      PottsModel potts_model;
      double kernel; // sum_{i}sum_{i'~i}I(z_i'==z_i)
      double beta; // interaction parameter of the Potts model
      double beta_est;
      bool beta_converged;
      field<sp_mat> V; // neighborhood graph
      // SW algorithm relevant variables
      field<sp_mat> u; // bond variables
      field<vec> visited;
      field<vec> groups;
      int ngroup;
    } potts_paras;
    
    struct MCMCSample
    {
      umat c;
      umat z;
      cube pi;
      vec beta;
      cube mu;
      cube Sigma;
      mat lambda;
      mat d;
    } paras_sample;
    
    struct GibbsSamplerProfile
    {
      vec loglik;
    } gibbs_profile;
    
    struct GibbsSamplerControl
    {
      int total_iter;
      int burnin;
      int nsample;
      String method; // fix/SW for estimating beta
      int potts_burnin;
      int potts_nsample;
      double step_size;
      double tol;
    } gibbs_control;
    
    int iter;
    
  public:
    void load_data(const mat &X, const mat &xy, const vec &ns)
    {
      dat.X = X;
      dat.xy = xy;
      dat.ns = ns;
      dat.N = dat.X.n_cols;
      dat.J = dat.X.n_rows;
      dat.L = dat.ns.n_elem;
      dat.R2 = max(dat.X, 1) - min(dat.X, 1);
      dat.R2 = dat.R2 % dat.R2;
      dat.section_idx = vec(dat.L + 1, fill::zeros);
      dat.section_idx.subvec(1, dat.L) = cumsum(dat.ns);
    }
    
    void set_hyper_paras(const int C, const int R, const int k,
      const double alpha0, const mat Psi0, const int n0, 
      const double nu1 = 0.5, const double nu2 = 0.5, 
      const double beta_min = 0, const double beta_max = 4)
    {
      hyper_paras.C = C;
      hyper_paras.R = R;
      hyper_paras.alpha0 = alpha0;
      hyper_paras.Psi0 = Psi0;
      hyper_paras.n0 = n0;
      hyper_paras.nu1 = nu1;
      hyper_paras.nu2 = nu2;
      hyper_paras.beta_min = beta_min;
      hyper_paras.beta_max = beta_max;
      hyper_paras.k = k;
    }
    
    void set_gibbs_control(int burnin, int nsample, String method,
      int potts_burnin, int potts_nsample, double step_size, double tol)
    {
      gibbs_control.burnin = burnin;
      gibbs_control.nsample = nsample;
      gibbs_control.total_iter = gibbs_control.burnin + 
        gibbs_control.nsample;
      gibbs_control.method = method;
      gibbs_control.potts_burnin = potts_burnin;
      gibbs_control.potts_nsample = potts_nsample;
      gibbs_control.step_size = step_size;
      gibbs_control.tol = tol;
    }
    
    void construct_graph()
    {
      for(int l = 0; l < dat.L; l++)
      {
        if(dat.ns(l) < hyper_paras.k)
        {
          cout << "ERROR: Number of cells/spots is smaller than the ";
          cout << "specified minimum number of neighbors (k)" << endl;
        }
        sp_mat graph(dat.ns(l), dat.ns(l));
        mat xy_l = dat.xy.rows(dat.section_idx(l), dat.section_idx(l+1)-1);
        vec one = ones<vec>(dat.ns(l));
        for(int n = 0; n < dat.ns(l); n++)
        {
          mat distance = xy_l - kron(xy_l.row(n), one);
          distance = sum(distance % distance , 1);
          uvec neighbor = sort_index(distance);
          for(int m = 1; m <= hyper_paras.k; m++)
          {
            graph(n, neighbor(m)) = 1;
            graph(neighbor(m), n) = 1;
          }
        }
        potts_paras.V(l) = graph;
      }
    }
    
    void initialize_paras(const uvec &init_c, const uvec &init_z,
      const mat &init_mu, const double init_beta)
    {
      paras.c = init_c;
      paras.z = init_z;
      paras.mu = init_mu;
      paras.lambda.ones(dat.J);
      paras.d.zeros(dat.J);
      paras.pi.zeros(hyper_paras.C, hyper_paras.R);
      paras.Sigma.ones(dat.J, dat.J);
      update_pi();
      update_Sigma();
      // Potts model relevant variables
      potts_paras.beta = init_beta;
      potts_paras.beta_est = init_beta;
      potts_paras.beta_converged = false;
      potts_paras.V = field<sp_mat>(dat.L);
      construct_graph();
      potts_paras.u = field<sp_mat>(dat.L);
      potts_paras.visited = field<vec>(dat.L);
      potts_paras.groups = field<vec>(dat.L);
      for(int l = 0; l < dat.L; l++)
      {
        potts_paras.u(l) = sp_mat(dat.ns(l), dat.ns(l));
        potts_paras.visited(l) = vec(dat.ns(l), fill::zeros);
        potts_paras.groups(l) = vec(dat.ns(l), fill::zeros);
      }
      potts_paras.ngroup = 0;
      if(gibbs_control.method == "SW")
      {
        potts_paras.potts_model.initialize_model(
          potts_paras.V(0), hyper_paras.R);
      }
      // profile
      gibbs_profile.loglik.zeros(gibbs_control.total_iter);
    }
    
    // sample from a multivariate normal distribution in
    // the form of N(Omega_inv * mu, Omega_inv)
    void update_mu()
    {
      for(int c = 0; c < hyper_paras.C; c++)
      {
        uvec idx = find(paras.c == c);
        mat Omega = idx.n_elem * paras.Sigma_inv;
        Omega.diag() += 1.0 / paras.lambda / dat.R2;
        vec mu = paras.Sigma_inv * sum(dat.X.cols(idx), 1) + 
          paras.d / paras.lambda / dat.R2;
        
        // sampling from MVN using Cholesky decomposition
        mat R = chol(Omega);
        vec temp = solve(R.t(), mu);
        vec z = randn<vec>(dat.J);
        paras.mu.col(c) = solve(R, z + temp);
      }
    }
  
    void update_pi()
    {
      for(int r = 0; r < hyper_paras.R; r++)
      {
        uvec idx_r = find(paras.z == r);
        uvec c_r = paras.c.elem(idx_r);
        vec alpha(hyper_paras.C);
        for(int c = 0; c < hyper_paras.C; c++)
        {
          uvec idx_cr = find(c_r == c);
          alpha(c) = idx_cr.n_elem + hyper_paras.alpha0;
        }
        // sample from a Dirichlet distribution
        paras.pi.col(r) = rdirichlet(alpha);
      }
    }
    
    void update_Sigma()
    {
      int fc_n0 = dat.N + hyper_paras.n0;
      mat temp = dat.X - paras.mu.cols(paras.c);
      mat fc_Psi0 = temp * temp.t() + hyper_paras.Psi0;
      paras.Sigma = riwish(fc_n0, fc_Psi0);
      paras.Sigma_inv = inv(paras.Sigma);
    }
    
    void update_lambda()
    {
      Function rgig("rgig");
      double a = 2.0 * hyper_paras.nu2;
      double bj;
      double p = hyper_paras.nu1 - hyper_paras.C / 2.0;
      
      for(int j = 0; j < dat.J; j++)
      {
        rowvec temp = paras.mu.row(j) - paras.d(j);
        bj = accu(temp % temp) / dat.R2(j);
        SEXP gig = rgig(_["n"] = 1, _["lambda"] = p, _["chi"] = bj, 
          _["psi"] = a);
        paras.lambda(j) = *REAL(gig);
      }
    }
    
    void update_d()
    {
      vec fc_mu = sum(paras.mu, 1) / hyper_paras.C;
      vec fc_sd = sqrt(paras.lambda % dat.R2 / hyper_paras.C);
      vec z = randn(dat.J);
      paras.d = z % fc_sd + fc_mu;
    }
    
    void update_c()
    {
      mat log_prob(dat.N, hyper_paras.C);
      mat XtSigmaInvmu = dat.X.t() * paras.Sigma_inv * paras.mu;
      vec XtSigmaInvX(dat.N);
      for(int n = 0; n < dat.N; n++)
      {
        XtSigmaInvX(n) = -0.5 * as_scalar(dat.X.col(n).t() * 
          paras.Sigma_inv * dat.X.col(n));
      }
      for(int c = 0; c < hyper_paras.C; c++)
      {
        double mutSigmaInvmu = -0.5 * as_scalar(paras.mu.col(c).t() * 
          paras.Sigma_inv * paras.mu.col(c));
        log_prob.col(c) = XtSigmaInvX + mutSigmaInvmu + XtSigmaInvmu.col(c);
      }
      log_prob += log(paras.pi.cols(paras.z).t());
      
      // convert log_prob to prob
      mat prob(dat.N, hyper_paras.C);
      prob = log_prob - kron(ones<rowvec>(hyper_paras.C), max(log_prob, 1));
      prob = exp(prob);
      prob = normalise(prob, 1, 1);
      
      // draw samples
      vec u = randu(dat.N);
      mat cumprob = cumsum(prob, 1);
      for(int n = 0; n < dat.N; n++)
      {
        paras.c(n) = as_scalar(find(cumprob.row(n) > u(n), 1, "first"));
      }
    }
    
    void update_u(int l)
    {
      potts_paras.u(l).zeros();
      for(int n = 0; n < dat.ns(l); n++)
      {
        for(sp_mat::col_iterator it = potts_paras.V(l).begin_col(n);
          it != potts_paras.V(l).end_col(n); ++it)
        {
          if(it.row() > n && (paras.z(dat.section_idx(l)+n) == 
            paras.z(dat.section_idx(l)+it.row())))
          {
            potts_paras.u(l)(it.row(), n) = (randu<double>() * 
              exp(potts_paras.beta)) > 1; 
          }
        }
      }
      potts_paras.u(l) += potts_paras.u(l).t();
    }
    
    void depth_first_search(int l, int n)
    {
      potts_paras.visited(l)(n) = 1;
      potts_paras.groups(l)(n) = potts_paras.ngroup - 1;
      for(sp_mat::col_iterator it = potts_paras.V(l).begin_col(n);
        it != potts_paras.V(l).end_col(n); ++it)
      {
        if(potts_paras.visited(l)(it.row()) == 0 && 
          potts_paras.u(l)(it.row(), n) == 1)
        {
          depth_first_search(l, it.row());
        }
      }
    }
    
    void update_groups(int l)
    {
      potts_paras.groups(l).zeros();
      potts_paras.visited(l).zeros();
      potts_paras.ngroup = 0;
      for(int n = 0; n < dat.ns(l); n++)
      {
        if(potts_paras.visited(l)(n) == 0)
        {
          potts_paras.ngroup += 1;
          depth_first_search(l, n);
        }
      }
    }
    
    void update_z()
    {
      for(int l = 0; l < dat.L; l++)
      {
        update_u(l);
        update_groups(l);
        for(int k = 0; k < potts_paras.ngroup; k++)
        {
          uvec idx = find(potts_paras.groups(l) == k);
          mat temp = paras.pi.rows(paras.c(dat.section_idx(l)+idx));
          rowvec prob = sum(log(temp), 0);
          prob = exp(prob - max(prob));
          prob = prob / sum(prob);
          rowvec cumprob = cumsum(prob);
          int label = as_scalar(find(cumprob > randu<double>(), 1, "first"));
          paras.z(dat.section_idx(l)+idx).fill(label);
        }
      }
    }
    
    void update_kernel()
    {
      potts_paras.kernel = 0;
      for(int n = 0; n < dat.ns(0); n++)
      {
        for(sp_mat::col_iterator it = potts_paras.V(0).begin_col(n);
          it != potts_paras.V(0).end_col(n); ++it)
        {
          if(it.row() > n && paras.z(n) == paras.z(it.row()))
          {
            potts_paras.kernel += 1;
          }
        }
      }
    }
    
    void update_beta()
    {
      update_kernel();
      double beta_proposed = potts_paras.beta +
        (2 * randu<double>() - 1) * gibbs_control.step_size;
      beta_proposed = max(hyper_paras.beta_min, beta_proposed);
      beta_proposed = min(hyper_paras.beta_max, beta_proposed);
      // generate Potts kernel values based on the proposed beta
      vec kernel_proposed = potts_paras.potts_model.run(beta_proposed,
          gibbs_control.potts_burnin, gibbs_control.potts_nsample);
      kernel_proposed = kernel_proposed * (potts_paras.beta - beta_proposed);
      double offset = max(kernel_proposed);
      double log_partition = -log(gibbs_control.potts_nsample) + offset +
        log(sum(exp(kernel_proposed - offset)));
      double log_accept_ratio = log_partition + 
        (beta_proposed - potts_paras.beta) * potts_paras.kernel;
      if(log(randu<double>()) < log_accept_ratio)
      {
        potts_paras.beta = beta_proposed;
      }
    }
    
    void check_beta_converge()
    {
      int size = 100; // check every 100 iterations
      if(iter > 0 && iter < gibbs_control.burnin && 
        (iter+1) % size == 0 && !potts_paras.beta_converged)
      {
        double beta_avg = mean(paras_sample.beta.subvec(iter-size+1, iter));
        if(abs(beta_avg - potts_paras.beta_est) < gibbs_control.tol)
        {
          potts_paras.beta_converged = true;
          gibbs_control.method = "fix";
          potts_paras.beta = beta_avg;
        }
        potts_paras.beta_est = beta_avg;
      }
      else if(iter == gibbs_control.burnin && !potts_paras.beta_converged)
      {
        cout << "\nbeta didn't converge within " << gibbs_control.burnin;
        cout << " burnin steps, consider increasing the number of\n";
        cout << "burnin after checking the trace plot of beta ";
        cout << "(plot(1:BASS@burnin, BASS@samples$beta))\n";
        gibbs_control.method = "fix";
        potts_paras.beta = potts_paras.beta_est;
      }
    }
    
    void update_loglik()
    {
      paras.loglik = -0.5 * dat.N * log(det(paras.Sigma));
      for(int n = 0; n < dat.N; n++)
      {
        vec temp = dat.X.col(n) - paras.mu.col(paras.c(n));
        paras.loglik -= 0.5 * as_scalar(temp.t() * paras.Sigma_inv * temp);
      }
    }
    
    void initialize_paras_sample()
    {
      paras_sample.c.zeros(dat.N, gibbs_control.nsample);
      paras_sample.z.zeros(dat.N, gibbs_control.nsample);
      paras_sample.pi.zeros(hyper_paras.C, hyper_paras.R, 
        gibbs_control.nsample);
      paras_sample.beta.zeros(gibbs_control.burnin);
      paras_sample.mu.zeros(dat.J, hyper_paras.C, gibbs_control.nsample);
      paras_sample.Sigma.zeros(dat.J, dat.J, gibbs_control.nsample);
      paras_sample.lambda.zeros(dat.J, gibbs_control.nsample);
      paras_sample.d.zeros(dat.J, gibbs_control.nsample);
    }
    
    void save_paras_sample()
    {
      if(iter < gibbs_control.burnin)
      {
        paras_sample.beta(iter) = potts_paras.beta;
      } 
      else
      {
        int mcmc_iter = iter - gibbs_control.burnin;
        paras_sample.c.col(mcmc_iter) = paras.c;
        paras_sample.z.col(mcmc_iter) = paras.z;
        paras_sample.pi.slice(mcmc_iter) = paras.pi;
        paras_sample.mu.slice(mcmc_iter) = paras.mu;
        paras_sample.Sigma.slice(mcmc_iter) = paras.Sigma;
        paras_sample.lambda.col(mcmc_iter) = paras.lambda;
        paras_sample.d.col(mcmc_iter) = paras.d;
      }
    }
    
    void save_gibbs_profile()
    {
      update_loglik();
      gibbs_profile.loglik(iter) = paras.loglik;
    }
    
    void monitor_mcmc()
    {
      if((iter+1) % 100 == 0)
      {
        if(iter < gibbs_control.burnin)
        {
          cout << "\rBurnin iteration: " << iter + 1 << 
            "/" << gibbs_control.burnin << flush;
        }
        else
        {
          cout << "\rSampling iteration: " << 
            iter - gibbs_control.burnin + 1 << 
            "/" << gibbs_control.nsample << flush;
        }
      }
    }
    
    void run()
    {
      initialize_paras_sample();
      for(iter = 0; iter < gibbs_control.total_iter; iter++)
      {
        update_c();
        update_z();
        update_mu();
        update_d();
        update_lambda();
        update_pi();
        update_Sigma();
        if(gibbs_control.method == "SW")
        {
          update_beta();
        }
        save_paras_sample();
        save_gibbs_profile();
        check_beta_converge();
        monitor_mcmc();
      }
      cout << endl << "done" << endl;
    }
    
    List get_output()
    {
      List output = List::create(
        _["c"] = paras_sample.c,
        _["z"] = paras_sample.z,
        _["pi"] = paras_sample.pi,
        _["beta"] = paras_sample.beta,
        _["mu"] = paras_sample.mu,
        _["Sigma"] = paras_sample.Sigma,
        _["lambda"] = paras_sample.lambda,
        _["d"] = paras_sample.d,
        _["beta_est"] = potts_paras.beta_est,
        _["loglik"] = gibbs_profile.loglik
      );
      return output;
    }
};


// [[Rcpp::export]]
List BASSFit(const arma::mat X, const arma::mat xy, const arma::vec ns, 
  const int C, const int R, const int k, const double alpha0, 
  const arma::mat Psi0, const int n0, const int burnin, const int nsample, 
  const String method, const int potts_burnin, const int potts_nsample, 
  const double step_size, const double tol, const arma::uvec init_c, 
  const arma::uvec init_z, const arma::mat init_mu, const double init_beta
)
{
  wall_clock timer;
  timer.tic();
  BASSModel model;
  
  model.load_data(X, xy, ns);
  model.set_hyper_paras(C, R, k, alpha0, Psi0, n0);
  model.set_gibbs_control(burnin, nsample, method,
    potts_burnin, potts_nsample, step_size, tol);
  model.initialize_paras(init_c, init_z, init_mu, init_beta);
  model.run();
  
  double elapsed = timer.toc();
  List output = model.get_output();
  output.push_back(elapsed, "elapsed_time");
  return output;
}

