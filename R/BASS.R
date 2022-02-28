# Authors: Zheng Li, Xiang Zhou

#' BASS object stores all the data and parameters used to run the algorithm
#'
#' @slot X L-list of raw gene expression matrix for L samples
#' @slot X_run N x J matrix of processed expression matrix for running BASS
#' @slot xy L-list of location matrix for L samples
#' @slot L Number of samples
#' @slot P Number of genes
#' @slot Ns L-list of number of cells/spots for L samples
#' @slot C Number of cell types
#' @slot R Number of spatial domains
#' @slot init_method Initialization method for both cell types and 
#'  tissue structures
#' @slot cov_struc Covariance structure for gene expression features
#' @slot kappa Prior parameter of gene expression
#' @slot alpha0 Concentration parameter of Dirichlet distribution
#' @slot a Shape parameter of the inverse gamma prior
#' @slot b Scale parameter of the inverse gamma prior
#' @slot W0 W0*I represents the scale matrix of Wishart prior
#' @slot n0 Degrees of freedom of Wishart prior
#' @slot k Number of neighbors for kNN graph
#' @slot burn_in Number of burn-in iterations
#' @slot samples Number of posterior samples
#' @slot beta_est_approach Approach to estimate beta ("FIXED", "ACCUR_EST")
#' @slot beta Fixed value of beta if the betaEstApproach is set to be "FIXED"
#' @slot beta_max upper bound of beta
#' @slot epsilon Uniform random walk step size
#' @slot beta_tol tolerance for checking beta convergence
#' @slot B Number of burn-in steps in Potts sampling
#' @slot M Number of Potts samples to approximate Potts partition ratio
#' @slot res Result list
#' @slot res_postprocess Result list after post-processing
#'
setClass("BASS", slots = list(
  X = "list",
  X_run = "array",
  xy = "list",
  L = "numeric",
  P = "numeric",
  Ns = "numeric",
  C = "numeric",
  R = "numeric",
  init_method = "character",
  cov_struc = "character",
  kappa = "numeric",
  alpha0 = "numeric",
  a = "numeric",
  b = "numeric",
  W0 = "numeric",
  n0 = "numeric",
  k = "numeric",
  burn_in = "numeric",
  samples = "numeric",
  beta_est_approach = "character",
  beta = "numeric",
  beta_max = "numeric",
  epsilon = "numeric",
  beta_tol = "numeric",
  B = "numeric",
  M = "numeric",
  res = "list",
  res_postprocess = "list"
  ))


#' Create BASS object
#' @export
createBASSObject <- function(
  X, xy, C, R, 
  init_method = c("kmeans", "mclust"), 
  cov_struc = c("EEE", "EII"),
  kappa = 0.01, 
  alpha0 = 1, 
  a = 2, 
  b = 1, 
  W0 = 1,
  n0 = 1,
  k = 4, 
  burn_in = 10000, 
  samples = 10000, 
  beta_est_approach = c("FIXED", "ACCUR_EST"),
  beta = 1,
  beta_max = 4,
  epsilon = 0.1,
  beta_tol = 0.1,
  B = 100,
  M = 100
  )
{
  L <- length(X)
  P <- nrow(X[[1]])
  gene_names <- rownames(X[[1]])

  # input checking
  if(length(X) != length(xy))
    stop("Number of samples (L) in expression list and location list", 
      " do not match")
  for(l in 1:L)
  {
    if(ncol(X[[l]]) != nrow(xy[[l]])){
      stop("Number of cells in expression matrix and location matrix",
        " do not match")
    }
    if(!identical(colnames(X[[l]]), rownames(xy[[l]]))){
      stop("Order of cells in expression matrix and location matrix",
        " do not match")
    }
    if(nrow(X[[l]]) != P){
      stop("Number of genes is not consistent across samples.",
        " Please use a common set of genes.")
    }
    if(!identical(gene_names, rownames(X[[l]]))){
      stop("Order of genes in expression matrices do not match")
    }

    if(!is.matrix(X[[l]])){
      X[[l]] <- as.matrix(X[[l]])
      cat("Expression data coerced to a matrix\n")
    }
    if(!is.matrix(xy[[l]])){
      xy[[l]] <- as.matrix(xy[[l]])
      cat("Location data coerced to a matrix\n")
    }
  }

  # create object
  BASS <- new(Class = "BASS", 
    X = X, xy = xy, L = L, P = P, Ns = sapply(X, ncol), C = C, R = R,
    init_method = init_method[1], cov_struc = cov_struc[1], 
    kappa = kappa, alpha0 = alpha0, a = a, b = b, W0 = W0, n0 = n0, 
    k = k, burn_in = burn_in, samples = samples, 
    beta_est_approach = beta_est_approach[1], beta = beta, 
    beta_max = beta_max, epsilon = epsilon, beta_tol = beta_tol,
    B = B, M = M)
  rm(X)
  rm(xy)
  showWelcomeMessage(BASS)

  return(BASS)
}


#' Show message
showWelcomeMessage <- function(BASS)
{
  cat("***************************************\n")
  cat("  Bayesian Analytics for Spatial Segmentation (BASS)\n")
  cat("  Authors: Zheng Li, Xiang Zhou\n")
  cat("  Affiliate: Department of Biostatistics, University of Michigan\n")
  cat("  INPUT INFO:\n")
  cat("    - Number of samples:", BASS@L, "\n")
  cat("    - Number of spots/cells:", BASS@Ns, "\n")
  cat("    - Number of genes:", BASS@P, "\n")
  cat("    - Potts interaction parameter estimation approach:", 
    BASS@beta_est_approach, "\n")
  if(BASS@beta_est_approach == "FIXED"){
    cat("    - Potts interaction parameter:", BASS@beta, "\n")
  }
  cat("    - Variance-covariance structure of gene expression features:",
    BASS@cov_struc, "\n")
  cat("  To list all hyper-parameters, Type listAllHyper(BASS_object)\n")
  cat("***************************************\n")
}


#' List all hyper-paramters of BASS
#' @export
listAllHyper <- function(BASS)
{
  cat("***************************************\n")
  cat("  Bayesian Analytics for Spatial Segmentation (BASS)\n")
  cat("  Please refer to the paper for details of paramters\n")
  cat("  ALL HYPER-PARAMETERS:\n")
  cat("    - Number of cell types C: ", BASS@C, "\n", sep = "")
  cat("    - Number of spatial domains R: ", BASS@R, "\n", sep = "")
  cat("    - Initialization method: ", BASS@init_method, "\n", sep = "")
  cat("    - Covariance structure: ", BASS@cov_struc, "\n", sep = "")
  if(BASS@cov_struc == "EII"){
    cat("    - Shape parameter of the inverse gamma prior a: ", 
      BASS@a, "\n", sep = "")
    cat("    - Scale parameter of the inverse gamma prior b: ", 
      BASS@b, "\n", sep = "")
  } else if(BASS@cov_struc == "EEE"){
    cat("    - Scale matrix of the Wishart prior W0: ", 
      BASS@W0, "I", "\n", sep = "")
    cat("    - Degrees of freedom of the Wishart prior n0: ", 
      BASS@n0, "\n", sep = "")
  }
  cat("    - MCMC burn-in samples: ", BASS@burn_in, "\n", sep = "")
  cat("    - MCMC posterior samples: ", BASS@samples, "\n", sep = "")
  cat("    - Potts interaction parameter estimation approach: ", 
    BASS@beta_est_approach, "\n", sep = "")
  if(BASS@beta_est_approach == "FIXED"){
    cat("    - Potts interaction parameter: ", BASS@beta, "\n", sep = "")
  } else if(BASS@beta_est_approach == "ACCUR_EST"){
    cat("    - Number of burn-in steps for Potts sampling: ", 
      BASS@B, "\n", sep = "")
    cat("    - Number of Potts samples to approximate Potts partition ratio: ", 
      BASS@M, "\n", sep = "")
    cat("    - Uniform random walk step size epsilon: ", 
      BASS@epsilon, "\n", sep = "")
    cat("    - Convergence tolerance for beta: ", 
      BASS@beta_tol, "\n", sep = "")
  }
  cat("    - Upper bound for beta: ", BASS@beta_max, "\n", sep = "")
  cat("    - Concentration parameter of Dirichlet distribution alpha0: ", 
    BASS@alpha0, "\n", sep = "")
  cat("    - Number of neighbors for the kNN graph: ", 
    BASS@k, "\n", sep = "")
  cat("***************************************\n")
}


#' customize printing message
setMethod(
  f = "show",
  signature = "BASS",
  definition = 
    function(object)
    {
      showWelcomeMessage(object)
    }
  )


#' Preprocess expression data
#' @export
BASS.preprocess <- function(
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
{
  geneSelect <- match.arg(geneSelect)
  X_run <- do.call(cbind, BASS@X)
  # 1.Library size normalization + log transformation
  if(doLogNormalize){
    cat("***** Log-normalize gene expression data *****\n")
    X_run <- scater::normalizeCounts(X_run, log = TRUE)
  }

  # 2.Gene selection
  if(geneSelect == "sparkx" & BASS@P > nSE){
    cat("***** Select spatially expressed genes with sparkx *****\n")
    genes <- lapply(1:BASS@L, function(l){
      capture.output(sparkx.l <- SPARK::sparkx(BASS@X[[l]], BASS@xy[[l]]))
      sparkx.l <- sparkx.l$res_mtest[order(sparkx.l$res_mtest$adjustedPval), ]
      genes.l <- head(rownames(sparkx.l), n = nSE)
    })
    genes <- unique(unlist(genes))
    X_run <- X_run[genes, ]
  } else if(geneSelect == "hvgs" & BASS@P > nHVG){
    cat("***** Select highly variable genes *****\n")
    dec <- scran::modelGeneVar(X_run)
    genes <- scran::getTopHVGs(dec, n = nHVG)
    X_run <- X_run[genes, ]
  }
  cat("***** Exclude genes with 0 expression *****\n")
  idx_rm <- apply(X_run, 1, sum) == 0
  X_run <- X_run[!idx_rm, ]

  # 3.Dimension reduction
  if(doPCA){
    cat("***** Reduce data dimension with PCA *****\n")
    X_run <- apply(X_run, MARGIN = 1, scale, 
      center = T, scale = scaleFeature)
    Q <- prcomp(X_run, scale. = F)$rotation
    X_run <- (X_run %*% Q)[, 1:nPC]
  } else{
    X_run <- t(X_run)
  }

  # 4.Batch effect correction
  if(doBatchCorrect & BASS@L > 1){
    cat("***** Correct batch effect with Harmony *****\n")
    X_run <- harmony::HarmonyMatrix(
      data_mat = X_run,
      meta_data = rep(1:BASS@L, BASS@Ns), 
      do_pca = F, verbose = F)
  }
  BASS@X_run <- X_run
  return(BASS)
}


#' Run BASS algorithm
#' @export
BASS.run <- function(BASS)
{
  xy <- do.call(rbind, BASS@xy)
  res <- Gibbs(Xin = BASS@X_run, xy = xy, Ns = BASS@Ns, C = BASS@C, R = BASS@R,
    initMethod = BASS@init_method, covStruc = BASS@cov_struc, 
    kappa = BASS@kappa, alpha0 = BASS@alpha0, a = BASS@a, b = BASS@b,
    W0in = BASS@W0 * diag(ncol(BASS@X_run)), n0 = BASS@n0, k = BASS@k, 
    warmUp = BASS@burn_in, numSamples = BASS@samples, 
    betaEstApproach = BASS@beta_est_approach, betaIn = BASS@beta, 
    betaMax = BASS@beta_max, epsilon = BASS@epsilon, betaTol = BASS@beta_tol,
    M = BASS@M, B = BASS@B)
  BASS@res <- res
  return(BASS)
}


#' Postprocess BASS results
#' @export
BASS.postprocess <- function(BASS)
{
  # adjust for label switching
  c_ls <- label.switching::label.switching(method = c("ECR-ITERATIVE-1"),
    z = t(BASS@res$c) + 1, K = BASS@C)
  z_ls <- label.switching::label.switching(method = c("ECR-ITERATIVE-1"),
    z = t(BASS@res$z) + 1, K = BASS@R)

  # 1.cell type labels
  idx_l <- c(0, cumsum(BASS@Ns))
  c_est_ls <- c_ls$clusters[1, ]
  c_est_ls <- lapply(1:BASS@L, function(l){
      c_est_ls[(idx_l[l]+1):idx_l[l+1]]
    })

  # 2.tissue structure labels
  z_est_ls <- z_ls$clusters[1, ]
  z_est_ls <- lapply(1:BASS@L, function(l){
      z_est_ls[(idx_l[l]+1):idx_l[l+1]]
    })

  # 3.cell type proportion matrix
  ncell_cr <- table(unlist(c_est_ls), unlist(z_est_ls))
  # certain tissue domains or cell types may not contain any cells
  r_wcell <- which(1:BASS@R %in% colnames(ncell_cr))
  c_wcell <- which(1:BASS@C %in% rownames(ncell_cr))
  pi_est_ls <- matrix(0, BASS@C, BASS@R)
  pi_est_ls[c_wcell, r_wcell] <- ncell_cr
  ncell_r <- apply(pi_est_ls, 2, sum)
  ncell_r[ncell_r == 0] <- 1
  pi_est_ls <- pi_est_ls %*% diag(1 / ncell_r)

  postprocess <- list(
    c_ls = c_est_ls,
    z_ls = z_est_ls,
    pi_ls = pi_est_ls
    )
  BASS@res_postprocess <- postprocess
  return(BASS)
}

