# Authors: Zheng Li, Xiang Zhou

#' BASS object
#' 
#' BASS object stores data, settings, hyper-parameters, MCMC controlling 
#' parameters, and results for the spatial transcriptomic analysis.
#'
#' @slot X An L-list of raw gene expression count matrices for L tissue 
#'   sections. For the lth tissue section, the matrix has dimension PxNl 
#'   representing P common genes shared across tissue sections and Nl 
#'   cells/spots in the lth tissue section.
#' @slot X_run An NxJ gene expression feature matrix obtained from 
#'   \code{BASS.preprocess} for running the BASS algorithm. Rows are cells from 
#'   all tissue sections and columns are low-dimensional gene expression 
#'   features. Modify this slot if you have pre-processed the data by yourself.
#' @slot xy An L-list of spatial coordinates matrices for L tissue sections.
#'  For the lth tissue section, the matrix has dimension Nlx2 representing the 
#'  same Nl cells/spots as in the expression count matrix and their x and y 
#'  spatial coordinates.
#' @slot L Numeric. Number of tissue sections.
#' @slot P Numeric. Number of genes.
#' @slot Ns Numeric Vector. Number of cells/spots in each of the L tissue 
#'   sections.
#' @slot C Numeric. Number of cell types.
#' @slot R Numeric. Number of spatial domains.
#' @slot init_method Character. Initialize the cell type clusters and spatial 
#'   domains with either k-means clustering (\code{kmeans}) or Gaussian mixture 
#'   model (\code{mclust}).
#' @slot cov_struc Character. Variance-covariance structure for the gene 
#'   expression features. \code{EEE} assumes an unstructured and common 
#'   variance-covariance structure across cell types while \code{EII} assumes 
#'   the same independent variance-covariance structure across cell types.
#' @slot kappa Numeric. Prior parameter in the gene expression model when 
#'   \code{cov_struc} is set to be \code{EII}.
#'   That is x_ij | c_i = c ~ N(mu_c, sigma2 / kappa).
#' @slot a Numeric. Shape parameter of the inverse-gamma prior on the variance 
#'   parameter when \code{cov_struc} is set to be \code{EII}.
#' @slot b Numeric. Scale parameter of the inverse-gamma prior on the variance 
#'   parameter when \code{cov_struc} is set to be \code{EII}.
#' @slot w0 Numeric. w0*I is the scale matrix of the inverse-Wishart prior on 
#'   the variance-covariance matrix when \code{cov_struc} is set to be 
#'   \code{EEE}, where I is an identity matrix.
#' @slot n0 Numeric. Degrees of freedom of the inverse-Wishart prior on the 
#'   variance-covariance matrix when \code{cov_struc} is set to be \code{EEE}.
#' @slot alpha0 Numeric. Concentration parameter of the Dirichlet prior  
#'   specified on the cell type composition vector in each spatial domain.
#' @slot  k Numeric. Minimum number of neighbors for each cell/spot based on 
#'   the Euclidean distance.
#' @slot burn_in Numeric. Number of burn-in interactions in the MCMC algorithm.
#' @slot samples Numeric. Number of posterior samples in the MCMC algorithm.
#' @slot beta_est_appraoch Character. Fix the cell-cell interaction parameter 
#'   to be \code{beta} (\code{FIXED}) or estimate the parameter based on the 
#'   data (\code{ACCUR_EST}) at hand.
#' @slot  beta Numeric. Pre-specified cell-cell interaction parameter if 
#'   \code{beta_est_appraoch} is set to be \code{FIXED}.
#' @slot  beta_max Numeric. Upper bound of the cell-cell interaction parameter.
#' @slot epsilon Numeric. Step size of a uniform random walk.
#' @slot beta_tol Numeric. Threshold of convergence for the cell-cell 
#'   interaction parameter.
#' @slot B Numeric. Number of Potts samples to approximate the partition ratio.
#' @slot M Numeric. Number of burn-in iterations in Potts sampling. 
#' @slot res List. Posterior samples of all the parameters. The list contains:
#'   \itemize{
#'     \item mu: JxCxS array with each element (j, c, s) representing the mean 
#'       of the jth feature in cell type c in the posterior sample s. 
#'     \item sigma2: JxCxS array with each element (j, c, s) representing the 
#'       variance of the jth feature in cell type c in the posterior sample s 
#'       when \code{cov_struc} is set to be \code{EII}.
#'     \item Sigma: JxJxS array with each slice s representing the JxJ 
#'       variance-covariance matrix in the gene expression model in the 
#'       posterior sample s when \code{cov_struc} is set to be \code{EEE}.
#'     \item pi: CxRxS array with each element (c, r, s) representing the 
#'       proportion of cell type c in spatial domain r in the posterior 
#'       sample s.
#'     \item z: NxS matrix of spatial domain labels with cells/spots as rows 
#'       and posterior samples as columns. Cells/spots across all tissue 
#'       sections have been stacked into a vector.
#'     \item init_z: N-vector of initial spatial domain labels obtained using 
#'       the method specified by \code{init_method}. Cells/spots across all 
#'       tissue sections have been stacked into a vector.
#'     \item c: NxS matrix of cell type labels with cells/spots as rows and 
#'       posterior samples as columns. Cells/spots across all tissue sections 
#'       have been stacked into a vector.
#'     \item init_c: N-vector of initial cell type labels obtained using the 
#'       method specified by \code{init_method}. Cells/spots across all tissue 
#'       sections have been stacked into a vector.
#'     \item burninBeta: Numeric. MCMC samples of the cell-cell interaction 
#'       parameter \code{beta}.
#'     \item beta: Numeric. Estimate of the cell-cell interaction parameter 
#'       \code{beta}.
#'     \item acceptBetaP: Numeric. Acceptance rate of the Metropolis-Hastings 
#'       algorithm for sampling the cell-cell interaction parameter \code{beta}.
#'     \item lambda: JxS matrix with each element (j, s) representing the jth 
#'       feature-specific scaling factor of mu_jc in the posterior sample s.
#'     \item d: JxS matrix with each element (j, s) representing the jth 
#'       feature-specific prior mean of mu_jc in the posterior sample s.
#'     \item R2: Jx1 matrix of squared ranges for each expression feature.
#'   }
#'  Refer to the supplementary file of the paper for the details of all the 
#'  parameters.
#' @slot res_postprocess List. Estimates of cell type labels, spatial domain 
#'   labels, and cell type composition in each spatial domain after 
#'   post-processing the posterior samples in \code{res}. The list contains:
#'   \itemize{
#'      \item c_ls: List. Cell type labels for each of the L tissue sections.
#'      \item z_ls: List. Spatial domain labels for each of the L tissue 
#'        sections.
#'      \item p_ls: Matrix. CxR matrix with the cth row and rth column 
#'        representing the proportion of cell type c in the spatial domain r.
#'   }
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
  a = "numeric",
  b = "numeric",
  w0 = "numeric",
  n0 = "numeric",
  alpha0 = "numeric",
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


#' Create a BASS object
#'
#' Create a BASS object to store data, settings, hyper-parameters, MCMC 
#' controlling parameters, and results for the spatial transcriptomic analysis.
#' 
#' @param X An L-list of gene expression count matrices for L tissue sections.
#'   For the lth tissue section, the matrix has dimension PxNl representing P 
#'   common genes shared across tissue sections and Nl cells/spots in the
#'   lth tissue section.
#' @param xy An L-list of spatial coordinates matrices for L tissue sections.
#'  For the lth tissue section, the matrix has dimension Nlx2 representing the 
#'  same Nl cells/spots as in the expression count matrix and their x and y 
#'  spatial coordinates.
#' @param C Numeric. Number of cell types.
#' @param R Numeric. Number of spatial domains.
#' @param init_method Character. Initialize the cell type clusters and spatial 
#'   domains with either k-means clustering (\code{kmeans}) or Gaussian mixture 
#'   model (\code{mclust}).
#' @param cov_struc Character. Variance-covariance structure for the gene 
#'   expression features. \code{EEE} assumes an unstructured and common 
#'   variance-covariance structure across cell types while \code{EII} assumes 
#'   the same independent variance-covariance structure across cell types.
#' @param kappa Numeric. Prior parameter in the gene expression model when 
#'   \code{cov_struc} is set to be \code{EII}.
#'   That is x_ij | c_i = c ~ N(mu_c, sigma2 / kappa).
#' @param a Numeric. Shape parameter of the inverse-gamma prior on the variance 
#'   parameter when \code{cov_struc} is set to be \code{EII}.
#' @param b Numeric. Scale parameter of the inverse-gamma prior on the variance 
#'   parameter when \code{cov_struc} is set to be \code{EII}.
#' @param w0 Numeric. w0*I is the scale matrix of the inverse-Wishart prior on 
#'   the variance-covariance matrix when \code{cov_struc} is set to be 
#'   \code{EEE}, where I is an identity matrix.
#' @param n0 Numeric. Degrees of freedom of the inverse-Wishart prior on the 
#'   variance-covariance matrix when \code{cov_struc} is set to be \code{EEE}.
#' @param alpha0 Numeric. Concentration parameter of the Dirichlet prior  
#'   specified on the cell type composition vector in each spatial domain.
#' @param  k Numeric. Minimum number of neighbors for each cell/spot based on 
#'   the Euclidean distance.
#' @param burn_in Numeric. Number of burn-in interactions in the MCMC algorithm.
#' @param samples Numeric. Number of posterior samples in the MCMC algorithm.
#' @param beta_est_appraoch Character. Fix the cell-cell interaction parameter 
#'   to be \code{beta} (\code{FIXED}) or estimate the parameter based on the 
#'   data (\code{ACCUR_EST}) at hand.
#' @param  beta Numeric. Pre-specified cell-cell interaction parameter if 
#'   \code{beta_est_appraoch} is set to be \code{FIXED}.
#' @param  beta_max Numeric. Upper bound of the cell-cell interaction parameter.
#' @param epsilon Numeric. Step size of a uniform random walk.
#' @param beta_tol Numeric. Threshold of convergence for the cell-cell 
#'   interaction parameter.
#' @param B Numeric. Number of Potts samples to approximate the partition ratio.
#' @param M Numeric. Number of burn-in iterations in Potts sampling. 
#'
#' @return A BASS object. Refer to \linkS4class{BASS} for a complete list of
#'   slots in the object.
#' @export
#'
createBASSObject <- function(
  X, xy, C, R, 
  init_method = c("kmeans", "mclust"), 
  cov_struc = c("EEE", "EII"),
  kappa = 0.01, 
  a = 0.01, 
  b = 0.01, 
  w0 = 1,
  n0 = 1,
  alpha0 = 1, 
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
  for(l in 1:L)
  {
    if(!is.matrix(X[[l]])){
      X[[l]] <- as.matrix(X[[l]])
      cat("Expression data coerced to a matrix\n")
    }
    if(!is.matrix(xy[[l]])){
      xy[[l]] <- as.matrix(xy[[l]])
      cat("Location data coerced to a matrix\n")
    }
  }
  
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
  }

  # create object
  BASS <- new(Class = "BASS", 
    X = X, xy = xy, L = L, P = P, Ns = sapply(X, ncol), C = C, R = R,
    init_method = init_method[1], cov_struc = cov_struc[1], 
    kappa = kappa, alpha0 = alpha0, a = a, b = b, w0 = w0, n0 = n0, 
    k = k, burn_in = burn_in, samples = samples, 
    beta_est_approach = beta_est_approach[1], beta = beta, 
    beta_max = beta_max, epsilon = epsilon, beta_tol = beta_tol,
    B = B, M = M)
  rm(X)
  rm(xy)
  showWelcomeMessage(BASS)

  return(BASS)
}


#' Show BASS object
#' @noRd
showWelcomeMessage <- function(BASS)
{
  cat("***************************************\n")
  cat("  INPUT INFO:\n")
  cat("    - Number of tissue sections:", BASS@L, "\n")
  cat("    - Number of cells/spots:", BASS@Ns, "\n")
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


#' Show BASS settings
#' 
#' Show a complete list of settings, hyper-parameters, and MCMC controlling
#' parameters used by BASS.
#' 
#' @param BASS A BASS object created by \code{createBASSObject}.
#' 
#' @export
listAllHyper <- function(BASS)
{
  cat("***************************************\n")
  cat("  Please refer to the paper for details of the paramters\n")
  cat("  ALL HYPER-PARAMETERS:\n")
  cat("    - Number of cell types C: ", BASS@C, "\n", sep = "")
  cat("    - Number of spatial domains R: ", BASS@R, "\n", sep = "")
  cat("    - Initialization method: ", BASS@init_method, "\n", sep = "")
  cat("    - Covariance structure: ", BASS@cov_struc, "\n", sep = "")
  if(BASS@cov_struc == "EII"){
    cat("    - Shape parameter of the inverse-gamma prior a: ", 
      BASS@a, "\n", sep = "")
    cat("    - Scale parameter of the inverse-gamma prior b: ", 
      BASS@b, "\n", sep = "")
  } else if(BASS@cov_struc == "EEE"){
    cat("    - Scale matrix of the inverse-Wishart prior W0: ", 
      BASS@w0, "I", "\n", sep = "")
    cat("    - Degrees of freedom of the inverse-Wishart prior n0: ", 
      BASS@n0, "\n", sep = "")
  }
  cat("    - Number of MCMC burn-in iterations: ", BASS@burn_in, "\n", sep = "")
  cat("    - Number of MCMC posterior samples: ", BASS@samples, "\n", sep = "")
  cat("    - Potts interaction parameter estimation approach: ", 
    BASS@beta_est_approach, "\n", sep = "")
  if(BASS@beta_est_approach == "FIXED"){
    cat("    - Potts interaction parameter: ", BASS@beta, "\n", sep = "")
  } else if(BASS@beta_est_approach == "ACCUR_EST"){
    cat("    - Number of burn-in interations in Potts sampling: ", 
      BASS@B, "\n", sep = "")
    cat("    - Number of Potts samples to approximate the partition ratio: ", 
      BASS@M, "\n", sep = "")
    cat("    - Step size of a uniform random walk epsilon: ", 
      BASS@epsilon, "\n", sep = "")
    cat("    - Threshold of convergence for beta: ", 
      BASS@beta_tol, "\n", sep = "")
  }
  cat("    - Upper bound for beta: ", BASS@beta_max, "\n", sep = "")
  cat("    - Concentration parameter of the Dirichlet prior alpha0: ", 
    BASS@alpha0, "\n", sep = "")
  cat("    - Minimum number of neighbors for each cell/spot based on the ",
    "Euclidean distance: ", BASS@k, "\n", sep = "")
  cat("***************************************\n")
}


#' Overload print for BASS object
#' @noRd
setMethod(
  f = "show",
  signature = "BASS",
  definition = 
    function(object)
    {
      showWelcomeMessage(object)
    }
  )


#' Pre-process gene expression data
#' 
#' Pre-process the gene expression data by first combining cells/spots across 
#' all L tissue sections, conduct library size normalization followed by a 
#' log2-transformation (after adding a pseudo-count of 1), select feature genes,  
#' perform dimension reduction with PCA on the normalized expression matrix to 
#' extract J low-dimensional expression features, and perform batch effect 
#' adjustment with Harmony (Korsunsky et al., 2019) to align expression data 
#' from different tissue sections.
#' 
#' @param BASS A BASS object created by \code{createBASSObject}.
#' @param doLogNormalize Logical. Whether to perform the library size 
#'   normalization followed by a log2-transformation (after adding a 
#'   pseudo-count of 1).
#' @param geneSelect Character. Perform feature selection by either selecting 
#'   spatially variable genes with SPARK-X (Zhu et al., 2021) (\code{sparkx}) 
#'   or highly variable genes with the scran package (\code{hvgs}).
#' @param doPCA Logical. Whether to perform PCA on the normalized expression 
#'   matrix to extract \code{nPC} low-dimensional expression features.
#' @param scaleFeature Logical. Whether to center and standardize each gene 
#'  to have a mean of 0 and standard deviation of 1 before performing PCA.
#' @param doBatchCorrect Logical. Whether to perform batch effect adjustment 
#'   with Harmony (Korsunsky et al., 2019) to align expression data from 
#'   different tissue sections.
#' @param nHVG Numeric. Number of highly variable genes to select if 
#'   \code{geneSelect} is set to be \code{hvgs}.
#' @param nSE Numeric. Number of spatially variable genes to select if 
#'   \code{geneSelect} is set to be \code{sparkx}.
#' @param nPC Numeric. Number of top PCs to extract if \code{doPCA} is set to 
#'   be \code{TRUE}.
#' 
#' @return A BASS object with the processed gene expression data stored in the 
#'   \code{X_run} slot. Refer to \linkS4class{BASS} for a complete list of slots 
#'   in the object.
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


#' Run BASS inference algorithm
#' 
#' Run the Gibbs sampling algorithm in combination with the Metroplis-Hastings 
#' algorithm to perform parameter inference.
#' 
#' @param BASS A BASS object containing the pre-processed gene expression data 
#'   in the \code{X_run} slot. The pre-processing can be conducted either with 
#'   the provided \code{BASS.preprocess} function or by the users.
#' 
#' @return A BASS object with posterior samples of all the parameters stored in 
#'   the \code{res} slot. Refer to \linkS4class{BASS} for details of the 
#'   \code{res} slot and a complete list of slots in the object.
#' @export
BASS.run <- function(BASS)
{
  xy <- do.call(rbind, BASS@xy)
  res <- Gibbs(Xin = BASS@X_run, xy = xy, Ns = BASS@Ns, C = BASS@C, R = BASS@R,
    initMethod = BASS@init_method, covStruc = BASS@cov_struc, 
    kappa = BASS@kappa, alpha0 = BASS@alpha0, a = BASS@a, b = BASS@b,
    W0in = BASS@w0 * diag(ncol(BASS@X_run)), n0 = BASS@n0, k = BASS@k, 
    warmUp = BASS@burn_in, numSamples = BASS@samples, 
    betaEstApproach = BASS@beta_est_approach, betaIn = BASS@beta, 
    betaMax = BASS@beta_max, epsilon = BASS@epsilon, betaTol = BASS@beta_tol,
    M = BASS@M, B = BASS@B)
  BASS@res <- res
  return(BASS)
}


#' Post-process the posterior sampling results
#' 
#' Post-process the posterior sampling results to address the label switching 
#' issue associated with the sampling of cell type labels and spatial domain 
#' labels based on the ECR-1 algorithm, estimate the cell type labels and 
#' spatial domain labels as the mode of all their posterior samples, and 
#' estimate the cell type composition in each spatial domain based on the final 
#' estimates of cell type and spatial domain labels.
#'   
#' @param BASS A BASS object output from \code{BASS.run}.
#'
#' @return A BASS object with the final estimates of cell type labels, spatial 
#'   domain labels, and cell type composition in each spatial domain stored in 
#'   the \code{res_postprocess} slot. Refer to \linkS4class{BASS} for details of 
#'   the \code{res_postprocess} slot and a complete list of slots in the object.
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

