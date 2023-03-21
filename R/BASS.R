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
#' @slot X_run A JxN gene expression feature matrix obtained from 
#'   \code{BASS.preprocess} for running the BASS algorithm. Columns are cells 
#'   from all tissue sections and rows are low-dimensional gene expression 
#'   features. Modify this slot if you have pre-processed the data by yourself.
#' @slot xy An L-list of spatial coordinates matrices for L tissue sections.
#'  For the lth tissue section, the matrix has dimension Nlx2 representing the 
#'  same Nl cells/spots as in the expression count matrix and their x and y 
#'  spatial coordinates.
#' @slot L Numeric. Number of tissue sections.
#' @slot P Numeric. Number of genes.
#' @slot Ns Numeric vector. Number of cells/spots in each of the L tissue 
#'   sections.
#' @slot C Numeric. Number of cell types.
#' @slot R Numeric. Number of spatial domains.
#' @slot init_method Character. Initialize the cell type clusters and spatial 
#'   domains with either k-means clustering (\code{kmeans}) or Gaussian mixture 
#'   model (\code{mclust}).
#' @slot psi0 Numeric. psi0*I is the scale matrix of the inv-Wishart prior on 
#'   the variance-covariance matrix, where I is an identity matrix.
#' @slot n0 Numeric. Degrees of freedom of the inv-Wishart prior on the 
#'   variance-covariance matrix.
#' @slot alpha0 Numeric. Concentration parameter of the Dirichlet prior  
#'   specified on the cell type composition vector in each spatial domain.
#' @slot k Numeric. Minimum number of neighbors for each cell/spot based on 
#'   the Euclidean distance.
#' @slot burnin Numeric. Number of burn-in iterations in the MCMC.
#' @slot nsample Numeric. Number of posterior samples in the MCMC.
#' @slot beta_method Character. Fix the cell-cell interaction parameter 
#'   to be \code{beta} (\code{fix}) or estimate the parameter based on the 
#'   data at hand (\code{SW}).
#' @slot  beta Numeric. Pre-specified cell-cell interaction parameter if 
#'   \code{beta_method} is set to be \code{fix} or otherwise its initial value.
#' @slot step_size Numeric. Step size of a uniform random walk.
#' @slot tol Numeric. Threshold of convergence for the cell-cell 
#'   interaction parameter.
#' @slot potts_burnin Numeric. Number of burn-in iterations in Potts sampling.
#' @slot potts_nsample Numeric. Number of Potts samples to approximate the 
#'   partition ratio.
#' @slot samples List. Samples from the MCMC. The list contains:
#'   \itemize{
#'     \item c: NxS matrix of cell type labels with cells/spots as rows and 
#'       posterior samples as columns. Cells/spots across all tissue sections 
#'       have been stacked into a vector.
#'     \item z: NxS matrix of spatial domain labels with cells/spots as rows 
#'       and posterior samples as columns. Cells/spots across all tissue 
#'       sections have been stacked into a vector.
#'     \item beta: Numeric. Burn-in samples of the cell-cell interaction 
#'       parameter \code{beta}.
#'     \item pi: CxRxS array with each element (c, r, s) representing the 
#'       proportion of cell type c in spatial domain r in the posterior 
#'       sample s.
#'     \item mu: JxCxS array with each element (j, c, s) representing the mean 
#'       of the jth feature in cell type c in the posterior sample s. 
#'     \item Sigma: JxJxS array with each slice s representing the JxJ 
#'       variance-covariance matrix in the gene expression model in the 
#'       posterior sample s.
#'     \item d: JxS matrix with each element (j, s) representing the jth 
#'       feature-specific prior mean of mu_jc in the posterior sample s.
#'     \item lambda: JxS matrix with each element (j, s) representing the jth 
#'       feature-specific scaling factor of mu_jc in the posterior sample s.
#'     \item loglik: log-likelihood in each MCMC iteration.
#'   }
#'  Refer to the supplementary file of the paper for details of all the 
#'  parameters.
#' @slot results List. Estimates of cell type labels, spatial domain 
#'   labels, and cell type composition in each spatial domain after 
#'   post-processing the posterior samples in \code{samples}.
#'   The list contains:
#'   \itemize{
#'     \item init_c: N-vector of initial cell type labels obtained using the 
#'       method specified by \code{init_method}. Cells/spots across all tissue 
#'       sections have been stacked into a vector.
#'     \item init_z: N-vector of initial spatial domain labels obtained using 
#'       the method specified by \code{init_method}. Cells/spots across all 
#'       tissue sections have been stacked into a vector.
#'      \item c: List. Cell type labels for each of the L tissue sections.
#'      \item z: List. Spatial domain labels for each of the L tissue sections.
#'      \item p: Numeric matrix. CxR matrix with the cth row and rth column 
#'        representing the proportion of cell type c in the spatial domain r.
#'     \item beta: Numeric. Estimate of the cell-cell interaction parameter 
#'       \code{beta}.
#'     \item loglik: log-likelihood.
#'   }
#' @slot elapsed_time: Elapsed time for running MCMC in seconds.
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
  psi0 = "numeric",
  n0 = "numeric",
  alpha0 = "numeric",
  k = "numeric",
  burnin = "numeric",
  nsample = "numeric",
  beta_method = "character",
  beta = "numeric",
  step_size = "numeric",
  tol = "numeric",
  potts_nsample = "numeric",
  potts_burnin = "numeric",
  samples = "list",
  results = "list",
  elapsed_time = "numeric"
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
#' @param psi0 Numeric. psi0*I is the scale matrix of the inv-Wishart prior on 
#'   the variance-covariance matrix, where I is an identity matrix.
#' @param n0 Numeric. Degrees of freedom of the inv-Wishart prior on the 
#'   variance-covariance matrix.
#' @param alpha0 Numeric. Concentration parameter of the Dirichlet prior  
#'   specified on the cell type composition vector in each spatial domain.
#' @param k Numeric. Minimum number of neighbors for each cell/spot based on 
#'   the Euclidean distance.
#' @param burnin Numeric. Number of burn-in iterations in the MCMC.
#' @param nsample Numeric. Number of posterior samples in the MCMC.
#' @param beta_method Character. Fix the cell-cell interaction parameter 
#'   to be \code{beta} (\code{fix}) or estimate the parameter based on the 
#'   data at hand (\code{SW}) .
#' @param  beta Numeric. Pre-specified cell-cell interaction parameter if 
#'   \code{beta_method} is set to be \code{fix} or otherwise its initial value.
#' @param step_size Numeric. Step size of a uniform random walk.
#' @param tol Numeric. Threshold of convergence for the cell-cell interaction 
#'   parameter.
#' @param potts_burnin Numeric. Number of burn-in iterations in Potts sampling. 
#' @param potts_nsample Numeric. Number of Potts samples to approximate the 
#'   partition ratio.
#'
#' @return A BASS object. Refer to \linkS4class{BASS} for a complete list of
#'   slots in the object.
#' @export
#'
createBASSObject <- function(
  X, xy, C, R, 
  init_method = c("kmeans", "mclust"), 
  psi0 = 1,
  n0 = 1,
  alpha0 = 1, 
  k = 4, 
  burnin = 2000, 
  nsample = 5000,
  beta_method = c("SW", "fix"),
  beta = 1,
  step_size = 0.1,
  tol = 1e-3,
  potts_burnin = 10,
  potts_nsample = 10
  )
{
  L <- length(X)
  for(l in 1:L)
  {
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
    init_method = match.arg(init_method), alpha0 = alpha0, psi0 = psi0, 
    n0 = n0, k = k, burnin = burnin, nsample = nsample, beta = beta, 
    beta_method = match.arg(beta_method), step_size = step_size, tol = tol,
    potts_burnin = potts_burnin, potts_nsample = potts_nsample)
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
  cat("    - Potts interaction parameter estimation method:", 
    BASS@beta_method, "\n")
  if(BASS@beta_method == "fix"){
    cat("      - Fix Potts interaction parameter to be:", BASS@beta, "\n")
  } else if(BASS@beta_method == "SW"){
    cat("      - Estimate Potts interaction parameter with SW algorithm\n")
  }
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
  cat("    - Scale matrix of the inverse-Wishart prior Psi0: ", 
    BASS@psi0, "I", "\n", sep = "")
  cat("    - Degrees of freedom of the inverse-Wishart prior n0: ", 
    BASS@n0, "\n", sep = "")
  cat("    - Number of MCMC burn-in iterations: ", BASS@burnin, "\n", sep = "")
  cat("    - Number of MCMC posterior samples: ", BASS@nsample, "\n", sep = "")
  cat("    - Potts interaction parameter estimation approach: ", 
    BASS@beta_method, "\n", sep = "")
  if(BASS@beta_method == "fix"){
    cat("      - Fix Potts interaction parameter to be: ", 
      BASS@beta, "\n", sep = "")
  } else if(BASS@beta_method == "SW"){
    cat("    - Number of burn-in interations in Potts sampling: ", 
      BASS@potts_burnin, "\n", sep = "")
    cat("    - Number of Potts samples to approximate the partition ratio: ", 
      BASS@potts_nsample, "\n", sep = "")
    cat("    - Step size of a uniform random walk: ", 
      BASS@step_size, "\n", sep = "")
    cat("    - Threshold of convergence for beta: ", 
      BASS@tol, "\n", sep = "")
  }
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
      meta_data = as.character(rep(1:BASS@L, BASS@Ns)), 
      do_pca = F, verbose = F)
  }
  BASS@X_run <- t(X_run)
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
#'   the \code{samples} slot. Refer to \linkS4class{BASS} for details of the 
#'   \code{samples} slot and a complete list of slots in the object.
#' @importFrom mclust Mclust mclustBIC
#' @export
BASS.run <- function(BASS)
{
  xy <- do.call(rbind, BASS@xy)
  # initialize cell type labels, spatial domain labels, and
  # mean parameters of gene expression features
  if(BASS@init_method == "kmeans"){
    km_c <- stats::kmeans(x = t(BASS@X_run), centers = BASS@C, 
      iter.max = 100, nstart = 1000)
    km_z <- stats::kmeans(x = t(BASS@X_run), centers = BASS@R,
      iter.max = 100, nstart = 1000)
    init_c <- km_c$cluster - 1
    init_z <- km_z$cluster - 1
    init_mu <- t(km_c$centers)
  } else if(BASS@init_method == "mclust"){
    mclust_c <- mclust::Mclust(data = t(BASS@X_run), G = BASS@C,
      modelNames = "EEE", verbose = FALSE)
    mclust_z <- mclust::Mclust(data = t(BASS@X_run), G = BASS@R,
      modelNames = "EEE", verbose = FALSE)
    init_c <- mclust_c$classification - 1
    init_z <- mclust_z$classification - 1
    init_mu <- mclust_c$parameters$mean
  }
  
  output <- BASSFit(X = BASS@X_run, xy = xy, ns = BASS@Ns, 
    C = BASS@C, R = BASS@R, k = BASS@k, alpha0 = BASS@alpha0, 
    Psi0 = BASS@psi0 * diag(nrow(BASS@X_run)), n0 = BASS@n0,
    burnin = BASS@burnin, nsample = BASS@nsample, method = BASS@beta_method,
    potts_burnin = BASS@potts_burnin, potts_nsample = BASS@potts_nsample,
    step_size = BASS@step_size, tol = BASS@tol, init_c = init_c, 
    init_z = init_z, init_mu = init_mu, init_beta = BASS@beta)
  
  BASS@samples <- output[c("c", "z", "beta", "pi", "mu", "Sigma", "d", 
    "lambda", "loglik")]
  BASS@results[c("init_c", "init_z", "beta")] <- list(init_c, init_z, 
    output$beta_est)
  BASS@results["loglik"] <- mean(BASS@samples$loglik)
  BASS@elapsed_time <- output$elapsed_time
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
#' @param adjustLS Logical. Adjust for label switching using the ECR-1 
#'   algorithm (TRUE) or directly use the posterior mode of cell type or 
#'   spatial domain labels as the final estimates without adjusting for label 
#'   switching (FALSE).
#'
#' @return A BASS object with the final estimates of cell type labels, spatial 
#'   domain labels, and cell type composition in each spatial domain stored in 
#'   the \code{results} slot. Refer to \linkS4class{BASS} for details of the 
#'   \code{results} slot and a complete list of slots in the object.
#' @export
BASS.postprocess <- function(BASS, adjustLS = TRUE)
{
  cat("Post-processing...\n")
  if(adjustLS){ # adjust for label switching
    nTypes <- max(BASS@samples$c) + 1
    nDomains <- max(BASS@samples$z) + 1
    if(nTypes > 1){
      invisible(capture.output(c_ls <- label.switching::label.switching(
        method = c("ECR-ITERATIVE-1"), z = t(BASS@samples$c) + 1, 
        K = nTypes))) 
    }
    if(nDomains > 1){
      invisible(capture.output(z_ls <- label.switching::label.switching(
        method = c("ECR-ITERATIVE-1"), z = t(BASS@samples$z) + 1, 
        K = nDomains)))
    }

    # 1.cell type labels
    section_idx <- c(0, cumsum(BASS@Ns))
    if(nTypes > 1){
      c_est_ls <- c_ls$clusters[1, ]
    } else{
      c_est_ls <- rep(1, sum(BASS@Ns))
    }
    c_est_ls <- lapply(1:BASS@L, function(l){
        c_est_ls[(section_idx[l]+1):section_idx[l+1]]
      })

    # 2.spatial domain labels
    if(nDomains > 1){
      z_est_ls <- z_ls$clusters[1, ]
    } else{
      z_est_ls <- rep(1, sum(BASS@Ns))
    }
    z_est_ls <- lapply(1:BASS@L, function(l){
        z_est_ls[(section_idx[l]+1):section_idx[l+1]]
      })
  } else{ # use posterior mode
    getmode <- function(v)
    {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    c_est_ls <- apply(BASS@samples$c, MARGIN = 1, getmode) + 1
    z_est_ls <- apply(BASS@samples$z, MARGIN = 1, getmode) + 1
    section_idx <- c(0, cumsum(BASS@Ns))
    c_est_ls <- lapply(1:BASS@L, function(l){
      c_est_ls[(section_idx[l]+1):section_idx[l+1]]
    })
    z_est_ls <- lapply(1:BASS@L, function(l){
      z_est_ls[(section_idx[l]+1):section_idx[l+1]]
    })
  }
  
  # 3.cell type proportion matrix
  ncell_cr <- table(unlist(c_est_ls), unlist(z_est_ls))
  # certain spatial domains or cell types may not contain any cells
  r_wcell <- which(1:BASS@R %in% colnames(ncell_cr))
  c_wcell <- which(1:BASS@C %in% rownames(ncell_cr))
  pi_est_ls <- matrix(0, BASS@C, BASS@R)
  pi_est_ls[c_wcell, r_wcell] <- ncell_cr
  ncell_r <- apply(pi_est_ls, 2, sum)
  ncell_r[ncell_r == 0] <- 1
  pi_est_ls <- pi_est_ls %*% diag(1 / ncell_r)

  results <- list(
    c = c_est_ls,
    z = z_est_ls,
    pi = pi_est_ls
    )
  BASS@results[c("c", "z", "pi")] <- results
  cat("done\n")
  return(BASS)
}

