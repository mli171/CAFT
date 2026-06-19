#' Unrestricted CAFT Estimation
#'
#' @description
#' Fits the unrestricted rank-based compositional accelerated failure time
#' model for each taxon in microbiome count data with zero cells. Zero counts
#' are treated as censored observations below sample-specific detection limits.
#' This function performs only taxon-level unrestricted estimation and does not
#' conduct hypothesis testing, multiple-testing adjustment, or bootstrap
#' calibration.
#'
#' @param otu.table The community OTU table (or taxa count table). Each row
#'  corresponds to a sample and each column corresponds to one OTU (taxa).
#' @param x Design matrix containing all covariates.
#'  All K-level categorical variables should be converted to K-1 indicator
#'  variables before input.
#' @param filter.thresh a real value between 0 and 1 for OTU table sample
#'  presence filtering. Any OTUs present in fewer than \code{filter.thresh}
#'  proportion of samples are filtered out. We set the default to be 0.05.
#' @param regularize A logical value. If \code{TRUE}, adds a small penalty to
#'   the rank-based score to improve numerical stability. This can be useful
#'   when the unpenalized rank equations have flat directions, heavy censoring,
#'   many ties, or when the optimizer has convergence difficulty. Default is
#'   \code{TRUE}.
#' @param n.cores Integer. Number of CPU cores to use for parallel computation.
#'  Default is \code{1} (no parallelism). If \code{n.cores > 1}, the function
#'  runs tasks in parallel using \code{foreach}/\code{doParallel} with a PSOCK
#'  cluster. On typical desktops/laptops, a good choice is
#'  \code{max(1L, parallel::detectCores() - 1L)}.
#'
#' @details
#' CAFT fits a rank-based accelerated failure time model to the negative
#' log-relative abundance of each taxon. Zero counts are treated as censored
#' observations below the sample-specific detection limit, avoiding the use of
#' pseudocounts.
#'
#' This function performs only the unrestricted taxon-level estimation step. It
#' returns the unrestricted coefficient estimates and the transformed
#' survival-type data needed by downstream CAFT testing functions. It does not
#' compute the median-based null value, restricted estimates, score test
#' statistics, p-values, q-values, or bootstrap-calibrated p-values.
#'
#' @return An object of class \code{"caft_est"} containing unrestricted
#'   taxon-level coefficient estimates, censoring variables, the centered design
#'   matrix, and filtering/fitting status indicators.
#' \describe{
#' \item{call}{The matched function call.}
#' \item{otu.table}{The OTU count table used in the unrestricted CAFT
#'  estimation after input checking.}
#' \item{lib.size}{The sequencing library size for each sample.}
#' \item{x.raw}{The original, uncentered covariate matrix used in the model.}
#' \item{x}{The centered covariate matrix used in the unrestricted rank-based
#'  AFT estimation.}
#' \item{x.name}{Names of the covariates in the design matrix.}
#' \item{filter.thresh}{The prevalence filtering threshold used in the
#'  analysis.}
#' \item{regularize}{Logical value indicating whether regularized rank-based
#'  AFT estimation was used.}
#' \item{n.data}{Number of samples used in the unrestricted estimation.}
#' \item{n.taxa}{Number of taxa in the OTU table.}
#' \item{n.param}{Number of covariates in the design matrix.}
#' \item{taxa.name}{Names of the taxa in the OTU table.}
#' \item{taxa.used}{Names of taxa that passed the prevalence filter and had
#'  successful unrestricted estimation.}
#' \item{t.star.all}{Transformed censored abundance outcomes used in the
#'  rank-based AFT estimation.}
#' \item{delta.all}{Censoring indicators for \code{t.star.all}, where one
#'  indicates an observed nonzero count and zero indicates a censored zero
#'  count.}
#' \item{beta.est}{A matrix that includes the unrestricted estimated
#'  coefficients obtained by fitting each OTU count to the CAFT log-linear
#'  model. The number of columns matches the number of covariates from
#'  \code{x}.}
#' \item{skip.otu}{Indicator for OTUs skipped during taxa presence filtering.}
#' \item{skip.rare}{Indicator for OTUs skipped during taxa presence filtering.}
#' \item{skip.fail.rank.fit.pen}{Indicator for OTUs whose unrestricted
#'  regularized rank-based AFT estimation failed.}
#' \item{fit.error.messages}{Error messages from failed unrestricted
#'  taxon-level fits, if any.}
#' }
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#'
#' @export
#'
#' @examples
#' data(Colon)
#'
#' count.tab <- Colon$otu
#' sample.tab <- Colon$meta
#'
#' keep.sample <- !is.na(sample.tab$age)
#' count.tab <- count.tab[keep.sample, , drop = FALSE]
#' sample.tab <- sample.tab[keep.sample, , drop = FALSE]
#'
#' keep.otu <- colSums(count.tab > 0) > 1
#' count.tab <- count.tab[, keep.otu, drop = FALSE]
#'
#' ## Use a small subset of taxa for a fast example.
#' prev <- colSums(count.tab > 0)
#' keep.top <- names(sort(prev, decreasing = TRUE))[seq_len(min(10L, length(prev)))]
#' count.tab <- count.tab[, keep.top, drop = FALSE]
#'
#' ## Remove samples with zero library size after subsetting taxa.
#' keep.lib <- rowSums(count.tab) > 0
#' count.tab <- count.tab[keep.lib, , drop = FALSE]
#' sample.tab <- sample.tab[keep.lib, , drop = FALSE]
#'
#' Disease1 <- Disease2 <- rep(0, nrow(sample.tab))
#' Disease1[sample.tab$disease == "CRC"] <- 1
#' Disease2[sample.tab$disease == "adenoma"] <- 1
#'
#' Age <- as.numeric(sample.tab$age)
#' Gender <- as.numeric(factor(sample.tab$gender)) - 1
#'
#' x.test <- cbind(CRC = Disease1, adenoma = Disease2)
#' x.adj <- cbind(Age = Age, Gender = Gender)
#' x <- cbind(x.test, x.adj)
#'
#' est <- caft_estimate(
#'   otu.table = count.tab,
#'   x = x,
#'   regularize = TRUE,
#'   n.cores = 1L
#' )
#'
#' print(est)
caft_estimate <- function(otu.table,
                          x,
                          filter.thresh = 0.05,
                          regularize = TRUE,
                          n.cores = 1L){

  if (!is.matrix(otu.table)) {otu.table <- as.matrix(otu.table)}

  n.data  <- NROW(otu.table)
  n.taxa  <- NCOL(otu.table)

  if (missing(x) || is.null(x)) {
    stop("The full design matrix 'x' must be provided.", call. = FALSE)
  }
  x <- as.matrix(x)

  if (NROW(otu.table) != NROW(x)) {
    stop(" Number of samples not match between OTU table and covairates matrix!")
  }
  if (anyNA(otu.table)) {
    stop("otu.table contains missing values. Please remove or impute missing values before calling caft_estimate().")
  }
  if (anyNA(x)) {
    stop("The covariate matrix contains missing values. Please remove or impute missing values before calling caft_estimate().")
  }

  lib.size <- rowSums(otu.table)
  if (any(lib.size <= 0)) {
    stop("All samples must have positive total library size.")
  }
  # center the covariates
  x.raw <- x
  x <- scale(x.raw, center = TRUE, scale = FALSE)
  x <- as.matrix(x)

  n.param <- NCOL(x)
  n.cores <- as.integer(n.cores)

  taxa.name <- colnames(otu.table)
  if (is.null(taxa.name)) {
    taxa.name <- paste0("taxon", seq_len(n.taxa))
    colnames(otu.table) <- taxa.name
  }
  x.name <- colnames(x.raw)
  if (is.null(x.name)) {
    x.name <- paste0("x", seq_len(n.param))
    colnames(x.raw) <- x.name
    colnames(x) <- x.name
  }

  # Relative abundance
  ra.all <- otu.table / lib.size
  # detection limits: Censored relative abundance: every row should same
  lim.ra <- matrix(1 / lib.size, nrow = n.data, ncol = n.taxa)

  #--------------------------
  # create survival data
  #--------------------------
  log_neg <- -log10(lim.ra)
  log_pos <- -log10(ra.all)

  t.star.all <- matrix(log_neg, nrow=n.data, ncol=n.taxa)
  idx <- (otu.table > 0)
  t.star.all[idx] <- log_pos[idx]

  t.star.all <- as.data.frame(t.star.all)
  colnames(t.star.all) <- taxa.name

  delta.all <- as.data.frame((otu.table > 0) * 1L)
  colnames(delta.all) <- taxa.name

  # fit_one_taxon(): one function for seq and parallel computing
  fit_one_taxon <- function(ii) {
    tstar <- t.star.all[[ii]]
    delta.1 <- delta.all[[ii]]

    if (sum(delta.1) <= n.data * filter.thresh) {
      return(list(
        beta = rep(NA_real_, n.param),
        skip_rare = 1L,
        skip_fail_fit = 0L,
        error_message = NA_character_
      ))
    }

    fit0.pen <- try(
      estimate.rank.aft(
        y = tstar,
        delta = delta.1,
        x = x,
        Gamma = NULL,
        Lambda = diag(n.param),
        Gamma.ginv = NULL,
        Lambda.ginv = diag(n.param),
        b = NULL,
        beta = NULL,
        test = TRUE,
        regularize = regularize,
        tol = 1e-12
      ),
      silent = TRUE
    )

    if (inherits(fit0.pen, "try-error")) {
      return(list(
        beta = rep(NA_real_, n.param),
        skip_rare = 0L,
        skip_fail_fit = 1L,
        error_message = as.character(fit0.pen)
      ))
    }

    list(
      beta = as.numeric(fit0.pen$beta),
      skip_rare = 0L,
      skip_fail_fit = 0L,
      error_message = NA_character_
    )
  }

  n.cores <- min(n.cores, n.taxa)

  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores, type = "PSOCK")
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(cl)

    res_phase1 <- foreach::foreach(
      ii = seq_len(n.taxa),
      .errorhandling = "pass"
    ) %dopar% {
      fit_one_taxon(ii)
    }
  } else {
    res_phase1 <- lapply(seq_len(n.taxa), fit_one_taxon)
  }

  # deal with errors possiblely from foreach
  res_phase1 <- lapply(res_phase1, function(z) {
    if (inherits(z, "error")) {
      list(
        beta = rep(NA_real_, n.param),
        skip_rare = 0L,
        skip_fail_fit = 1L,
        error_message = "Parallel unrestricted taxon-level fit failed."
      )
    } else {
      z
    }
  })

  beta.est <- do.call(rbind, lapply(res_phase1, `[[`, "beta"))
  beta.est <- as.data.frame(beta.est)
  colnames(beta.est) <- paste0("b", seq_len(n.param), ".est")
  rownames(beta.est) <- taxa.name

  skip.rare <- vapply(res_phase1, function(z) z$skip_rare, integer(1))
  skip.fail.rank.fit.pen <- vapply(res_phase1, function(z) z$skip_fail_fit, integer(1))
  fit.error.messages <- vapply(res_phase1, function(z) z$error_message, character(1))

  names(skip.rare) <- taxa.name
  names(skip.fail.rank.fit.pen) <- taxa.name
  names(fit.error.messages) <- taxa.name

  taxa.used <- taxa.name[skip.rare == 0L & skip.fail.rank.fit.pen == 0L & stats::complete.cases(beta.est)]

  out <- list(
    call = match.call(),
    otu.table = otu.table,
    lib.size = lib.size,
    x.raw = x.raw,
    x = x,
    x.name = x.name,
    filter.thresh = filter.thresh,
    regularize = regularize,
    n.data = n.data,
    n.taxa = n.taxa,
    n.param = n.param,
    taxa.name = taxa.name,
    taxa.used = taxa.used,
    t.star.all = t.star.all,
    delta.all = delta.all,
    beta.est = beta.est,
    skip.otu = skip.rare,
    skip.rare = skip.rare,
    skip.fail.rank.fit.pen = skip.fail.rank.fit.pen,
    fit.error.messages = fit.error.messages
  )

  class(out) <- "caft_est"
  out

}

#' Print Method for Unrestricted CAFT Estimation
#'
#' @param x An object of class \code{"caft_est"}.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print caft_est
#' @export
print.caft_est <- function(x, ...) {
  cat("CAFT unrestricted estimation\n")
  cat("  Samples:", x$n.data, "\n")
  cat("  Taxa:", x$n.taxa, "\n")
  cat("  Covariates:", x$n.param, "\n")
  cat("  Taxa used:", length(x$taxa.used), "\n")
  cat("  Taxa skipped by rarity filter:", sum(x$skip.rare == 1L), "\n")
  cat(
    "  Taxa failed in unrestricted fit:",
    sum(x$skip.fail.rank.fit.pen == 1L),
    "\n"
  )

  invisible(x)
}
