#' Rank-Based Compositional Analysis using Log-Linear Models for Microbiome Data with Zero Cells
#'
#' @description This is a novel statistical framework for differential abundance
#'  analysis of microbiome data, termed the Compositional Accelerated Failure
#'  Time (CAFT) model. The CAFT model addresses zero read counts by treating
#'  them as censored observations below the detection limit, similar to
#'  censoring mechanisms employed in survival analysis. This approach is
#'  inherently resistant to multiplicative bias, eliminates the need for
#'  pseudocounts, and addresses compositional bias through the
#'  establishment of appropriate score test procedures.
#'
#' In addition to the default asymptotic CAFT inference, the function optionally
#' provides bootstrap calibration of taxon-level p-values. Bootstrap calibration
#' is mainly intended as a sensitivity analysis or preferred calibration method
#' when the number of taxa is small in the data set.
#'
#' @param otu.table The community OTU table (or taxa count table). Each row
#'  corresponds to a sample and each column corresponds to one OTU (taxa).
#' @param x.test Covariate(s) of interest. This can be a vector, matrix, or data
#'  frame. Multi-level categorical variables should be coded as indicator
#'  variables before input. Use \code{x.test} together with \code{x.adj} when
#'  the default hypothesis matrix \code{Gamma} is desired.
#' @param x.adj Covariates that need to be adjusted. Requirements are the
#'  same as \code{x.test} above.  Use if default Gamma is desired.
#' @param x Optional design matrix containing all covariates.
#'  Default value is \code{NULL}. Only used
#'  if matrix \code{Gamma} is provided as input. It can be a vector, matrix,
#'  or dataframe. All K-level categorical variables should be converted to K-1
#'  indicator variables before input.
#' @param Gamma Optional hypothesis matrix specifying the linear combinations
#'  of regression coefficients to be tested.  Not used if data are entered as
#'  \code{x.test} and \code{x.adj}.
#' @param b Optional numeric vector specifying the right-hand side of the null
#'  hypothesis \eqn{H_0: \Gamma \beta_j = b}. Its length must equal the number
#'  of rows of \code{Gamma}. If \code{b = NULL}, CAFT uses the median-based
#'  null value estimated across taxa.
#' @param filter.thresh a real value between 0 and 1 for OTU table sample
#'  presence filtering. Any OTUs present in fewer than \code{filter.thresh}
#'  proportion of samples are filtered out. We set the default to be 0.05.
#' @param fdr.nominal the nominal false discover rate (FDR). The default is 0.2.
#' @param adjust.method a character string. Use multiple comparison/testing
#'  adjustment methods to control the family-wise error rate/false discover
#'  rate. Default to "\code{BH}". See \code{\link{p.adjust}} for the details.
#' @param regularize a logical value. If TRUE, adds a small penalty to the
#'   rank-based score to stabilize estimation and reduce small-sample
#'   bias.  Use this when the unpenalized rank equations have flat directions,
#'   there is near-separation, heavy censoring/ties, or the optimizer fails
#'   to converge. Default is TRUE.
#' @param test.method a character string. The variance estimator used to estimate
#'   the variance of the score equation. The methods of 'Cox', 'rank', and
#'   'martingale' are available. Default is 'rank'.
#' @param n.cores Integer. Number of CPU cores to use for parallel computation.
#'  Default is \code{1} (no parallelism). If \code{n.cores > 1}, the function
#'  runs tasks in parallel using \code{foreach}/\code{doParallel} with a PSOCK
#'  cluster. On typical desktops/laptops, a good choice is
#'  \code{max(1L, parallel::detectCores() - 1L)}.
#' @param boot.B Integer. Number of bootstrap replicates used for empirical
#'  calibration of taxon-level p-values. The default is \code{0L}, which
#'  disables bootstrap calibration. At least two bootstrap replicates are
#'  required for bootstrap calibration.
#' @param boot.parallel Logical. If \code{TRUE} and \code{boot.B >= 2}, the
#'  bootstrap loop is parallelized over bootstrap replicates. The default is
#'  \code{FALSE}.
#' @param boot.n.cores Integer. Number of CPU cores used for outer bootstrap
#'  parallelization when \code{boot.parallel = TRUE}. The default is
#'  \code{1L}.
#' @param boot.seed Optional integer seed for reproducible bootstrap resampling.
#'  The default is \code{NULL}.
#' @param boot.return.dist Logical. If \code{TRUE}, return the bootstrap-centered
#'  beta distribution and the names of taxa used in the bootstrap calibration.
#'  The default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, progress messages are shown during
#'   bootstrap computation. The default is \code{FALSE}.
#' @param return.mr.resid Logical; if \code{TRUE}, return subject-level
#'   residual contribution matrices from the observed CAFT fit. Default is
#'   \code{FALSE}. The returned residuals are the compensated Cox-type score
#'   contributions computed by \code{mySi.no.surv()} and evaluated at the
#'   restricted estimate used in the score test.
#'
#' @details
#' CAFT fits a rank-based accelerated failure time model to the negative
#' log-relative abundance of each taxon, treating zero counts as censored values
#' below the sample-specific detection limit. The resulting taxon-specific
#' regression coefficients are compared through compositional contrasts.
#'
#' With the default \code{x.test}/\code{x.adj} interface, the function tests the
#' covariates in \code{x.test} while adjusting for the covariates in
#' \code{x.adj}. Internally, this corresponds to a default hypothesis matrix
#' \code{Gamma} selecting the columns of \code{x.test}. For more general
#' hypotheses, users may provide a full design matrix \code{x}, a hypothesis
#' matrix \code{Gamma}, and optionally a null value \code{b}, corresponding to
#' \eqn{H_0: \Gamma \beta_j = b}. If \code{b = NULL}, the null value is estimated
#' as the median-based reference across taxa.
#'
#' When \code{boot.B >= 2}, CAFT performs fixed-design bootstrap calibration.
#' The observed model is first fit, and the bootstrap is carried out on the taxa
#' that pass the original filtering step. In each bootstrap replicate, samples
#' are resampled with replacement from the OTU table while the covariates are
#' kept fixed. The model is then refit with \code{filter.thresh = 0} on the same
#' set of taxa. Bootstrap p-values are computed from the centered bootstrap
#' distribution of the tested coefficient contrasts. For a single testing
#' variable, both a two-sided empirical beta-tilde p-value and a studentized
#' chi-square/Wald bootstrap p-value are returned. For multiple testing
#' variables, the bootstrap p-value is based on a multivariate statistic;
#' in this case \code{p.boot} and \code{p.boot.chi} are identical.
#'
#' If \code{boot.parallel = TRUE}, the outer bootstrap loop is parallelized over
#' bootstrap replicates. To avoid nested parallelism, the inner CAFT fitting
#' routine is run with \code{n.cores = 1L} inside each bootstrap worker.
#'
#' @return Return a list consisting of
#' \describe{
#' \item{beta.est}{A matrix that include the estimated coefficients by fitting
#'  each OTU count to the log-linear model. The number of columns should match
#'  the number of covariates from \code{x.test} and \code{x.adj}.}
#' \item{betahat.median}{The selected median of all betas, which was used in the
#'  restricted score test of the significance of differential abundant OTU. It can
#'  also be used to calculated the effect size of each OTU, defined as the difference
#'  from the median of all betas.}
#' \item{b.null}{The null value actually used in the restricted score test. This
#'  equals \code{betahat.median} when \code{b = NULL}, and equals the
#'  user-specified \code{b} otherwise.}
#' \item{rank.teststat}{test statistics from the proposed score test}
#' \item{rank.teststat.norm}{Numeric matrix with the normalized (z-score)
#'  rank-based test statistics for each taxon (rows) and each tested contrast
#'   (columns, in the order of \code{Gamma}'s rows). Under the null, entries are
#'   approximately \eqn{N(0,1)}; larger absolute values indicate stronger
#'   evidence against the null. Though the per-taxon p-values used for FDR control
#'   already be provided through \code{rank.teststat} and \code{p.otu}.}
#' \item{skip.otu}{the names of skipped OTU during taxa presence filtering
#'  above}
#' \item{p.otu}{p-values for individual OTU association tests}
#' \item{df.test}{Degrees of freedom per taxon in the test.}
#' \item{beta.est.r}{Matrix of constrained (restricted) estimates from the
#'  restricted score test.}
#' \item{p.detected.otu}{detected significantly differential abundant taxa
#'  (denoted by the column names of the OTU table) at the nominal FDR based on
#'  \code{p.otu}}
#' \item{q.otu}{q-values comes from adjusting p-values by the Bonferroni
#'  correction or the adjustment method specified in \code{adjust.method}. For
#'  individual OTU association tests}
#' \item{q.detected.otu}{detected significantly differential abundant taxa
#'  (denoted by the column names of the OTU table) at the nominal FDR
#'  based on \code{q.otu}}
#' \item{p.boot}{Bootstrap-calibrated empirical beta-tilde p-values. If
#'  \code{boot.B < 2}, this is a named vector of \code{NA} values.}
#' \item{q.boot}{Multiplicity-adjusted values from \code{p.boot}.}
#' \item{p.boot.detected.otu}{Taxa detected using \code{p.boot < fdr.nominal}.
#'  Empty if bootstrap calibration is not run.}
#' \item{q.boot.detected.otu}{Taxa detected using
#' \code{q.boot < fdr.nominal}. Empty if bootstrap calibration is not run.}
#' \item{p.boot.chi}{Studentized chi-square/Wald bootstrap p-values. For
#'  multiple testing variables, these are the same as \code{p.boot}. If
#'  \code{boot.B < 2}, this is a named vector of \code{NA} values.}
#' \item{q.boot.chi}{Multiplicity-adjusted values from \code{p.boot.chi}.}
#' \item{p.boot.chi.detected.otu}{Taxa detected using
#'  \code{p.boot.chi < fdr.nominal}. Empty if bootstrap calibration is not run.}
#' \item{q.boot.chi.detected.otu}{Taxa detected using
#' \code{q.boot.chi < fdr.nominal}. Empty if bootstrap calibration is not run.}
#' \item{boot.median}{Matrix of bootstrap null values used to center the
#'  bootstrap coefficient estimates. When \code{b} is user-specified, these
#'  values equal the specified null value across bootstrap replicates.}
#' \item{boot.n.success}{Number of bootstrap replicates that produced complete
#'  bootstrap null values.}
#' \item{boot.n.fail}{Number of bootstrap replicates skipped or failed.}
#' \item{boot.beta.tilde}{Returned only when \code{boot.return.dist = TRUE}.
#'  Array of bootstrap-centered tested coefficient contrasts.}
#' \item{boot.taxa}{Returned only when \code{boot.return.dist = TRUE}. Names of
#'  taxa used in the bootstrap calibration.}
#' \item{mr.resid}{Returned only when \code{return.mr.resid = TRUE}. A list of
#'  length equal to the number of taxa. Each element is an \eqn{n \times p}
#'  matrix of subject-level compensated Cox-type residual (Martingale residual)
#'  contributions, evaluated at the restricted estimate used in the score test.}
#' }
#'
#' @import stats
#' @import graphics
#' @import phyloseq
#' @import MASS
#' @importFrom MASS ginv Null
#' @import ICSNP
#' @import foreach
#' @import doParallel
#' @importFrom doRNG %dorng%
#' @export
#'
#' @examples
#' \donttest{
#' library(CAFT)
#' data(Colon)
#' library(phyloseq)
#'
#' count.tab <- t(as.data.frame(as.matrix(otu_table(Colon))))
#' sample.tab <- as.data.frame(as.matrix(sample_data(Colon)))
#'
#' pNA <- which(is.na(sample.tab$age))
#' if (length(pNA) > 0) {
#'   count.tab <- count.tab[-pNA, , drop = FALSE]
#'   sample.tab <- sample.tab[-pNA, , drop = FALSE]
#' }
#'
#' p_otu <- which(rowSums(t(count.tab) > 0) > 1)
#' count.tab <- count.tab[, p_otu, drop = FALSE]
#'
#' Disease1 <- Disease2 <- rep(0, NROW(sample.tab))
#' Disease1[sample.tab$disease == "CRC"] <- 1
#' Disease2[sample.tab$disease == "adenoma"] <- 1
#'
#' Age <- as.numeric(sample.tab$age)
#' Gender <- as.numeric(factor(sample.tab$gender)) - 1
#'
#' x.test <- cbind(CRC = Disease1, adenoma = Disease2)
#' x.adj <- cbind(Age = Age, Gender = Gender)
#'
#' ## Default CAFT analysis
#' res.CAFT <- caft(
#'   otu.table = count.tab,
#'   x.test = x.test,
#'   x.adj = x.adj
#' )
#'
#' ## Bootstrap-calibrated CAFT analysis; use a small B for examples only.
#' res.CAFT.boot <- caft(
#'   otu.table = count.tab,
#'   x.test = x.test,
#'   x.adj = x.adj,
#'   boot.B = 10,
#'   boot.parallel = FALSE,
#'   boot.seed = 1
#' )
#'
#' ## User-specified linear hypothesis using x, Gamma, and b.
#' x.all <- cbind(CRC = Disease1, adenoma = Disease2, Age = Age, Gender = Gender)
#' Gamma <- matrix(c(1, -1, 0, 0), nrow = 1)
#'
#' res.CAFT.gamma <- caft(
#'   otu.table = count.tab,
#'   x = x.all,
#'   Gamma = Gamma,
#'   b = 0
#' )
#' }
caft <- function(otu.table, x.test = NULL, x.adj = NULL, x = NULL, Gamma = NULL,
                 b = NULL,
                 filter.thresh = 0.05, fdr.nominal = 0.20, adjust.method = "BH",
                 regularize=TRUE, test.method="rank", n.cores=1L,
                 boot.B = 0L,
                 boot.parallel = FALSE,
                 boot.n.cores = 1L,
                 boot.seed = NULL,
                 boot.return.dist = FALSE,
                 verbose = FALSE,
                 return.mr.resid = FALSE) {

  if (!is.null(boot.seed)) set.seed(boot.seed)

  if (!is.matrix(otu.table)) {otu.table <- as.matrix(otu.table)}

  using.xtest <- is.null(x)
  if (using.xtest && is.null(x.test)) {
    stop("No data provided: either x.test, or x together with Gamma, is required.")
  }
  if (!using.xtest && !is.null(x.test)) {
    stop("Only one of x or x.test can be specified, not both.")
  }
  if (using.xtest && !is.null(Gamma)) {
    stop("Gamma should not be specified when using the x.test/x.adj interface.")
  }
  if (!using.xtest && is.null(Gamma)) {
    stop("Gamma must be specified when x is supplied directly.")
  }
  if (using.xtest) {x <- cbind(x.test, x.adj)}
  if (!is.matrix(x)) {x <- as.matrix(x)}
  if (NROW(otu.table) != NROW(x)) {
    stop(" Number of samples not match between OTU table and covairates matrix!")
  }
  if (!is.null(Gamma)) {
    if (!is.matrix(Gamma)) {
      Gamma <- matrix(Gamma, nrow = 1)
    }
    Gamma <- as.matrix(Gamma)
    storage.mode(Gamma) <- "numeric"
    if (ncol(Gamma) != NCOL(x)) {
      stop("Number of columns in Gamma must match the number of columns in x.")
    }
    Gamma.rank <- qr(Gamma)$rank
    if (Gamma.rank != nrow(Gamma)) {
      stop("Gamma matrix provided does not have full row rank.")
    }
  }

  # missing values
  pNA.x <- which(apply(x, 1, function(x) any(is.na(x)))) # covariates
  pNA.otu <- which(apply(otu.table, 1, function(x) any(is.na(x)))) # covariates
  pNA <- unique(c(pNA.x, pNA.otu))
  if (length(pNA) > 0) {
    warning("Missing values deleted!")
    otu.table <- otu.table[-pNA, ]
    x <- x[-pNA, ]
  }

  n.param <- NCOL(x)
  # set up Gamma and Lambda matrices
  if (is.null(Gamma)) {
    n.test <- NCOL(x.test)
    Gamma <- diag(n.param)[1:n.test, , drop = FALSE]
    if (n.test < n.param) {
      Lambda <- diag(n.param)[(n.test + 1):n.param, , drop = FALSE]
    } else {
      Lambda <- NULL
    }
  } else {
    n.test <- nrow(Gamma)
    if (n.test < n.param) {
      Lambda <- t(MASS::Null(t(Gamma)))
    } else {
      Lambda <- NULL
    }
  }
  Gamma.ginv <- MASS::ginv(Gamma)
  Lambda.ginv <- if (is.null(Lambda)) NULL else MASS::ginv(Lambda)

  if (is.null(b)) {
    b.user <- NULL
  } else {
    b.user <- as.numeric(b)
    if (length(b.user) != n.test) {
      stop("Length of b must equal the number of rows of Gamma!")
    }
  }

  x.test.use <- x %*% t(Gamma)
  x.test.use <- as.matrix(x.test.use)

  if (is.null(Lambda)) {
    x.adj.use <- NULL
  } else {
    x.adj.use <- x %*% t(Lambda)
    x.adj.use <- as.matrix(x.adj.use)
  }

  if (boot.B >=2L && !is.null(x.adj.use)) {
    X.adj <- cbind(1, x.adj.use)
    x.test.use <- residuals(lm.fit(x = X.adj, y = x.test.use))
    x.test.use <- as.matrix(x.test.use)
    if (NCOL(x.test.use) != n.test) {
      x.test.use <- matrix(x.test.use, ncol = n.test)
    }
  }

  ## Matrix to rebuild x from tested + nuisance coordinates
  if (is.null(Lambda)) {
    A <- Gamma
  } else {
    A <- rbind(Gamma, Lambda)
  }

  outer.parallel <- boot.B >= 2L && boot.parallel && boot.n.cores > 1L

  if (outer.parallel && n.cores > 1L) {
    warning("Nested parallelism detected: setting caft_fit(n.cores = 1) inside outer bootstrap workers.")
    fit.n.cores <- 1L
  } else {
    fit.n.cores <- n.cores
  }

  if (is.null(x.adj.use)) {
    x.obs <- x.test.use
  } else {
    x.obs <- cbind(x.test.use, x.adj.use)
  }
  x.obs <- x.obs %*% solve(t(A))

  fit.obs <- caft_fit(
    otu.table      = otu.table,
    x              = x.obs,
    Gamma          = Gamma,
    Gamma.ginv     = Gamma.ginv,
    Lambda         = Lambda,
    Lambda.ginv    = Lambda.ginv,
    b              = b.user,
    filter.thresh  = filter.thresh,
    fdr.nominal    = fdr.nominal,
    adjust.method  = adjust.method,
    regularize     = regularize,
    test.method    = test.method,
    n.cores        = fit.n.cores,
    return.mr.resid = return.mr.resid
  )

  empty_boot_results <- function(fit.obs, boot.B) {
    n.otu <- length(fit.obs$p.otu)

    otu.names <- names(fit.obs$p.otu)
    if (is.null(otu.names) || length(otu.names) != n.otu) {
      otu.names <- colnames(otu.table)
    }

    fit.obs$p.boot <- setNames(rep(NA_real_, n.otu), otu.names)
    fit.obs$q.boot <- setNames(rep(NA_real_, n.otu), otu.names)
    fit.obs$p.boot.chi <- setNames(rep(NA_real_, n.otu), otu.names)
    fit.obs$q.boot.chi <- setNames(rep(NA_real_, n.otu), otu.names)

    fit.obs$p.boot.detected.otu <- character(0)
    fit.obs$q.boot.detected.otu <- character(0)
    fit.obs$p.boot.chi.detected.otu <- character(0)
    fit.obs$q.boot.chi.detected.otu <- character(0)

    fit.obs$boot.median <- matrix(
      NA_real_,
      nrow = 0L,
      ncol = n.test,
      dimnames = list(NULL, paste0("test", seq_len(n.test)))
    )

    fit.obs$boot.n.success <- 0L
    fit.obs$boot.n.fail <- as.integer(max(boot.B, 0L))

    fit.obs
  }

  if (boot.B < 2L) {
    if (boot.B == 1L) {
      warning("At least two bootstrap replicates are required; bootstrap calibration skipped.")
    }
    return(empty_boot_results(fit.obs, boot.B))
  }

  ## Null bootstrap (using filtered taxa)
  use <- !is.na(fit.obs$p.otu)

  if (sum(use) == 0L) {
    warning("No taxa passed filtering; bootstrap calibration skipped.")
    return(empty_boot_results(fit.obs, boot.B))
  }

  n.taxa.boot <- sum(use)

  Gamma.use <- as.matrix(Gamma)
  storage.mode(Gamma.use) <- "numeric"

  beta.obs.all <- as.matrix(fit.obs$beta.est)
  storage.mode(beta.obs.all) <- "numeric"

  beta.obs <- beta.obs.all[use, , drop = FALSE] %*% t(Gamma.use)

  beta.null.obs <- as.numeric(fit.obs$b.null[seq_len(n.test)])
  beta.tilde.obs <- sweep(beta.obs, 2, beta.null.obs, "-")

  beta.boot <- array(NA_real_, dim = c(boot.B, n.taxa.boot, n.test))
  boot.median <- matrix(NA_real_, nrow = boot.B, ncol = n.test)
  colnames(boot.median) <- paste0("test", seq_len(n.test))

  otu.boot.base <- otu.table[, use, drop = FALSE]

  if (outer.parallel) {
    # cl <- parallel::makeCluster(boot.n.cores)
    # doParallel::registerDoParallel(cl)
    # on.exit({
    #   try(parallel::stopCluster(cl), silent = TRUE)
    #   foreach::registerDoSEQ()
    # }, add = TRUE)

    # doParallel::registerDoParallel(cores = boot.n.cores)
    # on.exit(foreach::registerDoSEQ(), add = TRUE)

    if (.Platform$OS.type == "windows") {

      ## Windows: PSOCK cluster
      cl <- parallel::makeCluster(boot.n.cores)
      doParallel::registerDoParallel(cl)

      on.exit({
        try(parallel::stopCluster(cl), silent = TRUE)
        foreach::registerDoSEQ()
      }, add = TRUE)

    } else {

      ## macOS / Linux: fork-style backend
      doParallel::registerDoParallel(cores = boot.n.cores)
      on.exit(foreach::registerDoSEQ(), add = TRUE)

    }

    rng.seed <- if (is.null(boot.seed)) {
      sample.int(.Machine$integer.max, 1L)
    } else {
      boot.seed
    }

    boot.res <- foreach::foreach(
      b = seq_len(boot.B),
      .combine = "c",
      .inorder = TRUE,
      .errorhandling = "pass",
      .packages = c("CAFT", "MASS"),
      .export = c("caft_fit"),
      .options.RNG = rng.seed
    ) %dorng% {

      idx <- sample.int(nrow(otu.table), size = nrow(otu.table), replace = TRUE)
      otu.boot <- otu.boot.base[idx, , drop = FALSE]

      fit.b <- tryCatch(
        caft_fit(
          otu.table       = otu.boot,
          x               = x.obs,
          Gamma           = Gamma,
          Gamma.ginv      = Gamma.ginv,
          Lambda          = Lambda,
          Lambda.ginv     = Lambda.ginv,
          b               = b.user,
          filter.thresh   = 0,
          fdr.nominal     = fdr.nominal,
          adjust.method   = adjust.method,
          regularize      = regularize,
          test.method     = test.method,
          n.cores         = 1L,
          return.mr.resid = FALSE
        ),
        error = function(e) e
      )

      if (inherits(fit.b, "error")) {
        list(list(
          ok = FALSE,
          beta = matrix(NA_real_, nrow = n.taxa.boot, ncol = n.test),
          median = rep(NA_real_, n.test),
          msg = conditionMessage(fit.b)
        ))
      } else {
        beta.b.all <- as.matrix(fit.b$beta.est)
        storage.mode(beta.b.all) <- "numeric"

        beta.b <- beta.b.all[, , drop = FALSE] %*% t(Gamma.use)

        list(list(
          ok = TRUE,
          beta = beta.b,
          median = as.numeric(fit.b$b.null[seq_len(n.test)]),
          msg = NA_character_
        ))
      }
    }

  } else {

    boot.res <- vector("list", boot.B)

    for (b in seq_len(boot.B)) {
      if (isTRUE(verbose) && boot.B >= 20L && b %% 20L == 0L) {
        message("Bootstrap replicate ", b, " of ", boot.B)
      }

      idx <- sample.int(nrow(otu.table), size = nrow(otu.table), replace = TRUE)
      otu.boot <- otu.boot.base[idx, , drop = FALSE]

      fit.b <- tryCatch(
        caft_fit(
          otu.table       = otu.boot,
          x               = x.obs,
          Gamma           = Gamma,
          Gamma.ginv      = Gamma.ginv,
          Lambda          = Lambda,
          Lambda.ginv     = Lambda.ginv,
          b               = b.user,
          filter.thresh   = 0,
          fdr.nominal     = fdr.nominal,
          adjust.method   = adjust.method,
          regularize      = regularize,
          test.method     = test.method,
          n.cores         = 1L,
          return.mr.resid = FALSE
        ),
        error = function(e) e
      )

      if (inherits(fit.b, "error")) {
        boot.res[[b]] <- list(
          ok = FALSE,
          beta = matrix(NA_real_, nrow = n.taxa.boot, ncol = n.test),
          median = rep(NA_real_, n.test),
          msg = conditionMessage(fit.b)
        )
      } else {
        beta.b.all <- as.matrix(fit.b$beta.est)
        storage.mode(beta.b.all) <- "numeric"

        beta.b <- beta.b.all[, , drop = FALSE] %*% t(Gamma.use)

        boot.res[[b]] <- list(
          ok = TRUE,
          beta = beta.b,
          median = as.numeric(fit.b$b.null[seq_len(n.test)]),
          msg = NA_character_
        )
      }
    }
  }

  for (b in seq_len(boot.B)) {
    if (inherits(boot.res[[b]], "error")) {
      warning("Bootstrap replicate ", b, " failed and was skipped: ", conditionMessage(boot.res[[b]]))
      next
    }
    if (isTRUE(boot.res[[b]]$ok)) {
      beta.boot[b, , ] <- boot.res[[b]]$beta
      boot.median[b, ] <- boot.res[[b]]$median
    } else {
      warning( "Bootstrap replicate ", b, " failed and was skipped: ", boot.res[[b]]$msg)
    }
  }

  ## Center each bootstrap beta by its null median
  beta.boot.tilde <- array(NA_real_, dim = c(boot.B, n.taxa.boot, n.test))
  for (k in seq_len(n.test)) {
    beta.boot.tilde[, , k] <- sweep(beta.boot[, , k, drop = TRUE], 1, boot.median[, k], "-")
  }

  ########################
  ## Bootstrap p-values ##
  ########################

  #   compute bootstrap p-value based on CI
  p.boot.short <- rep(NA_real_, n.taxa.boot)
  if (n.test == 1L) {
    beta.boot.tilde.mat <- beta.boot.tilde[, , 1, drop = TRUE]
    beta.tilde.obs.vec <- as.numeric(beta.tilde.obs[, 1])

    for (j in seq_len(n.taxa.boot)) {
      good <- is.finite(beta.boot.tilde.mat[, j]) & is.finite(beta.tilde.obs.vec[j])
      if (sum(good) > 0L) {
        p.right <- sum(beta.boot.tilde.mat[good, j] >= beta.tilde.obs.vec[j])
        p.left  <- sum(beta.boot.tilde.mat[good, j] <= beta.tilde.obs.vec[j])

        p.boot.short[j] <- min(1, (1 + 2 * min(p.right, p.left)) / (1 + sum(good))
        )
      }
    }
  } else {
    ## multiple variables in x.test
    for (j in seq_len(n.taxa.boot)) {

      boot.mat <- matrix(beta.boot.tilde[, j, ], nrow = boot.B, ncol = n.test)
      obs.vec <- as.numeric(beta.tilde.obs[j, ])

      good <- stats::complete.cases(boot.mat) & all(is.finite(obs.vec))

      if (sum(good) > n.test) {

        boot.good <- boot.mat[good, , drop = FALSE]
        mu.j <- colMeans(boot.good)
        boot.centered <- sweep(boot.good, 2, mu.j, "-")
        obs.centered <- obs.vec - mu.j

        Sigma.j <- stats::cov(boot.centered)
        Sigma.inv.j <- MASS::ginv(Sigma.j)

        T.obs.j <- as.numeric(t(obs.centered) %*% Sigma.inv.j %*% obs.centered)
        T.boot.j <- rowSums((boot.centered %*% Sigma.inv.j) * boot.centered)

        p.boot.short[j] <- (1 + sum(T.boot.j >= T.obs.j, na.rm = TRUE)) / (1 + sum(is.finite(T.boot.j)))
      }
    }
  }

  #   compute bootstrap Wald test
  p.boot.chi.short <- rep(NA_real_, n.taxa.boot)
  if (n.test == 1L) {
    beta.boot.tilde.mat <- beta.boot.tilde[, , 1, drop = TRUE]
    beta.tilde.obs.vec <- as.numeric(beta.tilde.obs[, 1])
    mu <- colMeans(beta.boot.tilde.mat, na.rm = TRUE)
    boot.centered <- sweep(beta.boot.tilde.mat, 2, mu, "-")
    obs.centered <- beta.tilde.obs.vec - mu
    sigma <- sqrt(colMeans(boot.centered^2, na.rm = TRUE))
    null.Z <- sweep(boot.centered, 2, sigma, "/")
    Z <- obs.centered / sigma

    for (j in seq_len(n.taxa.boot)) {
      good <- is.finite(null.Z[, j]) & is.finite(Z[j])
      if (sum(good) > 0L) {
        p.boot.chi.short[j] <- (1 + sum(null.Z[good, j]^2 >= Z[j]^2)) / (1 + sum(good))
      }
    }
  } else {
    ## For multiple testing variables, p.boot.short is already based on
    ## the multivariate studentized quadratic/Wald statistic.
    p.boot.chi.short <- p.boot.short
  }

  p.boot <- rep(NA_real_, length(fit.obs$p.otu))
  p.boot[use] <- p.boot.short
  q.boot <- rep(NA_real_, length(fit.obs$p.otu))
  q.boot[use] <- p.adjust(p.boot[use], method = adjust.method)
  names(p.boot) <- colnames(otu.table)
  names(q.boot) <- colnames(otu.table)
  fit.obs$p.boot <- p.boot
  fit.obs$q.boot <- q.boot

  p.boot.chi <- rep(NA_real_, length(fit.obs$p.otu))
  p.boot.chi[use] <- p.boot.chi.short
  q.boot.chi <- rep(NA_real_, length(fit.obs$p.otu))
  q.boot.chi[use] <- p.adjust(p.boot.chi[use], method = adjust.method)
  names(p.boot.chi) <- colnames(otu.table)
  names(q.boot.chi) <- colnames(otu.table)
  fit.obs$p.boot.chi <- p.boot.chi
  fit.obs$q.boot.chi <- q.boot.chi

  fit.obs$p.boot.detected.otu <- colnames(otu.table)[which(p.boot < fdr.nominal)]
  fit.obs$q.boot.detected.otu <- colnames(otu.table)[which(q.boot < fdr.nominal)]
  fit.obs$p.boot.chi.detected.otu <- colnames(otu.table)[which(p.boot.chi < fdr.nominal)]
  fit.obs$q.boot.chi.detected.otu <- colnames(otu.table)[which(q.boot.chi < fdr.nominal)]

  fit.obs$boot.median <- boot.median
  fit.obs$boot.n.success <- sum(stats::complete.cases(boot.median))
  fit.obs$boot.n.fail <- boot.B - fit.obs$boot.n.success

  if (boot.return.dist) {
    fit.obs$boot.beta.tilde <- beta.boot.tilde
    fit.obs$boot.taxa <- colnames(otu.table)[use]
  }

  return(fit.obs)

}
