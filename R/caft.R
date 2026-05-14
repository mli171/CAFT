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
#' In addition to the original asymptotic CAFT inference, this version optionally
#' provides permutation-based calibration of taxon-level p-values. When
#' permutation is requested, the tested covariate is first residualized against
#' the adjustment covariates, and empirical p-values are computed from permuted
#' test statistics.
#'
#' @param otu.table the community OTU table (or taxa count table). Each row
#'  corresponds to a sample and each column corresponds to one OTU (taxa).
#' @param x.test the covariates are of interest. It can be a vector, matrix, or
#'  data frame. All K-level categorical variables should be converted to K-1
#'  indicator variables before input.  Use if default Gamma is desired.
#' @param x.adj other covariates that need to be adjusted. Requirements are the
#'  same as \code{x.test} above.  Use if default Gamma is desired.
#' @param Gamma matrix to specify what hypothesis to test.  Not used if data are
#' entered as \code{x.test} and \code{x.adj}.
#' @param x data for all covariates.  Default value is \code{NULL}. Only used
#'  if matrix \code{Gamma} is provided as input. It can be a vector, matrix,
#'  or dataframe. All K-level categorical variables should be converted to K-1
#'  indicator variables before input.
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
#' @param perm.B Integer. Number of permutations used for empirical calibration
#'   of taxon-level p-values. Default is \code{0L}, which disables permutation
#'   and returns the original CAFT asymptotic results only.
#' @param perm.parallel Logical; if \code{TRUE} and \code{perm.B > 0}, the
#'   permutation loop is parallelized across outer workers. Default is
#'   \code{FALSE}.
#' @param perm.n.cores Integer. Number of CPU cores used for outer permutation
#'   parallelization when \code{perm.parallel = TRUE}. Default is \code{1L}.
#' @param perm.seed Optional integer random seed for reproducible permutations.
#'   Default is \code{NULL}.
#' @param return.mr.resid Logical; if \code{TRUE}, return subject-level
#'   residual contribution matrices from the observed CAFT fit. Default is
#'   \code{FALSE}. The returned residuals are the compensated Cox-type score
#'   contributions computed by \code{mySi.no.surv()} and evaluated at the
#'   restricted estimate used in the score test.
#'
#' @details
#' If \code{perm.B > 0}, the function first fits the observed CAFT model using
#' the residualized tested covariate when adjustment covariates are present, and
#' otherwise using the original tested covariate. It then permutes the rows of
#' the tested covariate, refits the model for each permutation, and computes
#' empirical p-values by comparing the observed test statistic to its
#' permutation distribution. Because permutation calibration can be
#' computationally expensive, it is recommended primarily when the number of
#' OTUs is not too large, for example fewer than 50.
#'
#' When \code{perm.parallel = TRUE}, outer permutation workers are used. To avoid
#' nested parallelism, the internal CAFT fitting routine is forced to run in
#' single-thread mode inside each permutation worker.
#'
#' Permutation mode is intended for the \code{x.test}/\code{x.adj} interface.
#' If \code{x} and \code{Gamma} are supplied directly, users should ensure that
#' the tested component being permuted is well defined.
#' @return Return a list consisting of
#' \describe{
#' \item{beta.est}{A matrix that include the estimated coefficients by fitting
#'  each OTU count to the log-linear model. The number of columns should match
#'  the number of covariates from \code{x.test} and \code{x.adj}.}
#' \item{betahat.median}{The selected median of all betas, which was used in the
#'  restricted score test of the significance of differential abundant OTU. It can
#'  also be used to calculated the effect size of each OTU, defined as the difference
#'  from the median of all betas.}
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
#' \item{MR.resid}{Returned only when \code{return.mr.resid = TRUE}. A list of
#'  length equal to the number of taxa. Each element is an \eqn{n \times p}
#'  matrix of subject-level compensated Cox-type residual (Martingale residual)
#'  contributions, evaluated at the restricted estimate used in the score test.}
#' \item{p.perm}{Empirical taxon-level p-values from the permutation
#'   distribution. Returned only when \code{perm.B > 0}.}
#' \item{p.perm.detected.otu}{detected significantly differential abundant taxa
#'  (denoted by the column names of the OTU table) at the nominal FDR based on
#'  \code{p.perm}}
#' \item{q.perm}{Multiplicity-adjusted empirical p-values obtained from
#'   \code{p.adjust(p.perm, method = adjust.method)}. Returned only when
#'   \code{perm.B > 0}.}
#' \item{q.perm.detected.otu}{detected significantly differential abundant taxa
#'  (denoted by the column names of the OTU table) at the nominal FDR based on
#'  \code{q.perm}.}
#' }
#' @import stats
#' @import graphics
#' @import phyloseq
#' @import MASS
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
#'
#' library(phyloseq)
#'
#' count.tab <- t(as.data.frame(as.matrix(otu_table(Colon))))
#' sample.tab <- as.data.frame(as.matrix(sample_data(Colon)))
#'
#' pNA <- which(is.na(sample.tab$age))
#' if (length(pNA) > 0) {
#'   count.tab <- count.tab[-pNA, ]
#'   sample.tab <- sample.tab[-pNA, ]
#' }
#'
#' p_otu <- which(rowSums(t(count.tab) > 0) > 1)
#' count.tab <- count.tab[, p_otu]
#'
#' Disease1 <- Disease2 <- rep(0, NROW(sample.tab))
#' Disease1[sample.tab$disease == "CRC"] <- 1
#' Disease2[sample.tab$disease == "adenoma"] <- 1
#'
#' Age <- as.numeric(sample.tab$age)
#' Gender <- as.numeric(factor(sample.tab$gender)) - 1
#'
#' x.test <- cbind(Disease1, Disease2)
#' x.adj  <- cbind(Age, Gender)
#'
#' ## CAFT analysis
#' res.CAFT <- caft(otu.table = count.tab, x.test = x.test, x.adj = x.adj)
#'
#' ## permutation-calibrated CAFT analysis
#'
#' res.CAFT.perm <- caft(
#'   otu.table = count.tab,
#'   x.test = x.test,
#'   x.adj = x.adj,
#'   perm.B = 1000,
#'   perm.parallel = FALSE,
#'   perm.seed = 1
#' )
#' }
caft <- function(otu.table, x.test = NULL, x.adj = NULL, x = NULL, Gamma = NULL,
                 filter.thresh = 0.05, fdr.nominal = 0.20, adjust.method = "BH",
                 regularize=TRUE, test.method="rank", n.cores=1L,
                 perm.B = 0L,
                 perm.parallel = FALSE,
                 perm.n.cores = 1L,
                 perm.seed = NULL,
                 return.mr.resid = FALSE) {

  if (!is.null(perm.seed)) set.seed(perm.seed)

  if (!is.matrix(otu.table)) {
    otu.table <- as.matrix(otu.table)
  }
  if (is.null(x) & is.null(x.test)) {
    stop("No data provided: Either x or x.test (and possibly x.adj) required")
  }
  if (!is.null(x) & !is.null(x.test)) {
    stop("Only one of x or x.test can be specified, not both")
  }
  if (is.null(x)) {
    x <- cbind(x.test, x.adj)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (NROW(otu.table) != NROW(x)) {
    stop(" Number of samples not match between OTU table and covairates matrix!")
  }
  if (!is.null(Gamma)) {
    if (is.null(x)) {
      stop("Gamma should not be specified if x.test (and possibly x.adj) is specified")
    }
    if (!("matrix" %in% class(Gamma))) Gamma <- matrix(Gamma, nrow = 1)
    Gamma.rank <- qr(Gamma)$rank
    if (Gamma.rank != nrow(Gamma)) {
      stop("Gamma matrix provided does not have full row rank")
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

  x.test.use <- x %*% t(Gamma)
  x.test.use <- as.matrix(x.test.use)

  if (is.null(Lambda)) {
    x.adj.use <- NULL
  } else {
    x.adj.use <- x %*% t(Lambda)
    x.adj.use <- as.matrix(x.adj.use)
  }

  if (perm.B > 0L && !is.null(x.adj.use)) {
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

  outer.parallel <- perm.B > 0L && perm.parallel && perm.n.cores > 1L

  if (outer.parallel && n.cores > 1L) {
    warning("Nested parallelism detected: \n Setting caft_fit(n.cores = 1) inside outer permutation workers.")
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
    filter.thresh  = filter.thresh,
    fdr.nominal    = fdr.nominal,
    adjust.method  = adjust.method,
    regularize     = regularize,
    test.method    = test.method,
    n.cores        = fit.n.cores,
    return.mr.resid = return.mr.resid
  )

  if (perm.B <= 0L) {
    fit.obs$p.perm <- NA_real_
    fit.obs$q.perm <- NA_real_
    fit.obs$p.perm.detected.otu <- NA_character_
    fit.obs$q.perm.detected.otu <- NA_character_
    return(fit.obs)
  }

  ## Permutation need

  tested <- !is.na(fit.obs$p.otu)

  if (n.test == 1L) {
    T.obs <- drop(as.matrix(fit.obs$rank.teststat.norm))
  } else {
    T.obs <- fit.obs$rank.teststat
  }

  if (outer.parallel) {
    doParallel::registerDoParallel(cores = perm.n.cores)
    on.exit(foreach::registerDoSEQ(), add = TRUE)

    T.perm <- foreach::foreach(
      b = seq_len(perm.B),
      .combine = "rbind",
      .errorhandling = "pass",
      .packages = "CAFT",
      .export = c("caft_fit"),
      .options.RNG = perm.seed
    ) %dorng% {

      idx <- sample.int(nrow(x.test.use))
      x.test.perm <- x.test.use[idx, , drop = FALSE]

      if (is.null(x.adj.use)) {
        x.perm <- x.test.perm
      } else {
        x.perm <- cbind(x.test.perm, x.adj.use)
      }
      x.perm <- x.perm %*% solve(t(A))

      fit.b <- caft_fit(
        otu.table      = otu.table,
        x              = x.perm,
        Gamma          = Gamma,
        Gamma.ginv     = Gamma.ginv,
        Lambda         = Lambda,
        Lambda.ginv    = Lambda.ginv,
        filter.thresh  = filter.thresh,
        fdr.nominal    = fdr.nominal,
        adjust.method  = adjust.method,
        regularize     = regularize,
        test.method    = test.method,
        n.cores        = 1L,
        return.mr.resid = FALSE
      )

      if (n.test == 1L) {
        drop(as.matrix(fit.b$rank.teststat.norm))
      } else {
        fit.b$rank.teststat
      }
    }

  } else {
    T.perm <- matrix(NA_real_, nrow = perm.B, ncol = length(T.obs))

    for (b in seq_len(perm.B)) {
      idx <- sample.int(nrow(x.test.use))
      x.test.perm <- x.test.use[idx, , drop = FALSE]

      if (is.null(x.adj.use)) {
        x.perm <- x.test.perm
      } else {
        x.perm <- cbind(x.test.perm, x.adj.use)
      }
      x.perm <- x.perm %*% solve(t(A))

      fit.b <- caft_fit(
        otu.table      = otu.table,
        x              = x.perm,
        Gamma          = Gamma,
        Gamma.ginv     = Gamma.ginv,
        Lambda         = Lambda,
        Lambda.ginv    = Lambda.ginv,
        filter.thresh  = filter.thresh,
        fdr.nominal    = fdr.nominal,
        adjust.method  = adjust.method,
        regularize     = regularize,
        test.method    = test.method,
        n.cores        = 1L,
        return.mr.resid = FALSE
      )

      if (n.test == 1L) {
        T.perm[b, ] <- drop(as.matrix(fit.b$rank.teststat.norm))
      } else {
        T.perm[b, ] <- fit.b$rank.teststat
      }
    }
  }

  p.perm <- rep(NA_real_, length(T.obs))

  if (n.test == 1L) {
    p.perm[tested] <- sapply(which(tested), function(j) {
      good <- !is.na(T.perm[, j]) & !is.na(T.obs[j])
      (1 + sum(abs(T.perm[good, j]) >= abs(T.obs[j]))) / (1 + sum(good))
    })
  } else {
    p.perm[tested] <- sapply(which(tested), function(j) {
      good <- !is.na(T.perm[, j]) & !is.na(T.obs[j])
      (1 + sum(T.perm[good, j] >= T.obs[j])) / (1 + sum(good))
    })
  }

  q.perm <- rep(NA_real_, length(T.obs))
  q.perm[tested] <- p.adjust(p.perm[tested], method = adjust.method)

  fit.obs$p.perm <- p.perm
  fit.obs$p.perm.detected.otu <- colnames(otu.table)[which(p.perm < fdr.nominal)]
  fit.obs$q.perm <- q.perm
  fit.obs$q.perm.detected.otu <- colnames(otu.table)[which(q.perm < fdr.nominal)]

  return(fit.obs)

}
