#' CAFT Restricted Score Test
#'
#' @description
#' Performs the CAFT restricted score test using an existing unrestricted
#' estimation object returned by \code{caft_estimate()}. This function separates
#' testing from unrestricted estimation: the unrestricted taxon-level
#' coefficient estimates are reused from \code{caft_estimate()}, and only the restricted
#' estimation and score testing steps are carried out for the specified
#' hypothesis.
#'
#' @param est An object of class \code{"caft_est"} returned by
#'   \code{caft_estimate()}.
#' @param Gamma Hypothesis matrix specifying the linear combinations
#'   of regression coefficients to be tested. Each row of \code{Gamma}
#'   corresponds to one tested contrast, and the number of columns must equal
#'   the number of covariates in \code{est$x}. This argument is required.
#' @param b Optional numeric vector specifying the right-hand side of the null
#'   hypothesis \eqn{H_{0j}: \Gamma \beta_j = b} for each taxon \eqn{j}. Its
#'   length must equal the number of rows of \code{Gamma}. If \code{b = NULL},
#'   CAFT uses a median-based reference value estimated from the unrestricted
#'   taxon-level contrasts \eqn{\Gamma \hat\beta_j}.
#' @param fdr.nominal The nominal threshold used to report detected taxa. The
#'   default is \code{0.20}.
#' @param adjust.method A character string specifying the p-value adjustment
#'   method used to compute \code{q.otu}. The default is \code{"BH"}. See
#'   \code{\link[stats]{p.adjust}}.
#' @param regularize A logical value. If \code{TRUE}, adds a small penalty to
#'   the rank-based score to improve numerical stability in the restricted
#'   estimation step. The default uses the value stored in \code{est$regularize}.
#' @param test.method A character string specifying the variance estimator used
#'   in the rank-based score test. The available options are \code{"rank"},
#'   \code{"Cox"}, and \code{"martingale"}. The default is \code{"rank"}.
#' @param n.cores Integer. Number of CPU cores to use for taxon-level parallel
#'   computation in the restricted testing step. The default is \code{1}.
#' @param return.mr.resid Logical. If \code{TRUE}, return subject-level residual
#'   contribution matrices evaluated at the restricted estimates. The default is
#'   \code{FALSE}.
#'
#' @details
#' Let \eqn{\beta_j} denote the regression coefficient vector for taxon
#' \eqn{j}. CAFT tests taxon-level departures from a compositional reference
#' value through
#' \deqn{H_{0j}: \Gamma \beta_j = b,}
#' where \code{Gamma} specifies the tested linear combinations of regression
#' coefficients and \code{b} is the corresponding null reference value.
#'
#' This function is intended for use after \code{caft_estimate()}. The
#' unrestricted estimates \eqn{\hat\beta_j} stored in \code{est$beta.est} are
#' not recomputed. This allows users to test multiple hypotheses, or multiple
#' choices of \code{b}, using the same unrestricted CAFT fit.
#'
#' The hypothesis matrix \code{Gamma} must be supplied to \code{caft_test()}.
#' This keeps unrestricted estimation separate from hypothesis specification:
#' \code{caft_estimate()} computes and stores unrestricted estimates, while
#' \code{caft_test()} computes the hypothesis-specific restricted estimates
#' and score tests for the supplied \code{Gamma} and \code{b}.
#'
#' If \code{b = NULL}, the null reference value is estimated from the
#' unrestricted taxon-level contrasts \eqn{\Gamma \hat\beta_j}. For a
#' one-dimensional contrast, the ordinary median across taxa is used. For
#' multiple contrasts, a robust multivariate location estimate from
#' \code{ICSNP::HR.Mest()} is used. If \code{b} is supplied by the user, the
#' supplied value is used directly as the null reference value.
#'
#' The output records whether the null value was user-specified through
#' \code{b.user.specified}. This is used by \code{caft_bootstrap()} to determine
#' whether the same fixed null value should be used in every bootstrap
#' replicate, or whether the median-based null value should be recomputed within
#' each bootstrap replicate.
#'
#' For each taxon, CAFT obtains the restricted estimate under
#' \eqn{H_{0j}: \Gamma \beta_j = b} and then computes the rank-based score test
#' statistic. Raw taxon-level p-values are returned as \code{p.otu}, and
#' multiplicity-adjusted values are returned as \code{q.otu}. Bootstrap
#' calibration is not performed by this function; it is handled separately by
#' \code{caft_bootstrap()}.
#'
#' @return An object of class \code{"caft_test"}, which also inherits from
#'   \code{"caft_est"}. The object contains the original unrestricted
#'   estimation object together with the restricted testing results.
#' \describe{
#' \item{beta.est}{Unrestricted taxon-level coefficient estimates inherited
#'  from \code{caft_estimate()}.}
#' \item{Gamma}{The hypothesis matrix used in the restricted score test.}
#' \item{Lambda}{The nuisance-parameter contrast matrix used in the restricted
#'  estimation step.}
#' \item{Gamma.ginv}{Generalized inverse of \code{Gamma}.}
#' \item{Lambda.ginv}{Generalized inverse of \code{Lambda}, when applicable.}
#' \item{n.test}{Number of tested contrasts, equal to the number of rows of
#'  \code{Gamma}.}
#' \item{betahat.median}{The median-based null reference value estimated from
#'  the unrestricted taxon-level contrasts \eqn{\Gamma \hat\beta_j} when
#'  \code{b = NULL}. If \code{b} is supplied by the user, this is returned as
#'  \code{NA} because the median-based null value is not computed.}
#' \item{b.null}{The null value actually used in the restricted score test. This
#'  equals \code{betahat.median} when \code{b = NULL}, and equals the
#'  user-specified \code{b} otherwise.}
#' \item{b.user.specified}{Logical indicator of whether \code{b} was supplied
#'  by the user. If \code{FALSE}, \code{b.null} was computed as the
#'  median-based reference value. If \code{TRUE}, \code{b.null} equals the
#'  user-specified null value.}
#' \item{b.user}{The user-specified null value, returned only when \code{b} is
#'  supplied. Otherwise \code{NULL}.}
#' \item{rank.teststat}{Rank-based score test statistics for individual taxa.}
#' \item{rank.teststat.norm}{Numeric matrix with normalized score statistics
#'  for each taxon and each tested contrast. Rows correspond to taxa and columns
#'  correspond to the rows of \code{Gamma}.}
#' \item{p.otu}{Raw p-values for individual taxon-level association tests.}
#' \item{df.test}{Degrees of freedom used in the taxon-level score tests.}
#' \item{beta.est.r}{Restricted coefficient estimates obtained under
#'  \eqn{H_{0j}: \Gamma \beta_j = b}.}
#' \item{p.marginal.otu}{Taxa with raw p-values smaller than
#'  \code{fdr.nominal} based on \code{p.otu}.}
#' \item{q.otu}{Multiplicity-adjusted p-values computed from \code{p.otu} using
#'  the method specified by \code{adjust.method}.}
#' \item{q.detected.otu}{Taxa with adjusted p-values smaller than
#'  \code{fdr.nominal}.}
#' \item{skip.fail.rank.test.pen}{Indicator for taxa whose restricted
#'  estimation or rank-based score test failed.}
#' \item{fit.error.messages.test}{Error messages from failed restricted
#'  taxon-level tests, if any.}
#' \item{fdr.nominal}{The nominal threshold used to report detected taxa.}
#' \item{adjust.method}{The p-value adjustment method used to compute
#'  \code{q.otu}.}
#' \item{test.method}{The variance estimator used in the score test.}
#' \item{mr.resid}{Returned only when \code{return.mr.resid = TRUE}. A list of
#'  subject-level residual contribution matrices evaluated at the restricted
#'  estimates.}
#' }
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom MASS ginv Null
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
#'
#' x <- cbind(x.test, x.adj)
#' Gamma <- cbind(
#'   diag(ncol(x.test)),
#'   matrix(0, nrow = ncol(x.test), ncol = ncol(x.adj))
#' )
#'
#' est <- caft_estimate(
#'   otu.table = count.tab,
#'   x = x,
#'   regularize = TRUE,
#'   n.cores = 1L
#' )
#'
#' res <- caft_test(est, Gamma = Gamma, b = c(0, 0))
#' print(res)
#'
#' ## Reuse the same unrestricted estimates with a different user-specified
#' ## null value.
#' res.b <- caft_test(est, Gamma = Gamma, b = c(0.1, 0.1))
#'
#' ## Example with a full design matrix and an explicit Gamma.
#' x.all <- cbind(CRC = Disease1, adenoma = Disease2, Age = Age, Gender = Gender)
#' est2 <- caft_estimate(otu.table = count.tab, x = x.all)
#' Gamma2 <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
#' res2 <- caft_test(est2, Gamma = Gamma2, b = c(0, 0))
caft_test <- function(est,
                      Gamma,
                      b = NULL,
                      fdr.nominal = 0.20,
                      adjust.method = "BH",
                      regularize = est$regularize,
                      test.method = "rank",
                      n.cores = 1L,
                      return.mr.resid = FALSE) {

  if (!inherits(est, "caft_est")) {
    stop("est must be an object returned by caft_estimate().")
  }
  if (!(test.method %in% c("Cox", "rank", "martingale"))) {
    stop("test.method must be one of 'Cox', 'rank', or 'martingale'.")
  }

  n.cores <- as.integer(n.cores)

  # extract the required information from caft_est class object
  x <- est$x
  otu.table <- est$otu.table
  t.star.all <- est$t.star.all
  delta.all <- est$delta.all
  beta.est <- est$beta.est
  n.data <- est$n.data
  n.taxa <- est$n.taxa
  n.param <- est$n.param
  taxa.name <- est$taxa.name
  filter.thresh <- est$filter.thresh

  # set up Gamma and Lambda matrices
  if (missing(Gamma) || is.null(Gamma)) {
    stop("'Gamma' must be provided to caft_test().", call. = FALSE)
  }
  Gamma <- as.matrix(Gamma)
  if (NCOL(Gamma) != n.param) {
    stop("ncol(Gamma) must equal the number of covariates in est$x.")
  }
  n.test <- NROW(Gamma)
  if (n.test < 1L) {stop("Gamma must have at least one row.")}
  if (n.test < n.param) {
    Lambda <- t(MASS::Null(t(Gamma)))
  } else {
    Lambda <- NULL
  }

  Gamma.ginv <- MASS::ginv(Gamma)
  Lambda.ginv <- if (is.null(Lambda)) NULL else MASS::ginv(Lambda)

  beta.gamma <- as.matrix(beta.est) %*% t(Gamma)

  if (is.null(b)) {
    if (n.test == 1L) {
      betahat.median <- stats::median(beta.gamma[, 1], na.rm = TRUE)
    } else {
      beta.gamma.use <- beta.gamma[stats::complete.cases(beta.gamma), , drop = FALSE]
      if (nrow(beta.gamma.use) < 2L) {
        stop("Not enough complete taxon-level estimates to compute multivariate median.")
      }
      med.fit <- ICSNP::HR.Mest(beta.gamma.use, na.action = stats::na.omit)
      if (!is.null(med.fit$center)) {
        betahat.median <- med.fit$center
      } else if (!is.null(med.fit$mu)) {
        betahat.median <- med.fit$mu
      } else {
        stop("Could not extract the multivariate median from ICSNP::HR.Mest().")
      }
    }
    betahat.median <- as.numeric(betahat.median)
    b.null <- betahat.median
  } else {
    b.null <- as.numeric(b)
    if (length(b.null) != n.test) {
      stop("Length of b must equal the number of rows of Gamma.")
    }
    betahat.median <- rep(NA_real_, n.test)
  }

  fit_one_taxon_test <- function(ii) {
    tstar <- t.star.all[[ii]]
    delta.1 <- delta.all[[ii]]

    if (sum(delta.1) <= n.data * filter.thresh) {
      return(list(
        p = NA_real_,
        test = NA_real_,
        z = rep(NA_real_, n.test),
        df = NA_real_,
        beta_r = rep(NA_real_, n.param),
        mr.resid = NULL,
        skip_fail_test = 0L,
        error_message = NA_character_
      ))
    }

    res.pen <- try(
      estimate.rank.aft(
        y = tstar,
        delta = delta.1,
        x = x,
        Gamma = Gamma,
        Lambda = Lambda,
        Gamma.ginv = Gamma.ginv,
        Lambda.ginv = Lambda.ginv,
        b = b.null,
        beta = NULL,
        test = TRUE,
        regularize = regularize,
        tol = 1e-12
      ),
      silent = TRUE
    )

    if (inherits(res.pen, "try-error")) {
      return(list(
        p = NA_real_,
        test = NA_real_,
        z = rep(NA_real_, n.test),
        df = NA_real_,
        beta_r = rep(NA_real_, n.param),
        mr.resid = NULL,
        skip_fail_test = 1L,
        error_message = as.character(res.pen)
      ))
    }

    temp.rank <- try(
      test.rank.aft(res.pen, score = test.method),
      silent = TRUE
    )

    if (inherits(temp.rank, "try-error")) {
      return(list(
        p = NA_real_,
        test = NA_real_,
        z = rep(NA_real_, n.test),
        df = NA_real_,
        beta_r = rep(NA_real_, n.param),
        mr.resid = NULL,
        skip_fail_test = 1L,
        error_message = as.character(temp.rank)
      ))
    }

    mr.resid <- NULL

    if (return.mr.resid) {
      beta.cur <- if (is.null(res.pen$beta)) res.pen$beta.r else res.pen$beta

      mr.resid <- mySi.no.surv.resid(
        beta = beta.cur,
        y = tstar,
        x = x,
        delta = delta.1
      )$s
    }

    list(
      p = as.numeric(temp.rank$p.value),
      test = as.numeric(temp.rank$test),
      z = as.numeric(temp.rank$z.score),
      df = as.numeric(temp.rank$df),
      beta_r = as.numeric(res.pen$beta.r),
      mr.resid = mr.resid,
      skip_fail_test = 0L,
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

    res_phase2 <- foreach::foreach(
      ii = seq_len(n.taxa),
      .errorhandling = "pass"
    ) %dopar% {
      fit_one_taxon_test(ii)
    }
  } else {
    res_phase2 <- lapply(seq_len(n.taxa), fit_one_taxon_test)
  }

  res_phase2 <- lapply(res_phase2, function(z) {
    if (inherits(z, "error")) {
      list(
        p = NA_real_,
        test = NA_real_,
        z = rep(NA_real_, n.test),
        df = NA_real_,
        beta_r = rep(NA_real_, n.param),
        mr.resid = NULL,
        skip_fail_test = 1L,
        error_message = "Parallel restricted taxon-level test failed."
      )
    } else {
      z
    }
  })

  p.otu <- vapply(res_phase2, function(z) z$p, numeric(1))
  rank.teststat <- vapply(res_phase2, function(z) z$test, numeric(1))
  df.test <- vapply(res_phase2, function(z) z$df, numeric(1))
  skip.fail.rank.test.pen <- vapply(
    res_phase2,
    function(z) z$skip_fail_test,
    integer(1)
  )

  rank.teststat.norm <- do.call(rbind, lapply(res_phase2, `[[`, "z"))
  rank.teststat.norm <- as.matrix(rank.teststat.norm)

  if (ncol(rank.teststat.norm) != n.test) {
    rank.teststat.norm <- matrix(
      rank.teststat.norm,
      nrow = n.taxa,
      ncol = n.test,
      byrow = TRUE
    )
  }

  beta.est.r <- do.call(rbind, lapply(res_phase2, `[[`, "beta_r"))
  beta.est.r <- as.data.frame(beta.est.r)
  colnames(beta.est.r) <- paste0("b", seq_len(n.param), ".est.r")
  rownames(beta.est.r) <- taxa.name

  names(p.otu) <- taxa.name
  names(rank.teststat) <- taxa.name
  names(df.test) <- taxa.name
  names(skip.fail.rank.test.pen) <- taxa.name
  rownames(rank.teststat.norm) <- taxa.name
  colnames(rank.teststat.norm) <- paste0("test", seq_len(n.test))

  q.otu <- stats::p.adjust(p.otu, method = adjust.method)
  p.marginal.otu <- taxa.name[!is.na(p.otu) & p.otu < fdr.nominal]
  q.detected.otu <- taxa.name[!is.na(q.otu) & q.otu < fdr.nominal]

  fit.error.messages.test <- vapply(
    res_phase2,
    function(z) z$error_message,
    character(1)
  )
  names(fit.error.messages.test) <- taxa.name

  if (return.mr.resid) {
    mr.resid <- lapply(res_phase2, `[[`, "mr.resid")
    names(mr.resid) <- taxa.name
  } else {
    mr.resid <- NULL
  }

  out <- est

  out$Gamma <- Gamma
  out$Lambda <- Lambda
  out$Gamma.ginv <- Gamma.ginv
  out$Lambda.ginv <- Lambda.ginv
  out$n.test <- n.test
  out$betahat.median <- betahat.median
  out$b.null <- b.null
  out$b.user.specified <- !is.null(b)
  out$b.user <- if (is.null(b)) NULL else b.null
  out$rank.teststat <- rank.teststat
  out$rank.teststat.norm <- rank.teststat.norm
  out$p.otu <- p.otu
  out$df.test <- df.test
  out$beta.est.r <- beta.est.r
  out$p.marginal.otu <- p.marginal.otu
  out$q.otu <- q.otu
  out$q.detected.otu <- q.detected.otu
  out$skip.fail.rank.test.pen <- skip.fail.rank.test.pen
  out$fit.error.messages.test <- fit.error.messages.test
  out$fdr.nominal <- fdr.nominal
  out$adjust.method <- adjust.method
  out$test.method <- test.method

  if (return.mr.resid) {
    out$mr.resid <- mr.resid
  }

  class(out) <- c("caft_test", "caft_est")
  out
}

#' Print Method for CAFT Restricted Score Test
#'
#' @param x An object of class \code{"caft_test"}.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print caft_test
#' @export
print.caft_test <- function(x, ...) {
  cat("CAFT restricted score test\n")
  cat("  Samples:", x$n.data, "\n")
  cat("  Taxa:", x$n.taxa, "\n")
  cat("  Covariates:", x$n.param, "\n")
  cat("  Tested contrasts:", x$n.test, "\n")
  cat("  FDR level:", x$fdr.nominal, "\n")
  cat("  Adjustment method:", x$adjust.method, "\n")
  cat("  Taxa with p <", x$fdr.nominal, ":", length(x$p.marginal.otu), "\n")
  cat("  Taxa with q <", x$fdr.nominal, ":", length(x$q.detected.otu), "\n")
  cat(
    "  Taxa failed in restricted test:",
    sum(x$skip.fail.rank.test.pen == 1L, na.rm = TRUE),
    "\n"
  )

  invisible(x)
}
