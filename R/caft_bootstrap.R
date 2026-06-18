#' Bootstrap Calibration for CAFT Tests
#'
#' @description
#' Provides bootstrap calibration of taxon-level p-values for an existing CAFT
#' restricted score test object returned by \code{caft_test()}. Bootstrap
#' calibration is mainly intended as a sensitivity analysis or preferred
#' calibration method when the number of taxa is small in the data set.
#'
#' @param res An object of class \code{"caft_test"} returned by
#'  \code{caft_test()}.
#' @param boot.B Integer. Number of bootstrap replicates used for empirical
#'  calibration of taxon-level p-values. The default is \code{1000L}. At least
#'  two bootstrap replicates are required for bootstrap calibration.
#' @param boot.parallel Logical. If \code{TRUE}, the bootstrap loop is
#'  parallelized over bootstrap replicates. The default is \code{FALSE}.
#' @param boot.n.cores Integer. Number of CPU cores used for outer bootstrap
#'  parallelization when \code{boot.parallel = TRUE}. The default is
#'  \code{1L}.
#' @param boot.seed Optional integer seed for reproducible bootstrap resampling.
#'  The default is \code{NULL}.
#' @param boot.return.dist Logical. If \code{TRUE}, return the bootstrap-centered
#'  beta distribution and the names of taxa used in the bootstrap calibration.
#'  The default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, progress messages are shown during
#'  sequential bootstrap computation. The default is \code{FALSE}.
#'
#' @details
#' This function is intended to be used after \code{caft_estimate()} and
#' \code{caft_test()}. The observed restricted score test is not recomputed.
#' Instead, \code{caft_bootstrap()} uses the hypothesis matrix \code{Gamma}, the
#' null value \code{b.null}, the nominal FDR threshold, and the adjustment
#' method stored in the \code{caft_test} object.
#'
#' CAFT bootstrap calibration uses a fixed-design bootstrap procedure. Samples
#' are resampled with replacement from the observed OTU table while the
#' covariates are kept fixed. The bootstrap is carried out on the taxa that
#' pass the original filtering step. In each bootstrap replicate, the
#' unrestricted CAFT model is refit using \code{caft_estimate()} with
#' \code{filter.thresh = 0} on the same set of taxa. The inner CAFT fitting
#' routine is run with \code{n.cores = 1L} inside each bootstrap replicate to
#' avoid nested parallelism.
#'
#' Bootstrap p-values are computed from the centered bootstrap distribution of
#' the tested coefficient contrasts. If the original test used the default
#' median-based null value, the bootstrap null value is recomputed within each
#' bootstrap replicate from the bootstrap taxon-level contrasts. If the original
#' test used a user-specified \code{b}, the same fixed null value is used in
#' every bootstrap replicate.
#'
#' For a single testing variable, both a two-sided empirical beta-tilde
#' bootstrap p-value and a studentized chi-square/Wald bootstrap p-value are
#' returned. For multiple testing variables, the bootstrap p-value is based on a
#' multivariate statistic; in this case \code{p.boot} and \code{p.boot.chi} are
#' identical.
#'
#' If \code{boot.parallel = TRUE}, the outer bootstrap loop is parallelized over
#' bootstrap replicates using \code{foreach}/\code{doParallel}. Reproducible
#' parallel bootstrap sampling is handled through \code{doRNG}.
#'
#' @return The input \code{caft_test} object with additional bootstrap
#'  calibration results:
#' \describe{
#' \item{p.boot}{Bootstrap-calibrated empirical beta-tilde p-values. If
#'  \code{boot.B < 2}, this is a named vector of \code{NA} values.}
#' \item{q.boot}{Multiplicity-adjusted values from \code{p.boot}, computed using
#'  the adjustment method stored in \code{res$adjust.method}.}
#' \item{p.boot.marginal.otu}{Taxa detected using \code{p.boot < fdr.nominal}.
#'  Empty if bootstrap calibration is not run.}
#' \item{q.boot.detected.otu}{Taxa detected using
#'  \code{q.boot < fdr.nominal}. Empty if bootstrap calibration is not run.}
#' \item{p.boot.chi}{Studentized chi-square/Wald bootstrap p-values. For
#'  multiple testing variables, these are the same as \code{p.boot}. If
#'  \code{boot.B < 2}, this is a named vector of \code{NA} values.}
#' \item{q.boot.chi}{Multiplicity-adjusted values from \code{p.boot.chi}.}
#' \item{p.boot.chi.marginal.otu}{Taxa detected using
#'  \code{p.boot.chi < fdr.nominal}. Empty if bootstrap calibration is not run.}
#' \item{q.boot.chi.detected.otu}{Taxa detected using
#'  \code{q.boot.chi < fdr.nominal}. Empty if bootstrap calibration is not run.}
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
#' }
#'
#' @importFrom foreach foreach registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doRNG %dorng%
#' @importFrom MASS ginv
#' @export
#'
#' @examples
#' data(Colon)
#'
#' count.tab <- t(as(phyloseq::otu_table(Colon), "matrix"))
#' sample.tab <- as(phyloseq::sample_data(Colon), "data.frame")
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
#' keep.top <- names(sort(prev, decreasing = TRUE))[seq_len(min(5L, length(prev)))]
#' count.tab <- count.tab[, keep.top, drop = FALSE]
#'
#' ## Remove samples with zero library size after subsetting taxa.
#' keep.lib <- rowSums(count.tab) > 0
#' count.tab <- count.tab[keep.lib, , drop = FALSE]
#' sample.tab <- sample.tab[keep.lib, , drop = FALSE]
#'
#' Disease1 <- rep(0, nrow(sample.tab))
#' Disease1[sample.tab$disease == "CRC"] <- 1
#'
#' Age <- as.numeric(sample.tab$age)
#' Gender <- as.numeric(factor(sample.tab$gender)) - 1
#'
#' x.test <- cbind(CRC = Disease1)
#' x.adj <- cbind(Age = Age, Gender = Gender)
#'
#' x <- cbind(x.test, x.adj)
#'
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
#' res <- caft_test(est, Gamma = Gamma, b = 0)
#'
#' \donttest{
#' ## Use a very small number of replicates only for illustration.
#' res.boot <- caft_bootstrap(
#'   res,
#'   boot.B = 2L,
#'   boot.seed = 1L
#' )
#'
#' res.boot$boot.n.success
#' head(res.boot$p.boot)
#' }
caft_bootstrap <- function(res,
                           boot.B = 1000L,
                           boot.parallel = FALSE,
                           boot.n.cores = 1L,
                           boot.seed = NULL,
                           boot.return.dist = FALSE,
                           verbose = FALSE) {

  if (!inherits(res, "caft_test")) {
    stop("res must be an object returned by caft_test().")
  }

  boot.B <- as.integer(boot.B)
  boot.n.cores <- as.integer(boot.n.cores)

  # extract the required information from caft_test class object
  otu.table <- res$otu.table
  x <- res$x
  Gamma <- res$Gamma
  n.test <- res$n.test
  fdr.nominal <- res$fdr.nominal
  adjust.method <- res$adjust.method
  regularize <- res$regularize

  if (is.null(otu.table) || is.null(x) || is.null(Gamma)) {
    stop("res must contain otu.table, x, and Gamma.")
  }

  if (!is.matrix(otu.table)) {otu.table <- as.matrix(otu.table)}
  if (!is.matrix(x)) {x <- as.matrix(x)}

  Gamma <- as.matrix(Gamma)
  n.otu <- length(res$p.otu)
  otu.names <- names(res$p.otu)
  if (is.null(otu.names) || length(otu.names) != n.otu) {
    otu.names <- colnames(otu.table)
  }

  empty_boot_results <- function(out, boot.B) {
    out$p.boot <- setNames(rep(NA_real_, n.otu), otu.names)
    out$q.boot <- setNames(rep(NA_real_, n.otu), otu.names)
    out$p.boot.chi <- setNames(rep(NA_real_, n.otu), otu.names)
    out$q.boot.chi <- setNames(rep(NA_real_, n.otu), otu.names)

    out$p.boot.marginal.otu <- character(0)
    out$q.boot.detected.otu <- character(0)
    out$p.boot.chi.marginal.otu <- character(0)
    out$q.boot.chi.detected.otu <- character(0)

    out$boot.median <- matrix(
      NA_real_,
      nrow = 0L,
      ncol = n.test,
      dimnames = list(NULL, paste0("test", seq_len(n.test)))
    )

    out$boot.n.success <- 0L
    out$boot.n.fail <- as.integer(max(boot.B, 0L))

    class(out) <- unique(c("caft_boot", class(out)))
    out
  }

  if (boot.B < 2L) {
    if (boot.B == 1L) {
      warning("At least two bootstrap replicates are required; bootstrap calibration skipped.")
    }
    return(empty_boot_results(res, boot.B))
  }

  # Null bootstrap (using filtered taxa)
  use <- !is.na(res$p.otu)

  if (sum(use) == 0L) {
    warning("No taxa passed filtering; bootstrap calibration skipped.")
    return(empty_boot_results(res, boot.B))
  }

  n.taxa.boot <- sum(use)
  taxa.boot <- otu.names[use]

  beta.obs.all <- as.matrix(res$beta.est)
  storage.mode(beta.obs.all) <- "numeric"

  beta.obs <- beta.obs.all[use, , drop = FALSE] %*% t(Gamma)
  beta.null.obs <- as.numeric(res$b.null[seq_len(n.test)])
  beta.tilde.obs <- sweep(beta.obs, 2, beta.null.obs, "-")

  otu.boot.base <- otu.table[, use, drop = FALSE]

  b.user.specified <- isTRUE(res$b.user.specified)
  b.fixed <- as.numeric(res$b.null[seq_len(n.test)])

  compute_boot_null <- function(beta.gamma.b) {
    beta.gamma.b <- as.matrix(beta.gamma.b)
    storage.mode(beta.gamma.b) <- "numeric"

    if (b.user.specified) {
      return(b.fixed)
    }

    if (n.test == 1L) {
      return(stats::median(beta.gamma.b[, 1], na.rm = TRUE))
    }

    beta.gamma.use <- beta.gamma.b[
      stats::complete.cases(beta.gamma.b),
      ,
      drop = FALSE
    ]

    if (nrow(beta.gamma.use) < 2L) {
      return(rep(NA_real_, n.test))
    }

    med.fit <- tryCatch(
      ICSNP::HR.Mest(beta.gamma.use, na.action = stats::na.omit),
      error = function(e) e
    )

    if (inherits(med.fit, "error")) {
      return(rep(NA_real_, n.test))
    }

    if (!is.null(med.fit$center)) {
      return(as.numeric(med.fit$center[seq_len(n.test)]))
    }

    if (!is.null(med.fit$mu)) {
      return(as.numeric(med.fit$mu[seq_len(n.test)]))
    }

    rep(NA_real_, n.test)
  }

  fit_one_boot <- function(bb) {
    idx <- sample.int(nrow(otu.table), size = nrow(otu.table), replace = TRUE)
    otu.boot <- otu.boot.base[idx, , drop = FALSE]

    est.b <- try(
      caft_estimate(
        otu.table = otu.boot,
        x = x,
        filter.thresh = 0,
        regularize = regularize,
        n.cores = 1L
      ),
      silent = TRUE
    )

    if (inherits(est.b, "try-error")) {
      return(list(
        ok = FALSE,
        beta = matrix(NA_real_, nrow = n.taxa.boot, ncol = n.test),
        median = rep(NA_real_, n.test),
        msg = "caft_estimate failed in this bootstrap replicate."
      ))
    }

    beta.b.all <- as.matrix(est.b$beta.est)

    if (nrow(beta.b.all) != n.taxa.boot) {
      return(list(
        ok = FALSE,
        beta = matrix(NA_real_, nrow = n.taxa.boot, ncol = n.test),
        median = rep(NA_real_, n.test),
        msg = "Bootstrap fit returned an unexpected number of taxa."
      ))
    }

    beta.b <- beta.b.all %*% t(Gamma)
    boot.null <- compute_boot_null(beta.b)

    if (any(!is.finite(boot.null))) {
      return(list(
        ok = FALSE,
        beta = matrix(NA_real_, nrow = n.taxa.boot, ncol = n.test),
        median = rep(NA_real_, n.test),
        msg = "Bootstrap null value could not be computed."
      ))
    }

    list(
      ok = TRUE,
      beta = beta.b,
      median = boot.null,
      msg = NA_character_
    )
  }

  beta.boot <- array(
    NA_real_,
    dim = c(boot.B, n.taxa.boot, n.test),
    dimnames = list(NULL, taxa.boot, paste0("test", seq_len(n.test)))
  )

  boot.median <- matrix(
    NA_real_,
    nrow = boot.B,
    ncol = n.test,
    dimnames = list(NULL, paste0("test", seq_len(n.test)))
  )

  outer.parallel <- isTRUE(boot.parallel) && boot.n.cores > 1L

  if (outer.parallel) {
    rng.seed <- if (is.null(boot.seed)) {
      sample.int(.Machine$integer.max, 1L)
    } else {
      boot.seed
    }
    cl <- parallel::makeCluster(boot.n.cores, type = "PSOCK")
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(cl)

    boot.res <- foreach::foreach(
      bb = seq_len(boot.B),
      .inorder = TRUE,
      .errorhandling = "pass",
      .packages = c("CAFT", "ICSNP"),
      .export = c("caft_estimate"),
      .options.RNG = rng.seed
    ) %dorng% {
      fit_one_boot(bb)
    }
  } else {
    if (!is.null(boot.seed)) {set.seed(boot.seed)}
    boot.res <- vector("list", boot.B)
    for (bb in seq_len(boot.B)) {
      if (isTRUE(verbose) && boot.B >= 20L && bb %% 20L == 0L) {
        message("Bootstrap replicate ", bb, " of ", boot.B)
      }
      boot.res[[bb]] <- fit_one_boot(bb)
    }
  }

  for (bb in seq_len(boot.B)) {
    z <- boot.res[[bb]]
    if (!is.list(z) || is.null(z$ok)) {
      if (isTRUE(verbose)) {
        message("Bootstrap replicate ", bb, " failed and was skipped: unexpected bootstrap return object.")
      }
      next
    }
    if (isTRUE(z$ok)) {
      beta.boot[bb, , ] <- z$beta
      boot.median[bb, ] <- z$median
    } else {
      if (isTRUE(verbose)) {
        message("Bootstrap replicate ", bb, " failed and was skipped: ", z$msg)
      }
    }
  }

  # Center each bootstrap beta by its null median
  beta.boot.tilde <- array(NA_real_, dim = c(boot.B, n.taxa.boot, n.test), dimnames = list(NULL, taxa.boot, paste0("test", seq_len(n.test))))
  for (kk in seq_len(n.test)) {
    beta.boot.tilde[, , kk] <- sweep(beta.boot[, , kk, drop = TRUE], 1, boot.median[, kk], "-")
  }

  ########################
  ## Bootstrap p-values ##
  ########################

  # compute bootstrap p-value based on CI
  p.boot.short <- rep(NA_real_, n.taxa.boot)
  if (n.test == 1L) {
    beta.boot.tilde.mat <- beta.boot.tilde[, , 1, drop = TRUE]
    beta.tilde.obs.vec <- as.numeric(beta.tilde.obs[, 1])

    for (jj in seq_len(n.taxa.boot)) {
      boot.vec <- beta.boot.tilde.mat[, jj]
      obs.val <- beta.tilde.obs.vec[jj]
      good <- is.finite(boot.vec) & is.finite(obs.val)
      if (sum(good) > 0L) {
        p.right <- sum(boot.vec[good] >= obs.val)
        p.left <- sum(boot.vec[good] <= obs.val)
        p.boot.short[jj] <- min(1, (1 + 2 * min(p.right, p.left)) / (1 + sum(good)))
      }
    }
  } else {
    for (jj in seq_len(n.taxa.boot)) {
      boot.mat <- matrix(beta.boot.tilde[, jj, ], nrow = boot.B, ncol = n.test)
      obs.vec <- as.numeric(beta.tilde.obs[jj, ])
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
        good.T <- is.finite(T.boot.j) & is.finite(T.obs.j)

        if (sum(good.T) > 0L) {
          p.boot.short[jj] <- (1 + sum(T.boot.j[good.T] >= T.obs.j)) / (1 + sum(good.T))
        }
      }
    }
  }

  # compute bootstrap Wald test
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
    for (jj in seq_len(n.taxa.boot)) {
      good <- is.finite(null.Z[, jj]) & is.finite(Z[jj])
      if (sum(good) > 0L) {
        p.boot.chi.short[jj] <- (1 + sum(null.Z[good, jj]^2 >= Z[jj]^2)) / (1 + sum(good))
      }
    }
  } else {
    # For multiple testing variables, p.boot.short is already based on
    # the multivariate studentized quadratic/Wald statistic.
    p.boot.chi.short <- p.boot.short
  }

  p.boot <- setNames(rep(NA_real_, n.otu), otu.names)
  p.boot[use] <- p.boot.short
  q.boot <- setNames(rep(NA_real_, n.otu), otu.names)
  q.boot[use] <- stats::p.adjust(p.boot[use], method = adjust.method)
  p.boot.chi <- setNames(rep(NA_real_, n.otu), otu.names)
  p.boot.chi[use] <- p.boot.chi.short
  q.boot.chi <- setNames(rep(NA_real_, n.otu), otu.names)
  q.boot.chi[use] <- stats::p.adjust(p.boot.chi[use], method = adjust.method)

  out <- res

  out$p.boot <- p.boot
  out$q.boot <- q.boot
  out$p.boot.chi <- p.boot.chi
  out$q.boot.chi <- q.boot.chi

  out$p.boot.marginal.otu <- otu.names[!is.na(p.boot) & p.boot < fdr.nominal]
  out$q.boot.detected.otu <- otu.names[!is.na(q.boot) & q.boot < fdr.nominal]
  out$p.boot.chi.marginal.otu <- otu.names[!is.na(p.boot.chi) & p.boot.chi < fdr.nominal]
  out$q.boot.chi.detected.otu <- otu.names[!is.na(q.boot.chi) & q.boot.chi < fdr.nominal]

  out$boot.median <- boot.median
  out$boot.n.success <- sum(stats::complete.cases(boot.median))
  out$boot.n.fail <- boot.B - out$boot.n.success

  if (isTRUE(boot.return.dist)) {
    out$boot.beta.tilde <- beta.boot.tilde
    out$boot.taxa <- taxa.boot
  }

  class(out) <- unique(c("caft_boot", class(out)))
  out
}
