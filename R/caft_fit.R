#' Internal single-fit engine for CAFT
#'
#' Low-level worker used by `caft()` after construction of the full design
#' matrix `x`, hypothesis matrix `Gamma`, and associated generalized inverses.
#'
#' @keywords internal
#' @noRd
caft_fit <- function(otu.table, x, Gamma, Gamma.ginv, Lambda, Lambda.ginv,
                     filter.thresh = 0.05, fdr.nominal = 0.20, adjust.method = "BH",
                     regularize=TRUE, test.method="rank", n.cores=1L,
                     return.mr.resid = FALSE) {

  # center the covariates
  x <- as.matrix(x - t(replicate(NROW(x), colMeans(x))))

  n.data  <- NROW(otu.table)
  n.taxa  <- NCOL(otu.table)
  n.param <- NCOL(x)
  n.test  <- nrow(Gamma)

  taxa.name <- colnames(otu.table)
  if (is.null(taxa.name)) {
    taxa.name <- paste("tstar", 1:n.taxa)
  }
  ra.all <- otu.table / rowSums(otu.table)
  lib.size <- rowSums(otu.table)
  # Censored relative abundance: every row should same
  lim.ra <- matrix(1 / lib.size, nrow = n.data, ncol = n.taxa)

  #--------------------------
  # create survival data
  #--------------------------
  log_neg <- -log10(lim.ra)
  log_pos <- -log10(ra.all)
  t.star.all <- matrix(log_neg,
                       nrow = nrow(otu.table),
                       ncol = ncol(otu.table)
  )
  idx <- (otu.table > 0)
  t.star.all[idx] <- log_pos[idx]
  t.star.all <- as.data.frame(t.star.all)
  colnames(t.star.all) <- paste("tstar", 1:n.taxa)

  delta.all <- ((otu.table > 0)^2)
  colnames(delta.all) <- paste("delta", 1:n.taxa)

  data <- data.frame(t.star.all, delta.all)

  if(n.cores > 1L){

    doParallel::registerDoParallel(cores=n.cores)

    col_tstar_idx  <- grep("^tstar", colnames(data))
    col_delta_idx  <- grep("^delta", colnames(data))
    stopifnot(length(col_tstar_idx) >= n.taxa, length(col_delta_idx) >= n.taxa)

    ncoef <- NCOL(x)

    res_phase1 <- foreach(
      ii = 1:n.taxa,
      .errorhandling = "pass"
    ) %dopar% {
      tstar   <- data[, col_tstar_idx[ii]]
      delta.1 <- data[, col_delta_idx[ii]]

      if (sum(delta.1) <= n.data * filter.thresh) {
        return(list(
          beta = rep(NA_real_, ncoef),
          skip_rare = 1L,
          skip_fail_fit = 0L
        ))
      }

      fit0.pen <- try(estimate.rank.aft(
        y = tstar, delta = delta.1, x = x,
        Gamma = NULL, Lambda = diag(n.param),
        Gamma.ginv = NULL, Lambda.ginv = diag(n.param),
        b = NULL, beta = NULL,
        test = TRUE, regularize = regularize, tol = 1e-12
      ), silent = TRUE)

      if (inherits(fit0.pen, "try-error")) {
        list(beta = rep(NA_real_, ncoef), skip_rare = 0L, skip_fail_fit = 1L)
      } else {
        list(beta = as.numeric(fit0.pen$beta), skip_rare = 0L, skip_fail_fit = 0L)
      }
    }

    beta.est <- do.call(rbind, lapply(res_phase1, `[[`, "beta"))
    beta.est <- as.data.frame(beta.est)
    colnames(beta.est) <- paste0("b", 1:NCOL(x), ".est")
    skip.rare <- vapply(res_phase1, function(z) z$skip_rare, integer(1))
    skip.fail.rank.fit.pen <- vapply(res_phase1, function(z) z$skip_fail_fit, integer(1))

    if (n.test==1) {
      betahat.median = median( as.matrix(beta.est) %*% t(Gamma), na.rm = T)
    } else {
      betahat.median = ICSNP::HR.Mest(as.matrix(beta.est) %*% t(Gamma), na.action=na.omit)$center
    }
    if (n.test==n.param) Lambda=NULL

    res_phase2 <- foreach(
      ii = 1:n.taxa,
      .errorhandling = "pass"
    ) %dopar% {
      tstar   <- data[, col_tstar_idx[ii]]
      delta.1 <- data[, col_delta_idx[ii]]

      if (sum(delta.1) <= n.data * filter.thresh) {
        return(list(
          p = NA_real_, test = NA_real_,
          z = rep(NA_real_, n.test), df = NA_real_,
          beta_r = rep(NA_real_, ncoef),
          skip_fail_test = 0L
        ))
      }

      res.pen <- try(estimate.rank.aft(
        y = tstar, delta = delta.1, x = x,
        Gamma = Gamma, Lambda = Lambda,
        Gamma.ginv = Gamma.ginv, Lambda.ginv = Lambda.ginv,
        b = betahat.median, beta = NULL,
        test = TRUE, regularize = regularize, tol = 1e-12
      ), silent = TRUE)

      if (inherits(res.pen, "try-error")) {
        return(list(
          p = NA_real_, test = NA_real_,
          z = rep(NA_real_, n.test), df = NA_real_,
          beta_r = rep(NA_real_, ncoef),
          skip_fail_test = 1L
        ))
      }

      # temp.rank <- try(test.rank.aft(res.pen, score = "rank"), silent = TRUE)
      temp.rank <- try(test.rank.aft(res.pen, score = test.method), silent = TRUE)
      if (inherits(temp.rank, "try-error")) {
        list(
          p = NA_real_, test = NA_real_,
          z = rep(NA_real_, n.test), df = NA_real_,
          beta_r = rep(NA_real_, ncoef),
          skip_fail_test = 1L
        )
      } else {
        mr.resid <- NULL

        if (return.mr.resid) {
          beta.cur <- if (is.null(res.pen$beta)) res.pen$beta.r else res.pen$beta

          mr.resid <- mySi.no.surv.resid(
            beta  = beta.cur,
            y     = tstar,
            x     = x,
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
          skip_fail_test = 0L
        )
      }
    }

    p.rank.pen            <- vapply(res_phase2, function(z) z$p,    numeric(1))
    rank.teststat             <- vapply(res_phase2, function(z) z$test, numeric(1))
    df.rank.pen           <- vapply(res_phase2, function(z) z$df,   numeric(1))
    skip.fail.rank.test.pen <- vapply(res_phase2, function(z) z$skip_fail_test, integer(1))
    rank.teststat.norm        <- do.call(rbind, lapply(res_phase2, `[[`, "z"))
    if (!is.matrix(rank.teststat.norm)) rank.teststat.norm <- matrix(rank.teststat.norm, nrow = n.taxa, ncol = n.test, byrow = TRUE)
    beta.est.r <- do.call(rbind, lapply(res_phase2, `[[`, "beta_r"))
    if (!is.matrix(beta.est.r)) beta.est.r <- matrix(beta.est.r, nrow = n.taxa, ncol = ncoef, byrow = TRUE)
    beta.est.r <- as.data.frame(beta.est.r)
    if (return.mr.resid) {
      mr.resid <- lapply(res_phase2, `[[`, "mr.resid")
    } else {
      mr.resid <- NULL
    }
  }else{

    #--------------------------
    # unconstrained parameter estimates
    #--------------------------

    beta.est <- beta.est.r <- as.data.frame(matrix(NA, n.taxa, NCOL(x)))
    colnames(beta.est) <- paste0("b", 1:NCOL(x), ".est")

    skip.rare <- skip.fail.rank.fit <- rep(0, n.taxa)
    skip.fail.rank.fit.pen <- rep(0, n.taxa)

    if (return.mr.resid) {
      mr.resid <- vector("list", n.taxa)
    } else {
      mr.resid <- NULL
    }

    for (ii in 1:n.taxa) {
      tstar <- data[, grep("tstar", colnames(data))[ii]]
      delta.1 <- data[, grep("delta", colnames(data))[ii]]
      if (sum(delta.1) <= n.data * filter.thresh) {
        out <- rep(NA, NCOL(x))
        skip.rare[ii] <- 1
      } else {
        fit0.pen <- try(estimate.rank.aft(
          y = tstar, delta = delta.1, x = x,
          Gamma = NULL, Lambda = diag(n.param),
          Gamma.ginv = NULL, Lambda.ginv = diag(n.param),
          b = NULL, beta = NULL,
          test = TRUE, regularize = regularize,
          tol = 10^-12
        ))
        if (inherits(fit0.pen, "try-error")) {
          beta.est[ii, ] <- rep(NA, NCOL(x))
          skip.fail.rank.fit.pen[ii] <- 1
        } else {
          beta.est[ii, ] <- fit0.pen$beta
        }
      }
    }

    #--------------------------
    # test equal to median
    #--------------------------

    if (n.test==1) {
      betahat.median = median( as.matrix(beta.est) %*% t(Gamma), na.rm = T)
    }else {
      betahat.median = ICSNP::HR.Mest(as.matrix(beta.est) %*% t(Gamma), na.action=na.omit)$center
    }
    if (n.test==n.param) Lambda=NULL

    rank.teststat <- df.rank.pen <- rep(NA, n.taxa)
    rank.teststat.norm <- matrix(NA, ncol = n.test, nrow = n.taxa)
    p.rank.pen <- rep(NA, n.taxa)
    skip.fail.rank.test.pen <- rep(NA, n.taxa)
    for (ii in 1:n.taxa) {
      tstar <- data[, grep("tstar", colnames(data))[ii]]
      delta.1 <- data[, grep("delta", colnames(data))[ii]]
      if (sum(delta.1) > n.data * filter.thresh) {
        res.pen <- try(estimate.rank.aft(
          y = tstar, delta = delta.1, x = x,
          Gamma = Gamma, Lambda = Lambda,
          Gamma.ginv = Gamma.ginv, Lambda.ginv = Lambda.ginv,
          b = betahat.median,
          beta = NULL, test = TRUE,
          regularize = regularize, tol = 10^-12
        ))
        if (inherits(res.pen, "try-error")) {
          p.rank.pen[ii] <- NA
          rank.teststat[ii] <- NA
          rank.teststat.norm[ii, ] <- NA
          skip.fail.rank.test.pen[ii] <- 1
        } else {
          # temp.rank <- try(test.rank.aft(res.pen, score = "rank"))
          temp.rank <- try(test.rank.aft(res.pen, score = test.method))
          if (inherits(temp.rank, "try-error")) {
            p.rank.pen[ii] <- NA
            rank.teststat.norm[ii, ] <- NA
            rank.teststat[ii] <- NA
            skip.fail.rank.test.pen[ii] <- 1
          } else {
            p.rank.pen[ii] <- temp.rank$p.value
            rank.teststat[ii] <- temp.rank$test
            rank.teststat.norm[ii, ] <- temp.rank$z.score
            df.rank.pen[ii] <- temp.rank$df
            beta.est.r[ii, ] <- res.pen$beta.r
            if (return.mr.resid) {
              beta.cur <- if (is.null(res.pen$beta)) res.pen$beta.r else res.pen$beta

              mr.resid[[ii]] <- mySi.no.surv.resid(
                beta  = beta.cur,
                y     = tstar,
                x     = x,
                delta = delta.1
              )$s
            }
          }
        }
      }
    }
  }

  p.otu <- p.rank.pen
  q.otu <- p.adjust(p.otu, method = adjust.method)

  out <- list(
    beta.est = beta.est,
    betahat.median = betahat.median,
    rank.teststat = rank.teststat,
    rank.teststat.norm = rank.teststat.norm,
    skip.otu = skip.rare,
    p.otu = p.otu,
    df.test = df.rank.pen,
    beta.est.r = beta.est.r,
    p.detected.otu = colnames(otu.table)[which(p.otu < fdr.nominal)],
    q.otu = q.otu,
    q.detected.otu = colnames(otu.table)[which(q.otu < fdr.nominal)]
  )

  if (return.mr.resid) {
    out$mr.resid <- mr.resid
  }

  return(out)

}
