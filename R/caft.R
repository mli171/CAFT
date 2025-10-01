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
#' @param otu.table the community OTU table (or taxa count table). Each row
#'  corresponds to a sample and each column corresponds to one OTU (taxa).
#' @param x.test the covariates are of interest. It can be a vector, matrix, or data
#'    frame. All K-level categorical variables should be converted to K-1
#'    indicator variables before input.  Use if default Gamma is desired.
#' @param x.cov other covariates that need to be adjusted. Requirements are the
#'    same as \code{X.test} above.  Use if default Gamma is desired.
#' @param Gamma matrix to specify what hypothesis to test.  Not used if data are entered as X.test and X.cov.
#' @param x data for all covariates.  Only used if matrix Gamma is provided as input. It can be a vector, matrix, or data
#'    frame. All K-level categorical variables should be converted to K-1 indicator variables before input.
#' @param filter.thresh a real value between 0 and 1 for OTU table sample presence
#'    filtering. Any OTUs present in fewer than \code{filter.thresh} proportion
#'    of samples are filtered out. We set the default to be 0.05.
#' @param fdr.nominal the nominal FDR. The default is 0.2.
#' @param adjust.method a character string. Use multiple comparison/testing
#'  adjustment methods to control the family-wise error rate/false discover
#'  rate. Default to "BH". See \code{\link{p.adjust}} for the details.
#' @param n.cores the number of cores to be used for parallel computing. The default is 1.
#' @return Return a list consisting of
#' \describe{
#' \item{est.rank.gs.pen}{A matrix that include the estimated coefficients by fitting
#' each OTU count to the log-linear model. The number of columns should match
#' the number of covairates from \code{Y} and \code{C}.}
#' \item{b.median.pen}{The selected median of all betas, which was used in the
#' restricted score test of the significance of differential abundant OTU. It can
#' also be used to calculated the effect size of each OTU, defined as the difference
#' from the median of all betas.}
#' \item{test.rank}{test statistics from the proposed score test}
#' \item{skip.otu}{the names of skipped OTU during taxa presence filtering
#'          above}
#' \item{p.otu}{p-values for individual OTU association tests}
#' \item{q.otu}{q-values (adjusted p-values by the Bonferroni correction or the
#'    adjustment method specified in \code{adjust.method}.)
#'          for individual OTU association tests}
#' \item{p.detected.otu}{detected significantly differential abundant taxa
#'          (denoted by the column names of the OTU table) at the nominal FDR
#'          based on \code{p.otu}}
#' \item{q.detected.otu}{detected significantly differential abundant taxa
#'          (denoted by the column names of the OTU table) at the nominal FDR
#'          based on \code{q.otu}}
#'}
#' @import stats
#' @import graphics
#' @import phyloseq
#' @import MASS
#' @import ICSNP
#' @import foreach
#' @import doParallel
#' @export
#' @examples
#' library(CAFT)
#' data(Colon)
#'
#' library(phyloseq)
#'
#' count.tab = t(as.data.frame(as.matrix(otu_table(Colon))))
#' sample.tab = as.data.frame(as.matrix(sample_data(Colon)))
#' tax.tab = as.data.frame(as.matrix(tax_table(Colon)))
#'
#' p = sample.tab$study_name %in% "WirbelJ_2018"
#' sample.tab = sample.tab[p,]
#' count.tab = count.tab[p, ]
#'
#' pNA = which(is.na(sample.tab$age))
#' if(length(pNA) > 0){
#' count.tab = count.tab[-pNA, ]
#'   sample.tab = sample.tab[-pNA,]
#' }
#'
#' ### otu presence filtering
#' p_otu = which(rowSums(t(count.tab) > 0) > 1)
#' count.tab = count.tab[,p_otu]
#' tax.tab = tax.tab[p_otu,]
#'
#' Disease = as.numeric(factor(sample.tab$disease, levels = c("healthy", "CRC"))) - 1
#' Age = as.numeric(sample.tab$age)
#' Gender = as.numeric(factor(sample.tab$gender)) - 1
#'
#' # CAFT
#' res.CAFT = caft(otu.table=count.tab, x.test=Disease,
#'                 x.cov=data.frame(Age=Age, Gender=Gender),
#'                 filter.thresh=0.06, adjust.method="BH")
#'
#' # CAFT (parallel version)
#' res.CAFT = caft(otu.table=count.tab, x.test=Disease,
#'                 x.cov=data.frame(Age=Age, Gender=Gender),
#'                 filter.thresh=0.06, adjust.method="BH",
#'                 n.cores=4)
caft <- function(otu.table, x.test = NULL, x.cov = NULL, x = NULL, Gamma = NULL, filter.thresh = 0.05, fdr.nominal = 0.20,
                   adjust.method = "BH", n.cores=1L) {
  if (!is.matrix(otu.table)) {
    otu.table <- as.matrix(otu.table)
  }
  if (is.null(x) & is.null(x.test)) {
    stop("No data provided: Either x or x.test (and possibly x.cov) required")
  }
  if (!is.null(x) & !is.null(x.test)) {
    stop("Only one of x or x.test can be specified, not both")
  }
  if (is.null(x)) {
    x <- cbind(x.test, x.cov)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (NROW(otu.table) != NROW(x)) {
    stop(" Number of samples not match between OTU table and covairates matrix!")
  }
  if (!is.null(Gamma)) {
    if (is.null(x)) {
      stop("Gamma should not be specified if x.test (and possibly x.cov) is specified")
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
    warnings("Missing values deleted!")
    otu.table <- otu.table[-pNA, ]
    x <- x[-pNA, ]
  }

  # center the covariates
  x <- as.matrix(x - t(replicate(NROW(x), colMeans(x))))

  n.data <- NROW(otu.table)
  n.taxa <- NCOL(otu.table)
  n.param <- NCOL(x)
  if (is.null(Gamma)) {
    n.test <- NCOL(x.test)
  } else {
    n.test <- nrow(Gamma)
  }

  taxa.name <- colnames(otu.table)
  if (is.null(taxa.name)) {
    taxa.name <- paste("tstar", 1:n.taxa)
  }
  ra.all <- otu.table / rowSums(otu.table)
  lib.size <- rowSums(otu.table)
  # Censored relative abundance: every row should same
  lim.ra <- matrix(1 / lib.size, nrow = n.data, ncol = n.taxa)

  # set up Gamma and Lambda matrices
  if (is.null(Gamma)) {
    Gamma <- diag(n.param)[1:n.test, , drop = FALSE]
    if (n.test < n.param) {
      Lambda <- diag(n.param)[(n.test + 1):n.param, , drop = FALSE]
    } else {
      Lambda <- NULL
    }
  } else {
    if (n.test < n.param) {
      Lambda <- t(MASS::Null(t(Gamma)))
    } else {
      Lambda <- NULL
    }
  }
  Gamma.ginv <- MASS::ginv(Gamma)
  Lambda.ginv <- MASS::ginv(Lambda)

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

  #--------------------------
  # unconstrained parameter estimates
  #--------------------------

  if(n.cores > 1L){

    doParallel::registerDoParallel(core=n.cores)

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
        test = TRUE, regularize = TRUE, tol = 1e-12
      ), silent = TRUE)

      if (inherits(fit0.pen, "try-error")) {
        list(beta = rep(NA_real_, ncoef), skip_rare = 0L, skip_fail_fit = 1L)
      } else {
        list(beta = as.numeric(fit0.pen$beta), skip_rare = 0L, skip_fail_fit = 0L)
      }
    }

    est.rank.gs.pen <- do.call(rbind, lapply(res_phase1, `[[`, "beta"))
    colnames(est.rank.gs.pen) <- paste0("b", 1:NCOL(x), ".est")
    skip.rare <- vapply(res_phase1, function(z) z$skip_rare, integer(1))
    skip.fail.rank.fit.pen <- vapply(res_phase1, function(z) z$skip_fail_fit, integer(1))

    if (n.test==1) {
      b.median.pen = median( as.matrix(est.rank.gs.pen) %*% t(Gamma), na.rm = T)
    } else {
      b.median.pen = ICSNP::HR.Mest(as.matrix(est.rank.gs.pen) %*% t(Gamma), na.action=na.omit)$center
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
        b = b.median.pen, beta = NULL,
        test = TRUE, regularize = TRUE, tol = 1e-12
      ), silent = TRUE)

      if (inherits(res.pen, "try-error")) {
        return(list(
          p = NA_real_, test = NA_real_,
          z = rep(NA_real_, n.test), df = NA_real_,
          beta_r = rep(NA_real_, ncoef),
          skip_fail_test = 1L
        ))
      }

      temp.rank <- try(test.rank.aft(res.pen, score = "rank"), silent = TRUE)
      if (inherits(temp.rank, "try-error")) {
        list(
          p = NA_real_, test = NA_real_,
          z = rep(NA_real_, n.test), df = NA_real_,
          beta_r = rep(NA_real_, ncoef),
          skip_fail_test = 1L
        )
      } else {
        list(
          p = as.numeric(temp.rank$p.value),
          test = as.numeric(temp.rank$test),
          z = as.numeric(temp.rank$z.score),
          df = as.numeric(temp.rank$df),
          beta_r = as.numeric(res.pen$beta.r),
          skip_fail_test = 0L
        )
      }
    }

    p.rank.pen            <- vapply(res_phase2, function(z) z$p,    numeric(1))
    test.rank.pen         <- vapply(res_phase2, function(z) z$test, numeric(1))
    df.rank.pen           <- vapply(res_phase2, function(z) z$df,   numeric(1))
    skip.fail.rank.test.pen <- vapply(res_phase2, function(z) z$skip_fail_test, integer(1))
    test.rank.pen.norm <- do.call(rbind, lapply(res_phase2, `[[`, "z"))
    if (!is.matrix(test.rank.pen.norm)) test.rank.pen.norm <- matrix(test.rank.pen.norm, nrow = n.taxa, ncol = n.test, byrow = TRUE)
    est.rank.gs.pen.r <- do.call(rbind, lapply(res_phase2, `[[`, "beta_r"))
    if (!is.matrix(est.rank.gs.pen.r)) est.rank.gs.pen.r <- matrix(est.rank.gs.pen.r, nrow = n.taxa, ncol = ncoef, byrow = TRUE)

  }else{

    est.rank.gs.pen <- est.rank.gs.pen.r <- as.data.frame(matrix(NA, n.taxa, NCOL(x)))
    colnames(est.rank.gs.pen) <- paste0("b", 1:NCOL(x), ".est")

    skip.rare <- skip.fail.rank.fit <- rep(0, n.taxa)
    skip.fail.rank.fit.pen <- rep(0, n.taxa)

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
          test = TRUE, regularize = T,
          tol = 10^-12
        ))
        if (inherits(fit0.pen, "try-error")) {
          est.rank.gs.pen[ii, ] <- rep(NA, NCOL(x))
          skip.fail.rank.fit.pen[ii] <- 1
        } else {
          est.rank.gs.pen[ii, ] <- fit0.pen$beta
        }
      }
    }

    #--------------------------
    # test equal to median
    #--------------------------

    if (n.test==1) {
      b.median.pen = median( as.matrix(est.rank.gs.pen) %*% t(Gamma), na.rm = T)
    }else {
      b.median.pen = ICSNP::HR.Mest(as.matrix(est.rank.gs.pen) %*% t(Gamma), na.action=na.omit)$center
    }
    if (n.test==n.param) Lambda=NULL

    test.rank.pen <- df.rank.pen <- rep(NA, n.taxa)
    test.rank.pen.norm <- matrix(NA, ncol = n.test, nrow = n.taxa)
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
          b = b.median.pen,
          beta = NULL, test = TRUE,
          regularize = T, tol = 10^-12
        ))
        if (inherits(res.pen, "try-error")) {
          p.rank.pen[ii] <- NA
          test.rank.pen[ii] <- NA
          test.rank.pen.norm[ii, ] <- NA
          skip.fail.rank.test.pen[ii] <- 1
        } else {
          temp.rank <- try(test.rank.aft(res.pen, score = "rank"))
          if (inherits(temp.rank, "try-error")) {
            p.rank.pen[ii] <- NA
            test.rank.pen.norm[ii, ] <- NA
            test.rank.pen[ii] <- NA
            skip.fail.rank.test.pen[ii] <- 1
          } else {
            p.rank.pen[ii] <- temp.rank$p.value
            test.rank.pen[ii] <- temp.rank$test
            test.rank.pen.norm[ii, ] <- temp.rank$z.score
            df.rank.pen[ii] <- temp.rank$df
            est.rank.gs.pen.r[ii, ] <- res.pen$beta.r
          }
        }
      }
    }
  }

  p.otu <- p.rank.pen
  q.otu <- p.adjust(p.otu, method = adjust.method)
  return(list(
    est.rank.gs.pen = est.rank.gs.pen,
    b.median.pen = b.median.pen,
    test.rank = test.rank.pen,
    skip.otu = skip.rare,
    p.otu = p.otu,
    df.test = df.rank.pen,
    est.rank.gs.pen.r = est.rank.gs.pen.r,
    p.detected.otu = colnames(otu.table)[which(p.otu < fdr.nominal)],
    q.otu = q.otu,
    q.detected.otu = colnames(otu.table)[which(q.otu < fdr.nominal)]
  ))
}
