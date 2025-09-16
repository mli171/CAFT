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
#' @param Y the covariates are of interest. It can be a vector, matrix, or data
#'    frame. A K-level categorical variable will need to be binarized to K-1
#'    covariates.
#' @param C other covariates that need to be adjusted. Requirements are the
#'    same as \code{Y} above.
#' @param filter.thresh a real value between 0 and 1 for OTU table sample presence
#'    filtering. Any OTUs present in fewer than \code{filter.thresh} proportion
#'    of samples are filtered out. We set the default to be 0.05.
#' @param fdr.nominal the nominal FDR. The default is 0.2.
#' @param adjust.method a character string. Use multiple comparison/testing
#'  adjustment methods to control the family-wise error rate/false discover
#'  rate. Default to "BH". See \code{\link{p.adjust}} for the details.
#'
#' @return Return a list consisting of
#' \describe{
#' \item{est.rank.gs.pen}{A matrix that include the estimated coefficients by fitting
#' each OTU count to the log-linear model. The number of columns should match
#' the number of covairates from \code{Y} and \code{C}.}
#' \item{b1.median.pen}{The selected median of all betas, which was used in the
#' restricted score test of the significance of differential abundant OTU. It can
#' also be used to calculated the effect size of each OTU, defined as the difference
#' from the median of all betas.}
#' \item{test.rank}{test statistics from the proposed score test}
#' \item{skip.otu}{the names of skipped OTU during taxa presence filtering
#'          above}
#' \item{p.otu}{p-values for individual OTU association tests}
#' \item{p.detected.otu}{detected significantly differential abundant taxa
#'          (denoted by the column names of the OTU table) at the nominal FDR
#'          based on \code{p.otu}}
#' \item{q.otu}{q-values (adjusted p-values by the Bonferroni correction or the
#'    adjustment method specified in \code{adjust.method}.)
#'          for individual OTU association tests}
#' \item{q.detected.otu}{detected significantly differential abundant taxa
#'          (denoted by the column names of the OTU table) at the nominal FDR
#'          based on \code{q.otu}}
#'}
#' @import stats
#' @import graphics
#' @import phyloseq
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
#' res.CAFT = caft(otu.table=count.tab, Y=Disease,
#'                C=data.frame(Age=Age, Gender=Gender),
#'                filter.thresh=0.06, adjust.method="BH")
#'
caft = function(otu.table,
                Y,
                C,
                filter.thresh=0.05,
                fdr.nominal=0.20,
                adjust.method="BH"){

  if (is.matrix(otu.table)){otu.table = as.matrix(otu.table)}

  x = cbind(Y, C)
  if (is.matrix(x)) {x = as.matrix(x)}
  if (NROW(otu.table) != NROW(x)) {stop("\n\n Number of samples not match
                                        between OTU table and covairates
                                        matrix!")}

  # missing values
  pNA.x = which(apply(x, 1, function(x) any(is.na(x)))) # covariates
  pNA.otu = which(apply(otu.table, 1, function(x) any(is.na(x)))) # covariates
  pNA = unique(c(pNA.x, pNA.otu))
  if(length(pNA) > 0){
    warnings("\n\n Missing values deleted!")
    otu.table = otu.table[-pNA, ]
    x = x[-pNA,]
  }

  # center the covariates
  x = as.matrix(x - t(replicate(NROW(x), colMeans(x))))

  n.data = NROW(otu.table)
  n.taxa = NCOL(otu.table)

  taxa.name = colnames(otu.table)
  if(is.null(taxa.name)){taxa.name = paste("tstar", 1:n.taxa)}

  ra.all = otu.table/rowSums(otu.table)
  lib.size = rowSums(otu.table)
  # Censored relative abundance: every row should same
  lim.ra = matrix(1/lib.size, nrow = n.data, ncol = n.taxa)

  #--------------------------
  # create survival data
  #--------------------------
  log_neg <- -log10(lim.ra)
  log_pos <- -log10(ra.all)
  t.star.all <- matrix(log_neg,
                       nrow = nrow(otu.table),
                       ncol = ncol(otu.table))
  idx <- (otu.table > 0)
  t.star.all[idx] <- log_pos[idx]
  t.star.all = as.data.frame(t.star.all)
  colnames(t.star.all) = paste("tstar", 1:n.taxa)

  delta.all = ((otu.table > 0)^2)
  colnames(delta.all) = paste("delta", 1:n.taxa)

  data = data.frame(t.star.all, delta.all)

  est.rank.gs.pen = as.data.frame(matrix(NA, n.taxa, NCOL(x)))
  colnames(est.rank.gs.pen) = paste0("b", 1:NCOL(x), ".est")

  skip.rare = skip.fail.rank.fit = rep(0, n.taxa)
  skip.fail.rank.fit.pen = rep(0, n.taxa)

  #--------------------------
  # parameter estimates
  #--------------------------

  for (ii in 1:n.taxa){
    tstar = data[, grep("tstar", colnames(data))[ii]]
    delta.1 = data[, grep("delta", colnames(data))[ii]]
    if (sum(delta.1) <= n.data*filter.thresh){
      out = rep(NA, NCOL(x))
      skip.rare[ii]  = 1
    }else{
      fit0.pen = try(estimate.rank.aft(y=tstar, delta=delta.1, x=x,
                                       Gamma=NULL, b=NULL, beta=NULL,
                                       test=TRUE, regularize = T,
                                       tol=10^-12))
      if (inherits(fit0.pen, "try-error")){
        est.rank.gs.pen[ii, ] = rep(NA, NCOL(x))
        skip.fail.rank.fit.pen[ii] = 1
      }else{
        est.rank.gs.pen[ii, ]  = fit0.pen$beta
      }
    }
  }

  #--------------------------
  # test equal to median
  #--------------------------

  b1.median.pen = median(est.rank.gs.pen$b1.est, na.rm = T)
  test.rank.pen = test.rank.pen.norm = rep(NA, n.taxa)
  p.rank.pen  = rep(NA, n.taxa)
  skip.fail.rank.test.pen = rep(NA, n.taxa)
  for (ii in 1:n.taxa){
    tstar = data[, grep("tstar", colnames(data))[ii]]
    delta.1 = data[, grep("delta", colnames(data))[ii]]
    if (sum(delta.1) > n.data*filter.thresh){
      # first covariate is the testing covariate
      res.pen = try(estimate.rank.aft(y=tstar, delta=delta.1, x=x,
                                      Gamma=c(1,rep(0,NCOL(x)-1)),
                                      b=b1.median.pen,
                                      beta=NULL, test=TRUE,
                                      regularize=T, tol=10^-12))
      if (inherits(res.pen, "try-error")){
        p.rank.pen[ii] = NA
        test.rank.pen[ii] = NA
        test.rank.pen.norm[ii] = NA
        skip.fail.rank.test.pen[ii] =  1
      }else{
        temp.rank = try(test.rank.aft(res.pen, score="rank"))
        if (inherits(temp.rank, "try-error")){
          p.rank.pen[ii] = NA
          test.rank.pen.norm[ii] = NA
          test.rank.pen[ii] = NA
          skip.fail.rank.test.pen[ii] =  1
        }else{
          p.rank.pen[ii] = temp.rank$p.value
          test.rank.pen[ii] = temp.rank$test
          test.rank.pen.norm[ii] = temp.rank$z.score
        }
      }
    }
  }

  return(list(est.rank.gs.pen=est.rank.gs.pen,
              b1.median.pen=b1.median.pen,
              test.rank=test.rank.pen,
              skip.otu=skip.rare,
              p.otu = p.rank.pen,
              p.detected.otu = colnames(otu.table)[which(p.otu < fdr.nominal)],
              q.otu = p.adjust(p.otu, method=adjust.method),
              q.detected.otu = colnames(otu.table)[which(q.otu < fdr.nominal)]))
}
