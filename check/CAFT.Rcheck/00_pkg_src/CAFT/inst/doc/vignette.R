## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----include=FALSE------------------------------------------------------------
library(CAFT)

## ----eval=FALSE---------------------------------------------------------------
# pak::pak("CAFT_1.0.tar.gz") # local directory

## ----eval=FALSE---------------------------------------------------------------
# remotes::install_github("mli171/CAFT", dependencies = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# browseVignettes("CAFT")

## ----eval=FALSE---------------------------------------------------------------
# vignette("vignette", package = "CAFT")

## ----eval=FALSE---------------------------------------------------------------
# caft(otu.table, x.test = NULL, x.adj = NULL, x = NULL, Gamma = NULL,
#   filter.thresh = 0.05, fdr.nominal = 0.2, adjust.method = "BH",
#   regularize = TRUE, n.cores = 1L)

## -----------------------------------------------------------------------------
data("URT")

throat.otu.table    <- URT$otu
throat.meta         <- URT$meta
throat.otu.taxonomy <- URT$tax

filter.out.sam = which(throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage != "None")
throat.otu.table.filter = throat.otu.table[-filter.out.sam,]
throat.meta.filter = throat.meta[-filter.out.sam,]

dim(throat.otu.table.filter)

cens.prop = colMeans(throat.otu.table.filter == 0, na.rm = T)
mean(cens.prop)

## -----------------------------------------------------------------------------
x.test = ifelse(throat.meta.filter$SmokingStatus == "NonSmoker", 0, 1)
x.adj  = ifelse(throat.meta.filter$Sex == "Male", 0, 1)
res.CAFT = caft(otu.table=throat.otu.table.filter, x.test=x.test, x.adj=x.adj, 
                filter.thresh = 0.10, fdr.nominal = 0.10)

## -----------------------------------------------------------------------------
res.CAFT$q.detected.otu

## ----fig.width=7, fig.height=5, out.width="80%", fig.align='center'-----------
groups = interaction(x.test, x.adj, sep = "_", drop = TRUE)
groups = factor(groups,
                labels = c("Non-Smoker\nMale", "Smoker\nMale",
                           "Non-Smoker\nFemale", "Smoker\nFemale"))

boxplot.otu = "2831"
throat.otu.taxonomy[which(colnames(throat.otu.table.filter) == boxplot.otu),]
otuboxplot(plot.otu=boxplot.otu, count.data=throat.otu.table.filter, 
           plot.title = "<i>Prevotellaceae Prevotella</i>",groups=groups)

## -----------------------------------------------------------------------------
library(phyloseq)

data(Colon)

count.tab = t(as.data.frame(as.matrix(otu_table(Colon))))
sample.tab = as.data.frame(as.matrix(sample_data(Colon)))
tax.tab = as.data.frame(as.matrix(tax_table(Colon)))

dim(count.tab)

## -----------------------------------------------------------------------------
pNA = which(is.na(sample.tab$age))
if(length(pNA) > 0){
  count.tab = count.tab[-pNA, ]
  sample.tab = sample.tab[-pNA,]
}
# No missing values from gender

## otu presence filtering
p_otu = which(rowSums(t(count.tab) > 0) > 1)
count.tab = count.tab[,p_otu]
tax.tab = tax.tab[p_otu,]

dim(count.tab)

cens.prop = colMeans(count.tab == 0, na.rm = T)
mean(cens.prop)

## -----------------------------------------------------------------------------
Disease1 = Disease2 = rep(0, NROW(sample.tab)) # healthy
Disease1[sample.tab$disease == "CRC"] = 1
Disease2[sample.tab$disease == "adenoma"] = 1

Age = as.numeric(sample.tab$age)
Gender = as.numeric(factor(sample.tab$gender)) - 1

x.test = cbind(Disease1, Disease2)
x.adj  = cbind(Age, Gender)

## -----------------------------------------------------------------------------
res.CAFT = caft(otu.table=count.tab, x.test=x.test, x.adj=x.adj)

## ----eval=FALSE---------------------------------------------------------------
# res.CAFT = caft(otu.table=count.tab, x.test=x.test, x.adj=x.adj, n.cores=4)

## -----------------------------------------------------------------------------
x = cbind(x.test, x.adj)
Gamma = matrix(c(1,-1,0,0), nrow=1, ncol=4)
Gamma
res.CAFT = caft(otu.table=count.tab, x=x, Gamma=Gamma)

## -----------------------------------------------------------------------------
Gamma = matrix(c(1,1,-1,0,0,0,0,0), nrow=2, ncol=4)
Gamma
res.CAFT = caft(otu.table=count.tab, x=x, Gamma=Gamma)

## -----------------------------------------------------------------------------
sessionInfo()

