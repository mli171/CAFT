# CAFT
Rank-Based Compositional Analysis using Log-Linear Models for Microbiome Data with Zero Cells.

## Overview
In this study, we introduce a novel statistical framework for differential
abundance analysis of microbiome data, termed the Compositional Accelerated
Failure Time (CAFT) model. The CAFT model addresses zero read counts by
treating them as censored observations below the detection limit, similar
to censoring mechanisms employed in survival analysis. This approach is inherently 
resistant to multiplicative bias, eliminates the need for pseudocounts, and
addresses compositional bias through the establishment of appropriate score
test procedures. For FDR control, we utilize and expand the idea from Efron’s
empirical null distribution to achieve better FDR control.

## Package download and installation
You can install the version of CAFT from Github:

```{r}
# install.packages("devtools")
devtools::install_github("mli171/CAFT")
```

## Open the Vignette in R

```{r}
browseVignettes("CAFT")
```

## CAFT: Compositional Rank-Based Analysis using AFT models

The main function in CAFT package is:

```{r}
caft()
```

## An example of using the caft function

Apply 'caft' to a dataset from the study of gut microbiome data set focusing on
the adult colorectal cancer using the stool samples [(Pasolli et al.,2017)](https://www.nature.com/articles/nmeth.4468).

```{r}
library(CAFT)
data(Colon)

library(phyloseq)

count.tab = t(as.data.frame(as.matrix(otu_table(Colon))))
sample.tab = as.data.frame(as.matrix(sample_data(Colon)))
tax.tab = as.data.frame(as.matrix(tax_table(Colon)))

p = sample.tab$study_name %in% "WirbelJ_2018"
sample.tab = sample.tab[p,]
count.tab = count.tab[p, ]

pNA = which(is.na(sample.tab$age))
if(length(pNA) > 0){
  count.tab = count.tab[-pNA, ]
  sample.tab = sample.tab[-pNA,]
}

### otu presence filtering
p_otu = which(rowSums(t(count.tab) > 0) > 1)
count.tab = count.tab[,p_otu]
tax.tab = tax.tab[p_otu,]

Disease = as.numeric(factor(sample.tab$disease, levels = c("healthy", "CRC"))) - 1
Age = as.numeric(sample.tab$age)
Gender = as.numeric(factor(sample.tab$gender)) - 1

# CAFT
res.CAFT = caft(otu.table=count.tab, Y=Disease,
                C=data.frame(Age=Age, Gender=Gender),
                filter.thresh=0.06, adjust.method="BH")
```

## References

Pasolli E, Schiffer L, Manghi P, Renson A, Obenchain V, Truong D, Beghini F, 
Malik F, Ramos M, Dowd J, Huttenhower C, Morgan M, Segata N, Waldron L (2017). 
“Accessible, curated metagenomic data through ExperimentHub.” *Nat. Methods*, 
14(11), 1023–1024. ISSN 1548-7091, 1548-7105, doi:10.1038/nmeth.4468.

