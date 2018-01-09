#!/usr/bin/env Rscript

library(adegenet)
library(data.table)

############
# FUNCTION #
############

###########
# GLOBALS #
###########

populations_genepop <- "output/stacks_populations/for_pca/populations.snps.gen"

########
# MAIN #
########

snps <- read.genepop(populations_genepop)
pop(snps) <- gsub("[[:digit:]]+", "", indNames(snps))

snp_tab <- tab(snps, NA.method = "mean")
snps_pca <- dudi.pca(snp_tab, scannf = FALSE, nf = Inf)
s.class(snps_pca$li, pop(snps), col=rainbow(nPop(snps)))
add.scatter.eig(snps_pca$eig[1:10], xax=1, yax=2)

var_exp <- snps_pca$eig / sum(snps_pca$eig)
pv_labs <- paste0(
    c("PC1", "PC2"),
    " (",
    signif(var_exp[c(1, 2)] * 100, 3),
    "%)")

pd_long <- melt(pd, id.vars = "rn")
pd_long[, PC := as.numeric(gsub("[[:alpha:]]+", "", variable))]
pd_long[, population := gsub("[[:digit:]]+", "", rn)]

ggplot(pd_long[PC < 11], aes(x = population, y = value, colour = population)) +
    facet_wrap(~ PC) +
    geom_point()


pd <- as.data.table(snps_pca$li, keep.rownames = TRUE)
pd[, population := gsub("[[:digit:]]+", "", rn)]

ggplot(pd,
       aes(x = Axis1, y = Axis2, colour = population)) +
    xlab(pv_labs[1]) + ylab(pv_labs[2]) +        
    geom_point()


pd
