#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(SNPRelate)

###########
# GLOBALS #
###########

vcf_file <- "output/stacks_populations/for_pca/populations.snps.vcf"
gds_file <- "temp/gstacks.gds"
threads <- 8

########
# MAIN #
########

if(!file.exists(gds_file)) {
    # convert the VCF to GDS (temporary)
    snpgdsVCF2GDS(vcf_file,
                  gds_file,
                  method = "biallelic.only",
                  verbose = TRUE)
}

gds <- snpgdsOpen(gds_file)

pca <- snpgdsPCA(gds,
                 autosome.only = FALSE,
                 num.thread = threads)

# set up plotting data
plotdata <- data.table(sample_id = pca$sample.id,
                       PCA1 = pca$eigenvect[, 1],
                       PCA2 = pca$eigenvect[, 2])
plotdata[, population := factor(gsub("[[:digit:]]", "", sample_id))]
pv_labs <- paste0(
    c("PCA1", "PCA2"),
    " (",
    signif(pca$varprop[c(1, 2)] * 100, 3),
    "%)")
ggplot(plotdata, aes(x = PCA1, y = PCA2, colour = population)) +
    xlab(pv_labs[1]) + ylab(pv_labs[2]) +        
    # scale_color_brewer(palette = "Set1",
    #                    guide = guide_legend(title = NULL)) +
    geom_point(size = 4, alpha = 0.8)


