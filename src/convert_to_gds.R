#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(SNPRelate)

snpgdsVCF2GDS(snakemake@input[[1]],
    snakemake@output[[1]],
    method = "biallelic.only",
    verbose = TRUE)

sessionInfo()
