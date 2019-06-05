#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(adegenet)
library(data.table)

###########
# GLOBALS #
###########

plink_file <- snakemake@input[["plink"]]
whitelist_file <- snakemake@output[["whitelist"]]
popmap_file <- snakemake@output[["popmap"]]

# dev
# plink_file <- "output/popgen/plink.raw"

########
# MAIN #
########

snp_data <- read.PLINK(plink_file)

# make popmap
ind_names <- indNames(snp_data)
popmap <- data.table(individual = ind_names,
           population = gsub("[[:digit:]]+", "", ind_names))
fwrite(popmap,
       popmap_file, sep = "\t",
       col.names = FALSE)


# make whitelist
loc_names <- data.table(locNames(snp_data))
whitelist <- loc_names[, tstrsplit(V1, "_")]
whitelist[, V3 := NULL]
fwrite(whitelist,
       whitelist_file, sep = "\t",
       col.names = FALSE)

# log
sessionInfo()

