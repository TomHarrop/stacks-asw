#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(adegenet)
library(data.table)


###########
# GLOBALS #
###########

vcf <- snakemake@input[["vcf"]]
imiss_file <- snakemake@input[["imiss"]]
whitelist_file <- snakemake@output[["whitelist"]]
popmap_file <- snakemake@output[["popmap"]]
indiv_missing_rate <- as.numeric(snakemake@params[["imiss_rate"]])

# dev
# vcf <- "output/popgen/mapped/locusfilter.vcf"
# imiss_file <- "output/popgen/mapped/stats_locusfilter.imiss"
# indiv_missing_rate <- 0.2

########
# MAIN #
########

snp_data <- VariantAnnotation::readVcf(vcf)
imiss <- fread(imiss_file)

# individual missingness and population size filters
keep_indiv <- imiss[F_MISS < indiv_missing_rate, unique(INDV)]
imiss[, population := gsub("[^[:alpha:]]+", "", INDV)]
pop_cut <- imiss[F_MISS < indiv_missing_rate,
                 .N,
                 by = population][, quantile(N, 0.05)]
keep_pop <- imiss[F_MISS < indiv_missing_rate,
                  .N,
                  by = population][N > pop_cut, unique(population)]

# make popmap
popmap <- imiss[
    INDV %in% keep_indiv & population %in% keep_pop,
    .(individual = INDV,
      population)]

fwrite(popmap,
       popmap_file,
       sep = "\t",
       col.names = FALSE)

# make whitelist
loc_names <- data.table(names(snp_data))
whitelist <- loc_names[, tstrsplit(V1, "_")]
whitelist[, V3 := NULL]
fwrite(whitelist,
       whitelist_file, sep = "\t",
       col.names = FALSE)

# log
sessionInfo()

