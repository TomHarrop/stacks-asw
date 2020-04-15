#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(VariantAnnotation)

###########
# GLOBALS #
###########

vcf <- snakemake@input[["vcf"]]
imiss_file <- snakemake@input[["imiss"]]
whitelist_file <- snakemake@output[["whitelist"]]
popmap_file <- snakemake@output[["popmap"]]
indiv_missing_rate <- as.numeric(snakemake@params[["imiss_rate"]])
fai_file <- snakemake@input[["fai"]]


# dev
# vcf <- "output/popgen/mapped/locusfilter.vcf"
# imiss_file <- "output/popgen/mapped/stats_locusfilter.imiss"
# indiv_missing_rate <- 0.2
# fai_file <- "output/map_to_genome/draft_genome.fasta.fai"

n_pops <- c("Coromandel", "Ruakura", "Taranaki", "Wellington", "Greymouth")
s_pops <- c("Lincoln", "O", "MararoaDowns", "Mossburn", "Fortrose")

########
# MAIN #
########

snp_data <- readVcf(vcf)
imiss <- fread(imiss_file)

# get the long boiz
fai_names <- c("Chr", "chr_length")
fai <- fread(fai_file, select = 1:2, col.names = fai_names)
keep_chr <- fai[chr_length > 1e5, Chr]

# subset by snps on those chr
snp_set <- snp_data[seqnames(rowRanges(snp_data)) %in% keep_chr, ]

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

popmap[population %in% n_pops, group := "North"]
popmap[population %in% s_pops, group := "South"]

fwrite(popmap,
       popmap_file,
       sep = "\t",
       col.names = FALSE)

# make whitelist
loc_names <- data.table(names(snp_set))
whitelist <- loc_names[, tstrsplit(V1, "_")]
fwrite(whitelist,
       whitelist_file, sep = "\t",
       col.names = FALSE)

# log
sessionInfo()

