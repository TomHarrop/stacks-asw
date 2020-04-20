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
para_data_file <- snakemake@input[["para_data"]]
fai_file <- snakemake@input[["fai"]]

indiv_missing_rate <- as.numeric(snakemake@params[["imiss_rate"]])

whitelist_file <- snakemake@output[["whitelist"]]
geo_file <- snakemake@output[["geo"]]
ns_file <- snakemake@output[["ns"]]
para_file <- snakemake@output[["para"]]

# dev
# vcf <- "output/060_popgen/populations_filtered.vcf"
# imiss_file <- "output/060_popgen/stats_locusfilter.imiss"
# indiv_missing_rate <- 0.2
# fai_file <- "output/005_ref/ref.fasta.fai"
# para_data_file <- "data/para_sample_info.tsv"

n_pops <- c("Coromandel", "Ruakura", "Taranaki", "Wellington", "Greymouth")
s_pops <- c("Lincoln", "O", "MararoaDowns", "Mossburn", "Fortrose")

########
# MAIN #
########

snp_data <- readVcf(vcf)
imiss <- fread(imiss_file)
para_data <- fread(para_data_file)

# get the long boiz
fai_names <- c("Chr", "chr_length")
fai <- fread(fai_file, select = 1:2, col.names = fai_names)
keep_chr <- fai[chr_length > 1e5, Chr]

# subset by snps on those chr
snp_set <- snp_data[seqnames(rowRanges(snp_data)) %in% keep_chr, ]

# make whitelist
loc_names <- data.table(names(snp_set))
whitelist <- loc_names[, tstrsplit(V1, ":")]
whitelist[, V3 := NULL]

# individual missingness and population size filters
keep_indiv <- imiss[F_MISS < indiv_missing_rate, unique(INDV)]
imiss[, population := gsub("[^[:alpha:]_]+", "", INDV)]
pop_cut <- imiss[F_MISS < indiv_missing_rate,
                 .N,
                 by = population][, quantile(N, 0.05)]
keep_pop <- imiss[F_MISS < indiv_missing_rate,
                  .N,
                  by = population][N > pop_cut, unique(population)]

# make popmaps
full_popmap <- imiss[
  INDV %in% keep_indiv & population %in% keep_pop,
  .(individual = INDV,
    population)]

geo_popmap <- full_popmap[startsWith(population, "geo_")]
geo_popmap[, population := sub("geo_", "", population)]
geo_popmap[population %in% n_pops, group := "North"]
geo_popmap[population %in% s_pops, group := "South"]

ns_popmap <- geo_popmap[, .(individual, population = group)]

para_popmap <- para_data[
  , .(individual = paste("para", Individual, sep = "_"),
      population = tolower(sub("-", "", pop_para_past)),
      group = tolower(sub("-", "", Parasitism)))][
        individual %in% keep_indiv]

# write outputs
fwrite(whitelist,
       whitelist_file,
       sep = "\t",
       col.names = FALSE)

fwrite(geo_popmap,
       geo_file,
       sep = "\t",
       col.names = FALSE)

fwrite(ns_popmap,
       ns_file,
       sep = "\t",
       col.names = FALSE)

fwrite(para_popmap,
       para_file,
       sep = "\t",
       col.names = FALSE)

# log
sessionInfo()
