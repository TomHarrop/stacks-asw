#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(SNPRelate)

###########
# GLOBALS #
###########

key_file <- 'data/SQ0003.txt'
vcf_file <- "output/stacks_populations/r0.8/populations.snps.vcf"
gds_file <- "temp/gstacks.gds"
threads <- 8

########
# MAIN #
########

# read the key file
key_data <- fread(key_file)
indiv_data <- key_data[, .(
    individual = Sample,
    Flowcell, Lane, Platename)]
indiv_data[, population := gsub("[[:digit:]]+", "", individual)]

if(!file.exists(gds_file)) {
    # convert the VCF to GDS (temporary)
    snpgdsVCF2GDS(vcf_file,
                  gds_file,
                  method = "biallelic.only",
                  verbose = TRUE)
}

# select samples to keep
keep_populations <- c("Coromandel",
                      "Taranaki",
                      "Wellington",
                      "Mossburn",
                      "Fortrose",
                      "Reefton",
                      "Greymouth",
                      "Lincoln",
                      "Ruakura")
keep_indiv <- indiv_data[population %in% keep_populations, unique(individual)]

# read the GDS
gds <- snpgdsOpen(gds_file)

# subset the SNPs
sample_ids <- keep_indiv[keep_indiv %in%
                             read.gdsn(index.gdsn(gds, "sample.id"))]
snpset <- snpgdsSelectSNP(
    gds,
    sample.id = sample_ids,
    maf = 0.05,
    autosome.only = FALSE,
    verbose = TRUE)
snpset_ids <- unlist(snpset)


# snpset <- snpgdsLDpruning(
#     gds,
#     sample.id = sample_ids,
#     autosome.only = FALSE)
# snpset_ids <- unlist(snpset)

pca <- snpgdsPCA(gds,
                 sample.id = sample_ids,
                 snp.id = snpset_ids,
                 autosome.only = FALSE,
                 num.thread = threads)

# set up plotting data
plot_data <- data.table(pca$eigenvect)
setnames(plot_data,
         names(plot_data),
         paste0("EV", c(1:length(names(plot_data)))))
plot_data[, sample_id := pca$sample.id]
plot_data[, population := factor(gsub("[[:digit:]]", "", sample_id))]

long_pd <- melt(plot_data, id.vars = c("sample_id", "population"))
long_pd[, eigenvect := as.numeric(gsub("[^[:digit:]]", "", variable))]

ggplot(long_pd[eigenvect < 10], aes(x = population, y = value, colour = population)) +
    facet_wrap(~ variable) +
    geom_point(position = position_jitter(width = 0.2))
