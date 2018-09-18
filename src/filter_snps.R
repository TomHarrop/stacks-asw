#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(SNPRelate)

###########
# GLOBALS #
###########

gds_file <- snakemake@input[[1]]
ped_file <- snakemake@params[["ped_file"]]
maf <- snakemake@params[["maf"]]
missing_rate <- snakemake@params[["missing_rate"]]
sample_missing_quantile <- snakemake@params[["sample_missing_quantile"]]


#dev
# gds_file <- 'output/060_pop_genet/snps.gds'
# maf <- 0.05
# missing_rate <- 0.2
# sample_missing_quantile <- 0.8


########
# MAIN #
########

# Load the data
gds_data <- snpgdsOpen(gds_file)

# Filter the SNPs
snp_set <- snpgdsSelectSNP(gds_data,
                           maf = maf,
                           autosome.only = FALSE,
                           verbose = TRUE,
                           missing.rate = missing_rate)
snp_set_ids <-  unlist(snp_set)

# Get samples where missing rate is higher than sample missing quantile
sample_missing_rates <- snpgdsSampMissRate(gdsobj = gds_data,
                                           snp.id = snp_set_ids,
                                           with.id = TRUE)

kept_indivs <- sample_missing_rates[sample_missing_rates < quantile(sample_missing_rates,
                                                                    sample_missing_quantile)]

# Write the PED
snpgdsGDS2PED(gds_data,
              ped.fn = ped_file,
              snp.id = snp_set_ids,
              sample.id = names(kept_indivs),
              verbose = TRUE)

# Log session info
sessionInfo()

