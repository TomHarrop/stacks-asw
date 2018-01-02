#!/usr/bin/env Rscript

library(data.table)

#############
# FUNCTIONS #
#############

GenerateMessage <- function(message.text) {
    message(paste0("[ ", date(), " ]: ", message.text))
}

ParseIndividualLoci <- function(tag_file) {
    # Number of assembled loci:
    # for i in *.tags.tsv.gz; do zcat $i | cut -f 3 | tail -n 1; done
    GenerateMessage(paste0("Reading ",
                           tag_file))
    my_tags <- fread(paste0("zgrep -v '^#' ", tag_file),
                     header = FALSE,
                     sep = "\t")
    my_tags[, length(unique(V3))]
}

ParseIndividualPolymorphicLoci <- function(allele_file) {
    # Number of polymorphic loci:
    # for i in *.alleles.tsv.gz; do zcat $i | grep -v "^#" | cut -f 3 | sort | uniq | wc -l; done
    GenerateMessage(paste0("Reading ",
                           allele_file))
    my_alleles <- fread(paste0("zgrep -v '^#' ", allele_file),
                        header = FALSE,
                        sep = "\t")
    my_alleles[, length(unique(V3))]
}

ParseIndividualSNPs <- function(snp_file) {
    # Number of SNPs:
    # for i in *.snps.tsv.gz; do zcat $i | grep -v "^#" | cut -f 5 | grep -c "E"; done
    GenerateMessage(paste0("Reading ",
                           snp_file))
    my_snps <- fread(paste0("zgrep -v '^#' ", snp_file),
                     header = FALSE,
                     sep = "\t")
    dim(my_snps[V5 == "E"])[1]
}

###########
# GLOBALS #
###########

tags_file <- snakemake@input[["tags_file"]]
snps_file <- snakemake@input[["snps_file"]]
alleles_file <- snakemake@input[["alleles_file"]]

output_sample_stats <- snakemake@output[["sample_stats"]]
log_file <- snakemake@log[["log"]]

# dev
# tags_file <- "output/stacks_denovo/Coromandel23.tags.tsv.gz"
# alleles_file <- "output/stacks_denovo/Coromandel23.alleles.tsv.gz"
# snps_file <- "output/stacks_denovo/Coromandel23.snps.tsv.gz"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# get individual name
individual <- gsub("^([^\\.]+).+", "\\1", basename(tags_file))

# run the counts
sample_stats <- data.table(
    individual = individual,
    assembled_loci = ParseIndividualLoci(
        tag_file = tags_file),
    polymorphic_loci = ParseIndividualPolymorphicLoci(
        allele_file = alleles_file),
    snps = ParseIndividualSNPs(
        snp_file = snps_file))

# write output
fwrite(sample_stats, output_sample_stats)

# write session info
sessionInfo()