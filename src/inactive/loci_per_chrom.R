#!/usr/bin/env Rscript

library(data.table)

vcf_file <- "output/060_popgen/populations.ns.all.vcf"
fai_file <- "output/005_ref/ref.fasta.fai"

# get a list of loci to match with bs results
vcf_names <- c("chrom", "pos", "id", "ref", "alt", "qual")
loci <- fread(cmd = paste('grep -v "^#"', vcf_file),
              header = FALSE,
              col.names = vcf_names,
              select = 1:6)
loci_per_chrom <- loci[, .(N = length(unique(pos))), by = chrom]
fwrite(loci_per_chrom, "data/loci_per_chrom.csv")


