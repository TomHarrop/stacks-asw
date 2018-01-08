#!/usr/bin/env Rscript

library(data.table)

#############
# FUNCTIONS #
#############

GenerateMessage <- function(message.text) {
    message(paste0("[ ", date(), " ]: ", message.text))
}

ParsePopulationsStats <- function(hapstats_file, sumstats_file){
    # Number of loci
    # cat batch_1.haplotypes.tsv | sed '1d' | wc -l 
    # Reads every locus regardless of population (and pop map specified) don't
    # have to provide 'defaultpop'
    # Polymorphic loci
    # cat batch_1.sumstats.tsv | grep -v "^#" | cut -f 2 | sort -n | uniq | wc -l 
    # Only takes the locus ID but by sorting them means no need to specify a
    # particular pop with a popmap
    # SNPs
    # cat batch_1.sumstats.tsv | grep -v "^#" | cut -f 2,5 | sort -n | uniq | wc -l 
    # Takes the locus ID and the SNP column position and works regardless of
    # which pop the locus was found in
    
    # parse hapstats
    GenerateMessage(paste0("Reading ",
                           hapstats_file))
    my_hapstats <- fread(paste0("grep -v '^#' ",
                                hapstats_file))
    my_loci <- my_hapstats[, length(unique(V1))]
    
    # parse sumstats
    GenerateMessage(paste0("Reading ",
                           sumstats_file))
    my_sumstats <- fread(paste0("grep -v '^#' ",
                                sumstats_file))
    my_polymorphic_loci <- my_sumstats[, length(unique(V1))]
    my_snps <- dim(unique(my_sumstats, by = c("V1", "V4")))[1]
    
    # return stats
    data.table(
        assembled_loci = my_loci,
        polymorphic_loci = my_polymorphic_loci,
        snps = my_snps
    )
}

###########
# GLOBALS #
###########

hapstats_file <- snakemake@input[["hapstats"]]
sumstats_file <- snakemake@input[["sumstats"]]
r <- snakemake@wildcards[["r"]]
output_pop_stats <- snakemake@output[["pop_stats"]]
log_file <- snakemake@log[["log"]]

# dev
# hapstats_file <- "output/stacks_populations/r0.8/populations.hapstats.tsv"
# sumstats_file <- "output/stacks_populations/r0.8/populations.sumstats.tsv"
# r <- "0.8"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# get the populations summary
population_stats <- ParsePopulationsStats(hapstats_file, sumstats_file)
population_stats[, r := r]

# write output
fwrite(population_stats, output_pop_stats)

# write session info
sessionInfo()
