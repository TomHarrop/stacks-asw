#!/usr/bin/env Rscript

library(data.table)

#############
# FUNCTIONS #
#############

GenerateMessage <- function(message.text) {
    message(paste0("[ ", date(), " ]: ", message.text))
}

ParsePopulationsStats <- function(stacks_dir){
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
    
    # check sumstats exists
    if (length(list.files(stacks_dir,
                          pattern = "populations.sumstats.tsv")) == 0) {
        stop(paste("populations.sumstats.tsv not found in", stacks_dir))
    }
    
    # check if haplotypes exists
    if (length(list.files(stacks_dir,
                          pattern = "populations.haplotypes.tsv")) == 0) {
        stop(paste("populations.haplotypes.tsv not found in", stacks_dir))
    }
    
    # parse hapstats
    GenerateMessage(paste0("Reading ",
                           stacks_dir,
                           "/populations.haplotypes.tsv"))
    my_hapstats <- fread(paste0("grep -v '^#' ",
                                stacks_dir,
                                "/populations.haplotypes.tsv"))
    my_loci <- my_hapstats[, length(unique(V1))]
    
    # parse sumstats
    GenerateMessage(paste0("Reading ",
                           stacks_dir,
                           "/populations.sumstats.tsv"))
    my_sumstats <- fread(paste0("grep -v '^#' ",
                                stacks_dir,
                                "/populations.sumstats.tsv"))
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

stacks_dir <- snakemake@params[["pop_dir"]]
output_pop_stats <- snakemake@output[["pop_stats"]]
log_file <- snakemake@log[["log"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


# get the populations summary
population_stats <- ParsePopulationsStats(stacks_dir)



# write session info
sessionInfo()
