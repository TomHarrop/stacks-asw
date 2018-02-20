#!/usr/bin/env Rscript

library(data.table)

#############
# FUNCTIONS #
#############

GenerateMessage <- function(message.text) {
    message(paste0("[ ", date(), " ]: ", message.text))
}

CalculateCovstats <- function(tags_file) {
    # read the file
    GenerateMessage(paste0("Reading ", tags_file))
    my_tags <- fread(paste0("zgrep -v '^#' ", tags_file))
    
    # get the individual name
    individual <- gsub("^([^\\.]+).+", "\\1", basename(tags_file))
    
    # get a list of non-blacklisted loci
    my_loci <- my_tags[V12 != 1 & V13 != 1, unique(V3)]
    
    # count the unique reads
    total_reads <- my_tags[V3 %in% my_loci & V9 != "",
                           length(unique(V9))]
    final_coverage <- my_tags[V3 %in% my_loci & V9 != "",
                              length(unique(V9)),
                              by = V3][, mean(V1)]
    
    # make the table
    return(
        data.table(
            individual = individual,
            final_coverage_mean = final_coverage,
            n_reads = total_reads))
}

###########
# GLOBALS #
###########

tags_file <- snakemake@input[["tags_file"]]
output_covstats <- snakemake@output[["covstats"]]
log_file <- snakemake@log[["log"]]

#dev
# tags_file <- "output/03_ustacks/Ruakura_9.tags.tsv.gz"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# run calc
covstats <- CalculateCovstats(tags_file)

# write output
fwrite(covstats, output_covstats)

# write session info
sessionInfo()


