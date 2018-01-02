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
    
    # count the unique reads
    total_reads <- my_tags[V9 != "",
                           length(unique(V9))]
    long_stats <- my_tags[V9 != "",
                          length(unique(V9)),
                          by = .(V3, V7)][, mean(V1), by = V7]
    
    # convert output
    covstats <- dcast(long_stats, . ~ V7, value.var = "V1")
    covstats[, `.` := NULL]
    covstats[, individual := individual]
    covstats[, total_reads := total_reads]
    
    return(covstats)
}

###########
# GLOBALS #
###########

tags_file <- snakemake@input[["tags_file"]]
output_covstats <- snakemake@output[["covstats"]]
log_file <- snakemake@log[["log"]]

#dev
# tags_file <- "output/stacks_denovo/Wellington30.tags.tsv.gz"

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


