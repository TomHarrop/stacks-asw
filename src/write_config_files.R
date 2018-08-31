#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

full_key_file <- snakemake@input[["key_file"]]
outdir <- snakemake@params[["outdir"]]

#############
# FUNCTIONS #
#############

WriteConfig <- function(flowcell_lane, full_key, outdir) {
    my_path <- paste0(outdir, "/", flowcell_lane, ".config")
    fwrite(full_key[fc_lane == flowcell_lane, .(sample_fullname, barcode)],
           my_path,
           sep = "\t",
           col.names = FALSE)
}

########
# MAIN #
########

full_key <- fread(full_key_file)

# fix spaces
full_key[, sample := gsub("[[:space:]]", "-", sample)]

# add sample details
full_key[, fc_lane := paste(flowcell, lane, sep = "_"),
         by = .(flowcell, lane)]
full_key[, sample_fullname := paste(sample,
                                    flowcell,
                                    lane,
                                    row,
                                    column,
                                    sep = "_"),
         by = .(sample,
                flowcell,
                lane,
                row,
                column)]

# write config files
split_key <- split(full_key, by = "fc_lane")
sapply(names(split_key), function(x)
    fwrite(split_key[[x]][, .(sample_fullname, barcode)],
           paste0(outdir, "/", x, ".config"),
           sep = "\t",
           col.names = FALSE))

# write popmap
fwrite(unique(full_key[control == "" & !grepl("\\.", sample)],
              by = "sample")[, .(
                  sample,
                  tolower(gsub("[[:digit:]]+", "", sample)))],
       paste0(outdir, "/population_map.txt"),
       sep = "\t",
       col.names = FALSE)

# save log
sessionInfo()

