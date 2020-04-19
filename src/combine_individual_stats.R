#!/usr/bin/env Rscript

# set log
log <- file(snakemake$log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)


#############
# FUNCTIONS #
#############


GetGcStats <- function(gc_file) {
    my_gc_stats <- fread(cmd = paste('grep "^#"', gc_file),
                         col.names = c("variable", "value"),
                         nrows = 4)
    my_gc_stats[, variable := sub("#", "", variable)]
    my_gc_stats_wide <- dcast(my_gc_stats, 1 ~ variable)[, `.` := NULL]
    return(my_gc_stats_wide)
}


GetReadCounts <- function(read_file) {
    my_stats <- fread(cmd = paste('grep "Input:"', read_file),
                      header = FALSE,
                      col.names = c("V1", "reads", "bases"))
    return(my_stats[, .(reads = as.numeric(gsub("[^[:digit:]]", "", reads)))])
}


###########
# GLOBALS #
###########

reads_files <- snakemake@input[["read_stats"]]
gc_files <- snakemake@input[["gc_stats"]]

read_stats_file <- snakemake@output[["read_stats"]]
gc_stats_file <- snakemake@output[["gc_stats"]]
gc_hist_file <- snakemake@output[["gc_hist"]]

# dev
# reads_files <- list.files("output/040_stats/reads",
#                           pattern = ".txt",
#                           full.names = TRUE)
# gc_files <- list.files("output/040_stats/gc_hist",
#                        pattern = ".txt",
#                        full.names = TRUE)

########
# MAIN #
########

# read counts
names(reads_files) <- sub(".txt", "", basename(reads_files))
reads_list <- lapply(reads_files, GetReadCounts)
read_counts <- rbindlist(reads_list, idcol = "individual")

# gc stats
names(gc_files) <- sub(".txt", "", basename(reads_files))
gc_stats_list <- lapply(gc_files, GetGcStats)
gc_stats <- rbindlist(gc_stats_list, idcol = "individual")

# gc hists
gc_hist_list <- lapply(gc_files, function(x)
    fread(cmd = paste('grep -v "#"', x),
          col.names = c("GC", "Count")))
gc_hists <- rbindlist(gc_hist_list, idcol = "individual")

# write output
fwrite(read_counts, read_stats_file)
fwrite(gc_stats, gc_stats_file)
fwrite(gc_hists, gc_hist_file)

# log
sessionInfo()
