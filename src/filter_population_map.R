#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

popmap_file <- snakemake@input[["popmap"]]
read_stats_file <- snakemake@input[["read_stats"]]
gc_stats_file <- snakemake@input[["gc_stats"]]
filtered_population_map <- snakemake@output[["map"]]
plot_file <- snakemake@output[["plot"]]

# dev
# popmap_file <- "output/stacks_config/population_map.txt"
# read_stats_file <- "output/combined_stats/reads.csv"
# gc_stats_file <- "output/combined_stats/gc_stats.csv"

########
# MAIN #
########

# read the popmap
popmap <- fread(popmap_file,
                header = FALSE,
                col.names = c("sample_name", "population"))
read_stats <- fread(read_stats_file)
gc_stats <- fread(gc_stats_file)

# filter by 0.5e6 < read count < 2e6
keep_reads <- read_stats[!is.na(reads) & 
                             0.5e6 < reads &
                             reads < 2e6,
                         unique(individual)]

# filter by median gc < 45
keep_gc <- gc_stats[Median <= 45, unique(individual)]

# produce a filtered popmap
filtered_popmap <- data.table(individual = intersect(keep_reads, keep_gc))
filtered_popmap[, population := gsub("[^[:alpha:]]+", "", individual)]

# do a log4 transform for plotting
hist_data <- read_stats[!is.na(reads),
                        .(l4reads = log(reads, 4)),
                        by = individual]
bins <- 100
qa <- hist_data[, seq(min(l4reads),
                      max(l4reads),
                      length.out = bins)]
bin_labels <- sapply(1:length(qa), function(i)
    mean(c(qa[i], qa[i + 1])))[c(1:bins - 1)]
hist_data[, lbin := as.numeric(
    as.character(
        cut(l4reads, breaks = qa,
            labels = bin_labels,
            include.lowest = TRUE)))]
lhist_log4 <- hist_data[, .(count = length(individual)), by = lbin]

# plot the histogram
filter_level <- as.numeric(cut(log(0.5e6, 4), breaks = qa, labels = bin_labels))
filter_plot <- mean(bin_labels[c(filter_level, filter_level + 1)])

gt <- paste0(
    "'",
    read_stats[, length(unique(individual))],
    " individuals. '*",
    "'Filter: ", 0.5e6, " reads. ",
    filtered_popmap[, length(unique(individual))],
    " individuals passed.'"
)
gp <- ggplot(lhist_log4, aes(x = lbin, y = count)) +
    scale_x_continuous(labels = function(x) 4^x) +
    ylab("Sample count") + xlab("Number of reads") +
    ggtitle(parse(text = gt)) +
    geom_vline(xintercept = filter_plot,
               linetype = 2) +
    geom_col()

# write output
fwrite(filtered_popmap,
       filtered_population_map,
       col.names = FALSE,
       sep = "\t")

ggsave(filename = plot_file,
       plot = gp,
       device = "pdf",
       width = 10,
       height = 7.5,
       units = "in")

# log
sessionInfo()
