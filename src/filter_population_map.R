#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

#############
# FUNCTIONS #
#############

TouchFlagFile <- function(individual, flag_dir) {
    flag_path <- paste0(flag_dir, "/", individual, ".tmp")
    file.create(flag_path)
}

###########
# GLOBALS #
###########

popmap_file <- snakemake@input[["popmap"]]
sample_dir <- snakemake@params[["sample_dir"]]
stats_file <- snakemake@input[["stats"]]
filtered_population_map <- snakemake@output[["map"]]
plot_file <- snakemake@output[["plot"]]
pop_counts <- snakemake@output[["pop_counts"]]

# popmap_file <- "output/stacks_config/population_map.txt"
# sample_dir <- "output/demux"
# stats_file <- "output/run_stats/read_stats.txt"

########
# MAIN #
########

# read the popmap
popmap <- fread(popmap_file,
                header = FALSE,
                col.names = c("sample_name", "population"))

# find paths for sample files
sample_files <- list.files(path = sample_dir, full.names = TRUE)
file_names <- data.table(filepath = sample_files,
                         sample_name = sapply(sample_files, function(x)
                             strsplit(basename(x), ".", fixed = TRUE)[[1]][1]))
all_samples <- merge(popmap, file_names, all.x = TRUE, all.y = FALSE)
all_samples[, bn := basename(filepath)]

# read the stats results
stats <- fread(stats_file)
stats[, bn := basename(filename)]

# merge the read counts
samples_with_readcount <- merge(all_samples,
                                stats[, .(bn, reads = n_scaffolds)],
                                by = "bn",
                                all.x = TRUE,
                                all.y = FALSE)

# filter samples by read count
filter <- 1857817 # 95% of samples with > 1857817 reads have mean coverage > 10
filtered_samples <- samples_with_readcount[!is.na(reads)][reads > filter]

# do a log4 transform for plotting
hist_data <- samples_with_readcount[!is.na(reads),
                                    .(l4reads = log(reads, 4)),
                                    by = sample_name]
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
lhist_log4 <- hist_data[, .(count = length(sample_name)), by = lbin]

# plot the histogram
filter_level <- as.numeric(cut(log(filter, 4), breaks = qa, labels = bin_labels))
filter_plot <- mean(bin_labels[c(filter_level, filter_level + 1)])

gt <- paste0(
    "'",
    samples_with_readcount[, length(unique(sample_name))],
    " individuals. '*",
    "'Filter: ", filter, " reads. ",
    filtered_samples[, length(unique(sample_name))],
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
fwrite(filtered_samples[, .(sample_name, population)],
       filtered_population_map,
       col.names = FALSE,
       sep = "\t")

fwrite(
    filtered_samples[, .(indivduals = length(unique(sample_name))),
                     by = population],
    pop_counts)

ggsave(filename = plot_file,
       plot = gp,
       device = "pdf",
       width = 10,
       height = 7.5,
       units = "in")

