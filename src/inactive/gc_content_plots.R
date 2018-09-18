#!/usr/bin/env Rscript

library(tidyverse)
library(ggridges)

# FUNCTIONS

ReadBBdukHist <- function(x) {
    read_tsv(x, comment = "#", col_names = c("GC", "count"))
}

ReadMeanGC <-function(x) {
    my_mean <- read_tsv(x, col_names = c("GC","count")) %>% 
        filter(GC == "#Mean") %>% 
        get("count",.) %>% 
        as.numeric()
    return(tibble(mean_gc=my_mean))
}

# get list of files
gc_file_list <- list.files("output/filtering/gc_hist", pattern = ".txt", full.names = TRUE)
names(gc_file_list) <- sub(".txt", "", basename(gc_file_list))

# read data
histogram_list <- lapply(gc_file_list, ReadBBdukHist)
gc_mean_list <- lapply(gc_file_list, ReadMeanGC)

# get GC means
gc_means <- bind_rows(gc_mean_list, .id = "individual") %>% 
    mutate(population = gsub("[[:digit:]]+", "", individual)) %>% 
    arrange(mean_gc)
individual_order <- gc_means$individual

# calculate a GC cutoff
cutoff <- 50
kept_indivs <- gc_means %>% filter(mean_gc < cutoff) %>% 
    select(individual) %>% unlist()

# read histograms
histogram_data <- bind_rows(histogram_list, .id = "individual") %>% 
    mutate(
        population = gsub("[[:digit:]]+", "", individual),
        individual = factor(individual, levels = individual_order),
        kept = individual %in% kept_indivs)

# plot GC content
ggplot(gc_means, aes(x = population,
                     y = mean_gc,
                     colour = population,
                     fill = population)) +
    geom_hline(yintercept = 50, linetype = 2, colour = alpha("black", 0.5)) +
    geom_boxplot(width = 0.5, outlier.size = - 1, fill = NA) +
    geom_point(position = position_jitter(width = 0.2))

# plot
ggplot(histogram_data, aes(x = GC, y = individual, height = count, fill = kept)) +
    facet_wrap(~ population, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    geom_density_ridges(stat = "identity", scale = 10, alpha = 0.5)
