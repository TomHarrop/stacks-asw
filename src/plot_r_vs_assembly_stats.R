#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########


# dev
pop_stats_file <- "output/run_stats/population_stats_combined.csv"

########
# MAIN #
########

pop_stats <- fread(pop_stats_file)
pd <- melt(pop_stats, id.vars = "r")

ggplot(pd, aes(x = r, y = value)) +
    facet_grid(variable ~ ., scales = "free_y") +
    geom_point()

