library(data.table)
library(ggplot2)

#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########


# dev
individual_stats_file <- "output/run_stats/individual_stats_combined.csv"
read_stats_file <- "output/run_stats/read_stats.txt"
covstats_file <- "output/run_stats/individual_covstats_combined.csv"

########
# MAIN #
########

# read data
individual_stats <- fread(individual_stats_file)
covstats <- fread(covstats_file)

# merge
setnames(individual_stats,
         c("assembled_loci", "polymorphic_loci", "snps"),
         c("Assembled loci (K)", "Polymorphic loci (K)", "SNPs (K)"))
merged_stats <- merge(
    covstats,
    individual_stats,
    by = "individual")

# melt
merged_stats[, population := gsub("[[:digit:]]+", "", individual)]
pd <- melt(merged_stats,
           id.vars = c("population",
                       "individual",
                       "primary",
                       "secondary",
                       "total_reads"))

# save output
saveRDS(pd, "report/cov_vs_reads.Rds")

# plot assembly stats
ggplot(pd, aes(x = primary, y = value/1e3, colour = population)) +
    theme(strip.placement = "outside",
          strip.background = element_blank()) +
    guides(colour = FALSE) +
        facet_grid(variable ~ population, scales = "free_y", switch = "y") +
    ylab(NULL) + xlab("Mean primary coverage") +
    geom_point(position = position_jitter(width = 0.3))

# plot coverage
ggplot(merged_stats, aes(x = total_reads/1e6,
                         y = primary,
                         colour = population,
                         fill = population)) +
    stat_ellipse(geom = "polygon",
                 alpha = 0.3,
                 level = 0.95,
                 colour = NA) +
    xlim(c(0, 7.5)) + ylim(c(0, 30)) +
    geom_point()


