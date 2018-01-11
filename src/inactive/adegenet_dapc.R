#!/usr/bin/env Rscript

library(adegenet)
library(data.table)
library(ggplot2)
library(ochRe)

############
# FUNCTION #
############

###########
# GLOBALS #
###########

key_file <- 'data/SQ0003.txt'
populations_genepop <- "output/stacks_populations/r0.8/populations.snps.gen"
adegenet_output <- "temp/adegenet.Rds"

# read the key file
key_data <- fread(key_file)
indiv_data <- key_data[, .(
    individual = Sample,
    Flowcell, Lane, Platename)]
indiv_data[, fc_lane := paste(Flowcell, Lane, sep = "_")]
indiv_data[, population := gsub("[[:digit:]]+", "", individual)]

# prepare snps
snps <- readRDS(adegenet_output)
keep_populations <- c("Coromandel",
                      "Taranaki",
                      "Wellington",
                      "Mossburn",
                      "Fortrose",
                      "Reefton",
                      "Greymouth",
                      "Lincoln",
                      "Ruakura")
keep_populations <- levels(pop(snps))
keep_populations <- indiv_data[, length(unique(individual)),
                               by = population][V1 > 19, unique(population)]

filtered_snps <- snps[as.character(pop(snps)) %in% keep_populations,
                      drop = TRUE]
maf <- minorAllele(filtered_snps)
keep_snps <- unique(names(maf[maf > 0.05]))
filtered_snps <- filtered_snps[loc = keep_snps, drop = TRUE]

# run dapc
dapc_results <- dapc(filtered_snps,
                     pop(filtered_snps),
                     n.pca = 100,
                     n.da = 9)

# what if we explicitly group by flowcell?
indiv_to_fc <- indiv_data[, structure(fc_lane,
                                      .Names = individual)]
fc_grp <- factor(indiv_to_fc[rownames(filtered_snps$tab)])
dapc_results <- dapc(filtered_snps,
                     fc_grp,
                     n.pca = 100,
                     n.da = 9)

# extract components
dapc_dt <- data.table(dapc_results$ind.coord, keep.rownames = TRUE)
setnames(dapc_dt, "rn", "individual")
dapc_dt[, population := gsub("[[:digit:]]+", "", individual)]
dapc_long <- melt(dapc_dt,
                  id.vars = c("individual", "population"))

# add contributions
pct <- data.table(
    variable = colnames(dapc_results$ind.coord),
    pct = round(100 * dapc_results$eig / sum(dapc_results$eig), 1))
pct_var <- pct[, structure(paste0(variable, " (", pct, "%)"),
                           .Names = variable)]
dapc_long[, pct_var := plyr::revalue(variable, pct_var)]

# add individual data
pd <- merge(dapc_long,
            indiv_data,
            by = c("individual", "population"),
            all.x = TRUE,
            all.y = FALSE)

# save plot data for RMD report
saveRDS(pd, "report/dapc_pd.Rds")
saveRDS(indiv_data, "report/dapc_indiv.Rds")

# 2d plots
pd_2d <- dcast(pd,
               individual + population + fc_lane ~ pct_var,
               value.var = "value")
ggplot(pd_2d, aes(x = get(pd[, levels(pct_var)][1]),
                  y = get(pd[, levels(pct_var)][2]),
                  colour = population)) +
    xlab(pd[, levels(pct_var)][1]) + ylab(pd[, levels(pct_var)][2]) +
    theme_minimal(base_size = 12) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = "Population")) +
    geom_point(alpha = 0.8,
               shape = 16,
               size = 2)
ggplot(pd_2d, aes(x = get(pd[, levels(pct_var)][1]),
                  y = get(pd[, levels(pct_var)][2]),
                  colour = fc_lane)) +
    xlab(pd[, levels(pct_var)][1]) + ylab(pd[, levels(pct_var)][2]) +
    theme_minimal(base_size = 12) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = "Flowcell")) +
    geom_point(alpha = 0.8,
               shape = 16,
               size = 2)

# components vs. population
ggplot(pd, aes(x = population, y = value, colour = population)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside") +
    facet_wrap(~ pct_var, ncol = 3, strip.position = "left") +
    ggtitle("Discriminants vs. populations") +
    xlab(NULL) + ylab(NULL) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.25),
               alpha = 0.6,
               shape = 16,
               size = 2)

# components vs fc_lane
ggplot(pd, aes(x = fc_lane, y = value, colour = population)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside") +
    facet_wrap(~ pct_var, ncol = 3, strip.position = "left") +
    ggtitle("Discriminants vs. flowcell") +
    xlab(NULL) + ylab(NULL) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.25),
               alpha = 0.6,
               shape = 16,
               size = 2)

# components vs plate
ggplot(pd, aes(x = Platename, y = value, colour = population)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside") +
    facet_wrap(~ pct_var, ncol = 3, strip.position = "left") +
    xlab(NULL) + ylab(NULL) +
    scale_colour_brewer(palette = "Set1",
                        guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.1),
               alpha = 0.5,
               shape = 16,
               size = 2)

# components vs population coloured by lane
ggplot(pd, aes(x = population, y = value, colour = fc_lane)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside") +
    facet_wrap(~ pct_var, ncol = 3, strip.position = "left") +
    ggtitle("Discriminants vs. flowcell") +
    xlab(NULL) + ylab(NULL) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.25),
               alpha = 0.6,
               shape = 16,
               size = 2)

# distribution of individuals over fc_lanes
ggplot(indiv_data, aes(x = population, y = fc_lane, colour = fc_lane)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_brewer(palette = "Set1", guide = FALSE) +
    xlab(NULL) + ylab(NULL) +
    geom_point(position = position_jitter(height = 0.2, width = 0.2))

