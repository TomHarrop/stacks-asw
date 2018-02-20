#!/usr/bin/env Rscript

library(adegenet)
library(data.table)

############
# FUNCTION #
############

###########
# GLOBALS #
###########

key_file <- 'data/SQ0003.txt'
populations_genepop <- "output/stacks_populations/r0.8/populations.snps.gen"
adegenet_output <- "temp/adegenet.Rds"

########
# MAIN #
########

# read the key file
key_data <- fread(key_file)
indiv_data <- key_data[, .(
    individual = Sample,
    Flowcell, Lane, Platename)]
indiv_data[, fc_lane := paste(Flowcell, Lane, sep = "_")]
indiv_data[, population := gsub("[[:digit:]]+", "", individual)]

# Plate vs. pop?
ggplot(indiv_data, aes(x = Flowcell, y = population)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.2))
ggplot(indiv_data, aes(x = Platename, y = population)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.2))
ggplot(indiv_data, aes(x = Platename, y = fc_lane)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.2))

# read SNPs
if(!file.exists(adegenet_output)){
    snps <- read.genepop(populations_genepop)
    pop(snps) <- gsub("[[:digit:]]+", "", indNames(snps))
    saveRDS(snps, adegenet_output)
} else {
    snps <- readRDS(adegenet_output)
}

# subset out weird popuations
keep_populations <- c("Coromandel",
                      "Taranaki",
                      "Wellington",
                      "Mossburn",
                      "Fortrose",
                      "Reefton",
                      "Greymouth",
                      "Lincoln",
                      "Ruakura")
filtered_snps <- snps[as.character(pop(snps)) %in% keep_populations,
                      drop = TRUE]

# calculate minor allele frequency
maf <- minorAllele(filtered_snps)
keep_snps <- unique(names(maf[maf > 0.05]))
filtered_snps <- filtered_snps[loc = keep_snps, drop = TRUE]

# impute missing SNPs by mean
snp_tab <- tab(filtered_snps, NA.method = "mean")

# run the PCA
snps_pca <- dudi.pca(snp_tab,
                     center = TRUE,
                     scale = FALSE,
                     scannf = FALSE,
                     nf = 9)
s.class(snps_pca$li, pop(filtered_snps), col=rainbow(nPop(filtered_snps)))
add.scatter.eig(snps_pca$eig[1:10], xax=1, yax=2)

# try a ggplot of the results
var_exp <- snps_pca$eig / sum(snps_pca$eig)
pv_labs <- paste0("PC", 1:length(var_exp), " (", round(var_exp*100, 1), "%)")
facet_labeller <- function(x) {
    list(pv_labs[x[, 1]])
}

pd <- data.table(snps_pca$li, keep.rownames = TRUE)
pd_long <- melt(pd, id.vars = "rn")
pd_long[, PC := as.numeric(gsub("[[:alpha:]]+", "", variable))]
pd_long[, population := gsub("[[:digit:]]+", "", rn)]
pd_data <- merge(pd_long,
                 indiv_data,
                 by.x = c("rn", "population"),
                 by.y = c("individual", "population"),
                 all.x = TRUE,
                 all.y = FALSE)

ggplot(pd_data, aes(x = population, y = value, colour = population)) +
    scale_colour_brewer(palette = "Paired") +
    facet_wrap(~ PC, ncol = 3, labeller = facet_labeller) +
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.6,
               shape = 16,
               size = 2)

ggplot(pd_data[PC < 11], aes(x = Platename, y = value, colour = population)) +
    facet_wrap(~ PC, ncol = 3, labeller = facet_labeller) +
    scale_colour_brewer(palette = "Paired") +
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.6,
               shape = 16,
               size = 2)

ggplot(pd_data[PC < 11], aes(x = fc_lane, y = value, colour = population)) +
    facet_wrap(~ PC) +
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.5,
               shape = 16,
               size = 2)

ggplot(pd_data[PC < 11], aes(x = population, y = value, colour = Platename)) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~ PC) +
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.5,
               shape = 16,
               size = 2)

ggplot(pd_data[PC < 11], aes(x = Flowcell, y = value, colour = population)) +
    facet_wrap(~ PC) +
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.5,
               shape = 16,
               size = 2)




pd <- as.data.table(snps_pca$li, keep.rownames = TRUE)
pd[, population := gsub("[[:digit:]]+", "", rn)]

ggplot(pd,
       aes(x = Axis1, y = Axis2, colour = population, fill = population)) +
    theme_minimal() +
    xlab(pv_labs[2]) + ylab(pv_labs[3]) +        
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    stat_ellipse(geom = "polygon",
                 alpha = 0.3,
                 colour = NA,
                 level = 0.95) +
    geom_point()
