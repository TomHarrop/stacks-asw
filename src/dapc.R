#!/usr/bin/env Rscript

library(adegenet)
library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

plink_file <- "output/popgen/plink.raw"

pop_order <- c(
    "Coromandel",
    "Ruakura",
    "Taranaki",
    "Wellington",
    "Reefton", 
    "Greymouth",
    "Lincoln",
    "O",
    "MararoaDowns",
    "Mossburn",
    "Fortrose")


########
# MAIN #
########

snp_data <- read.PLINK(plink_file)

# impute NAs on mean
na_means <- tab(snp_data, NA.method = "mean")
snps_imputed <- new("genlight", na_means)
ploidy(snps_imputed) <- 2

# try the pca
pca <- glPca(snps_imputed, nf = 9)
pct_var <- 100 * (pca$eig^2)/(sum(pca$eig^2))

# plot pca
pca_dt <- data.table(pca$scores, keep.rownames = TRUE)
setnames(pca_dt, "rn", "individual")
pca_long <- melt(pca_dt,
                 id.vars = "individual",
                 variable.name = "component",
                 value.name = "score")

pca_long[, population := gsub("[^[:alpha:]]+", "", individual)]

pca_pd <- merge(pca_long,
                data.table(component = paste0("PC", 1:length(pct_var)),
                           pct_var = pct_var),
                all.y = FALSE)
pca_pd[, population := factor(population, levels = pop_order)]
pca_pd[, facet_label := paste0(component, " (", round(pct_var, 1), "%)")]

ggplot(pca_pd, aes(x = population, y = score, colour = population)) +
    theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              panel.border = element_rect(fill = NA, colour = "black")) +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~facet_label) +
    scale_color_brewer(palette = "Set3", guide = FALSE) +
    geom_boxplot(fill = NA, colour = alpha("black", 0.5), outlier.colour = NA, width = 0.5) +
    geom_point(shape = 16, position = position_jitter(width = 0.2))


# dapc
pop(snps_imputed) <- factor(gsub("[^[:alpha:]]+", "", snps_imputed$ind.names),
                            levels = pop_order)
dapc_results <- dapc(snps_imputed, n.pca = 50, n.da = 10)
dapc_var <- 100 * (dapc_results$eig^2)/(sum(dapc_results$eig^2))

# cross validate
xvalDapc(tab(snps_imputed),
         pop(snps_imputed),
         n.pca.max = 50,
         training.set = 0.9,
         result = "groupMean",
         center = TRUE,
         scale = FALSE, 
         n.pca = 1:50,
         n.rep  = 50,
         xval.plot = TRUE)

# plot dapc
dapc_dt <- data.table(dapc_results$ind.coord, keep.rownames = TRUE)
setnames(dapc_dt, "rn", "individual")
dapc_long <- melt(dapc_dt,
                 id.vars = "individual",
                 variable.name = "component",
                 value.name = "score")

dapc_long[, population := gsub("[^[:alpha:]]+", "", individual)]


dapc_pd <- merge(dapc_long,
                data.table(component = paste0("LD", 1:length(dapc_var)),
                           dapc_var = dapc_var),
                all.y = FALSE)
dapc_pd[, population := factor(population, levels = pop_order)]
dapc_pd[, facet_label := paste0(component, " (", round(dapc_var, 1), "%)")]

ggplot(dapc_pd[dapc_var > 1], aes(x = population, y = score, colour = population)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.border = element_rect(fill = NA, colour = "black")) +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~facet_label) +
    scale_color_brewer(palette = "Set3", guide = FALSE) +
    geom_boxplot(fill = NA, colour = alpha("black", 0.5), outlier.colour = NA, width = 0.5) +
    geom_point(shape = 16, position = position_jitter(width = 0.2))
