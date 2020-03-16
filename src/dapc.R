#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(adegenet)
library(vcfR)
library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########


vcf_file <- snakemake@input[[1]]

dapc_plot_file <- snakemake@output[["dapc_plot"]]
pca_plot_file <- snakemake@output[["pca_plot"]]
dapc_xv_file <- snakemake@output[["dapc_xv"]]

# DEV
# vcf_file <- "output/popgen/mapped/stacks_populations/populations.snps.vcf"


# roughly north to south?
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

# construct genind
vcf <- read.vcfR(vcf_file)
snp_data <- vcfR2genlight(vcf)

# impute
na_means <- tab(snp_data, NA.method = "mean")
snps_imputed <- new("genlight",
                    na_means,
                    ploidy = 2)

pop(snps_imputed) <- factor(gsub("[^[:alpha:]]+", "", snps_imputed$ind.names),
                            levels = pop_order)

# run the pca
pca <- glPca(snps_imputed, nf = Inf)
pct_var <- 100 * (pca$eig^2)/(sum(pca$eig^2))

# set up a plot title
gt <- paste0(length(unique(snps_imputed$ind.names)),
             " individuals genotyped at ", 
             format(snps_imputed$n.loc, big.mark = ","),
             " loci")

# generate pca plot data
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

# plot the pca
pca_plot <- ggplot(pca_pd[component %in% paste0("PC", 1:9)],
                   aes(x = population, y = score, colour = population)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.border = element_rect(fill = NA, colour = "black")) +
    xlab(NULL) + ylab(NULL) +
    ggtitle(paste("PCA of", gt)) +
    facet_wrap(~facet_label) +
    scale_color_brewer(palette = "Set3", guide = FALSE) +
    geom_boxplot(fill = NA,
                 colour = alpha("black", 0.5),
                 outlier.colour = NA, 
                 width = 0.5) +
    geom_point(shape = 16,
               position = position_jitter(width = 0.2))

# run the dapc
dapc_results_opt <- dapc(snps_imputed,
                         n.pca = 10,
                         n.da = length(pop_order) - 1)

opt <- optim.a.score(dapc_results_opt,
                     n.pca = 1:50,
                     n.sim = 10)

dapc_results <- dapc(snps_imputed,
                     n.pca = opt$best,
                     n.da = length(pop_order) - 1)

dapc_var <- 100 * (dapc_results$eig^2)/(sum(dapc_results$eig^2))

# cross validate
xv <- xvalDapc(tab(snps_imputed),
               pop(snps_imputed),
               n.pca.max = 15,
               training.set = 0.9,
               result = "groupMean",
               center = TRUE,
               scale = FALSE, 
               n.pca = 1:15,
               n.rep  = 3,
               xval.plot = FALSE)

# generate dapc plot data
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

# plot the dapc
dapc_plot <- ggplot(dapc_pd[dapc_var > 1],
                    aes(x = population, y = score, colour = population)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.border = element_rect(fill = NA, colour = "black")) +
    xlab(NULL) + ylab(NULL) +
    ggtitle(paste("DAPC of", gt)) +
    facet_wrap(~facet_label) +
    scale_color_brewer(palette = "Set3", guide = FALSE) +
    geom_boxplot(fill = NA,
                 colour = alpha("black", 0.5),
                 outlier.colour = NA,
                 width = 0.5) +
    geom_point(shape = 16, position = position_jitter(width = 0.2))

# write output
ggsave(dapc_plot_file,
       dapc_plot,
       width = 10,
       height = 7.5,
       units = "in",
       device = cairo_pdf)
ggsave(pca_plot_file,
       pca_plot,
       width = 10,
       height = 7.5,
       units = "in",
       device = cairo_pdf)
saveRDS(xv, dapc_xv_file)

# write log
sessionInfo()
