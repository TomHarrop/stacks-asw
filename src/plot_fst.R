#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]])
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

fst_file <- snakemake@input[["fst"]]
fst_plot_file <- snakemake@output[["plot"]]

#fst_file <- "output/popgen/stacks_populations/populations.fst_summary.tsv"

# read data
fst_mat_dt <- fread(fst_file)

# make long
pairwise_fst <- melt(fst_mat_dt,
                     id.vars = "V1",
                     measure.vars = fst_mat_dt[, unique(V1)],
                     variable.name = "pop2",
                     value.name = "Fst")

# fill missing values
setnames(pairwise_fst, "V1", "pop1")
filled_fst <- rbind(pairwise_fst[!is.na(Fst)],
                    pairwise_fst[!is.na(Fst),
                                 .(pop1 = pop2, pop2 = pop1, Fst = Fst)])
fst_pd <- rbind(filled_fst,
      pairwise_fst[pop1 == pop2, .(pop1, pop2, Fst = 0)])

# get population order
fst_mat <- as.matrix(data.frame(dcast(fst_pd, pop1 ~ pop2), row.names = "pop1"))
hc <- hclust(as.dist(fst_mat), method = "ward.D2")
pop_order <- hc$labels[hc$order]

fst_pd[, pop1 := factor(pop1, levels = pop_order)]
fst_pd[, pop2 := factor(pop2, levels = pop_order)]

# draw the plot
gp <- ggplot(fst_pd, aes(x = pop1, y = pop2, fill = Fst)) +
    theme_minimal(base_family = "Lato") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_viridis_c(
        guide = guide_colorbar(
            title = expression(italic(F)["ST"])),
        na.value = NA) +
    coord_fixed () +
    geom_raster()

# write output
ggsave(fst_plot_file,
       gp,
       width = 10,
       height = 7.5,
       units = "in",
       device = cairo_pdf)

# log
sessionInfo()