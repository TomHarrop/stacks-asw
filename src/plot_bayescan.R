#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

# vcf_file <- "output/060_popgen/populations.ns.all.vcf"
# fai_file <- "output/005_ref/ref.fasta.fai"
# bs_file <- "output/080_bayescan/ns.all/bs/populations_fst.txt"

vcf_file <- snakemake@input[["vcf"]]
fai_file <- snakemake@input[["fai"]]
bs_file <- snakemake@input[["fst"]]

# get a list of loci to match with bs results
vcf_names <- c("chrom", "pos", "id", "ref", "alt", "qual")
loci <- fread(cmd = paste('grep -v "^#"', vcf_file),
              header = FALSE,
              col.names = vcf_names,
              select = 1:6)
loci[, snp_no := seq(1, .N)]

# read the FAI
fai_names <- c("Chr", "chr_length")
fai <- fread(fai_file, select = 1:2, col.names = fai_names)

# bayescan results
bs_names <- c("snp_no", "prob", "log10_PO", "qval", "alpha", "fst")
bs <- fread(bs_file,
            col.names = bs_names,
            skip = 1,
            header = FALSE)

# add locus info to bayescan results
bs_loci <- merge(bs, loci, by = "snp_no")
bs_loci[, nlog10_q := -log10(qval)]

# set up labelling for the contigs
n_to_label <- 20

# chromosome coordinates
chr_coords <- copy(fai)
setorder(chr_coords, -chr_length, Chr)
chr_coords[, chr_end := cumsum(chr_length)]
chr_coords[, chr_start := chr_end - chr_length + 1]
chr_coords[, lab_pos := chr_start + round(mean(c(chr_start, chr_end)), 0), by = Chr]
pos_to_label = chr_coords[, seq(1, max(chr_end), length.out = n_to_label)]

label_positions <- sapply(pos_to_label, function(x)
  chr_coords[, .I[which.min(abs(x - lab_pos))]])

chr_coords[label_positions, x_lab := Chr]
chr_coords[is.na(x_lab), x_lab := ""]

# set up manhattan plot
bs_coords <- merge(bs_loci,
                   chr_coords,
                   by.x = "chrom",
                   by.y = "Chr",
                   all.x = TRUE,
                   all.y = FALSE)
setorder(bs_coords, -chr_length, chrom, pos)
bs_coords[, bp_coord := pos + chr_start - 1]

# cycle through 6 colours to distinguish neighbouring contigs
n_cols <- 6
my_cols <- viridisLite::viridis(n_cols)
bs_coords[, chr_no := as.numeric(factor(chrom, levels = unique(chrom)))]
bs_coords[, pt_col := my_cols[[((chr_no - 1) %% n_cols) + 1]], by = chr_no]
col_scale <- unique(bs_coords, by = "chrom")[, structure(pt_col, names = chrom)]

# scale by shape to indicate q-values == 0
ylim_max <- bs_coords[is.finite(nlog10_q), max(nlog10_q)]
bs_coords[, plot_q := nlog10_q]
bs_coords[is.infinite(nlog10_q), plot_q := ylim_max]

# pick a qval cutoff
fdr_cutoff <- 1e-2

gp <- ggplot(bs_coords,
       aes(x = bp_coord,
           y = plot_q,
           colour = chrom,
           alpha = nlog10_q > -log10(fdr_cutoff),
           shape = is.infinite(nlog10_q))) +
  theme_minimal(base_size = 8 ) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   vjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab(NULL) + ylab(expression(-log[10](italic(Q)))) +
  scale_colour_manual(values = col_scale, guide = FALSE) +
  scale_x_continuous(breaks = chr_coords[, lab_pos],
                     labels = chr_coords[, x_lab]) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim_max + 0.5)) +
  scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1),
                     guide = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     guide = FALSE) +
  geom_point()

wo <- grid::convertUnit(unit(483, "pt"), "mm", valueOnly = TRUE)
ho <- grid::convertUnit(unit(483/2, "pt"), "mm", valueOnly = TRUE)
ggsave(snakemake@output[["pdf"]],
       gp,
       width = wo,
       height = ho,
       unit = "mm",
       device = cairo_pdf)

sessionInfo()

# all markers under fdr
# ggplot(bs_coords[qval < fdr_cutoff], aes(x = fst, y = qval, colour = alpha)) + geom_point()
# 
# # single contig
# plot_chr <- bs_coords[qval < fdr_cutoff, .N, by = chrom][which.max(N), chrom]
# x_max_chr <- bs_coords[chrom == plot_chr, unique(chr_length)]
# ggplot(bs_coords[chrom == plot_chr],
#        aes(x = pos, y = plot_q, colour = alpha)) +
#   xlim(c(1, x_max_chr)) +
#   geom_point()


