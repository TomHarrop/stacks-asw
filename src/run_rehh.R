#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(rehh)
library(ggplot2)


#############
# FUNCTIONS #
#############

RunIhs <- function(contig_vcf, min_pct){
  my_hh <- data2haplohh(hap_file = contig_vcf,
                        polarize_vcf = FALSE,
                        min_perc_geno.mrk = min_pct)
  my_scan <- scan_hh(my_hh,
                     polarized = FALSE,
                     discard_integration_at_border = FALSE)
  #return(ihh2ihs(my_scan, freqbin = 1))
  return(data.table(my_scan))
}

###########
# GLOBALS #
###########

pop_names <- names(snakemake@input)[names(snakemake@input) != "fai"]
print(pop_names)

fai_file <- snakemake@input[["fai"]]
pop1_files <- snakemake@input[[pop_names[[1]]]]
pop2_files <- snakemake@input[[pop_names[[2]]]]

# fai_file <- "output/005_ref/ref.fasta.fai"
# north_files <- list.files("output/100_ehh/ns.all",
#                           pattern = "North.*.phased.vcf.gz$",
#                           full.names = TRUE)
# south_files <- list.files("output/100_ehh/ns.all",
#                           pattern = "South.*.phased.vcf.gz$",
#                           full.names = TRUE)


########
# MAIN #
########

# run scans
pop1_scans <- rbindlist(lapply(pop1_files, RunIhs, min_pct = 0.9))
pop2_scans <- rbindlist(lapply(pop2_files, RunIhs, min_pct = 0.9))

# compare
xpehh <- data.table(ies2xpehh(scan_pop1 = data.frame(pop1_scans),
                              scan_pop2 = data.frame(pop2_scans),
                              popname1 = pop_names[[1]],
                              popname2 = pop_names[[2]],
                              p.adjust.method = "BH",
                              include_freq = TRUE))

# fix Chr coordinates
fai_names <- c("Chr", "chr_length")
chr_coords <- fread(fai_file, select = 1:2, col.names = fai_names)
setorder(chr_coords, -chr_length, Chr)
chr_coords[, chr_end := cumsum(chr_length)]
chr_coords[, chr_start := chr_end - chr_length + 1]

# set up manhattan plot
xpehh_coords <- merge(xpehh,
                      chr_coords,
                      by.x = "CHR",
                      by.y = "Chr",
                      all.x = TRUE,
                      all.y = FALSE)
setorder(xpehh_coords, -chr_length, CHR, POSITION)
xpehh_coords[, bp_coord := POSITION + chr_start - 1]
xpehh_coords[, CHR := factor(CHR, levels = gtools::mixedsort(unique(CHR)))]

# xlab_dt <- xpehh_coords[, .(bp_coord = mean(bp_coord)), by = CHR]

# interactive checks
# xpehh_coords[LOGPVALUE > 2] - check it out?
# cr <- calc_candidate_regions(xpehh, threshold = 4)
# manhattanplot(xpehh, pval=TRUE, cr = cr)

# quick plot, fixme
gp <- ggplot(xpehh_coords, aes(x = POSITION, y = LOGPVALUE)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(fill = NA)) +
  xlim(c(0, NA)) +
  xlab(NULL) + ylab(expression(-log[10](italic(p)))) +
  facet_grid(. ~ CHR, scales = "free_x", switch = "x") +
  geom_point()

wo <- grid::convertUnit(unit(483, "pt"), "mm", valueOnly = TRUE)
ho <- grid::convertUnit(unit(483/2, "pt"), "mm", valueOnly = TRUE)
ggsave(snakemake@output[["pdf"]],
       gp,
       width = wo,
       height = ho,
       unit = "mm",
       device = cairo_pdf)


fwrite(xpehh_coords, snakemake@output[["xpehh"]])

sessionInfo()

quit(save = "no")

###########
# NOT RUN #
###########


# haplotype plots?
# find the highest hit
sig_contig <- xpehh_coords[which.max(LOGPVALUE), CHR]
sig_pos <- xpehh_coords[which.max(LOGPVALUE), POSITION]

sig_n <- grep(sig_contig, north_files, value = TRUE)
sig_s <- grep(sig_contig, south_files, value = TRUE)

sig_n_hh <- data2haplohh(hap_file = sig_n,
                         polarize_vcf = FALSE,
                         min_perc_geno.mrk = 0.9)
sig_s_hh <- data2haplohh(hap_file = sig_s,
                         polarize_vcf = FALSE,
                         min_perc_geno.mrk = 0.9)


# get SNP position
marker_id_north <- which(sig_n_hh@positions == sig_pos)
marker_id_south <- which(sig_n_hh@positions == sig_pos)

north_furcation <- calc_furcation(sig_n_hh, mrk = marker_id_north)
south_furcation <- calc_furcation(sig_s_hh, mrk = marker_id_south)

plot(north_furcation)
plot(south_furcation)

north_haplen <- calc_haplen(north_furcation)
south_haplen <- calc_haplen(south_furcation)

plot(north_haplen)
plot(south_haplen)

haplen <- rbindlist(lapply(list(North = north_haplen,
                      South = south_haplen),
                 function(x) x$haplen),
          idcol = "pop")

haplen[, len := MAX - MIN]
haplen[, data.table(t(summary(len))), by = .(pop, DESCRIPTION)]

hl_s <- data.table(south_haplen$haplen)
hl_s[, hap := seq(1, .N)]




ggplot(hl_s, aes(x = MIN, xend = MAX, y = hap, yend = hap, colour = factor(ALLELE))) +
  scale_colour_viridis_d()   +
  geom_segment()

