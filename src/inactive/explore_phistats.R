library(data.table)
library(ggplot2)

fai_file <- "output/map_to_genome/draft_genome.fasta.fai"
phistats_file <- "output/popgen/mapped/stacks_populations/populations.phistats.tsv"

phistats <- fread(phistats_file, skip = 11)

fai_names <- c("Chr", "chr_length")
fai <- fread(fai_file, select = 1:2, col.names = fai_names)

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

# homemade manhattan plot
phistats_with_len <- merge(phistats, chr_coords, by = "Chr")
setorder(phistats_with_len, -chr_length, Chr, BP)
phistats_with_len[, bp_coord := BP + chr_start - 1]

# pick out the outliers
q99 <- phistats_with_len[, quantile(`phi_st`, 0.99)]
phistats_with_len[`phi_st` > q99, outlier := TRUE]
phistats_with_len[outlier == TRUE, point_colour := Chr]
phistats_with_len[is.na(outlier), point_colour := NA]

# order the contigs
phistats_with_len[
  , point_colour := factor(point_colour,
           levels = unique(gtools::mixedsort(point_colour, na.last = TRUE)))]

ggplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   vjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "mm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(breaks = chr_coords[, lab_pos],
                     labels = chr_coords[, x_lab]) +
  scale_colour_viridis_d() +
  geom_hline(yintercept = q99) +
  geom_point(mapping = aes(x = bp_coord,
                           y = `phi_st`),
             data = phistats_with_len[is.na(point_colour)]) +
  geom_point(mapping = aes(x = bp_coord,
                           y = `phi_st`,
                           colour = point_colour),
             data = phistats_with_len[!is.na(point_colour)])


d99 <- phistats_with_len[, quantile(`Smoothed D_est`, 0.99)]
ggplot(phistats_with_len, aes(x = bp_coord, y = `Smoothed D_est`)) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   vjust = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "mm")) +
  scale_x_continuous(breaks = chr_coords[, lab_pos],
                     labels = chr_coords[, x_lab]) +
  geom_hline(yintercept = d99) +
  geom_point()


# distance between markers
phistats_with_len[, prev_loc := c(0, bp_coord[-.N])]
phistats_with_len[, distance_from_prev := bp_coord - prev_loc]

phistats_with_len[, summary(distance_from_prev)]


plot_contig <- fai[3, Chr]
x_lim <- fai[Chr == plot_contig, c(0, chr_length)]

plot_dt <- phistats[Chr == plot_contig]

ggplot(plot_dt,
       aes(x = BP, y = phi_st)) +
  xlim(x_lim) +
  geom_path()

phi_with_len <- merge(phistats, fai, by = "Chr")
x <- phi_with_len[, .(mean(phi_st), .N), by = .(Chr, chr_length)]
x[, distance_bw_markers := chr_length / N]


ggplot(x, aes(y = distance_bw_markers, x = Chr)) +
  geom_point()


ggplot(x, aes(y = N, x = V2)) + geom_point()

ggplot(x[N>5], aes(y = distance_bw_markers, x = V2)) + geom_point()
x[N>5]
