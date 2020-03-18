library(data.table)
library(ggplot2)

fai_file <- "output/map_to_genome/draft_genome.fasta.fai"
phistats_file <- "output/popgen/mapped/stacks_populations/populations.phistats.tsv"

phistats <- fread(phistats_file, skip = 11)
setkey(phistats, Chr, BP)


fai <- fread(fai_file)
setorder(fai, -V2)

plot_contig <- fai[3, V1]
x_lim <- fai[V1 == plot_contig, c(0, V2)]

plot_dt <- phistats[plot_contig]

ggplot(plot_dt,
       aes(x = BP, y = phi_st)) +
  xlim(x_lim) +
  geom_path()

phi_with_len <- merge(phistats, fai, by.x = "Chr", by.y = "V1")
x <- phi_with_len[, .(mean(phi_st), .N), by = .(Chr, V2)]
x[, distance_bw_markers := V2 / N]


ggplot(x, aes(y = V1, x = V2)) +
  geom_point()


ggplot(x, aes(y = N, x = V2)) + geom_point()

ggplot(x[N>5], aes(y = distance_bw_markers, x = V2)) + geom_point()
x[N>5]
