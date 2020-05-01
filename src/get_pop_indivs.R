#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

popmap_file <- snakemake@input[["popmap"]]
outdir <- snakemake@output[["outdir"]]

popmap <- fread(popmap_file, header = FALSE)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

sapply(popmap[, unique(V2)], function(x)
  fwrite(popmap[V2 == x],
         file = paste0(outdir, "/", x, ".txt"),
         col.names = FALSE,
         sep = '\t'))

sessionInfo()
