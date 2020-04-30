#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

preview_file <- snakemake@input[["preview"]]
  
prv_wide <- fread(cmd = paste("grep",
                              "'('",
                              preview_file),
                  fill = TRUE,
                  sep = '\t')

prv_wide[, id := paste0("pop", 1:.N)]
prv <- melt(prv_wide,
            id.vars = "id")

prv[, id := factor(id, levels = gtools::mixedsort(unique(id)))]
prv[, c("n", "seg") := tstrsplit(gsub("[()]+", "", value), split = ", ")]
prv[, n := as.numeric(n)]
prv[, seg := as.numeric(seg)]

setorder(prv, id)

proj <- prv[prv[, .I[which.max(seg)], by = id][, V1], paste(n, collapse = ",")]

cat(proj, file = snakemake@output[["proj"]])

sessionInfo()
