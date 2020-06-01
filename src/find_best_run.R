#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

bestlhoods <-snakemake@input[["bestlhoods"]]

# bestlhoods <- list.files("output/095_fastsimcoal",
#            pattern = "*.bestlhoods",
#            recursive = TRUE,
#            full.names = TRUE)

########
# MAIN #
########

names(bestlhoods) <- sapply(strsplit(bestlhoods, "/"), function(x)
  paste(x[[3]], x[[4]], sep = "."))

# parse the names
runs <- rbindlist(lapply(bestlhoods, fread),
                  idcol = "filename",
                  fill = TRUE)
runs[, c("popset", "pruned", "model", "mig", "run") :=
       tstrsplit(filename, ".", fixed = TRUE)]
runs[, filename := NULL]
runs[, run := gsub(".*run([[:digit:]]{2,3}).*", "\\1", run)]

# calculate likelihood
runs[, Lhood := MaxObsLhood - MaxEstLhood]

# count the number of params for each model
info_cols <- c("MaxEstLhood",
               "MaxObsLhood",
               "popset",
               "pruned",
               "model", 
               "mig",
               "run",
               "Lhood")
param_cols <- names(runs)[!names(runs) %in% info_cols]

## AIC
# Count the number of parameters in the model
runs[, k := sum(!is.na(.SD)),
     .SDcols = param_cols,
     by = 1:nrow(runs)]
runs[, aic := (2 * k) - (2 * (MaxEstLhood/log10(exp(1))))]
setorder(runs, -Lhood)

# Max est lhood???
best_run <- runs[is.finite(MaxEstLhood)][which.max(MaxEstLhood)]
max_est <- runs[, max(MaxEstLhood)]

fwrite(runs,
       snakemake@output[["summary"]])
fwrite(best_run,
       snakemake@output[["bestrun"]])

# log
sessionInfo()
