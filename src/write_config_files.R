#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

outdir <- snakemake@params[["outdir"]]
para_key_file <- snakemake@input[["para_key_file"]]
geo_key_file <- snakemake@input[["geo_key_file"]]

# dev
# para_key_file <- "combined_key_data.csv"
# geo_key_file <- "data/reads/SQ0727.txt"


########
# MAIN #
########

geo_key <- fread(geo_key_file)
para_key <- fread(para_key_file)

# geo

# fix spaces
geo_key[, sample := gsub("[[:space:]]", "-", sample)]
geo_key[grepl("\\.", sample), sample := gsub("\\.", "-", sample)]


# add sample details
geo_key[, fc_lane := paste(flowcell, lane, sep = "_"),
         by = .(flowcell, lane)]
geo_key[, sample_fullname := paste(sample,
                                    flowcell,
                                    lane,
                                    row,
                                    column,
                                    sep = "_"),
         by = .(sample,
                flowcell,
                lane,
                row,
                column)]

# para 

# add sample details
para_key[, fc_lane := paste(key, lane, sep = "_"),
        by = .(key, lane)]
para_key[, sample_fullname := paste(sample_name,
                                   key,
                                   lane,
                                   sep = "_"),
        by = .(sample_name,
               key,
               lane)]


# write config files
split_geo_key <- split(geo_key, by = "fc_lane")
sapply(names(split_geo_key), function(x)
    fwrite(split_geo_key[[x]][, .(barcode, sample_fullname)],
           paste0(outdir, "/", x, ".config"),
           sep = "\t",
           col.names = FALSE))

split_para_key <- split(para_key, by = "fc_lane")
sapply(names(split_para_key), function(x)
  fwrite(split_para_key[[x]][, .(barcode, sample_fullname)],
         paste0(outdir, "/", x, ".config"),
         sep = "\t",
         col.names = FALSE))


# make popmap
popmap <- rbindlist(list(geo = geo_key[, data.table(sample)],
               para = para_key[, .(sample = sample_name)]),
          idcol = "expt")
popmap[, population := paste(tolower(gsub("[^[:alpha:]]+", "", sample)),
                             expt,
                             sep = "_")]

# write popmap
fwrite(unique(popmap[!grepl("NEG", sample), .(sample, population)]),
       sep = "\t",
       col.names = FALSE)

# save log
sessionInfo()

