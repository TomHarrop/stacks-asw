#!/usr/bin/env Rscript

library(data.table)

###########
# GLOBALS #
###########

input_csvs <- snakemake@input
output_csv <- snakemake@output[["combined"]]

########
# MAIN #
########

input_tables <- lapply(input_csvs, fread)
output_table <- rbindlist(input_tables)
fwrite(output_table, output_csv)
