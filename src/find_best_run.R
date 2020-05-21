#!/usr/bin/env Rscript

library(tidyverse)

###########
# GLOBALS #
###########
bestlhoods <-snakemake@input[["bestlhoods"]]
log_file <- snakemake@log[[1]]

########
# MAIN #
########
# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")



# Read in output files and summarise all runs
runs_summary <- purrr::map_dfr(bestlhoods, read_delim, delim = "\t") %>% 
				  select(MaxEstLhood, MaxObsLhood) %>%
  				  mutate(Lhood = MaxObsLhood - MaxEstLhood) %>%
  				  rownames_to_column("RunNumber") # Doing this assumes snakemake passes each run in order

#Pull out best run
best_run <- runs_summary %>% 
  			filter(Lhood == max(Lhood)) %>% 
  			select(RunNumber) %>% 
  			pull()

# Trim filepath to get model dir
model_dir <- str_remove(bestlhoods[1], "run[:digit:]+/[:graph:]+")

# Write out runs summary file, named after the best run
write_delim(runs_summary, paste0(model_dir, "best_run/best_is_", best_run, ".txt"), delim = " ")

# Write log
sessionInfo()

