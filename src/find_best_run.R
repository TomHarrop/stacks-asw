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




# Read in output files
runs_file <- purrr::map_dfr(bestlhoods, read_delim, delim = "\t") 

runs_summary <- runs_file %>% 
				select(MaxEstLhood, MaxObsLhood) %>%
  				mutate(Lhood = MaxObsLhood - MaxEstLhood) %>%
  				rownames_to_column("RunNumber") # Doing this assumes snakemake passes each run in order

message(runs_summary)

#Pull out best run
best_run <- runs_summary %>% 
  			filter(Lhood == max(Lhood)) %>% 
  			select(RunNumber) %>% 
  			pull()


# Trim filepath to get model dir
model_dir <- str_remove(bestlhoods[1], "run[:digit:]+/[:graph:]+")


# Create the best_run dir. Not sure if this is necessary, doing it just in case
ifelse(!dir.exists(file.path(model_dir, "best_run/")), 
  dir.create(file.path(model_dir, "best_run/")), FALSE)



message(paste("attempting to write file to ", paste0(model_dir, "best_run/best_is_", best_run, ".txt") ))
# Write out runs summary file, named after the best run
write_delim(runs_summary, paste0(model_dir, "best_run/best_is_", best_run, ".txt"), delim = " ")

## AIC
# Count the number of parameters in the model
k <- ncol(runs_file) -2

# Max est lhood
max_est <- max(runs_summary$MaxEstLhood)

# AIC
aic <- 2*k-2*(max_est/log10(exp(1)))

# Delta likelihood
delta_lhood <- runs_summary %>% 
  			   filter(Lhood == max(Lhood)) %>% 
  			   select(Lhood) %>% 
  			   pull()
write_lines(c(paste("AIC","DeltaL"), paste(aic, delta_lhood)), 
	paste0(model_dir, "best_run/AIC_", best_run, ".txt"))


# Write to log
sessionInfo()

