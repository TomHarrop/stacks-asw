library(data.table)
library(ggplot2)

adaptors_file <- "test/adaptors.txt"
result_file <- "test/result.txt"

MungBbmapResults <- function(my_raw_input) {
    raw_input <- copy(my_raw_input)
    raw_input[, file := gsub("\\..+", "", file)]
    raw_input[, c("individual", "flowcell", "lane", "row", "col", "file") :=
                  tstrsplit(file, "_")]
    raw_input[, reads := as.numeric(gsub("^([[:digit:]]+).*", "\\1", reads))]
    raw_input[, bases := as.numeric(gsub("^([[:digit:]]+).*", "\\1", bases))]
    return(raw_input)
}

raw_adaptors <- fread(adaptors_file)
raw_adaptors[, V3 := NULL]
setnames(raw_adaptors, c("V1", "V2", "V4"), c("file", "reads", "bases"))
adaptors <- MungBbmapResults(raw_adaptors)

raw_results <- fread(result_file, header = FALSE, col.names = c("file", "reads", "bases"))
results <- MungBbmapResults(raw_results)

input_output <- merge(adaptors,
      results,
      by = c("individual", "flowcell", "lane", "row", "col"),
      suffixes = c("_input", "_result"))

input_output[, pct_reads := 100 * reads_result / reads_input]
input_output[, pct_bases := 100 * bases_result / bases_input]
input_output[, population := gsub("[^[:alpha:]]+", "", individual)]
input_output[, fc_lane := paste(flowcell, lane, sep = "_")]

gp <- ggplot(input_output[population != "GBSNEG"],
       aes(x = population, y = pct_bases, colour = population)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.border = element_rect(fill = NA, colour = "black")) +
    facet_wrap(~fc_lane) +
    ylab("Bases kept (%)") + xlab(NULL) +
    ylim(c(0,100)) +
    scale_color_brewer(palette = "Set3", guide = FALSE) +
    geom_point(position = position_jitter())
gp
ggsave("test/adaptor_filtering.pdf", gp, width = 10, height = 7.5)

