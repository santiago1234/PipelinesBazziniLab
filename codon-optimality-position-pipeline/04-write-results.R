#### description
## this scripts reads the file: tmp_counts (tmp_counts is output by: 03-count-codons-by.py)
## makes a visualization plot for quality check
## and writtes the tmp_count to a tsv table
####

library(tidyverse)
library(parallel)
library(magrittr)

## usr arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  message("too few argumens, 1 needed")
  stop("usage: Rscript script.R  <out_file>")
  # by default the scripts reads as input a file called tmp.tab that was generated in the
  # previous step
}

input_table <- "tmp_counts"
out_file <- args[1]

codon_counts <- read.table("tmp_counts", stringsAsFactors = F, header = F,
                           col.names = c("position", "codon", "frequency")) %>% tbl_df

number_of_seqs <- filter(codon_counts, position == 1)$frequency %>% sum

codon_counts %<>% mutate(density = frequency / number_of_seqs, library = out_file)

## visualization

message("making quality check plot")
ggplot(codon_counts, aes(x = position, y = density)) +
  geom_smooth(method = "lm") +
  geom_line(alpha = 0.8) +
  facet_wrap(~ codon, nrow = 8) +
  geom_vline(xintercept = 50, color = "red", linetype = 2) +
  ylim(0, 0.09) + theme_bw() -> g

## save results
ggsave(g, filename = paste0(out_file, ".pdf"), width = 11, height = 10)
write_tsv(x = codon_counts, path = paste0(out_file, "-counts"))
message("task comppleated!!!")