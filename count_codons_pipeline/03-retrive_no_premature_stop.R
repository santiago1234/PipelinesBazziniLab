library(tidyverse)

###
# script to select transcits with no
# premature stop and length grater than
# user suplied argument
# usage: Rscript script.R <input_table> <output_table> <minimum_aa_length>
###


# usr arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
    message("too few argumens, at least 1 needed")
    stop("usage: Rscript script.R <input_table> <output_table> <minimum_aa_length>")
} else if (length(args) == 1) {
   args[2]<- "out_no_stop.tab"
    args[3] <- 100
} else if (length(args) == 2) {
    args[3] <- 100
}


# initialize variables
input_data <- args[1]
out_data <- args[2]
minimum_aa_length <- args[3]


message("getting transcripts with not premature stop")

input_data %>% read_tsv() %>% filter(stop_position == -1) %>%
    filter(transcrit_leng > minimum_aa_length) %>%
    write.table(file = out_data, sep = "\t", quote = FALSE,
      row.names = FALSE)

message("script completed >>>>")
