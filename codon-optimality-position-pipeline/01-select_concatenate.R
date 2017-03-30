library(tidyverse)
library(magrittr)
library(parallel)

## usr arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  message("too few argumens, 4 needed")
  stop("usage: Rscript script.R <library_table> <position_size> <out_file>")
}

## Initialize variable 
lib_df <- args[1]
size <- as.integer(args[2])
cores <- 30
out_file <- args[3]


get_begining_end <- function(seq, size = 50) {
  # Extract the first and last codons size of a dna string
  # Args:
  #   seq, str dna sequence
  #   size, int the size of the susbtring to extract begining and end
  # Retunrs: a str biging + end
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  length_seq <- floor(nchar(seq) / 3)
  seq <- substr(seq, 1, length_seq * 3)
  begining <- substr(seq, 1, size * 3)
  end <-  substrRight(seq, size * 3)
  paste(begining, end, sep = "")
}


message("extracting begining and end positions from: ", lib_df)
message(cores, " cores will be used")

lib_df %<>% read_tsv()
lib_df %<>% filter(transcrit_leng > size*3*2)
lib_df %<>% mutate(transcrit_sequence = mclapply(.$transcrit_sequence, 
                                                 FUN = function(x) get_begining_end(x, size),
                                                 mc.cores = cores) %>% unlist())
# save output table
write_tsv(lib_df, path = out_file)
