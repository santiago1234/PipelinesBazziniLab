library(tidyverse)
library(reshape2)
library(gridExtra)

###
# R script to make some visualization plots of the
# number of premature stop codons
# usage: Rscript <input_table> <output.pdf>
###

# usr arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
    # default output file
    args[2] = "plots.pdf"
}

data <- args[1] %>% read_tsv() %>% select(stop_position, transcrit_leng)

# Visualization

# Premature Stop Codom
(data$stop_position == -1) %>%
  table() %>% data.frame() %>%
  mutate(percent = Freq / (.$Freq[1] + .$Freq[2]),
         stop_codon = c("premature", "no_premature")) -> dat1
  

ggplot(dat1, aes(x = stop_codon, y = percent, alpha = stop_codon)) +
  geom_bar(stat="identity", width = 0.5, fill = "black") + theme_classic() +
  xlab("Stop Codon") + ylab("%") +
  scale_alpha_manual(values = c(0.7 , 0.4)) +
  ggtitle(paste(dat1$Freq[1] + dat1$Freq[2], "seqs")) +
  theme(legend.position = "none") -> plot1


# Transcript length
data %>% mutate(status = (data$stop_position != -1) %>%
                  map_chr(function(x) if(x) "premature stop" else "no premature stop")) %>%
  ggplot(aes(x = transcrit_leng, alpha = status)) +
  facet_grid(. ~ status) +
  geom_histogram(binwidth = 10) +
  theme_classic() +
  xlab("length aa") +
  scale_alpha_manual(values = c(0.9 , 0.6)) +
  theme(legend.position = "none") -> plot2

# Premature Stop Codom position distribution
data %>% filter(stop_position != -1) %>%
  ggplot(aes(x = stop_position)) +
  geom_histogram(binwidth = 2, alpha = 0.6) +
  ggtitle("First Premature Stop Codon Position") +
  theme_classic() -> plot3

# save plots
pdf(file = args[2], width = 11, height = 3)
grid.arrange(plot1, plot2, plot3, nrow = 1, widths = c(2, 6, 3))
dev.off()
