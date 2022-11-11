# load libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# load files with summary statistics
stat_files = list.files(path = "bamstat_plots", pattern = "*.summary$")

# iterate
initial_sequences_val = c()
properly_paired_and_mapped_val = c()
#
for(i in stat_files){
  # load summary statistics
  d1 = readLines(i)
  # extract numeric values:
  initial_sequences = grep(pattern = "^sequences:", x = d1) %>%
    (function(x){
      str_extract(pattern = "\\t[0-9]*?$", string = d1[x])
    }) %>%
    gsub(pattern = "\t", replacement = "") %>%
    as.integer()
  #
  initial_sequences_val = c(
    initial_sequences_val,
    initial_sequences
  )
  #
  properly_paired_and_mapped = grep(pattern = "reads mapped and paired", x = d1)  %>%
    (function(x){
      str_extract(pattern = "\\t[0-9]*?\\t#", string = d1[x])
    }) %>%
    gsub(pattern = "\t#{0,1}", replacement = "") %>%
    as.integer()
  #
  properly_paired_and_mapped_val = c(
    properly_paired_and_mapped_val,
    properly_paired_and_mapped
  )
}

# plot with samples
tiff(filename = "properly_paired_&_mapped.tif",
     res = 300, units = "in", width = 8, height = 8)
ggplot(data = data.frame(
  vals = properly_paired_and_mapped_val/initial_sequences_val
), aes(x = vals))+
  geom_histogram(color="black", fill="coral")+
  theme_bw()+
  labs(
    x = "Percentage of properly paired and mapped reads (%)",
    y = "Frequency"
  )
dev.off()

# plot with random values
set.seed(1)
range(rnorm(n=50,mean=.85,sd=.05))
#
tiff(filename = "properly_paired_&_mapped_random.tif",
     res = 300, units = "in", width = 8, height = 8)
ggplot(data = data.frame(
  vals = rnorm(n=50,mean=.9,sd=.1)
), aes(x = vals))+
  geom_histogram(color="black", fill="coral", bins=20)+
  theme_bw()+
  labs(
    x = "Percentage of properly paired and mapped reads (%)\n(random)",
    y = "Frequency"
  )
dev.off()
