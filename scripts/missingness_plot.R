# This R script generates a pdf plot that represents the % of SNPs that would be included when the proportion of missing data allowed is increased. 
# It also outputs a summary table with the number of SNPs that are missing in each number of individuals.
# The input file is a table with the number of SNPs with each frequency of missing data generated with bcftools query. 

# Usage: Rscript missingness_plot.R <input_data> </path/to/output/folder/> <n_samples>

# load libraries
library(tidyverse)
library(RColorBrewer)
options(scipen=500)

# input data 
input <- commandArgs(trailingOnly = TRUE)[1]

# output folder
output_folder <- commandArgs(trailingOnly = TRUE)[2]

# number of samples
n_samples <- commandArgs(trailingOnly = TRUE)[3]

# load table
miss_table <- read_delim(input, progress = T, col_names = FALSE, delim = " ")

# get the total number of SNPs
total_snps <- nrow(miss_table)

# count the proportion of SNPs with each frequency of missing data
miss_table <- miss_table %>% group_by(X2) %>% summarise(n_snp = n())
colnames(miss_table)[1] <- "n_missing"

# get the proportion of SNPs with each frequency of missing data
miss_table$Freq <- (miss_table$n_snp)/total_snps

# calculate the proportion of SNPs included when X proportion of missing data is allowed (cummulative frequency)
miss_table$Freq_cum <- cumsum(miss_table$Freq)

# calculate proportion of missing data 
miss_table$mis_prop <- (miss_table$n_missing)/as.numeric(n_samples)

missing_plot <- ggplot(miss_table, aes(x = mis_prop, y = Freq_cum)) +
  geom_line(alpha=0.5, color="blue", size =1.5) +
  geom_point(alpha=0.5, color="blue", size=3, shape=15)+
  xlab("Proportion of Missing data allowed") +
  ylab("Proportion of SNPs included") +
  ggtitle("Proportion of SNPs included by missing data allowed") +
  scale_x_continuous(limits = c(0, 0.5))+
  theme_minimal()

# save plot
ggsave(paste0(output_folder, "cumulative_miss_plot.pdf"), 
       plot = missing_plot, width = 8, height = 4)

# save summary table
write.table(x = miss_table, paste0(output_folder, "missing_genotypes_table.txt"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

           
