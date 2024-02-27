# This R script takes the mosdepth output bed files (mean depth per 10kbp window) of all individuals.
# It calculates the mean of the mean depth values of all individuals per row (as all bed files have the same number of rows - same number of windows).
# It then filters out windows with mean depth values greater than the mean + 0.5 * standard deviation.
# It also plots the distribution of mean depth values and saves the plot as a pdf.
# The filtered windows are saved as a bed file.

# Usage: Rscript mean_rd_filter_population_bed.R </path/to/input/data/> </path/to/sample_list.txt> </path/to/output/folder/>

# load libraries
library(tidyverse)
library(RColorBrewer)
options(scipen=500)

# input data folder
input_folder <- commandArgs(trailingOnly = TRUE)[1]

# output folder
output_folder <- commandArgs(trailingOnly = TRUE)[3]

# input file
input_file <- "_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_mosdepth.regions.bed.gz"

# list of individuals
samples_file <- commandArgs(trailingOnly = TRUE)[2]
samples <- readLines(samples_file)

# read bed file of the first individual into data frame
mean_rd_bed_pop <- read.table(paste0(input_folder, samples[1], input_file), header = FALSE)

# loop through the rest of the individuals
for (indiv in samples[-1]){
    # read bed file into data frame
    indiv_bed <- read.table(paste0(input_folder, indiv, input_file), header = FALSE)
    
    # add 4th column to global data frame
    mean_rd_bed_pop <- cbind(mean_rd_bed_pop, indiv_bed$V4)
}

# rename the columns
colnames(mean_rd_bed_pop) <- c("chromosome", "start", "end", samples)

# calculate mean of fourth column
mean_rd_bed_pop$mean_rd <- rowMeans(mean_rd_bed_pop[4:ncol(mean_rd_bed_pop)])

# add window number column
mean_rd_bed_pop$win_num <- seq_len(nrow(mean_rd_bed_pop))

# keep only first five columns
mean_rd_bed_pop <- mean_rd_bed_pop[c("chromosome", "start", "end", "mean_rd", "win_num")]

# remove the rows that contain "scaffold" in the chromosome column
mean_rd_bed_pop <- mean_rd_bed_pop[!grepl("scaffold", mean_rd_bed_pop$chromosome),]

# define cutoffs (median, median * 1.5, median * 1.8, 99th percentile, mean and mean+0.5*sd)
median <- median(mean_rd_bed_pop$mean_rd)
median1_5 <- median * 1.5
median1_8 <- median * 1.8
percent99 <- quantile(mean_rd_bed_pop$mean_rd, 0.99)
mean <- mean(mean_rd_bed_pop$mean_rd)
print(paste("Mean read depth is", mean))
mean_sd <- mean + 0.5*sd(mean_rd_bed_pop$mean_rd)

# plot the distribution of the mean depth
rd_plot_pop <- ggplot() +
    geom_histogram(data = mean_rd_bed_pop, aes(x = mean_rd), bins = 200) +
    scale_x_continuous(breaks = seq(0, 400, 10), limits = c(0, mean(mean_rd_bed_pop$mean_rd) * 2.5)) +
    geom_vline(aes(xintercept = median, color = "median"), linetype = "dashed", size = 0.3) +
    geom_vline(aes(xintercept = median1_5, color = "median1_5"), linetype = "dashed", size = 0.3) +
    geom_vline(aes(xintercept = median1_8, color = "median1_8"), linetype = "dashed", size = 0.3) +
    geom_vline(aes(xintercept = percent99, color = "percent99"), linetype = "dashed", size = 0.3) +
    geom_vline(aes(xintercept = mean, color = "mean"), linetype = "dashed", size = 0.3) +
    geom_vline(aes(xintercept = mean_sd, color = "mean_sd"), linetype = "dashed", size = 0.3) +
    scale_color_manual(values = c("median" = "black", "median1_5" = "blue", "median1_8" = "lightblue2", "percent99" = "green", "mean" = "red", "mean_sd" = "pink"),
                       labels = c("Median", "1.5*Median", "1.8*Median", "99th Percentile", "Mean", "Mean+0.5*sd"),
                       breaks = c("median", "median1_5", "median1_8", "percent99", "mean", "mean_sd"),
                       guide = guide_legend(override.aes = list(linetype = "solid", shape = 15))) +
    theme_classic()
          
# save plot
ggsave(filename = paste0(output_folder, "mean_rd_cutoffs_pop.pdf"),
       plot = rd_plot_pop,
       width = 4,
       height = 4,
       units = "in")

# calculate maximum depth based on median
maxdepth <- mean_sd

# see how many windows are excluded based on max depth (for summary table)
remaining_windows <- mean_rd_bed_pop %>% filter(mean_rd < maxdepth)
n_excluded_windows <- nrow(mean_rd_bed_pop) - nrow(remaining_windows)
print(paste("Number of windows filtered:", n_excluded_windows))

# generate bed of excluded windows
excluded_windows_bed <- mean_rd_bed_pop %>% filter(mean_rd >= maxdepth)

# save bed of windows to be filtered
write.table(x = excluded_windows_bed,
           file = paste0(output_folder, "excluded_windows_mean_rd_pop_filter.bed"),
           quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
           