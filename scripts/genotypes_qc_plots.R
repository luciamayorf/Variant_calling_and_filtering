# This R script generates a plots of the different statistics obtained with vcftools and plink.
# They serve as a diagnostic tool to evaluate the quality of the genotypes and the filtering steps.
# The input of the script is the prefix of the vcf analyzed with the script genotypes_qc_vcftools.sh

# Usage: Rscript genotypes_qc_plots.R <input_vcf_prefix> </path/to/input/tables/>

# input prefix 
input <- commandArgs(trailingOnly = TRUE)[1]

# input path
input_folder <- commandArgs(trailingOnly = TRUE)[2]

# load libraries
library(tidyverse)
library(ggplot2)
# library(plotly) --> these packages are not available at CESGA, I cannot create the interactive plots
# library(htmlwidgets)


#### Variant-based statistics ####

# plot variants quality
var_qual <- read_delim(paste0(input_folder, input, ".lqual"), delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
qual_plot <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
ggsave(paste0(input_folder, input, "_lqual_plot.pdf"), qual_plot)

# plot variants depth
var_depth <- read_delim(paste0(input_folder, input, ".ldepth.mean"), delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
depth_plot <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + xlim(0, 80)
summary(var_depth$mean_depth)
ggsave(paste0(input_folder, input, "_var_depth_plot.pdf"), depth_plot)

# plot variant missingness (proportion of missingness at each variant: how many individuals lack a genotype at a call site)
var_miss <- read_delim(paste0(input_folder, input, ".lmiss"), delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
miss_plot <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
summary(var_miss$fmiss)
ggsave(paste0(input_folder, input, "_var_miss_plot.pdf"), miss_plot)

# plot minor allele frequency
var_freq <- read_delim(paste0(input_folder, input, ".frq"), delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  
  # first we need to find minor allele frequencies
  var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  maf_plot <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

summary(var_freq$maf)
ggsave(paste0(input_folder, input, "_maf_plot.pdf"), maf_plot)


#### Individual based-statistics ####


# plot mean depth per individual
ind_depth <- read_delim(paste0(input_folder, input, ".idepth"), delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
ind_depth_plot <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, binwidth = 0.5) + theme_light()
ggsave(paste0(input_folder, input, "_ind_depth_plot.pdf"), ind_depth_plot)

# plot missing data per individual
ind_miss  <- read_delim(paste0(input_folder, input, ".imiss"), delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
ind_miss_plot <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
ggsave(paste0(input_folder, input, "_ind_miss_plot.pdf"), ind_miss_plot)

# plot heterozygosity and inbreeding coefficient per individual
ind_het <- read_delim(paste0(input_folder, input, ".het"), delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
  
  # add the population
  pop <- rep(NA, length(ind_het$ind))
  pop[grep("ca", ind_het$ind)] <- "captivity"
  pop[grep("sm", ind_het$ind)] <- "sierra morena"
  pop[grep("do", ind_het$ind)] <- "doñana"
  ind_het <- as_tibble(data.frame(ind_het, pop))

ind_het_plot <- ggplot(ind_het, aes(f, fill = pop, col = pop, text = ind)) +
                geom_histogram(alpha = 0.5, binwidth = 0.015) + 
                scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue")) + 
                scale_color_manual(values = c("darkorange", "darkgreen", "darkblue")) +
                theme_light()

ggsave(paste0(input_folder, input, "_ind_het_plot.pdf"), ind_het_plot)

# ind_het_plotly <- ggplotly(ind_het_plot, tooltip = c("ind"))
# htmlwidgets::saveWidget(ind_het_plotly, file = paste0(input_folder, input, "_ind_het_plot.html"))

# plot singletons per individual
ind_singletons <- read_delim(paste0(input_folder, input, ".singletons_per_ind"), delim = "\t", col_names = c("nsites", "indiv"), skip = 1)
ind_singletons_plot <- ggplot(ind_singletons, aes(nsites)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
ggsave(paste0(input_folder, input, "_ind_singletons_plot.pdf"), ind_singletons_plot)

# plot PCA
pca <- read_delim(paste0(input_folder, input, ".eigenvec"), delim = " ", col_names = F)
eigenval <- scan(paste0(input_folder, input, ".eigenval"))
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

  # add the population
  pop <- rep(NA, length(pca$ind))
  pop[grep("ca", pca$ind)] <- "captivity"
  pop[grep("sm", pca$ind)] <- "sierra morena"
  pop[grep("do", pca$ind)] <- "doñana"
  pca <- as_tibble(data.frame(pca, pop))
  
  # convert eigenvalues to percentage variance explained and plot
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  eigenval_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
                    ylab("Percentage variance explained") + theme_light()
  
  # plot PCA
  eigenvec_plot <- ggplot(pca, aes(x = PC1, y = PC2, col = pop, text = ind)) + 
                  geom_point(size = 2) +
                  scale_colour_manual(values = c("darkorange", "darkgreen", "darkblue")) + 
                  coord_equal() + theme_light() +
                  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
                  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

  # save plots
  ggsave(paste0(input_folder, input, "_eigenval_plot.pdf"), eigenval_plot)
  ggsave(paste0(input_folder, input, "_pca_plot.pdf"), eigenvec_plot)
  
  # indicate individuals name when clicking on each point (to be run in Rstudio)
  # pca_plotly <- ggplotly(eigenvec_plot, tooltip = c("ind"))
  # htmlwidgets::saveWidget(pca_plotly, file = paste0(input_folder, input, "_pca_plot.html"))
