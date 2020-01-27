############################################ 
## R script                               ##
## Project: IL17.TNF_Biologics_manuscript ##
## ITS prep                               ##
############################################

### Brief description:
### This script creates all of the necessary phyloseq objects that are used to generate manuscript figures/tables from ITS data.
### Please see the corresponding code for specific figures/tables.

### Imported files:
### 1) biom file - contains forward sequence data; generated with QIIME v.1.9.1; paired-end trimming of sequences 
### performed with Cutadapt using the following adapters: forward - TAGAGGAAGTAAAAGTCGTAA...TTACGACTTTTACTTCCTCTA; 
### reverse - TTYRCTRCGTTCTTCATC...GATGAAGAACGYAGYRAA
### 2) mapping file - contains metadata

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(Hmisc)
library(PMCMR)
library(PMCMRplus)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)

############################################################################
############################################################################
############################################################################

### Create phyloseq object ###

# import biom file
b = ".../IL17.TNF/ITS/outputs/8_biom_R/otu_table_R1.json"
biom = import_biom(b, taxaPrefix = F)

# import mapping file
# note: must leave A1 cell empty for R compatibility
m = ".../IL17.TNF/ITS/inputs/map/Map_IL17.TNF_ITS_all_v2_R.txt"
map = sample_data(read.table(m, header = TRUE, sep = "\t", row.names = 1))

# create phyloseq object
phy_ITS.R1 = phyloseq(otu_table(biom), tax_table(biom), map)

# provide column names to separate different taxonomic levels
colnames(tax_table(phy_ITS.R1)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# print row and column names
rownames(sample_data(phy_ITS.R1))
colnames(sample_data(phy_ITS.R1))

############################################################################
############################################################################
############################################################################

### Phyloseq object manipulation ###

# remove taxa with zero OTUs
phy_ITS.R1_z <- subset_taxa(phy_ITS.R1, rowSums(otu_table(phy_ITS.R1)) > 0)

# transform from counts to relative abundance
rel_abundance = function(x) { x/sum(x) }
phy_ITS.R1_zr <- transform_sample_counts(phy_ITS.R1_z, rel_abundance)

############################################################################
############################################################################
############################################################################

### Subset phyloseq objects ###

# filter samples
# note: this step is performed because the original biom file contains samples that are not relevant to the final analysis
phy_ITS.R1_human <- subset_samples(phy_ITS.R1_z, Host == "human")
phy_ITS.R1_human_clean <- subset_samples(phy_ITS.R1_human, Filter_all_analysis == "keep")

phy_ITS.R1_human_TNF <- subset_samples(phy_ITS.R1_human_clean, Treatment == "1_TNF")
phy_ITS.R1_human_IL17 <- subset_samples(phy_ITS.R1_human_clean, Treatment == "2_IL17")

phy_ITS.R1_human_TNF.B <- subset_samples(phy_ITS.R1_human_TNF, Timepoint_revised == "B")
phy_ITS.R1_human_TNF.C <- subset_samples(phy_ITS.R1_human_TNF, Timepoint_revised == "C")

phy_ITS.R1_human_IL17.B <- subset_samples(phy_ITS.R1_human_IL17, Timepoint_revised == "B")
phy_ITS.R1_human_IL17.C <- subset_samples(phy_ITS.R1_human_IL17, Timepoint_revised == "C")
phy_ITS.R1_human_IL17.D <- subset_samples(phy_ITS.R1_human_IL17, Timepoint_revised == "D")

# IL17 visits B and D do not have the same number of samples (some timepoints are missing from visit D)
# create phyloseq objects from timepoints B and D that have identical subjects (i.e. no missing timepoints)
IL17.D_subjects <- as.vector(sample_data(phy_ITS.R1_human_IL17.D)$Subject)
phy_ITS.R1_human_IL17.B.per.D <- subset_samples(phy_ITS.R1_human_IL17.B, (sample_data(phy_ITS.R1_human_IL17.B)$Subject) %in% IL17.D_subjects)

# check to make sure subjects are identical
rownames(sample_data(phy_ITS.R1_human_IL17.B.per.D))
rownames(sample_data(phy_ITS.R1_human_IL17.D))

# merge
phy_ITS.R1_human_TNF.B.C <- merge_phyloseq(phy_ITS.R1_human_TNF.B, phy_ITS.R1_human_TNF.C)
phy_ITS.R1_human_IL17.B.C <- merge_phyloseq(phy_ITS.R1_human_IL17.B, phy_ITS.R1_human_IL17.C)
phy_ITS.R1_human_IL17.B.D <- merge_phyloseq(phy_ITS.R1_human_IL17.B.per.D, phy_ITS.R1_human_IL17.D)
phy_ITS.R1_human_IL17.B.C.D <- merge_phyloseq(phy_ITS.R1_human_IL17.B, phy_ITS.R1_human_IL17.C, phy_ITS.R1_human_IL17.D)

phy_ITS.R1_human_TNF.B.C_IL17.B.C <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C, phy_ITS.R1_human_IL17.B.C)
phy_ITS.R1_human_TNF.B.C_IL17.B.D <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C, phy_ITS.R1_human_IL17.B.D)
phy_ITS.R1_human_TNF.B.C_IL17.B.C.D <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C, phy_ITS.R1_human_IL17.B.C.D)

########################################

# samples/subjects w/ low reads removed # 
phy_ITS.R1_human_clean_fR1 <- subset_samples(phy_ITS.R1_human_clean, Filter_alpha_beta != "filter_R1.R2")
phy_ITS.R1_human_TNF_fR1 <- subset_samples(phy_ITS.R1_human_clean_fR1, Treatment == "1_TNF")
phy_ITS.R1_human_IL17_fR1 <- subset_samples(phy_ITS.R1_human_clean_fR1, Treatment == "2_IL17")

phy_ITS.R1_human_TNF.B_fR1 <- subset_samples(phy_ITS.R1_human_TNF_fR1, Timepoint_revised == "B")
phy_ITS.R1_human_TNF.C_fR1 <- subset_samples(phy_ITS.R1_human_TNF_fR1, Timepoint_revised == "C")

phy_ITS.R1_human_IL17.B_fR1 <- subset_samples(phy_ITS.R1_human_IL17_fR1, Timepoint_revised == "B")
phy_ITS.R1_human_IL17.C_fR1 <- subset_samples(phy_ITS.R1_human_IL17_fR1, Timepoint_revised == "C")
phy_ITS.R1_human_IL17.D_fR1 <- subset_samples(phy_ITS.R1_human_IL17_fR1, Timepoint_revised == "D")

# IL17 visits B and D do not have the same number of samples (some timepoints are missing from visit D)
# create matching phyloseq objects from timepoints B and D (i.e. no missing timepoints)
IL17.D_subjects_fR1 <- as.vector(sample_data(phy_ITS.R1_human_IL17.D_fR1)$Subject)
phy_ITS.R1_human_IL17.B.per.D_fR1 <- subset_samples(phy_ITS.R1_human_IL17.B_fR1, (sample_data(phy_ITS.R1_human_IL17.B_fR1)$Subject) %in% IL17.D_subjects_fR1)

# check to make sure subjects are identical
rownames(sample_data(phy_ITS.R1_human_IL17.B.per.D_fR1))
rownames(sample_data(phy_ITS.R1_human_IL17.D_fR1))

# merge
phy_ITS.R1_human_TNF.B.C_fR1 <- merge_phyloseq(phy_ITS.R1_human_TNF.B_fR1, phy_ITS.R1_human_TNF.C_fR1)
phy_ITS.R1_human_IL17.B.C_fR1 <- merge_phyloseq(phy_ITS.R1_human_IL17.B_fR1, phy_ITS.R1_human_IL17.C_fR1)
phy_ITS.R1_human_IL17.B.D_fR1 <- merge_phyloseq(phy_ITS.R1_human_IL17.B.per.D_fR1, phy_ITS.R1_human_IL17.D_fR1)
phy_ITS.R1_human_IL17.B.C.D_fR1 <- merge_phyloseq(phy_ITS.R1_human_IL17.B_fR1, phy_ITS.R1_human_IL17.C_fR1, phy_ITS.R1_human_IL17.D_fR1)

phy_ITS.R1_human_TNF.B.C_IL17.B.C_fR1 <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C_fR1, phy_ITS.R1_human_IL17.B.C_fR1)
phy_ITS.R1_human_TNF.B.C_IL17.B.D_fR1 <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C_fR1, phy_ITS.R1_human_IL17.B.D_fR1)
phy_ITS.R1_human_TNF.B.C_IL17.B.C.D_fR1 <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C_fR1, phy_ITS.R1_human_IL17.B.C.D_fR1)

############################################################################
############################################################################
############################################################################

### Rarefy to even depth ###

# rarefaction performed on filtered phyloseq objects only (i.e. w/ suffix "fR1")
# depth 1000 to match 16S data
phy_ITS.R1_human_clean_fR1_even1000 <- rarefy_even_depth(phy_ITS.R1_human_clean_fR1, sample.size = 1000, rngseed = 711, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

############################################################################
############################################################################
############################################################################

### Subset rarefied phyloseq objects ###

phy_ITS.R1_human_TNF_fR1_even1000 <- subset_samples(phy_ITS.R1_human_clean_fR1_even1000, Treatment == "1_TNF")
phy_ITS.R1_human_IL17_fR1_even1000 <- subset_samples(phy_ITS.R1_human_clean_fR1_even1000, Treatment == "2_IL17")
phy_ITS.R1_human_TNF.B_fR1_even1000 <- subset_samples(phy_ITS.R1_human_TNF_fR1_even1000, Timepoint_revised == "B")
phy_ITS.R1_human_TNF.C_fR1_even1000 <- subset_samples(phy_ITS.R1_human_TNF_fR1_even1000, Timepoint_revised == "C")
phy_ITS.R1_human_IL17.B_fR1_even1000 <- subset_samples(phy_ITS.R1_human_IL17_fR1_even1000, Timepoint_revised == "B")
phy_ITS.R1_human_IL17.C_fR1_even1000 <- subset_samples(phy_ITS.R1_human_IL17_fR1_even1000, Timepoint_revised == "C")
phy_ITS.R1_human_IL17.D_fR1_even1000 <- subset_samples(phy_ITS.R1_human_IL17_fR1_even1000, Timepoint_revised == "D")
IL17.D_subjects_fR1_even1000 <- as.vector(sample_data(phy_ITS.R1_human_IL17.D_fR1_even1000)$Subject)
phy_ITS.R1_human_IL17.B.per.D_fR1_even1000 <- subset_samples(phy_ITS.R1_human_IL17.B_fR1_even1000, (sample_data(phy_ITS.R1_human_IL17.B_fR1_even1000)$Subject) %in% IL17.D_subjects_fR1_even1000)
phy_ITS.R1_human_TNF.B.C_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_TNF.B_fR1_even1000, phy_ITS.R1_human_TNF.C_fR1_even1000)
phy_ITS.R1_human_IL17.B.C_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_IL17.B_fR1_even1000, phy_ITS.R1_human_IL17.C_fR1_even1000)
phy_ITS.R1_human_IL17.B.D_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_IL17.B.per.D_fR1_even1000, phy_ITS.R1_human_IL17.D_fR1_even1000)
phy_ITS.R1_human_IL17.B.C.D_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_IL17.B_fR1_even1000, phy_ITS.R1_human_IL17.C_fR1_even1000, phy_ITS.R1_human_IL17.D_fR1_even1000)
phy_ITS.R1_human_TNF.B.C_IL17.B.C_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C_fR1_even1000, phy_ITS.R1_human_IL17.B.C_fR1_even1000)
phy_ITS.R1_human_TNF.B.C_IL17.B.D_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C_fR1_even1000, phy_ITS.R1_human_IL17.B.D_fR1_even1000)
phy_ITS.R1_human_TNF.B.C_IL17.B.C.D_fR1_even1000 <- merge_phyloseq(phy_ITS.R1_human_TNF.B.C_fR1_even1000, phy_ITS.R1_human_IL17.B.C.D_fR1_even1000)

############################################################################
############################################################################
############################################################################

### Statistics functions ###

# all plot statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75, per90 = per90))
}

# boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

# whiskers statistics: median, min value, max value
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

############################################################################
############################################################################
############################################################################

### Reads ###

# directory for storing files
dir = ".../IL17.TNF/ITS/outputs/10_reads/"

# background theme
bkg <- theme_few() +
  theme(axis.text.x = element_text(size = 16, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 11, face = "bold"))

# create dataset for plotting
d <- data.frame(sample_data(phy_ITS.R1_human_TNF.B.C_IL17.B.C.D))   

# create filenames
filename_plot = "ITS.R1_human_TNF.B.C_IL17.B.C.D_post.otu_reads_plot.pdf"  
filename_plot.stats = "ITS.R1_human_TNF.B.C_IL17.B.C.D_post.otu_reads_plot.stats.csv"
filename_stats = "ITS.R1_human_TNF.B.C_IL17.B.C.D_post.otu_reads_stats.txt"

# plot and save
p <- ggplot(data = d, aes(x = Treatment, y = Reads_R1_post_otu_picking, fill = Treatment)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
               color = "black", size = 0.5, width = 0.5, fill = "white") +
  geom_jitter(width = 0.1, size = 1.5) +
  scale_x_discrete(labels = c("TNFi", "IL-17i")) +
  expand_limits(y = 0) +
  xlab(NULL) + ylab("Reads post-OTU picking") +
  guides(fill = FALSE, color = FALSE) + # no legend
  bkg

fp = paste(dir, filename_plot, sep = "")
pdf(file = fp)
plot(p)
dev.off()

# obtain plot statistics and save
s <- aggregate(Reads_R1_post_otu_picking ~ Treatment, data = d, stats.all)
s <- t(s)
rownames(s) <- c("treatment", "mean", "sd", "median", "min", "max", 
                 "10%ile", "25%ile", "75%ile", "90%ile")

fps = paste(dir, filename_plot.stats, sep = "")
write.csv(s, file = fps)

# calculate statistics
mw <- wilcox.test(Reads_R1_post_otu_picking ~ Treatment, data = d, paired = FALSE)

# save calculations
fs = paste(dir, filename_stats, sep = "")
cat("Reads: ITS-R1 TNF.BC-IL17.BCD cohort\n\n", file = fs)
cat(reads[i], file = fs, append = TRUE)
cat("\n\n", file = fs, append = TRUE)
cat("Mann-Whitney U test:\n", file = fs, append = TRUE)
capture.output(mw, file = fs, append = TRUE)


