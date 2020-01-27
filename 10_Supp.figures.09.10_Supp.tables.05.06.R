############################################# 
## R script                                ##
## Project: IL17.TNF_Biologics_manuscript  ##
## Supplementary figures 09, 10            ##
## Supplementary tables 05, 06             ##
## Metagenomic data                        ##
#############################################

### Brief description:
### Supplementary figure 09: TNFi violin plots
### Supplementary figure 10: IL-17i violin plots
### Supplementary table 05: TNFi statistics
### Supplementary table 06: IL-17i statistics

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(ggplot2)
library(reshape2)
library(scales)
library(ggthemes)

############################################################################
############################################################################
############################################################################

### Metagenomic pathways ###

### Analysis was performed with QIIME v1.9.1 using the group_significance.py command.
### Comparisons were made pre vs post treatment in the TNFi and IL-17i cohorts.
### Only significant results were included in tables and plotted.
### QIIME analysis provided the Kruskal-Wallis results.
### Wilcoxon statistics were calculated using R commands (see below)

### Input files: 
### Column 1: sampleID
### Column 2: timepoint
### Remainder of columns: reletive abundance across samples of significant pathways for each sample

## Setup ##

# colors
col1 <- c("#c40018", "#ffb400")
col2 <- c("#8559a5", "#448ef6", "#53d397")

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title = element_text(size = 20, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(0, 0, 5, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  theme(strip.background = element_rect(fill = "black")) +
  theme(strip.text = element_text(color = "white", face = "bold", size = 10)) +
  theme(legend.position = "none") +
  theme(panel.spacing = unit(1.5, "lines"))

# directory
dir = ".../IL17.TNF/Metagenomics/jobs/4_pathways_violoin.plots_R/"

############################

### TNFi pathways ###

## Statistics
# includes all significant results from from QIIME group_significance.py script
# (i.e. outliers also included)

# read in data table
d.pathways <- read.table(file = paste(dir, "sig.pathways_TNF.B.C_all.txt", sep = ""),
                         header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                         na.strings = "NA")

# convert to data frame
d.pathways <- as.data.frame(d.pathways)

# create table for storing wilcoxon results
wilc.table <- matrix(data = NA, nrow = 15, ncol = 2)
colnames(wilc.table) <- c("pathways", "pre-maintenance")

# calculate wilcoxon statistic across all pathways
for (i in 1:15) {
  wilc.table[i,1] <- colnames(d.pathways)[i+1]
  wilc.table[i,2] <- wilcox.test(d.pathways[,i+1] ~ Pre_post_det, data = d.pathways, paired = TRUE)$p.value
}

# save
ft = paste(dir, "pathways_TNF.B.C_wilcoxon.csv", sep = "")
write.csv(file = ft, wilc.table)

##########

## Plotting
# excludes outliers

# read in data table
d.pathways <- read.table(file = paste(dir, "sig.pathways_TNF.B.C_figure.txt", sep = ""),
                         header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                         na.strings = "NA")

# convert to data frame
d.pathways <- as.data.frame(d.pathways)

# melt data for graphing
d <- melt(d.pathways, id.vars = "Pre_post_det")

# plot and save
p <- ggplot(d, aes(x = Pre_post_det, y = value, fill = Pre_post_det)) +
  facet_wrap(~ variable, scales="free_x", ncol = 2,
             labeller = label_wrap_gen(width = 40)) +
  geom_violin() + 
  geom_jitter(width = 0.1, size = 1) +
  scale_fill_manual(values = col1) +
  coord_flip() +
  scale_y_continuous(labels = scientific)  +
  scale_x_discrete(labels = c("TNFi maint", "TNFi pre")) +
  xlab("") + ylab("\nPathway relative abundance") +
  bkg

fp = paste(dir, "pathways_TNF.B.C_plot.pdf", sep = "")
pdf(file = fp, width = 8, height = 10)
plot(p)
dev.off()

############################

### IL-17i pathways ###

## Statistics
# includes all significant results from from QIIME group_significance.py script
# (i.e. outliers also included)

# read in data table
d.pathways <- read.table(file = paste(dir, "sig.pathways_IL17.B.C.D_all.txt", sep = ""),
                         header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                         na.strings = "NA")

# convert to data frame
d.pathways <- as.data.frame(d.pathways)

### pre-loading

# subset pre/loading samples
d.load <- d.pathways[d.pathways$Pre_post_det != "IL17.maint", ]

# create table for storing wilcoxon results
wilc.table.load <- matrix(data = NA, nrow = 8, ncol = 2)
colnames(wilc.table.load) <- c("pathway", "pre-loading")

# calculate wilcoxon statistic across all pathways
for (i in 1:8) {
  wilc.table.load[i,1] <- colnames(d.load)[i+1]
  wilc.table.load[i,2] <- wilcox.test(d.load[,i+1] ~ Pre_post_det, data = d.load, paired = TRUE)$p.value
}

# save
ft.load = paste(dir, "pathways_IL17.B.C.D_wilcoxon.loading.csv", sep = "")
write.csv(file = ft.load, wilc.table.load)

### pre-maintenance

# subset pre/maintenance samples
d.maint <- d.pathways[d.pathways$Pre_post_det != "IL17.load", ]
d.maint <- d.maint[-c(9, 12, 15, 20), ] # remove non-matching "B" visits

# create table for storing wilcoxon results
wilc.table.maint <- matrix(data = NA, nrow = 8, ncol = 2)
colnames(wilc.table.maint) <- c("pathway", "pre-maintenance")

# calculate wilcoxon statistic across all pathways
for (i in 1:8) {
  wilc.table.maint[i,1] <- colnames(d.maint)[i+1]
  wilc.table.maint[i,2] <- wilcox.test(d.maint[,i+1] ~ Pre_post_det, data = d.maint, paired = TRUE)$p.value
}

# save
ft.maint = paste(dir, "pathways_IL17.B.C.D_wilcoxon.maintenance.csv", sep = "")
write.csv(file = ft.maint, wilc.table.maint)

##########

## Plotting
# excludes outliers

# read in data table
d.pathways <- read.table(file = paste(dir, "sig.pathways_IL17.B.C.D_figure.txt", sep = ""),
                         header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                         na.strings = "NA")

# convert to data frame
d.pathways <- as.data.frame(d.pathways)

# melt data for graphing
d <- melt(d.pathways, id.vars = "Pre_post_det")

# plot and save
p <- ggplot(d, aes(x = Pre_post_det, y = value, fill = Pre_post_det)) +
  facet_wrap(~ variable, scales="free_x", ncol = 2,
             labeller = label_wrap_gen(width = 40)) +
  geom_violin() + 
  geom_jitter(width = 0.1, size = 1) +
  scale_fill_manual(values = col2) +
  coord_flip() +
  scale_x_discrete(labels = c("IL-17i maint", "IL-17i load", "IL-17i pre")) +
  scale_y_continuous(labels = scientific)  +
  xlab("") + ylab("\nPathway relative abundance") +
  bkg

fp = paste(dir, "pathways_IL17.B.C.D_plot.pdf", sep = "")
pdf(file = fp, height = 6, width = 8)
plot(p)
dev.off()

