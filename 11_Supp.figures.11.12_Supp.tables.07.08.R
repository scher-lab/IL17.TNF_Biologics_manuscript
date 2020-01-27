#############################################
## R script                                ##
## Project: IL17.TNF Biologics manuscript  ##
## Supplementary figures 11, 12            ##
## Supplementary tables 07, 08             ##
## 16S and ITS data                        ##
#############################################

### Brief description:

### Supplementary figures 11, 12:
### Panels A: Bacterial and fungal taxa correlations w/ local proteins, cytokines and metabolites for TNFi
### Panels B: Bacterial and fungal taxa correlations w/ local proteins, cytokines and metabolites for IL-17i

### Supplementary table 07: Statistics for proteins and cytokines
### Supplementary table 08: Statistics for metabolites

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(corrplot)
library(reshape2)

############################################################################
############################################################################
############################################################################

### Correlograms for proteins, cytokines and metabolites ###

### Only genera that had relative abundance of at least 0.5% in at least two samples were maintained.
### QIIME v1.9.1 was used to calculate Spearman correlations via the observation_metadata_correlation.py command.
### Correlation tables (input files) were created in Microsoft Excel.

### Input files:
### Column 1: significant genera
### Remainder of columns: Spearman's rank correlation coefficient's accross genera

# directory for storing files
dir = ".../IL17.TNF/16S.ITS/jobs/6_correlogram/"

############################

## TNFi: pre-maintenance ##

# FDR < 0.1
d <- read.table(file = paste(dir, "TNF.B.C_cytokines_metabolites_spearman_0.1fdr.txt", sep = ""),
                header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                na.strings = "NA")

pdf(file = paste(dir, "TNF.B.C_corrplot_spearman_0.1fdr.pdf", sep = ""))
corrplot(as.matrix(d), method = "circle", diag = TRUE, na.label = " ",
         is.corr = TRUE, cl.ratio = 0.18, cl.cex = 0.7,
         tl.cex = 0.7, tl.col = "black")
dev.off()

# FDR < 0.2 #
d <- read.table(file = paste(dir, "TNF.B.C_cytokines_metabolites_spearman_0.2fdr.txt", sep = ""),
                header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                na.strings = "NA")

pdf(file = paste(dir, "TNF.B.C_corrplot_spearman_0.2fdr.pdf", sep = ""))
corrplot(as.matrix(d), method = "circle", diag = TRUE, na.label = " ",
         is.corr = TRUE, cl.ratio = 0.25, cl.cex = 0.7,
         tl.cex = 0.7, tl.col = "black")
dev.off()

############################

## IL-17i: pre-loading-maintenance ##

# FDR < 0.1 #
d <- read.table(file = paste(dir, "IL17.B.C.D_cytokines_metabolites_spearman_0.1fdr.txt", sep = ""),
                header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                na.strings = "NA")

pdf(file = paste(dir, "IL17.B.C.D_corrplot_spearman_0.1fdr.pdf", sep = ""))
corrplot(as.matrix(d), method = "circle", diag = TRUE, na.label = " ",
         is.corr = TRUE, cl.ratio = 0.18, cl.cex = 0.7,
         tl.cex = 0.7, tl.col = "black")
dev.off()

# FDR < 0.2 #
d <- read.table(file = paste(dir, "IL17.B.C.D_cytokines_metabolites_spearman_0.2fdr.txt", sep = ""),
                header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                na.strings = "NA")

pdf(file = paste(dir, "IL17.B.C.D_corrplot_spearman_0.2fdr.pdf", sep = ""))
corrplot(as.matrix(d), method = "circle", diag = TRUE, na.label = " ",
         is.corr = TRUE, cl.ratio = 0.18, cl.cex = 0.7,
         tl.cex = 0.7, tl.col = "black")
dev.off()


############################################################################
############################################################################
############################################################################

### Statistics for proteins/cytokines ###

### Input files:
### Column 1: sample ID
### Column 2: timepoint
### Remainder of columns: protein/cytokine quantities across samples

## TNFi ##

# directory
dir = ".../IL17.TNF/16S/jobs/6_cytokines_R/"

# read in data table
d.cytokines <- read.table(file = paste(dir, "cytokines_TNF.B.C.txt", sep = ""),
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                          na.strings = "NA")

# convert to data frame
d.cytokines <- as.data.frame(d.cytokines)

# create table for storing wilcoxon results
wilc.table <- matrix(data = NA, nrow = 13, ncol = 2)
colnames(wilc.table) <- c("cytokine", "pre-maintenance")

# calculate wilcoxon between pre/maintenance visits for each cytokine
for (i in 1:13) {
  wilc.table[i,1] <- colnames(d.cytokines)[i+1]
  wilc.table[i,2] <- wilcox.test(d.cytokines[,i+1] ~ Pre_post_det, data = d.cytokines, paired = TRUE)$p.value
}

# save
ft = paste(dir, "cytokines_TNF.B.C_wilcoxon.csv", sep = "")
write.csv(file = ft, wilc.table)

############################

## IL17-i ###

# directory
dir = ".../IL17.TNF/16S/jobs/6_cytokines_R/"

# read in data table
d.cytokines <- read.table(file = paste(dir, "cytokines_IL17.B.C.D.txt", sep = ""),
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                          na.strings = "NA")

# convert to data frame
d.cytokines <- as.data.frame(d.cytokines)

##########

# subset pre/loading samples
d.load <- d.cytokines[d.cytokines$Pre_post_det != "3_maint", ]

# subset fecal samples
d.load.fecal <- d.load[-c(10,12:14)]
d.load.fecal <- na.omit(d.load.fecal)
d.load.fecal.25.33 <- d.load[c(1,10,12)]

# create table for storing wilcoxon results
wilc.table.load.fecal.a <- matrix(data = NA, nrow = 9, ncol = 2)
colnames(wilc.table.load.fecal.a) <- c("cytokine", "pre-loading")

wilc.table.load.fecal.b <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(wilc.table.load.fecal.b) <- c("cytokine", "pre-loading")

# calculate wilcoxon between pre/loading visits for fecal cytokines
for (i in 1:9) {
  wilc.table.load.fecal.a[i,1] <- colnames(d.load.fecal)[i+1]
  wilc.table.load.fecal.a[i,2] <- wilcox.test(d.load.fecal[,i+1] ~ Pre_post_det, data = d.load.fecal, paired = TRUE)$p.value
}

for (i in 1:2) {
  wilc.table.load.fecal.b[i,1] <- colnames(d.load.fecal.25.33)[i+1]
  wilc.table.load.fecal.b[i,2] <- wilcox.test(d.load.fecal.25.33[,i+1] ~ Pre_post_det, data = d.load.fecal.25.33, paired = TRUE)$p.value
}

# save
ft.load.fecal.a = paste(dir, "cytokines_IL17.B.C.D_wilcoxon.loading_25.33_fecal.a.csv", sep = "")
write.csv(file = ft.load.fecal.a, wilc.table.load.fecal.a)

ft.load.fecal.b = paste(dir, "cytokines_IL17.B.C.D_wilcoxon.loading_25.33_fecal.b.csv", sep = "")
write.csv(file = ft.load.fecal.b, wilc.table.load.fecal.b)

##########

# subset pre/maintenance samples
d.maint <- d.cytokines[d.cytokines$Pre_post_det != "2_load", ]

# subset fecal samples
d.maint.fecal <- d.maint[-c(10,12:14)]
d.maint.fecal <- na.omit(d.maint.fecal)
d.maint.fecal <- d.maint.fecal[-c(9, 12, 13, 18), ]
d.maint.fecal.25.33 <- d.maint[c(1,10,12)]
d.maint.fecal.25.33 <- d.maint.fecal.25.33[-c(9, 12, 15, 20), ]

# create table for storing wilcoxon results
wilc.table.maint.fecal.a <- matrix(data = NA, nrow = 9, ncol = 2)
colnames(wilc.table.maint.fecal.a) <- c("cytokine", "pre-maintenance")

wilc.table.maint.fecal.b <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(wilc.table.maint.fecal.b) <- c("cytokine", "pre-maintenance")

# calculate wilcoxon between pre/maintenance visits for fecal cytokines
for (i in 1:9) {
  wilc.table.maint.fecal.a[i,1] <- colnames(d.maint.fecal)[i+1]
  wilc.table.maint.fecal.a[i,2] <- wilcox.test(d.maint.fecal[,i+1] ~ Pre_post_det, data = d.maint.fecal, paired = TRUE)$p.value
}

for (i in 1:2) {
  wilc.table.maint.fecal.b[i,1] <- colnames(d.maint.fecal.25.33)[i+1]
  wilc.table.maint.fecal.b[i,2] <- wilcox.test(d.maint.fecal.25.33[,i+1] ~ Pre_post_det, data = d.maint.fecal.25.33, paired = TRUE)$p.value
}

# save
ft.maint.fecal.a = paste(dir, "cytokines_IL17.B.C.D_wilcoxon.maintenance_25.33.a_fecal.csv", sep = "")
write.csv(file = ft.maint.fecal.a, wilc.table.maint.fecal.a)

ft.maint.fecal.b = paste(dir, "cytokines_IL17.B.C.D_wilcoxon.maintenance_25.33.b_fecal.csv", sep = "")
write.csv(file = ft.maint.fecal.b, wilc.table.maint.fecal.b)

############################################################################
############################################################################
############################################################################

### Statistics for metabolites ###

### Input files:
### Column 1: sample ID
### Column 2: timepoint
### Remainder of columns: metabolite quantities across samples

## TNFi ##

# directory
dir = ".../IL17.TNF/16S/jobs/7_metabolites_R/"

# read in data table
d.metabolites <- read.table(file = paste(dir, "metabolites_TNF.B.C.txt", sep = ""),
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                          na.strings = "NA")

# convert to data frame
d.metabolites <- as.data.frame(d.metabolites)

# create table for storing wilcoxon results
wilc.table <- matrix(data = NA, nrow = 20, ncol = 2)
colnames(wilc.table) <- c("metabolite", "pre-maintenance")

# calculate wilcoxon between pre/maintenance visits for each cytokine
for (i in 1:20) {
  wilc.table[i,1] <- colnames(d.metabolites)[i+1]
  wilc.table[i,2] <- wilcox.test(d.metabolites[,i+1] ~ Pre_post_det, data = d.metabolites, paired = TRUE)$p.value
}

# save
ft = paste(dir, "metabolites_TNF.B.C_wilcoxon.csv", sep = "")
write.csv(file = ft, wilc.table)

############################

## IL-17i ##

# directory
dir = ".../IL17.TNF/16S/jobs/7_metabolites_R/"

# read in data table
d.metabolites <- read.table(file = paste(dir, "metabolites_IL17.B.C.D.txt", sep = ""),
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                          na.strings = "NA")

# convert to data frame
d.metabolites <- as.data.frame(d.metabolites)

##########

# subset pre/loading samples
d.load <- d.metabolites[d.metabolites$Pre_post_det != "3_maint", ]

# create table for storing wilcoxon results
wilc.table.load <- matrix(data = NA, nrow = 20, ncol = 2)
colnames(wilc.table.load) <- c("metabolite", "pre-loading")

# calculate wilcoxon between pre/loading visits for each cytokine
for (i in 1:20) {
  wilc.table.load[i,1] <- colnames(d.load)[i+1]
  wilc.table.load[i,2] <- wilcox.test(d.load[,i+1] ~ Pre_post_det, data = d.load, paired = TRUE)$p.value
}

# save
ft.load = paste(dir, "metabolites_IL17.B.C.D_wilcoxon.loading.csv", sep = "")
write.csv(file = ft.load, wilc.table.load)

##########

# subset pre/maintenance samples
d.maint <- d.metabolites[d.metabolites$Pre_post_det != "2_load", ]
d.maint <- d.maint[-c(9, 12, 15, 20), ]

# create table for storing wilcoxon results
wilc.table.maint <- matrix(data = NA, nrow = 20, ncol = 2)
colnames(wilc.table.maint) <- c("metabolite", "pre-maintenance")

# calculate wilcoxon between pre/maintenance visits for each metabolite
for (i in 1:20) {
  wilc.table.maint[i,1] <- colnames(d.maint)[i+1]
  wilc.table.maint[i,2] <- wilcox.test(d.maint[,i+1] ~ Pre_post_det, data = d.maint, paired = TRUE)$p.value
}

# save
ft.maint = paste(dir, "metabolites_IL17.B.C.D_wilcoxon.maintenance.csv", sep = "")
write.csv(file = ft.maint, wilc.table.maint)

