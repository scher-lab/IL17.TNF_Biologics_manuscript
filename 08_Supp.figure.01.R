############################################ 
## R script                               ##
## Project: IL17.TNF_Biologics_manuscript ##
## Supplementary figure 01                ##
## 16S data                               ##
############################################

### Brief description:

### Supplementary figure 01:
### Panel A: Taxa summary at the order level
### Panel B: LEfSe analysis

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(ggplot2)
library(ggthemes)
library(data.table)
library(scales)
library(RColorBrewer)
library(pals)
library(ggplot2)
library(reshape2)

############################################################################
############################################################################
############################################################################

### Create phyloseq object ###

# import biom file
b = ".../IL17.TNF/16S/outputs/13_biom_revision/otu_table_Hlt.vs.TNF.B.C.json"
biom = import_biom(b, taxaPrefix = F)

# import mapping file
# must leave A1 cell of mapping file empty for R compatibility
m = ".../IL17.TNF/16S/inputs/map/revision/Map_IL17.TNF_16S_human_Hlt.vs.TNF.B.C_R.txt"
map = sample_data(read.table(m, header = TRUE, sep = "\t", row.names = 1))

# create phyloseq object
ph = phyloseq(otu_table(biom), tax_table(biom), map)

# provide column names to separate different taxonomic levels
colnames(tax_table(ph)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# load the tree file (use 'unannotated.tree')
t = ".../greengenes/gg_13_8_otus/trees/97_otus_unannotated.tree"
tree = import_qiime(treefilename = t) 

# merge tree with phyloseq object
phy_16S = merge_phyloseq(ph, map, tree)

# print row and column names
rownames(sample_data(phy_16S))
colnames(sample_data(phy_16S))

# subset to only Healthy and TNF baseline samples
phy_16S_Hlt.TNF.B <- subset_samples(phy_16S, Timepoint_revised == "B")

############################################################################
############################################################################
############################################################################

### Taxa summary ###

## Setup ##

# background
bkg <- theme_classic() +
  theme(axis.text.x = element_text(size = 16, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(axis.title.y = element_text(size = 20, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 12, face = "bold", color = "black")) +
  theme(legend.spacing.x = unit(0.3, "cm")) +
  theme(legend.title = element_blank()) +
  theme(legend.justification = "top") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# directory for storing files
dir = ".../IL17.TNF/16S/jobs/3_taxa.summary_R/"

# list of phyloseq objects to process
pseqs <- list(physeq = phy_16S_Hlt.TNF.B)

# list of phyloseq names
pseq.names <- c("16S_human_Hlt.TNF.B")

# list of taxonomic ranks to process
tax.rank <- c("Order")

############################

## Plotting ##

# for each phyloseq object in pseqs, perform taxa summary by treatment-visit
for (i in seq_along(pseqs)) {
  
  # create new directory to store results
  dd = paste(pseq.names[i], "tx", sep = "_")
  dir.new = paste(dir, dd, "/", sep = "")
  dir.create(dir.new)
  
  # remove empty OTUs
  phylo <- subset_taxa(pseqs[[i]], rowSums(otu_table(pseqs[[i]])) > 0)
  
  # create filenames
  filename_table = paste(pseq.names[i], "taxa.summary", "tx", tax.rank[j], "table.csv", sep = "_")    
  filename_plot = paste(pseq.names[i], "taxa.summary", "tx", tax.rank[j], "plot.pdf", sep = "_")    
    
  # combine OTUs by taxanomic rank
  t <- tax_glom(phylo, taxrank = tax.rank[j])
    
  # merge samples according to mapping file category
  m <- merge_samples(t, group = "Treatment")
    
  # transform to relative abundance
  tm <- transform_sample_counts(m, rel_abundance)
        
  r <- merge(tax_table(tm)[,1:4], t(otu_table(tm)), by = "row.names", all = TRUE) # merge tax table and rel abundance counts from OTU table
  r$rel.abund.sum <- rowSums(r[,6:7]) # sum rel abundance counts across treatments
  r.sorted <- r[order(-r$rel.abund.sum),] # sort in descending order
  ft = paste(dir.new, filename_table, sep = "") # save
  write.csv(r.sorted, file = ft)
    
  # choose top taxa in alphabetical order for plot legend
  top.taxa <- sort(r.sorted$Order[1:12])
    
  # choose color palatte based on total number of taxa
  c <- sample(kovesi.rainbow_bgyr_35_85_c72(length(r.sorted[,$Order))
    
  # plot and save
  p <- plot_bar(tm, x = "Sample", fill = tax.rank[j]) + 
    scale_fill_manual(breaks = top.taxa, values = c) + # only display top taxa in legend
    scale_x_discrete(labels = c("Healthy", "TNFi\npre")) +
    xlab("") + ylab("Relative abundance") +
    bkg 
         
  fp = paste(dir.new, filename_plot, sep = "")
  pdf(file = fp, width = 6)
  plot(p)
  dev.off()
}

############################################################################
############################################################################
############################################################################

### LEfSe LDA ###

### Analysis was performed using LEfSe v1.0.7.
### LEfSe input file included taxa at the genus level.
### Only taxa that had relative abundance of at least 0.1% in at least one sample were maintained.
### Output file generated by run_lefse.py command without subclass analysis (.res) was used to create LDA graph (panel A).
### FDR significance was calculated using group_significance.py command in QIIME v1.9.1 via Kruksal-Wallis.

## Setup ##

# colors
col1 <- c("#929aab", "#ffb400")

# background theme
bkg <- theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 13, color = "black")) +
  theme(axis.title.x = element_text(size = 16, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 12, face = "bold", color = "black")) +
  theme(legend.spacing.x = unit(0.3, "cm")) +
  theme(legend.title = element_blank()) +
  theme(legend.justification = "top")

# directory for storing files
dir = ".../IL17.TNF/16S/jobs/4_lefse/"

# lefse files
lefse.files <- c("1_Hlt.vs.TNF_0.001relabund/human_Hlt.vs.TNF.baseline_genus_0.001relabund.res")
lefse.names <- c("16S_human_Hlt.vs.TNF.baseline_0.001relabund")

############################

## Plotting ##

# for each lefse file, create an LDA plot
for (i in seq_along(lefse.files)) {
  
  # read in LEfSe output file
  f = paste(dir, lefse.files[i], sep = "")
  lef <- read.table(f, sep = "\t", header = FALSE)
  
  # change column headings
  colnames(lef) = c("taxa", "log_highest_class_avg", "class", "LDA", "p_value")
  
  # remove all non-significant entries
  lef.subset <- subset(lef, LDA != "NA")
  
  # add column with abridged taxa names
  for (j in 1:nrow(lef.subset)) {
    
    # take the 2 lowest taxa classifications
    s <- lapply(strsplit(as.character(lef.subset[[1]]), '\\.'), tail, n = 2)
    
    # add them as separate columns
    t1 <- s[[j]][1]
    t2 <- s[[j]][2]
    
    lef.subset$taxa_1[j] <- paste(t1)
    if(is.na(t2)) {
      lef.subset$taxa_2[j] <- paste(t1)
      sr <- t1
    }
    else {
      lef.subset$taxa_2[j] <- paste(t2)
      sr <- t2
    }
    
    # use the lowest taxa classification as the shortened taxa name
    sr <- gsub("__c_", "](c)", sr)
    sr <- gsub("__o_", "](o)", sr)
    sr <- gsub("__f_", "](f)", sr)
    sr <- gsub("_c_", "(c)", sr)
    sr <- gsub("_o_", "(o)", sr)
    sr <- gsub("_f_", "(f)", sr)
    sr <- gsub("___", "__[", sr)
    sr <- gsub("UC__", "UC_[", sr)
    sr <- gsub("_$", "]", sr)
    lef.subset$taxa_s[j] <- paste(sr)
  }
  
  # find duplicates of shortened taxa names
  dt <- data.table(lef.subset$taxa_s)
  dup <- unlist(unique(dt[duplicated(dt) | duplicated(dt, fromLast=TRUE)]))
  
  # change names of duplicates to allow for appropriate data plotting
  if (length(dup) > 0) {
    for (j in seq_along(dup)) {
      for (k in 1:nrow(lef.subset)) {
        if (lef.subset$taxa_s[k] == dup[[j]]) {
          lef.subset$taxa_s[k] <- paste(lef.subset$taxa_1[k], lef.subset$taxa_2[k], sep = "|")
        }
      }
    }
  }
  
  # order factor levels for plotting
  lef.subset$taxa_s <- factor(lef.subset$taxa_s, levels = lef.subset$taxa_s[rev(order(lef.subset$class, lef.subset$LDA))])
  
  # plot
  p <- ggplot(data = lef.subset, aes(x = taxa_s, y = LDA, fill = class)) +
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = col1, labels = c("Healthy", "TNFi")) +
    xlab(NULL) + ylab("\nLDA score (log 10)") +
    bkg
  
  # create filename
  filename_plot = paste(lefse.names[i], "lefse", "LDA.plot.pdf", sep = "_")  
  
  # save plot   
  fp = paste(dir, filename_plot, sep = "")
  pdf(file = fp, height = 6, width = 7)
  plot(p)
  dev.off()
}

############################################################################
############################################################################
############################################################################

### LEfSe relative abundance ###

### Input file:
### Column 1: SampleID
### Column 2: Treatment (Healthy vs TNF)
### Remainder of columns: relative abundance across samples for each taxon 
### that was identified to be significant by LEfSe analysis

## Setup ##

# colors
col1 <- c("#929aab", "#ffb400")

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 9, color = "black")) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0, size = 11, color = "black")) +
  theme(axis.title.y = element_text(angle = 180, hjust = 0, size = 13, color = "black", face = "bold", margin = unit(c(0, 6, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8, face = "bold", color = "black")) +
  theme(legend.spacing.x = unit(0.2, "cm")) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")

# function for plot
stats.dot.median <- function(x) {
  m <- median(x)
  min <- min(x)
  max <- max(x)
  return(c(y = m, ymin = min, ymax = max))
}

############################

## Plotting ##

# read in table
d.rel.abund <- read.table(file = ".../IL17.TNF/16S/jobs/5_lefse.taxa_rel.abund/Hlt.vs.TNF_0.001relabund.txt",
                          header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                          na.strings = "NA")

# convert to data frame
d.rel.abund <- as.data.frame(d.rel.abund)

# melt data for plotting
d <- melt(d.rel.abund, id = "Treatment")

# plot and save
p <- ggplot(d, aes(x = variable, y = value, color = Treatment)) +
  geom_hline(yintercept=0, alpha = 0.1) +
  stat_summary(fun.data = stats.dot.median, geom = "pointrange", 
               position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_manual(values = col1, labels = c("Healthy", "TNFi")) +
  scale_x_discrete(position = "top") +
  # sqrt scale
  scale_y_continuous(trans = sqrt_trans(), position = "right", limits = c(0,0.2)) +
  xlab("") + ylab("") +
  bkg 

pdf(file = ".../IL17.TNF/16S/jobs/5_lefse.taxa_rel.abund/Hlt.vs.TNF_0.001relabund_plot.pdf",
    height = 5.5, width = 4)
plot(p)
dev.off()

# save plot stats
s <- aggregate(value ~ variable+Treatment, data = d, stats.dot.median)
s <- t(s)

write.csv(s, file = ".../IL17.TNF/16S/jobs/5_lefse.taxa_rel.abund/Hlt.vs.TNF_0.001relabund_stats.csv")




