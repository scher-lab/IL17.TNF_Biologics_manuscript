############################################ 
## R script                               ##
## Project: IL17.TNF Biologics manuscript ##
## Figure 05				  ##
## Supplementary figure 06		  ##
## Supplementary table 04 		  ##
## ITS data				  ##
############################################

### Brief description:

### Figure 05, Supplementary figure 06:
### Panels A and D: TNFi relative abundance line plots
### Panels B and E: IL-17i relative abundance line plots
### Panels C and F: Boxplots representing magnitude of relative abundance change in TNFi and IL-17i subsets

### Supplementary table 04:
### P-value calculations of changes in taxa relative abundance pre-post treatment with TNFi and IL-17i

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(PMCMR)
library(PMCMRplus)

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
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
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

### Line plots of relative abudnance over time ###

## Setup ##

# colors
# green: expansion
# orange: contraction
# grey: <10% change
col1 <- c("#53d397", "#fc8a15", "#bbbbbb")

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 21, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 22, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 13, face = "bold", color = "black")) +
  theme(legend.spacing.x = unit(0.3, "cm")) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top", legend.box = "horizontal")

# directory for storing files
dir = ".../IL17.TNF/ITS/jobs/5_specific.taxa_lineplots_R/"

############################

## Plotting - TNFi subset ##

# list of phyloseq objects to process
pseqs <- list(physeq = phy_ITS.R1_human_TNF.B.C)

# list of phyloseq names
pseq.names <- c("ITS.R1_human_TNF.B.C")

# list of specific taxa 
taxa.names <- c("o__Saccharomycetales",
                "g__Candida", 
		"s__Candida_albicans", "s__Saccharomyces_cerevisiae") 

# for each phyloseq object in pseqs, calculate the relative abudnance over time for specific taxa
for (i in seq_along(pseqs)) {
  
  # merge phyloseq objects at each taxnomic level
  order <- tax_glom(pseqs[[i]], taxrank = "Order", NArm = FALSE)
  genus <- tax_glom(pseqs[[i]], taxrank = "Genus", NArm = FALSE)
  species <- tax_glom(pseqs[[i]], taxrank = "Species", NArm = FALSE)
  
  # transform to relative abundance
  rel_abundance = function(x) { x/sum(x) } # function
  to <- transform_sample_counts(order, rel_abundance)
  tg <- transform_sample_counts(genus, rel_abundance)
  ts <- transform_sample_counts(species, rel_abundance)
  
  for (j in seq_along(taxa.names)) {
    
    # for each taxon, create new directory to store results
    dir.taxa = paste(dir, taxa.names[j], "/", sep = "")
    if (dir.exists(dir.taxa) == FALSE) {
      dir.create(dir.taxa)
    }  
    
    # create filenames
    filename_table = paste(pseq.names[i], taxa.names[j], "rel.abund", "table.csv", sep = "_")
    filename_plot = paste(pseq.names[i], taxa.names[j], "rel.abund", "plot.pdf", sep = "_")  
    
    # extract taxon relative abundance
    if (grepl("o__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(to, Order == taxa.names[j])
    } else if (grepl("g__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tg, Genus == taxa.names[j])
    } else if (grepl("s__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(ts, Species == taxa.names[j])
    }
    
    # create dataset from which to generate plot
    d <- as.data.frame(merge(sample_data(taxa.rel.abund), t(otu_table(taxa.rel.abund)), 
                             by = "row.names", all = TRUE))  
    
    # rename column with relative abundance values for easier plotting 
    colnames(d)[ncol(d)] <- "Taxa_rel_abundance" 
    
    # create categories by which to color lineplot: 1) >= 10% expansion, 2) >= 10% contraction, 3) <10% change
    # subset and spread dataset into pre/post columns
    # then calculate delta relative abundance
    d.rel <- d %>% 
      subset(select = c("Subject", "Pre_post", "Taxa_rel_abundance")) %>%
      spread(key = "Pre_post", value = "Taxa_rel_abundance") %>%
      mutate(Delta_rel_abundance = (`2_post`-`1_pre`)) %>%
      mutate(
        Change_type = case_when(
          abs(Delta_rel_abundance) >= 0.1 & sign(Delta_rel_abundance) == 1 ~ "1_expansion",
          abs(Delta_rel_abundance) >= 0.1 & sign(Delta_rel_abundance) == -1 ~ "2_contraction",
          TRUE ~ "3_no.change"
        )
      )
    
    # merge with original dataset
    d.final <- as.data.frame(merge(d, d.rel, by = "Subject", all = TRUE)) 
    
    # save
    ft = paste(dir.taxa, filename_table, sep = "")
    write.csv(d.final, file = ft)
    
    # create line plot and save
    p <- ggplot(data = d.final, aes(x = Timepoint_revised, y = Taxa_rel_abundance, group = Subject)) + 
      geom_line(aes(color = Change_type), size = 1) +
      geom_text_repel(data = subset(d, Timepoint_revised == "B"), 
                      aes(label = Subject), nudge_x = -0.3, size = 3, color = "black", segment.alpha = 0.3) +
      geom_point(shape = 20, size = 3, color = "black") +
      scale_color_manual(values = col1, labels = c(">=10% Expansion", ">=10% Contraction", "<10% Change")) +
      scale_x_discrete(labels=c("TNFi\npre", "TNFi\nmaint")) +
      scale_y_continuous(limits = c(0,1)) +
      xlab(NULL) + ylab("Relative abundance") + 
      bkg 
    
    fp = paste(dir.taxa, filename_plot, sep = "")
    pdf(file = fp)
    plot(p)
    dev.off()
  }
}

############################

## Plotting - IL-17i subset ##

# list of phyloseq objects to process
pseqs <- list(physeq = phy_ITS.R1_human_IL17.B.C.D)

# list of phyloseq names
pseq.names <- c("ITS.R1_human_IL17.B.C.D")

# list of specific taxa 
taxa.names <- c("o__Saccharomycetales",
                "g__Candida", 
                "s__Candida_albicans", "s__Saccharomyces_cerevisiae") 

# for each phyloseq object in pseqs, calculate the relative abudnance over time for specific taxa
for (i in seq_along(pseqs)) {
  
  # merge phyloseq objects to each taxnomic level
  order <- tax_glom(pseqs[[i]], taxrank = "Order", NArm = FALSE)
  genus <- tax_glom(pseqs[[i]], taxrank = "Genus", NArm = FALSE)
  species <- tax_glom(pseqs[[i]], taxrank = "Species", NArm = FALSE)
  
  # transform to relative abundance
  rel_abundance = function(x) { x/sum(x) } # function
  to <- transform_sample_counts(order, rel_abundance)
  tg <- transform_sample_counts(genus, rel_abundance)
  ts <- transform_sample_counts(species, rel_abundance)
  
  for (j in seq_along(taxa.names)) {
    
    # for each taxon, create new directory to store results
    dir.taxa = paste(dir, taxa.names[j], "/", sep = "")
    if (dir.exists(dir.taxa) == FALSE) {
      dir.create(dir.taxa)
    }  
    
    # create filenames
    filename_table = paste(pseq.names[i], taxa.names[j], "rel.abund", "table.csv", sep = "_")
    filename_plot = paste(pseq.names[i], taxa.names[j], "rel.abund", "plot.pdf", sep = "_")  
    
    # extract relative abundance of specific taxa
    if (grepl("o__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(to, Order == taxa.names[j])
    } else if (grepl("g__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tg, Genus == taxa.names[j])
    } else if (grepl("s__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(ts, Species == taxa.names[j])
    }
    
    # create dataset from which to generate plot
    d <- as.data.frame(merge(sample_data(taxa.rel.abund), t(otu_table(taxa.rel.abund)), 
                             by = "row.names", all = TRUE))  
    
    # rename column with relative abundance values for easier plotting 
    colnames(d)[ncol(d)] <- "Taxa_rel_abundance" 
    
    # create categories by which to color lineplot: 1) >= 10% expansion, 2) >= 10% contraction, 3) <10% change
    # subset and spread dataset into timepoint columns
    # then calculate delta relative abundance between B-C and C-D timepoints
    d.IL17 <- d %>%
      subset(select = c("Subject", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      mutate(Delta_rel_abundance.BC = (C-B)) %>%
      mutate(
        Change_type.BC = case_when(
          abs(Delta_rel_abundance.BC) >= 0.1 & sign(Delta_rel_abundance.BC) == 1 ~ "1_expansion",
          abs(Delta_rel_abundance.BC) >= 0.1 & sign(Delta_rel_abundance.BC) == -1 ~ "2_contraction",
          TRUE ~ "3_no.change"
        )
      ) %>%
      mutate(Delta_rel_abundance.CD = (D-C)) %>%
      mutate(
        Change_type.CD = case_when(
          abs(Delta_rel_abundance.CD) >= 0.1 & sign(Delta_rel_abundance.CD) == 1 ~ "1_expansion",
          abs(Delta_rel_abundance.CD) >= 0.1 & sign(Delta_rel_abundance.CD) == -1 ~ "2_contraction",
          is.na(Delta_rel_abundance.CD) ~ "4_na",
          TRUE ~ "3_no.change"
        )
      )
    
    # merge with original dataset
    d.final <- as.data.frame(merge(d, d.IL17, by = "Subject", all = TRUE))
    
    # save
    ft = paste(dir.taxa, filename_table, sep = "")
    write.csv(d.final, file = ft)
    
    # create line plot and save
    p <- ggplot(data = d.final, aes(x = Timepoint_revised, y = Taxa_rel_abundance, group = Subject)) + 
      geom_line(data = subset(d.final, Timepoint_revised != "D"), 
                aes(group = Subject, color = Change_type.BC), size = 1) +
      geom_line(data = subset(d.final, (Change_type.CD != "4_na" & Timepoint_revised != "B")), 
                aes(group = Subject, color = Change_type.CD), size = 1) +
      geom_text_repel(data = subset(d, Timepoint_revised == "B"), 
                      aes(label = Subject), nudge_x = -0.3, size = 3, color = "black", segment.alpha = 0.3) +
      geom_point(shape = 20, size = 3, color = "black") +
      scale_color_manual(values = col1, labels = c(">=10% Expansion", ">=10% Contraction", "<10% Change")) +
      scale_x_discrete(labels=c("IL-17i\npre", "IL-17i\nload", "IL-17i\nmaint")) +
      scale_y_continuous(limits = c(0,1)) +
      xlab(NULL) + ylab("Relative abundance") + 
      bkg
    
    fp = paste(dir.taxa, filename_plot, sep = "")
    pdf(file = fp)
    plot(p)
    dev.off()
  }
}

############################################################################
############################################################################
############################################################################

### Boxplots of delta relative abundance over time ###

## Setup ##

# colors
col1 <- c("#c40018", "#448ef6", "#8559a5")

# shapes
shape1 <- c(16,15,17)

# background theme
bkg <- theme_classic() +
  theme(axis.text.x = element_text(size = 20, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 22, face = "bold", color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.position = "none")

# directory for storing files
dir = ".../IL17.TNF/ITS/jobs/5_specific.taxa_delta.rel.abund_R/"

# list of phyloseq objects to process
pseqs <- list(physeq = phy_ITS.R1_human_TNF.B.C_IL17.B.C.D)

# list of phyloseq names
pseq.names <- c("ITS.R1_human_TNF.B.C_IL17.B.C.D")

# list of specific taxa 
taxa.names <- c("o__Saccharomycetales",
                "g__Candida", 
                "s__Candida_albicans", "s__Saccharomyces_cerevisiae") 

#############################

## Plotting ##

# for each phyloseq object in pseqs, calculate the relative abudnance over time for specific taxa
for (i in seq_along(pseqs)) {
  
  # merge phyloseq objects
  order <- tax_glom(pseqs[[i]], taxrank = "Order", NArm = FALSE)
  genus <- tax_glom(pseqs[[i]], taxrank = "Genus", NArm = FALSE)
  species <- tax_glom(pseqs[[i]], taxrank = "Species", NArm = FALSE)
  
  # transform to relative abundance
  to <- transform_sample_counts(order, rel_abundance)
  tg <- transform_sample_counts(genus, rel_abundance)
  ts <- transform_sample_counts(species, rel_abundance)
  
  for (j in seq_along(taxa.names)) {
    
    # for each taxon, create new directory to store results
    dir.taxa = paste(dir, taxa.names[j], "/", sep = "")
    if (dir.exists(dir.taxa) == FALSE) {
      dir.create(dir.taxa)
    }
    
    # create filenames
    filename_table = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "table_all.csv", sep = "_") 
    filename_table_TNF = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "table_TNF.csv", sep = "_") 
    filename_table_IL17.load = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "table_IL17.load.csv", sep = "_") 
    filename_table_IL17.maint = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "table_IL17.maint.csv", sep = "_") 
    filename_table_for.plot = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "table_for.plot.csv", sep = "_")
    filename_plot = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "plot.pdf", sep = "_")  
    filename_plot.stats = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "plot.stats.csv", sep = "_")
    filename_stats = paste(pseq.names[i], taxa.names[j], "delta.rel.abund", "stats.txt", sep = "_")     
    
    # relative abundance of specific taxa
    if (grepl("o__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(to, Order == taxa.names[j])
    } else if (grepl("g__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tg, Genus == taxa.names[j])
    } else if (grepl("s__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(ts, Species == taxa.names[j])
    } 
    
    # create dataset from which to generate plot
    d <- as.data.frame(merge(sample_data(taxa.rel.abund), t(otu_table(taxa.rel.abund)), 
                             by = "row.names", all = TRUE))  
    
    # rename column with relative abundance values for easier plotting 
    colnames(d)[ncol(d)] <- "Taxa_rel_abundance" 
    
    # save table
    ftab = paste(dir.taxa, filename_table, sep = "")
    write.csv(d, file = ftab) 
    
    ############################
    
    ## TNFi ##
	  
    # spread data across timepoints
    # calculate magnitude delta relative abundance
    d.TNF <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "1_TNF") %>% 
      filter(Timepoint_revised != "D") %>% 
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      mutate(Delta_rel_abundance = abs(C-B))
    
    # rename columns
    colnames(d.TNF)[colnames(d.TNF)=="B"] <- "1_pre"  
    colnames(d.TNF)[colnames(d.TNF)=="C"] <- "2_post"  
    
    # save table
    ftab.TNF = paste(dir.taxa, filename_table_TNF, sep = "")
    write.csv(d.TNF, file = ftab.TNF)
    
    ############################
    
    ## IL-17i loading ##
    
    # spread data across timepoints
    # calculate magnitude delta relative abundance
    d.IL17.load <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "2_IL17") %>% 
      filter(Timepoint_revised != "D") %>%
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      mutate(Delta_rel_abundance = abs(C-B))
    
    # rename treatment
    d.IL17.load$Treatment <- replace(as.character(d.IL17.load$Treatment), 
                                     as.character(d.IL17.load$Treatment) == "2_IL17", "2_IL17.load")
    # rename columns
    colnames(d.IL17.load)[colnames(d.IL17.load)=="B"] <- "1_pre"  
    colnames(d.IL17.load)[colnames(d.IL17.load)=="C"] <- "2_post"
    
    # save table
    ftab.IL17.load = paste(dir.taxa, filename_table_IL17.load, sep = "")
    write.csv(d.IL17.load, file = ftab.IL17.load)
    
    ############################
    
    ## IL-17i maintenance ##
    
    # spread data across timepoints
    # calculate magnitude delta relative abundance
    d.IL17.maint <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "2_IL17") %>% 
      filter(Timepoint_revised != "C") %>%
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      filter(!is.na(`D`)) %>%
      mutate(Delta_rel_abundance = abs(D-B))
    
    # rename treatment
    d.IL17.maint$Treatment <- replace(as.character(d.IL17.maint$Treatment), as.character(d.IL17.maint$Treatment) == "2_IL17", "3_IL17.maint")
    
    # rename columns
    colnames(d.IL17.maint)[colnames(d.IL17.maint)=="B"] <- "1_pre"  
    colnames(d.IL17.maint)[colnames(d.IL17.maint)=="D"] <- "2_post"
    
    # save table
    ftab.IL17.maint = paste(dir.taxa, filename_table_IL17.maint, sep = "")
    write.csv(d.IL17.maint, file = ftab.IL17.maint)
    
    ############################
    
    # join subsets
    d.final <- rbind(d.TNF, d.IL17.load, d.IL17.maint)
    
    # save dataset used for generating plot
    ftab.plot = paste(dir.taxa, filename_table_for.plot, sep = "")
    write.csv(d.final, file = ftab.plot)
    
    # plot TNFi vs IL-17i magnitude delta relative abundance and save
    p <- ggplot(data = d.final, aes(x = Treatment, y = Delta_rel_abundance, color = Treatment, shape = Treatment)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 1, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                   color = "black", size = 0.5, width = 0.5, fill = "white") +
      geom_jitter(width = 0.15, size = 3) +
      scale_color_manual(values = col1) +
      scale_x_discrete(labels = c("TNFi\npre-maint", "IL-17i\npre-load", "IL-17i\npre-maint")) +
      xlab(NULL) + ylab("Magnitude of relative abundance change") +
      ylim(0, 0.6) +
      bkg
    
    fp = paste(dir.taxa, filename_plot, sep = "")
    pdf(file = fp)
    plot(p)
    dev.off()
    
    # obtain plot statistics and save
    s <- aggregate(Delta_rel_abundance ~ Treatment, data = d.final, stats.all)
    s <- t(s)
    rownames(s) <- c("treatment", "mean", "sd", "median", "min", "max", 
                     "10%ile", "25%ile", "75%ile", "90%ile")
    
    fps = paste(dir.taxa, filename_plot.stats, sep="")
    write.csv(s, file = fps)
    
    #############################

    ## General Statistics ##
    
    # split IL-17i into loading and maintenance subcategories
    d.load <- d.final[d.final$Treatment != "3_IL17.maint", ]
    d.maint <- d.final[d.final$Treatment != "2_IL17.load", ]
    
    # mann-whitney
    mw.load <- wilcox.test(Delta_rel_abundance ~ Treatment, data = d.load, paired = FALSE)
    mw.maint <- wilcox.test(Delta_rel_abundance ~ Treatment, data = d.maint, paired = FALSE)
    
    # save calculations
    fs = paste(dir.taxa, filename_stats, sep = "")
    cat("Statistical Tests:\n\n", file = fs)
    cat("Delta relative abundance\n", file = fs, append = TRUE)
    cat(taxa.names[j], file = fs, append = TRUE) 
    cat("\n\n", file = fs, append = TRUE)
    cat(pseq.names[i], file = fs, append = TRUE)
    cat("\n\n==========================================\n\n", file = fs, append = TRUE)
    cat("Mann-Whitney test (non-parametric) - loading\n", file = fs, append = TRUE)
    capture.output(mw.load, file = fs, append = TRUE)
    cat("\n", file = fs, append = TRUE)
    cat("Mann-Whitney test (non-parametric) - maintenance\n", file = fs, append = TRUE)
    capture.output(mw.maint, file = fs, append = TRUE)
  }
}

###########################################################################
############################################################################
############################################################################

### P-value table of taxa relative abundance change between pre and post visits within TNFi and IL-17i cohorts ###

# list of phyloseq objects to process
pseqs <- list(physeq = phy_ITS.R1_human_TNF.B.C_IL17.B.C.D)

# list of specific taxa
taxa.names <- c("p__Ascomycota",
                "c__Saccharomycetes", "c__Tremellomycetes",
                "o__Saccharomycetales", "o__Tremellales", 
                "f__Saccharomycetaceae", "f__Mycosphaerellaceae", "f__Glomerellaceae",
                "f__Plectosphaerellaceae", "f__Sclerotiniaceae", "f__Mucoraceae",
                "f__Incertae_sedis",
                "g__Saccharomyces", "g__Candida", "g__Penicillium", 
                "g__Debaryomyces", "g__Glomerella", "g__Plectosphaerella",
                "g__Aspergillus",
                "s__Pichia_sp_1_TMS_2011", "s__Penicillium_charlesii", "s__Hypocrea_pseudokoningii",
                "s__Fusarium_sp_4_TMS_2011", "s__Candida_fermenticarens", "s__Trichoderma_longibrachiatum",
                "s__Rhodotorula_sp_RhSoW01", "s__Aspergillus_niger", "s__Saccharomyces_cerevisiae", 
                "s__Candida_albicans", "s__Candida_dubliniensis", "s__Candida_parapsilosis", 
                "s__Candida_tropicalis") 

# for each phyloseq object in pseqs, calculate the relative abudnance
# over time for specific taxa
for (i in seq_along(pseqs)) {
  
  # merge phyloseq objects to each taxnomic level
  phylum <- tax_glom(pseqs[[i]], taxrank = "Phylum", NArm = FALSE)
  class <- tax_glom(pseqs[[i]], taxrank = "Class", NArm = FALSE)
  order <- tax_glom(pseqs[[i]], taxrank = "Order", NArm = FALSE)
  family <- tax_glom(pseqs[[i]], taxrank = "Family", NArm = FALSE)
  genus <- tax_glom(pseqs[[i]], taxrank = "Genus", NArm = FALSE)
  species <- tax_glom(pseqs[[i]], taxrank = "Species", NArm = FALSE)
  
  # transform to relative abundance
  tp <- transform_sample_counts(phylum, rel_abundance)
  tc <- transform_sample_counts(class, rel_abundance)
  to <- transform_sample_counts(order, rel_abundance)
  tf <- transform_sample_counts(family, rel_abundance)
  tg <- transform_sample_counts(genus, rel_abundance)
  ts <- transform_sample_counts(species, rel_abundance)
  
  # create table for storing data
  p.val <- matrix(data = NA, nrow = 32, ncol = 4)
  colnames(p.val) <- c("Taxa", "p.value.TNF", "p.value.IL17.load", "p.value.IL17.maint")
  
  for (j in seq_along(taxa.names)) {
    
    # for each taxon, create new directory to store results
    dir.taxa = paste(".../IL17.TNF/ITS/jobs/5_specific.taxa_rel.abund.wilcox_R/", taxa.names[j], "/", sep = "")
    if (dir.exists(dir.taxa) == FALSE) {
      dir.create(dir.taxa)
    }
    
    # extract relative abundance of specific taxa
    if (grepl("p__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tp, Phylum == taxa.names[j])
    } else if (grepl("c__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tc, Class == taxa.names[j])
    } else if (grepl("o__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(to, Order == taxa.names[j])
    } else if (grepl("f__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tf, Family == taxa.names[j])
    } else if (grepl("g__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tg, Genus == taxa.names[j])
    } else if (grepl("s__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(ts, Species == taxa.names[j])
    }
    
    # create dataset
    d <- as.data.frame(merge(sample_data(taxa.rel.abund), t(otu_table(taxa.rel.abund)), 
                             by = "row.names", all = TRUE))  
    
    # rename column with relative abundance values for easier plotting 
    # specify the column number for Incertae_sedis as it has multiple OTU IDs; others have only one OTU ID
    if (taxa.names[j] == "f__Incertae_sedis") {
      colnames(d)[44] <- "Taxa_rel_abundance" # multiple OTUs; this is the column that contains the OTU of interest
    } else {
      colnames(d)[ncol(d)] <- "Taxa_rel_abundance" 
    }
    
    ####################
    
    ## TNFi ##
    
    # subset dataset
    d.TNF <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "1_TNF") %>% 
      filter(Timepoint_revised != "D")
    
    write.csv(d.TNF, file = paste(dir.taxa, "TNF.data.csv", sep = "/"))
    
    ####################
    
    ## IL-17i loading ##
    
    # subset dataset
    d.IL17.load <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "2_IL17") %>% 
      filter(Timepoint_revised != "D")
    
    write.csv(d.IL17.load, file = paste(dir.taxa, "IL17.load.data.csv", sep = "/"))
    
    ####################
    
    ## IL-17i maintenance ##
    
    # subset dataset
    d.IL17.maint <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "2_IL17") %>% 
      filter(Timepoint_revised != "C") %>%
      filter(Subject != "cos8" & Subject != "cos11" & Subject != "cos14" & Subject != "cos18")
    
    write.csv(d.IL17.maint, file = paste(dir.taxa, "IL17.maint.data.csv", sep = "/"))
    
    ####################
    
    # calculate wilcoxon between pre-post visits for each subset
    wilc.TNF <- wilcox.test(Taxa_rel_abundance ~ Timepoint_revised, data = d.TNF, paired = TRUE)
    wilc.IL17.load <- wilcox.test(Taxa_rel_abundance ~ Timepoint_revised, data = d.IL17.load, paired = TRUE)
    wilc.IL17.maint <- wilcox.test(Taxa_rel_abundance ~ Timepoint_revised, data = d.IL17.maint, paired = TRUE)
    
    # store stats
    p.val[j, 1] <- taxa.names[j]
    p.val[j, 2] <- wilc.TNF$p.value
    p.val[j, 3] <- wilc.IL17.load$p.value
    p.val[j, 4] <- wilc.IL17.maint$p.value
  }
  
  # save table of stats
  write.csv(p.val, file = ".../IL17.TNF/ITS/jobs/5_specific.taxa_rel.abund.wilcox_R/p.value.table.csv")
}

############################################################################
############################################################################
############################################################################

### P-value table of magnitude relative abundance change between TNfi and IL-17i cohorts ###

# list of phyloseq objects to process
pseqs <- list(physeq = phy_ITS.R1_human_TNF.B.C_IL17.B.C.D)

# list of specific taxa
taxa.names <- c("p__Ascomycota",
                "c__Saccharomycetes", "c__Tremellomycetes",
                "o__Saccharomycetales", "o__Tremellales",
                "f__Saccharomycetaceae", "f__Mycosphaerellaceae", "f__Glomerellaceae",
                "f__Plectosphaerellaceae", "f__Sclerotiniaceae", "f__Mucoraceae",
                "f__Incertae_sedis",
                "g__Saccharomyces", "g__Candida", "g__Penicillium", 
                "g__Debaryomyces", "g__Glomerella", "g__Plectosphaerella",
                "g__Aspergillus",
                "s__Pichia_sp_1_TMS_2011", "s__Penicillium_charlesii", "s__Hypocrea_pseudokoningii",
                "s__Fusarium_sp_4_TMS_2011", "s__Candida_fermenticarens", "s__Trichoderma_longibrachiatum",
                "s__Rhodotorula_sp_RhSoW01", "s__Aspergillus_niger", "s__Saccharomyces_cerevisiae", 
                "s__Candida_albicans", "s__Candida_dubliniensis", "s__Candida_parapsilosis", 
                "s__Candida_tropicalis") 

# for each phyloseq object in pseqs, calculate the relative abudnance
# over time for specific taxa
for (i in seq_along(pseqs)) {
  
  # merge phyloseq objects to each taxnomic level
  phylum <- tax_glom(pseqs[[i]], taxrank = "Phylum", NArm = FALSE)
  class <- tax_glom(pseqs[[i]], taxrank = "Class", NArm = FALSE)
  order <- tax_glom(pseqs[[i]], taxrank = "Order", NArm = FALSE)
  family <- tax_glom(pseqs[[i]], taxrank = "Family", NArm = FALSE)
  genus <- tax_glom(pseqs[[i]], taxrank = "Genus", NArm = FALSE)
  species <- tax_glom(pseqs[[i]], taxrank = "Species", NArm = FALSE)
  
  # transform to relative abundance
  tp <- transform_sample_counts(phylum, rel_abundance)
  tc <- transform_sample_counts(class, rel_abundance)
  to <- transform_sample_counts(order, rel_abundance)
  tf <- transform_sample_counts(family, rel_abundance)
  tg <- transform_sample_counts(genus, rel_abundance)
  ts <- transform_sample_counts(species, rel_abundance)
  
  # create table for storing data
  p.val <- matrix(data = NA, nrow = 32, ncol = 3)
  colnames(p.val) <- c("Taxa", "p.value.load", "p.value.maint")
  
  for (j in seq_along(taxa.names)) {
    
    # extract relative abundance of specific taxa
    if (grepl("p__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tp, Phylum == taxa.names[j])
    } else if (grepl("c__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tc, Class == taxa.names[j])
    } else if (grepl("o__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(to, Order == taxa.names[j])
    } else if (grepl("f__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tf, Family == taxa.names[j])
    } else if (grepl("g__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(tg, Genus == taxa.names[j])
    } else if (grepl("s__", taxa.names[j])) {
      taxa.rel.abund <- subset_taxa(ts, Species == taxa.names[j])
    }
    
    # create dataset
    d <- as.data.frame(merge(sample_data(taxa.rel.abund), t(otu_table(taxa.rel.abund)), 
                             by = "row.names", all = TRUE))  
    
    # rename column with relative abundance values for easier plotting 
    # specify the column number for Incertae_sedis as it has multiple OTU IDs; others have only one OTU ID
    if (taxa.names[j] == "f__Incertae_sedis") {
      colnames(d)[44] <- "Taxa_rel_abundance" # multiple OTUs; this is the column that contains the OTU of interest
    } else {
      colnames(d)[ncol(d)] <- "Taxa_rel_abundance" 
    }
    
    ####################
    
    ## TNFi ##
    
    # subset dataset
    # spread data across timepoints
    # calculate absolute delta relative abundance
    d.TNF <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "1_TNF") %>% 
      filter(Timepoint_revised != "D") %>% 
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      mutate(Delta_rel_abundance = abs(C-B))
    
    # rename columns
    colnames(d.TNF)[colnames(d.TNF)=="B"] <- "1_pre"  
    colnames(d.TNF)[colnames(d.TNF)=="C"] <- "2_post"  
    
    ####################
    
    ## IL-17i loading ##
    
    # subset dataset
    # spread data across timepoints
    # calculate absolute delta relative abundance
    d.IL17.load <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "2_IL17") %>% 
      filter(Timepoint_revised != "D") %>%
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      mutate(Delta_rel_abundance = abs(C-B))
    
    # rename treatment
    d.IL17.load$Treatment <- replace(as.character(d.IL17.load$Treatment), 
                                     as.character(d.IL17.load$Treatment) == "2_IL17", "2_IL17.load")
    # rename columns
    colnames(d.IL17.load)[colnames(d.IL17.load)=="B"] <- "1_pre"  
    colnames(d.IL17.load)[colnames(d.IL17.load)=="C"] <- "2_post"
    
    ####################
    
    ## IL-17i maintenance ##
    
    # subset dataset
    # spread data across timepoints
    # calculate absolute delta relative abundance
    d.IL17.maint <- d %>%
      subset(select = c("Subject", "Treatment", "Timepoint_revised", "Taxa_rel_abundance")) %>%
      filter(Treatment == "2_IL17") %>% 
      filter(Timepoint_revised != "C") %>%
      spread(key = "Timepoint_revised", value = "Taxa_rel_abundance") %>%
      filter(!is.na(`D`)) %>%
      mutate(Delta_rel_abundance = abs(D-B))
    
    # rename treatment
    d.IL17.maint$Treatment <- replace(as.character(d.IL17.maint$Treatment), as.character(d.IL17.maint$Treatment) == "2_IL17", "3_IL17.maint")
    
    # rename columns
    colnames(d.IL17.maint)[colnames(d.IL17.maint)=="B"] <- "1_pre"  
    colnames(d.IL17.maint)[colnames(d.IL17.maint)=="D"] <- "2_post"
    
    ####################
    
    # join subsets
    d.final <- rbind(d.TNF, d.IL17.load, d.IL17.maint)
    
    # split IL-17i into loading and maintenance subcategories
    d.load <- d.final[d.final$Treatment != "3_IL17.maint", ]
    d.maint <- d.final[d.final$Treatment != "2_IL17.load", ]
    
    # mann-whitney
    mw.load <- wilcox.test(Delta_rel_abundance ~ Treatment, data = d.load, paired = FALSE)
    mw.maint <- wilcox.test(Delta_rel_abundance ~ Treatment, data = d.maint, paired = FALSE)
    
    # store stats
    p.val[j, 1] <- taxa.names[j]
    p.val[j, 2] <- mw.load$p.value
    p.val[j, 3] <- mw.maint$p.value
  }
  
  # save table of stats
  write.csv(p.val, file = ".../IL17.TNF/ITS/jobs/5_specific.taxa_delta.rel.abund_R/p.value.table.csv")
}

