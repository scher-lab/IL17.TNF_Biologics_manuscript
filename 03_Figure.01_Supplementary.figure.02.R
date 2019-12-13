############################################ 
## R script                               ##
## Project: IL17.TNF_Biologics_manuscript ##
## Figure 01, Supplementary Figure 02     ##
## 16S data                               ##
## Author: JM                             ##
############################################

### Brief description:
### This script covers the code for Figure 01 and Supplementary Figure 02.

### Figure 01:
### Panel A: Alpha diversity, shannon
### Panel B: Beta diversity pcoa, bray-curtis
### Panel C: Beta diversity, distance variability between TNFi and IL-17i
### Panel D: Taxa summary

### Supplementary Figure 02:
### Panel A: Alpha diversity, observed OTUs
### Panel B: Alpha diversity, simpson

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(vegan)
library(ade4)
library(PMCMR)
library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(data.table)
library(scales)

############################################################################
############################################################################
############################################################################

### Statistics functions ###

# All plot statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
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

# Boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

# Whiskers statistics: median, min value, max value
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

############################################################################
############################################################################
############################################################################

### Alpha diversity boxplots ###

## Setup ##

# background theme
bkg <- theme_few() +
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 14, color = "black")) +
  theme(axis.title.y = element_text(size = 20, face = "bold", color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# directory for storing files
dir = ".../IL17.TNF/16S/jobs/1_alpha.div_boxplots_R/"

# list of phyloseq objects to process
pseqs <- list(physeq = phy_16S_human_TNF.B.C_IL17.B.C.D_even1000)

# list of phyloseq names
pseq.names <- c("16S_human_TNF.B.C_IL17.B.C.D_even1000")

# list of diversity calculations
divs <- c("Observed", "Shannon", "Simpson")

##############

## Plotting ##

# for each phyloseq object in pseq, calculate diversity indices listed in divs
for (i in seq_along(pseqs)) {
  
  # create new directory to store results
  dir.new = paste(dir, pseq.names[i], "/", sep = "")
  dir.create(dir.new)
  
  for (j in seq_along(divs)) {
    
    # create filenames
    filename_table.div = paste(pseq.names[i], "adiv", divs[j], "table_diversity.csv", sep = "_")
    filename_table.for.plot = paste(pseq.names[i], "adiv", divs[j], "table_for.plot.csv", sep = "_")
    filename_plot = paste(pseq.names[i], "adiv", divs[j], "plot.pdf", sep = "_")  
    filename_plot.stats = paste(pseq.names[i], "adiv", divs[j], "plot.stats.csv", sep = "_")
    filename_stats = paste(pseq.names[i], "adiv", divs[j], "stats.txt", sep = "_")
    
    # remove empty OTUs
    phylo <- subset_taxa(pseqs[[i]], rowSums(otu_table(pseqs[[i]])) > 0)
    
    # calculate diversity
    rich <- estimate_richness(phylo, split = TRUE, measures = divs[j])
    
    # remove "X" from sample names, which is added to numeric IDs (i.e. TNF samples)
    rownames(rich) <- gsub('^X','',rownames(rich))
    
    # save diversity table
    ftd = paste(dir.new, filename_table.div, sep = "")
    write.csv(rich, file = ftd)
    
    # create dataset from which to generate plot by merging mapping file with calculated diversity
    d <- as.data.frame(merge(sample_data(phylo), rich, by = "row.names", all = TRUE))
    
    # rename diversity calculation column for easier plotting 
    colnames(d)[colnames(d) == divs[j]] <- "Diversity" 
    
    # save dataset
    ftp = paste(dir.new, filename_table.for.plot, sep = "")
    write.csv(d, file = ftp)
    
    # choose colors for plot
    col1 <- c("#ffb400", "#c40018", "#53d397", "#448ef6", "#8559a5") # TNF.B.C_IL17.B.C.D
    
    # plot and save alpha diversity
    p <- ggplot(data = d, aes(x = Pre_post_treatment_det, y = Diversity, fill = Pre_post_treatment_det)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 0.8, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                   color = "black", size = 0.5, width = 0.5) +
      geom_jitter(width = 0.1, size = 1.5) +
      scale_x_discrete(labels = c("TNFi\npre", "TNFi\nmaint", "IL-17i\npre", "IL-17i\nload", "IL-17i\nmaint")) +
      scale_fill_manual(values = col1) +      
      xlab(NULL) +
      bkg
    
    if(divs[j] == "Observed") {
      p <- p + ylab("Observed OTUs") + expand_limits(y = c(0, max(d$Diversity)+50))
    } else if(divs[j] == "Shannon") {
      p <- p + ylab("Shannon index") + expand_limits(y = c(0, max(d$Diversity)+0.5))
    } else if(divs[j] == "Simpson") {
      p <- p + ylab("Simpson index") + expand_limits(y = c(0, max(d$Diversity)+0.2))
    }
    
    fp = paste(dir.new, filename_plot, sep = "")
    pdf(file = fp)
    plot(p)
    dev.off()
    
    # obtain plot statistics and save
    s <- aggregate(Diversity ~ Pre_post_treatment_det, data = d, stats.all)
    s <- t(s)
    rownames(s) <- c("tx.timepoint", "mean", "sd", "median", "min", "max", 
                     "10%ile", "25%ile", "75%ile", "90%ile")
    
    fps = paste(dir.new, filename_plot.stats, sep = "")
    write.csv(s, file = fps)
    
    # create file to save general statistics
    fs = paste(dir.new, filename_stats, sep = "")
    cat("Statistical Tests:\n\n", file = fs)
    cat("Alpha diversity - ", file = fs, append = TRUE)
    cat(divs[j], file = fs, append = TRUE) 
    cat("\n", file = fs, append = TRUE)
    cat(pseq.names[i], file = fs, append = TRUE)
    cat("\n\n==========================================\n\n", file = fs, append = TRUE)
    
    # calculate statistics
    # use wilcoxon paired non-parametric test w/o accounting for multiple comparisons testing
    
    # subset TNF and IL17
    d.TNF <- d[d$Treatment == "1_TNF", ]
    d.IL17 <- d[(d$Treatment == "2_IL17"), ]
    
    # wilcoxon for TNF
    wilc.TNF <- wilcox.test(Diversity ~ Pre_post_treatment_det, data = d.TNF, paired = TRUE)
    
    # save calculations
    cat("Wilcoxon test on human main cohort - TNF:\n", file = fs, append = TRUE)
    capture.output(wilc.TNF, file = fs, append = TRUE)
    cat("\n==========================================\n\n", file = fs, append = TRUE)
    
    # subset loading
    d.load <- d.IL17[d.IL17$Pre_post_det != "3_maint", ]
    
    # subset maintenance
    # need to trim baseline visit samples to match maintenance samples (not equal in number)
    d.IL17.visit.D_subject <- as.vector(d.IL17[(d.IL17$Timepoint_revised == "D"), ]$Subject)
    d.maint <- d.IL17[d.IL17$Subject %in% d.IL17.visit.D_subject, ]
    d.maint <- d.maint[d.maint$Pre_post_det != "2_load", ]
    
    # wilcoxon for IL17 loading and maintenance
    wilc.load <- wilcox.test(Diversity ~ Pre_post_treatment_det, data = d.load, paired = TRUE)
    wilc.maint <- wilcox.test(Diversity ~ Pre_post_treatment_det, data = d.maint, paired = TRUE)
    
    cat("Wilcoxon test on human main cohort - IL17 loading:\n", file = fs, append = TRUE)
    capture.output(wilc.load, file = fs, append = TRUE)
    cat("\n==========================================\n\n", file = fs, append = TRUE)
    cat("Wilcoxon test on human main cohort - IL17 maintenance:\n", file = fs, append = TRUE)
    capture.output(wilc.maint, file = fs, append = TRUE)
  }
}

############################################################################
############################################################################
############################################################################

### Beta diversity pcoa ###

## Setup ##

# colors
col1 <- c("#ffb400", "#c40018", "#53d397", "#448ef6", "#8559a5")

# background theme
bkg <- theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 12, color = "black", face = "bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.justification = "top") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))
  
# function to specify that axis labels have 2 decimal places
f.dec <- function(x){
  format(round(x, 2), nsmall = 2)
}

# directory for storing files
dir = ".../IL17.TNF/16S/jobs/2_beta.div_pcoa_R/"

# list of phyloseq objects
pseqs <- list(physeq = phy_16S_human_TNF.B.C_IL17.B.C.D_even1000)

# list of phyloseq names
pseq.names <- c("16S_human_TNF.B.C_IL17.B.C.D_even1000")

# list of distance methods
dists <- c("bray")

##############

## Plotting ##

# for each phyloseq object in pseqs, perform ordination using distance methods in dists
for (i in seq_along(pseqs)) {
  
  # create new directory to store results
  dir.new = paste(dir, pseq.names[i], "/", sep = "")
  dir.create(dir.new)
  
  for (j in seq_along(dists)) {
    
    # create filenames
    filename_matrix = paste(pseq.names[i], "bdiv", dists[j], "matrix.csv", sep = "_")  
    filename_plot_pre.post.tx = paste(pseq.names[i], "bdiv", dists[j], "plot_pre.post.tx.pdf", sep = "_")
    filename_plot_msq = paste(pseq.names[i], "bdiv", dists[j], "plot_msq.pdf", sep = "_")  
    filename_stats = paste(pseq.names[i], "bdiv", dists[j], "stats.txt", sep = "_")  
    
    # remove empty OTUs
    phylo <- subset_taxa(pseqs[[i]], rowSums(otu_table(pseqs[[i]])) > 0)
    
    # calculate distance matrix
    d <- distance(phylo, method = dists[j])
      
    # perform ordination
    ord <- ordinate(phylo, method = "PCoA", distance = d)    
    
    # save distance matrix
    fm = paste(dir.new, filename_matrix, sep = "")
    write.csv(data.matrix(d), file = fm)
    
    # choose colors
    my.cols = col1 
    
    # plot beta diversity
    p <- plot_ordination(phylo, ordination = ord, color = "Pre_post_treatment_det") +
      geom_vline(xintercept = 0, color = "#919190", size = 0.5) +
      geom_hline(yintercept = 0, color = "#919190", size = 0.5) +
      geom_point(size = 3) +
      scale_color_manual(values = my.cols, labels = c("TNFi pre", "TNFi maint", "IL-17i pre", "IL-17i load", "IL-17i maint")) +
      bkg +
      scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
      scale_y_continuous(labels = f.dec)   # 2 decimal places on y-axis
    
    # save plot
    fp = paste(dir.new, filename_plot_pre.post.tx, sep = "")
    pdf(file = fp)
    plot(p)
    dev.off()
  } 
}

############################################################################
############################################################################
############################################################################

### Beta diversity distance variability pre/post treatment ###

## Setup ##

# colors
col1 <- c("#c40018", "#448ef6", "#8559a5")

# shapes
shape1 <- c(16,15,17)

# background themes
bkg <- theme_classic() +
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 14, color = "black")) +
  theme(axis.title.y = element_text(size = 20, face = "bold", color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.position = "none")

bkg_qqplot <- theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 11, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8)) +
  theme(legend.title = element_text(size = 10)) +
  theme(plot.title = element_text(size = 12, color = "black", face = "bold"))

# directory for storing files
dir = ".../IL17.TNF/16S/jobs/2_beta.div_pre.post.tx_R/"

# list of phyloseq objects
pseqs <- list(physeq = phy_16S_human_TNF.B.C_IL17.B.C.D,
              physeq = phy_16S_human_TNF.B.C_IL17.B.C.D_even1000)

# list of phyloseq names
pseq.names <- c("16S_human_TNF.B.C_IL17.B.C.D",
                "16S_human_TNF.B.C_IL17.B.C.D_even1000")

# list of distance methods
dists <- c("bray")

##############

## Plotting ##

# for each phyloseq object in pseqs, calculate all unique pre-post distances 
# and compare distance values in TNF vs IL17 subjects
for (i in seq_along(pseqs)) {
  
  # create new directory to store results
  dir.new = paste(dir, pseq.names[i], "/", sep = "")
  dir.create(dir.new)
  
  for (j in seq_along(dists)) {
    
    # create filenames
    filename_matrix = paste(pseq.names[i], "bdiv_var.pre.post.tx", dists[j], "matrix.csv", sep = "_")  
    filename_dist = paste(pseq.names[i], "bdiv_var.pre.post.tx", dists[j], "dist.table.csv", sep = "_") 
    filename_plot = paste(pseq.names[i], "bdiv_var.pre.post.tx", dists[j], "plot.pdf", sep = "_")  
    filename_plot.stats = paste(pseq.names[i], "bdiv_var.pre.post.tx", dists[j], "plot.stats.csv", sep = "_")
    filename_stats_TNF.IL17.load = paste(pseq.names[i], "bdiv_var.pre.post.tx", dists[j], "stats_TNF.IL17.load.txt", sep = "_")
    filename_stats_TNF.IL17.maint = paste(pseq.names[i], "bdiv_var.pre.post.tx", dists[j], "stats_TNF.IL17.maint.txt", sep = "_")
    
    # remove empty OTUs
    phylo <- subset_taxa(pseqs[[i]], rowSums(otu_table(pseqs[[i]])) > 0)
    
    d <- distance(phylo, method = dists[j])    
    
    dm <- data.matrix(d, rownames.force = TRUE)
    
    # save distance matrix
    fm = paste(dir.new, filename_matrix, sep = "")
    write.csv(dm, file = fm)
    
    # store row and column names of distance matrix
    row = rownames(dm)
    col = colnames(dm)
    
    # store indices of upper triangle of distance matrix excluding diagnoal 
    # only upper triangle is used as matrix is symmetrical
    ind <- which(upper.tri(dm, diag = FALSE), arr.ind = TRUE)
    
    # designate boxplot categories
    cat <- c("1_TNF", "2_IL17_load", "3_IL17_maint")
    
    # create lists for storing sample names, treatment timepoints, distance values, and category names
    # each list contains a separate vector for each boxplot category
    sample.names.1 <- vector(mode = "list", length = length(cat))
    names(sample.names.1) <- c(paste(cat, sep=","))
    
    sample.names.2 <- vector(mode = "list", length = length(cat))
    names(sample.names.2) <- c(paste(cat, sep=","))
    
    tx.timepoint.1 <- vector(mode = "list", length = length(cat))
    names(tx.timepoint.1) <- c(paste(cat, sep=","))
    
    tx.timepoint.2 <- vector(mode = "list", length = length(cat))
    names(tx.timepoint.2) <- c(paste(cat, sep=","))
    
    distance.values <- vector(mode = "list", length = length(cat))
    names(distance.values) <- c(paste(cat, sep=","))
    
    category <- vector(mode = "list", length = length(cat))
    names(category) <- c(paste(cat, sep=","))
    
    for (k in 1:length(cat)) {
      sample.names.1[[k]] <- character(length(ind))
      sample.names.2[[k]] <- character(length(ind))
      tx.timepoint.1[[k]] <- character(length(ind))
      tx.timepoint.2[[k]] <- character(length(ind))
      distance.values[[k]] <- rep(-1, length(ind))
      category[[k]] <- character(length(ind))
    }
    
    # create counter
    c = 1
    
    # loop through indices of upper triangle of distance matrix 
    # store appropriate unique distance values according to criteria below
    for (q in 1:nrow(ind)) {
      
      # indices of row and column in distance matrix
      index.row <- ind[q,1]
      index.col <- ind[q,2]
      
      # process row sample
      row.tx <- toString(get_variable(phylo)[row[index.row], "Treatment"]) 
      row.subject <- toString(get_variable(phylo)[row[index.row], "Subject"])
      row.pre.post <- toString(get_variable(phylo)[row[index.row], "Pre_post"])
      row.tx.timepoint <- toString(get_variable(phylo)[row[index.row], "Treatment_timepoint_revised"])
        
      # process column sample
      col.tx <- toString(get_variable(phylo)[col[index.col], "Treatment"]) 
      col.subject <- toString(get_variable(phylo)[col[index.col], "Subject"])
      col.pre.post <- toString(get_variable(phylo)[col[index.col], "Pre_post"])
      col.tx.timepoint <- toString(get_variable(phylo)[col[index.col], "Treatment_timepoint_revised"])
      
      # compare row and column samples
      
      # if row and column samples are in the same treatment category
      if (row.tx == col.tx) {
        
        # and if row and column subjects are identical
        if (row.subject == col.subject) {
          
          # and if row visit is "pre" and column visit is "post"
          if (row.pre.post == "1_pre" & col.pre.post == "2_post") {
            
            # assign distance values to appropriate boxplot category
            if (row.tx.timepoint == "1_TNF.B" & col.tx.timepoint == "1_TNF.C") {
              graph.cat = "1_TNF"
            } else if (row.tx.timepoint == "2_IL17.B" & col.tx.timepoint == "2_IL17.C"){
              graph.cat = "2_IL17_load"
            } else if (row.tx.timepoint == "2_IL17.B" & col.tx.timepoint == "2_IL17.D") {
              graph.cat = "3_IL17_maint"
            }
            
            # assign names of samples, treatment-timepoints, boxplot category, and distance value
            sample.names.1[[graph.cat]][c] <- row[index.row]
            sample.names.2[[graph.cat]][c] <- col[index.col]
            tx.timepoint.1[[graph.cat]][c] <- row.tx.timepoint
            tx.timepoint.2[[graph.cat]][c] <- col.tx.timepoint
            distance.values[[graph.cat]][c] <- dm[index.row, index.col]
            category[[graph.cat]][c] <- graph.cat
            
            # update counter
            c = c + 1
          }
        }
      }
    } 
    
    # combine vectors from each list
    comb_sample.names.1 <- stack(sample.names.1)
    comb_sample.names.2 <- stack(sample.names.2)
    comb_tx.timepoint.1 <- stack(tx.timepoint.1)
    comb_tx.timepoint.2 <- stack(tx.timepoint.2)
    comb_distance.values <- stack(distance.values)
    comb_category <- stack(category)
    
    # create combined data frame
    comb.data <- data.frame(Sample1 = comb_sample.names.1$values, TxTimepoint1 = comb_tx.timepoint.1$values, 
                            Sample2 = comb_sample.names.2$values, TxTimepoint2 = comb_tx.timepoint.2$values,
                            Distance = as.numeric(comb_distance.values$values),
                            Category = comb_category$values)
    
    # remove empty values
    comb.data.s <- comb.data[(comb.data$Distance != -1), ]
    
    # save distance data
    fd = paste(dir.new, filename_dist, sep = "")
    write.csv(comb.data.s, file = fd)
    
    # plot and save
    p <- ggplot(data = comb.data.s, aes(x = Category, y = Distance, color = Category, shape = Category)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 1, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                   color = "black", size = 0.5, width = 0.5, fill = "white") +
      geom_jitter(width = 0.15, size = 3) +
      scale_color_manual(values = col1) +
      scale_shape_manual(values = shape1) +
      scale_x_discrete(labels = c("TNFi\npre-maint", "IL-17i\npre-load", "IL-17i\npre-maint")) +
      expand_limits(y = c(0, 1)) +
      xlab(NULL) + ylab("Bray-curtis distance") +
      bkg
    
    fp = paste(dir.new, filename_plot, sep = "")
    pdf(file = fp) 
    plot(p)
    dev.off() 
    
    # obtain plot statistics and save
    s <- aggregate(Distance ~ Category, data = comb.data.s, stats.all)
    s <- t(s)
    rownames(s) <- c("category", "mean", "sd", "median", "min", "max", 
                     "10%ile", "25%ile", "75%ile", "90%ile")
    
    fps = paste(dir.new, filename_plot.stats, sep = "")
    write.csv(s, file = fps)
    
    # calculate general statistics
    
    # keep TNF but subset IL17 into loading and maintenance subcategories
    # comparison is done between TNF and IL17
    d.IL17.load <- comb.data.s[comb.data.s$Category != "3_IL17_maint", ]
    d.IL17.maint <- comb.data.s[comb.data.s$Category != "2_IL17_load", ]
    
    # TNF-IL17 loading #
    
    # Mann-Whitney
    mw.IL17.load <- wilcox.test(Distance ~ Category, data = d.IL17.load, paired = FALSE)
    
    # save calculations
    fs.IL17.load = paste(dir.new, filename_stats_TNF.IL17.load, sep = "")
    cat("Statistical Tests - TNF-IL17 loading:\n\n", file = fs.IL17.load)
    cat("Beta diversity variability pre-post treatment\n", file = fs.IL17.load, append = TRUE)
    cat(dists[j], file = fs.IL17.load, append = TRUE) 
    cat("\n\n", file = fs.IL17.load, append = TRUE)
    cat(pseq.names[i], file = fs.IL17.load, append = TRUE)
    cat("\n\n==========================================\n\n", file = fs.IL17.load, append = TRUE)
    cat("Mann-Whitney test (non-parametric)\n", file = fs.IL17.load, append = TRUE)
    capture.output(mw.IL17.load, file = fs.IL17.load, append = TRUE)
    
    # TNF-IL17 maintenance #
    
    # Mann-Whitney
    mw.IL17.maint <- wilcox.test(Distance ~ Category, data = d.IL17.maint, paired = FALSE)
    
    # save calculations
    fs.IL17.maint = paste(dir.new, filename_stats_TNF.IL17.maint, sep = "")
    cat("Statistical Tests - IL17 maintenance:\n\n", file = fs.IL17.maint)
    cat("Beta diversity variability with treatment vs between subjects\n", file = fs.IL17.maint, append = TRUE)
    cat(dists[j], file = fs.IL17.maint, append = TRUE) 
    cat("\n\n", file = fs.IL17.maint, append = TRUE)
    cat(pseq.names[i], file = fs.IL17.maint, append = TRUE)
    cat("\n\n==========================================\n\n", file = fs.IL17.maint, append = TRUE)
    cat("Mann-Whitney test (non-parametric)\n", file = fs.IL17.maint, append = TRUE)
    capture.output(mw.IL17.maint, file = fs.IL17.maint, append = TRUE)
  }
} 

############################################################################
############################################################################
############################################################################

### Summarize taxa ###

## Setup ##

# colors
col <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
         "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
         "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
         "#771122", "#AA4455", "#DD7788", "#21e445", "#a43eef", "#26e27f", 
         "#b465ff", "#6d58f0", "#01c66d", "#f9016a", "#00d8b4", "#e50133",
         "#c23897", "#93d95c", "#1c6fd8", "#c093ff", "#2e821a", "#ef97ff",
         "#f19d00", "#7e9fff", "#e8c348", "#845fb9", "#ff87d1", "#008d69",
         "#ff667c", "#0192d3", "#ffe47e", "#1475be", "#ff864c", "#0088b9",
         "#ff6c56", "#aebbff", "#f3b35b", "#d79fd3", "#ff81b1", "#9bd49c",
         "#c83d67", "#a2588b", "#ff8da0", "#b5626e")

# background
bkg <- theme_classic() +
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 14, color = "black")) +
  theme(axis.title.y = element_text(size = 20, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 10, face = "bold", color = "black")) +
  theme(legend.spacing.x = unit(0.3, "cm")) +
  theme(legend.title = element_blank()) +
  theme(legend.justification = "top") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# directory for storing files
dir = ".../IL17.TNF/16S/jobs/3_taxa.summary_R/"

# list of phyloseq objects to process
pseqs <- list(physeq = phy_16S_human_TNF.B.C_IL17.B.C.D)

# list of phyloseq names
pseq.names <- c("16S_human_TNF.B.C_IL17.B.C.D")

# taxanomic rank
tax.rank <- c("Order")

##############

## Plotting ##

# for each phyloseq object in pseqs, perform taxa summary by treatment-visit
for (i in seq_along(pseqs)) {
  
  # create new directory to store results
  dd = paste(pseq.names[i], "pre.post.tx.det", sep = "_")
  dir.new = paste(dir, dd, "/", sep = "")
  dir.create(dir.new)

  # create filenames
  filename_table = paste(pseq.names[i], "taxa.summary", "pre.post.tx.det", "Order", "table.csv", sep = "_")    
  filename_plot = paste(pseq.names[i], "taxa.summary", "pre.post.tx.det", "Order", "plot.pdf", sep = "_")    
    
  # remove empty OTUs
  phylo <- subset_taxa(pseqs[[i]], rowSums(otu_table(pseqs[[i]])) > 0)
    
  # combine OTUs by taxanomic rank
  t <- tax_glom(phylo, taxrank = tax.rank)
    
  # merge samples according to mapping file category
  m <- merge_samples(t, group = "Pre_post_treatment_det")
    
  # transform to relative abundance
  tm <- transform_sample_counts(m, rel_abundance)
    
  # create table of relative abundance values and save
  r <- merge(tax_table(tm)[,1:4], t(otu_table(tm)), by = "row.names", all = TRUE) # merge tax table and rel abundance counts from OTU table
  r$rel.abund.sum <- rowSums(r[,6:10]) # sum rel abundance counts across treatment visits for each taxon
  r.sorted <- r[order(-r$rel.abund.sum),] # sort in descending order
  ft = paste(dir.new, filename_table, sep = "") # save
  write.csv(r.sorted, file = ft)
    
  # choose top taxa in alphabetical order for plot legend
  # criteria: sum across all categories > 0.01
  top.taxa <- sort(r.sorted$Order[1:13])
  
  # choose color palatte based on total number of taxa
  c <- sample(col, length(r.sorted$Order))
    
  # create plot
  p <- plot_bar(tm, x = "Sample", fill = "Order") + 
      scale_fill_manual(breaks = top.taxa, values = c) + # only display top taxa in legend
      scale_x_discrete(labels = c("TNFi\npre", "TNFi\nmaint", "IL-17i\npre", "IL-17i\nload", "IL-17i\nmaint")) +
      #scale_y_continuous(labels = percent) +
      xlab("") + ylab("Relative abundance") +
      bkg 
    
  # save plot      
  fp = paste(dir.new, filename_plot, sep = "")
  pdf(file = fp, width = 8)
  plot(p)
  dev.off()
}


