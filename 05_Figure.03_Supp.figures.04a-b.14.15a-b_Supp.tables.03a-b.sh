#############################################
## Bash script                             ##
## Project: IL17.TNF Biologics manuscript  ##
## Figure 03                               ##
## Supplementary figures 04a-b, 14, 15a-b  ##
## Supplementary tables 03a-b              ##
## 16S data                                ##
#############################################

### Brief description:

### Figure 03:
### Panel A: Bacterial co-occurrence network pre-IL-17i treatment
### Panel B: Bacterial co-occurrence network post-IL-17i treatment (maintenance visit)

### Supplementary figure 04a: Bacterial co-occurrence network pre-TNFi treatment
### Supplementary figure 04b: Bacterial co-occurrence network post-TNFi treatment (maintenance visit)

### Supplementary figure 14:  Expanded bacterial co-occurrence network post-TNFi treatment (maintenance visit)

### Supplementary figure 15a: Expanded bacterial co-occurrence network pre-IL-17i treatment
### Supplementary figure 15b: Expanded bacterial co-occurrence network post-IL-17i treatment (maintenance visit)

### Supplementary table 03a: TNFi bacterial co-occurrence network nodes
### Supplementary table 03b: IL-17i bacterial co-occurrence network nodes

### Below is a script that runs SparCC on the relevant OTU tables. 
### SparCC results are imported into Cytoscape to generage networks.
### Inkscape is used to edit the network graphs. 

#/bin/bash
conda activate sparcc

# Install SparCC in a python 2.6 environment

# These steps were repeated for each of the otu tables
# The tables were converted from biom to tsv format prior to running sparcc

OTU_TABLE=/some/data/table.txt

# Run SparCC on the original data
SparCC.py $OTU_TABLE -i 20 --cor_file="otu_table_cor.out"

# Generate 100 shuffled datasets
MakeBootstraps.py $OTU_TABLE -n 100 -t permutation_#.txt -p pvals/

# Run SparCC on each of the permutations in parallel
for i in {0..99}; do
    SparCC.py pvals/permutation_${i}.txt -i 20 --cor_file="pvals/perm_cor_${i}.txt"&
done

# Wait for them all to finish
wait;

# Compute the one and two sided psuedo P values
PseudoPvals.py otu_table_cor.out pvals/perm_cor_#.txt 5 -o pvals.one_sided.txt -t one_sided
PseudoPvals.py otu_table_cor.out pvals/perm_cor_#.txt 5 -o pvals.two_sided.txt -t two_sided

# Convert the files into a format that Cytoscape can import (this script requires pandas)
# Python script provided as a separate file
python convert_to_cytoscape.py pvals.one_sided.txt pvals.one_sided_cyto.txt

# Process input file to shorten taxa names as necessary

### Cytoscape
# File -> Import -> Network from file
# Edge-weighted sping embeded layout
# Small node size
# Dashed lines for negative correlations
# Export the plots as PDFs

# Additional Processing
# Futher modification was done using the program Inkscape on the exported PDF
# Changes included increasing the font size and darkening the edges connecting nodes
