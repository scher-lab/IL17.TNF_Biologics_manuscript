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

# Convert the files into a format SparCC can import (this script requires pandas)
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
