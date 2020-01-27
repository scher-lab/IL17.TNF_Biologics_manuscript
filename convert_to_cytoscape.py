#############################################
## Python script                           ##
## Project: IL17.TNF Biologics manuscript  ##
## Figure 03                               ##
## Supplementary figures 04a-b, 14, 15a-b  ##
## Supplementary tables 03a-b              ##
## 16S data                                ##
#############################################

### Brief description:

### Figure 03
### Panel A: Bacterial co-occurrence network pre-IL-17i treatment
### Panel B: Bacterial co-occurrence network post-IL-17i treatment (maintenance visit)

### Supplementary figure 04a: Bacterial co-occurrence network pre-TNFi treatment
### Supplementary figure 04b: Bacterial co-occurrence network post-TNFi treatment (maintenance visit)

### Supplementary figure 14:  Expanded bacterial co-occurrence network post-TNFi treatment (maintenance visit)

### Supplementary figure 15a: Expanded bacterial co-occurrence network pre-IL-17i treatment
### Supplementary figure 15b: Expanded bacterial co-occurrence network post-IL-17i treatment (maintenance visit)


#!/usr/env python

import pandas as pd
from sys import argv

# Read in the passed path
df = pd.read_csv(argv[1], sep='\t')

# Convert it to long format and rename
df = df.melt(id_vars='OTU_id')
df = df.rename(columns={'OTU_id': 'source_node',
                        'variable': 'target_node',
                        'value': 'edge'})
# Write to a new file
df.to_csv(argv[2], sep='\t', index=False)
