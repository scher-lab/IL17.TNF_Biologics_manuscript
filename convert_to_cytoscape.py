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
