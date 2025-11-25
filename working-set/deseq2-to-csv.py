# The purpose of this program is to process the DESeq2 datasets from our data processing step, filter the data based on significance, 
# and output the result in CSV format.

import pandas as pd

## Process salt DESeq2 results
salt = pd.read_csv("DESeq2_datasets/salt.tabular", sep="\t", header=None)
salt.columns = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# Filter for significant genes
sig_salt = salt[(salt['padj'] < 0.05) & (abs(salt['log2FoldChange']) >= 1)]

# Label up/downregulated genes
sig_salt.loc[:, 'regulation'] = sig_salt['log2FoldChange'].apply(lambda x: 'up' if x >= 1 else 'down')

# Save filtered table
sig_salt.to_csv("CSV_datasets/salt_significant_genes.csv", index=False)


## Process heat DESeq2 results
heat = pd.read_csv("DESeq2_datasets/heat.tabular", sep="\t", header=None)
heat.columns = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# Filter for significant genes
sig_heat = heat[(heat['padj'] < 0.05) & (abs(heat['log2FoldChange']) >= 1)]

# Label up/downregulated genes
sig_heat.loc[:, 'regulation'] = sig_heat['log2FoldChange'].apply(lambda x: 'up' if x >= 1 else 'down')

# Save filtered table
sig_heat.to_csv("CSV_datasets/heat_significant_genes.csv", index=False)