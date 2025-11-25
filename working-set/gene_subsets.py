# The purpose of this program is to sort the data into meaningful sets: 
# - Significant genes only belonging to heat stress
# - Significant genes only belonging to salt stress
# - Significant genes that are common between both stresses

import pandas as pd

# Load filtered CSVs
sig_salt = pd.read_csv("CSV_datasets/salt_significant_genes.csv")
sig_heat = pd.read_csv("CSV_datasets/heat_significant_genes.csv")

# Merge on gene_id
merged = pd.merge(sig_salt, sig_heat, on="gene_id", how="outer", suffixes=("_salt", "_heat"))

# Genes significant in both conditions
both = merged.dropna(subset=["log2FoldChange_salt", "log2FoldChange_heat"])

# Genes unique to salt
only_salt = merged[merged["log2FoldChange_heat"].isna()]

# Genes unique to heat
only_heat = merged[merged["log2FoldChange_salt"].isna()]

both.to_csv("sorted_gene_subsets/both_stress_genes.csv", index=False)
only_salt.to_csv("sorted_gene_subsets/only_salt_genes.csv", index=False)
only_heat.to_csv("sorted_gene_subsets/only_heat_genes.csv", index=False)
