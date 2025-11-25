# This program creates a CSV file that contains the union of the results from each stress

import pandas as pd

# Load filtered CSVs
sig_salt = pd.read_csv("CSV_datasets/salt_significant_genes.csv")
sig_heat = pd.read_csv("CSV_datasets/heat_significant_genes.csv")

# Merge on gene_id
merged = pd.merge(sig_salt, sig_heat, on="gene_id", how="outer", suffixes=("_salt", "_heat"))

# Save merged file
merged.to_csv("merged_dataset/merged_salt_heat_genes.csv", index=False)