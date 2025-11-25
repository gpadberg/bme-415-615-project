# This program outputs a venn diagram to visualize the counts of significant data between the two stresses: heat and salt.

import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Load filtered significant genes
sig_salt = pd.read_csv("CSV_datasets/salt_significant_genes.csv")
sig_heat = pd.read_csv("CSV_datasets/heat_significant_genes.csv")

# Get sets of gene IDs
salt_genes = set(sig_salt['gene_id'])
heat_genes = set(sig_heat['gene_id'])

# Create Venn diagram
plt.figure(figsize=(6,6))
venn2([salt_genes, heat_genes], set_labels=("Salt Stress", "Heat Stress"))
plt.title("Overlap of Significant Genes in Salt and Heat Stress")
plt.show()
