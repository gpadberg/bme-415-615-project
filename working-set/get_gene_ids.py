# The purpose of this program is to extract the gene IDs from the CSV files containing the processed and filtered data.

import pandas as pd

# Salt only
df = pd.read_csv("final_gene_subsets/only_salt_genes.csv")
genes = "\n".join(df['gene_id'].astype(str))
with open("significant_gene_ids/salt_gene_ids.txt", "w") as f:
    f.write(genes)

# Heat only
df = pd.read_csv("final_gene_subsets/only_heat_genes.csv")
genes = "\n".join(df['gene_id'].astype(str))
with open("significant_gene_ids/heat_gene_ids.txt", "w") as f:
    f.write(genes)

# Intersection of both stresses
df = pd.read_csv("final_gene_subsets/both_stress_genes.csv")
genes = "\n".join(df['gene_id'].astype(str))
with open("significant_gene_ids/both_gene_ids.txt", "w") as f:
    f.write(genes)
