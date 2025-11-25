# Extract gene IDs into upregulated and downregulated files for each dataset

import pandas as pd
import os

# Make sure output directory exists
os.makedirs("significant_gene_ids", exist_ok=True)

def save_up_down(input_csv, prefix):
    """
    Reads a CSV, splits genes into up/down based on 'regulation_*' column,
    and saves each list to a .txt file.
    """
    df = pd.read_csv(input_csv)

    # Identify the regulation column automatically
    reg_col = [c for c in df.columns if c.startswith("regulation")][0]

    # UP
    up = df[df[reg_col] == "up"]['gene_id'].astype(str)
    up_out = f"significant_gene_ids/{prefix}_up.txt"
    up.to_csv(up_out, index=False, header=False)

    # DOWN
    down = df[df[reg_col] == "down"]['gene_id'].astype(str)
    down_out = f"significant_gene_ids/{prefix}_down.txt"
    down.to_csv(down_out, index=False, header=False)

    print(f"Saved: {up_out}, {down_out}")

# Salt only
save_up_down("sorted_gene_subsets/only_salt_genes.csv", "salt")

# Heat only
save_up_down("sorted_gene_subsets/only_heat_genes.csv", "heat")

# Both stresses
save_up_down("sorted_gene_subsets/both_stress_genes.csv", "both")
