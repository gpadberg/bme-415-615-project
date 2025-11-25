# Extract gene IDs into upregulated and downregulated files for each dataset

import pandas as pd
import os

os.makedirs("significant_gene_ids", exist_ok=True)

def save_up_down(input_csv, prefix):
    """
    Reads a CSV, splits genes into up/down based on the correct 'regulation_*' column,
    and saves each list to a .txt file.
    """
    df = pd.read_csv(input_csv)

    # Pick the regulation column that exists
    reg_cols = [c for c in df.columns if c.startswith("regulation")]
    if len(reg_cols) == 0:
        raise ValueError(f"No regulation column found in {input_csv}")

    # If multiple regulation columns, pick the one that matches the prefix
    reg_col = reg_cols[0]  # default
    for c in reg_cols:
        if c.endswith(prefix):
            reg_col = c
            break

    # UP genes
    up = df[df[reg_col] == "up"]['gene_id'].astype(str)
    up_out = f"significant_gene_ids/{prefix}_up.txt"
    up.to_csv(up_out, index=False, header=False)

    # DOWN genes
    down = df[df[reg_col] == "down"]['gene_id'].astype(str)
    down_out = f"significant_gene_ids/{prefix}_down.txt"
    down.to_csv(down_out, index=False, header=False)

    print(f"Saved: {up_out} ({len(up)} genes), {down_out} ({len(down)} genes)")

# Salt only
save_up_down("sorted_gene_subsets/only_salt_genes.csv", "salt")

# Heat only
save_up_down("sorted_gene_subsets/only_heat_genes.csv", "heat")

# Both stresses
save_up_down("sorted_gene_subsets/both_stress_genes.csv", "both")
