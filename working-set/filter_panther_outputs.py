"""
This script returns the top 15 significant GO terms from PANTHER overrepresentation analyses
for four gene sets: heat_up, heat_down, salt_up, salt_down. These results tell us which biological
processes are most enriched in each gene set.

"""

import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parent # path to the folder holding this file

panther_file_root = root.with_stem("433_common_genes") # switch from working_set folder to 433_common_genes folder

panther_file = panther_file_root / "panther_output" # folder with the panther output files

def load_panther(filename):
    """
    Loads a panther overrepresentation output file and renames the FDR column to a generic name 'FDR'.

    """
    df = pd.read_csv(panther_file / filename, sep="\t", skiprows=11)

    # Gets default FDR column name and renames it to 'FDR'
    fdr_col = [c for c in df.columns if c.endswith("(FDR)")]
    if not fdr_col:
        raise ValueError(f"No FDR column found in {filename}")
    df = df.rename(columns={fdr_col[0]: "FDR"})
    return df

# Load all four outputs
heat_up   = load_panther("heat_up_analysis.txt")
heat_down = load_panther("heat_down_analysis.txt")
salt_up   = load_panther("salt_up_analysis.txt")
salt_down = load_panther("salt_down_analysis.txt")

# Filter to significant GO terms
alpha = 0.05 # significance threshold
heat_up_sig   = heat_up[heat_up["FDR"] <= alpha]
heat_down_sig = heat_down[heat_down["FDR"] <= alpha]
salt_up_sig   = salt_up[salt_up["FDR"] <= alpha]
salt_down_sig = salt_down[salt_down["FDR"] <= alpha]

# Take top N by FDR
n = 15 # number of top terms to return
heat_up_top   = heat_up_sig.nsmallest(n, "FDR")
heat_down_top = heat_down_sig.nsmallest(n, "FDR")
salt_up_top   = salt_up_sig.nsmallest(n, "FDR")
salt_down_top = salt_down_sig.nsmallest(n, "FDR")

print(heat_up_top)
print(heat_down_top)
print(salt_up_top)
print(salt_down_top)