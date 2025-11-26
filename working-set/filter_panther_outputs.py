"""
This script returns the top 15 significant GO terms from PANTHER overrepresentation analyses
for four gene sets: heat_up, heat_down, salt_up, salt_down. These results tell us which biological
processes are most enriched in each gene set.

"""

import pandas as pd
from pathlib import Path
import csv
import re

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

export_path = panther_file_root / "top_GO_terms"
export_path.mkdir(exist_ok=True)

# Choose a term column that usually exists in BP outputs
TERM_COL = "GO biological process complete"

# Pick the columns you care about (only keep ones that exist)
def pick_cols(df):
    wanted = [TERM_COL, "Fold Enrichment", "FDR", "P-value", "# in Query", "# in Reference", "Expected"]
    return [c for c in wanted if c in df.columns]

results = {
    "heat_up": (heat_up_sig, heat_up_top),
    "heat_down": (heat_down_sig, heat_down_top),
    "salt_up": (salt_up_sig, salt_up_top),
    "salt_down": (salt_down_sig, salt_down_top),
}

# Export per-condition files (full significant + top 15)
for name, (sig_df, top_df) in results.items():
    cols_sig = pick_cols(sig_df)
    cols_top = pick_cols(top_df)

    # export all significant as TSV
    # TSV is tab separated while CSV is comma separated
    # TSV may be good for processing since some GO terms have commas in them
    sig_df[cols_sig].to_csv(export_path / f"{name}_GO_significant.tsv", sep="\t", index=False)

    # export top 15 as CSV
    # these will be good for putting into tables for the actual report paper
    # CSVs are easier to open in Excel
    top_df[cols_top].to_csv(export_path / f"{name}_GO_top15.csv", index=False, quoting=csv.QUOTE_MINIMAL)

# Creates a master file with all top 15 GO terms from each condition
# This will be useful for plotting comparisons
combined = []
for name, (_, top_df) in results.items():
    df = top_df.copy()
    df = df[pick_cols(df)]
    df.insert(0, "condition", name)
    combined.append(df)

combined_df = pd.concat(combined, ignore_index=True)
combined_df.to_csv(export_path / "all_conditions_GO_top15.csv", index=False)