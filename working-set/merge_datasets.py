# Create a merged CSV containing the union of significant salt + heat results.

import pandas as pd
import os


def load_sig_file(path):
    """Load a significant-gene CSV and ensure required columns exist."""
    df = pd.read_csv(path)

    required = {"gene_id", "log2FoldChange", "regulation", "padj"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    return df


def merge_sig_sets(salt_path, heat_path, output_path):
    """Merge salt + heat significant gene sets using an outer join."""
    sig_salt = load_sig_file(salt_path)
    sig_heat = load_sig_file(heat_path)

    merged = pd.merge(
        sig_salt,
        sig_heat,
        on="gene_id",
        how="outer",
        suffixes=("_salt", "_heat")
    )

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    merged.to_csv(output_path, index=False)
    print(f"Merged dataset saved: {output_path} (n={len(merged)})")


# Run merge

merge_sig_sets(
    salt_path="CSV_datasets/salt_significant_genes.csv",
    heat_path="CSV_datasets/heat_significant_genes.csv",
    output_path="merged_dataset/merged_salt_heat_genes.csv"
)
