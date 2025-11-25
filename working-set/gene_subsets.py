# Sort DESeq2 significant gene datasets into:
#   - genes unique to salt
#   - genes unique to heat
#   - genes shared by both stresses

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


def sort_gene_sets(salt_path, heat_path, out_dir):
    """Merge salt + heat sets and split into 'both', 'salt only', and 'heat only'."""
    sig_salt = load_sig_file(salt_path)
    sig_heat = load_sig_file(heat_path)

    # Merge by gene_id
    merged = pd.merge(
        sig_salt,
        sig_heat,
        on="gene_id",
        how="outer",
        suffixes=("_salt", "_heat")
    )

    # Genes significant in both
    both = merged.dropna(subset=["log2FoldChange_salt", "log2FoldChange_heat"])

    # Unique to salt
    only_salt = merged[
        merged["log2FoldChange_heat"].isna() & 
        merged["log2FoldChange_salt"].notna()
    ]

    # Unique to heat
    only_heat = merged[
        merged["log2FoldChange_salt"].isna() &
        merged["log2FoldChange_heat"].notna()
    ]

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Save files
    both_path = os.path.join(out_dir, "both_stress_genes.csv")
    salt_only_path = os.path.join(out_dir, "only_salt_genes.csv")
    heat_only_path = os.path.join(out_dir, "only_heat_genes.csv")

    both.to_csv(both_path, index=False)
    only_salt.to_csv(salt_only_path, index=False)
    only_heat.to_csv(heat_only_path, index=False)

    print(f"Saved:\n  - {both_path} (common genes: {len(both)})\n"
          f"  - {salt_only_path} (salt-only: {len(only_salt)})\n"
          f"  - {heat_only_path} (heat-only: {len(only_heat)})")


# Run
sort_gene_sets(
    salt_path="CSV_datasets/salt_significant_genes.csv",
    heat_path="CSV_datasets/heat_significant_genes.csv",
    out_dir="sorted_gene_subsets"
)
