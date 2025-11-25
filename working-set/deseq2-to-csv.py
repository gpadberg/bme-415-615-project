# Process DESeq2 datasets, filter significant genes, and save cleaned CSV outputs.

import pandas as pd


def process_deseq_file(input_path, output_path, lfc_threshold=1.0, padj_threshold=0.05):
    """
    Process a DESeq2 .tabular file:
      - Assign column names
      - Filter significant genes
      - Label as up/down regulated
      - Save to CSV
    """
    # Read the DESeq2 output
    df = pd.read_csv(input_path, sep="\t", header=None)
    df.columns = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

    # Ensure columns are numeric
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")

    # Filter by significance thresholds
    sig = df[
        (df["padj"] < padj_threshold) &
        (df["log2FoldChange"].abs() >= lfc_threshold)
    ].copy()

    # Add regulation column
    sig["regulation"] = sig["log2FoldChange"].apply(
        lambda x: "up" if x >= lfc_threshold else "down"
    )

    # Save the cleaned dataset
    sig.to_csv(output_path, index=False)
    print(f"Saved: {output_path} (n={len(sig)})")


# Run the processing

# Salt
process_deseq_file(
    input_path="DESeq2_datasets/salt.tabular",
    output_path="CSV_datasets/salt_significant_genes.csv"
)

# Heat
process_deseq_file(
    input_path="DESeq2_datasets/heat.tabular",
    output_path="CSV_datasets/heat_significant_genes.csv"
)