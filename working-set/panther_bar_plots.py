"""
Reads exported PANTHER GO-term results (top-go-terms + significant terms) and produces plots.
Generates:
- 4 bar charts (one per condition) of top 15 enriched GO terms
- 1 bar chart of number of significant GO terms per condition
- 1 bar chart of overlap of significant GO terms between heat and salt conditions
- 1 heatmap of top GO terms across all conditions
Saves plots to go_plots/ folder inside the 433_common_genes/ directory.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    """
    Main function to load data and generate plots.

    """
    export_dir, plots_dir = resolve_project_paths()

    # Load top15 files
    top_files = {
        "heat_up": export_dir / "heat_up_GO_top15.csv",
        "heat_down": export_dir / "heat_down_GO_top15.csv",
        "salt_up": export_dir / "salt_up_GO_top15.csv",
        "salt_down": export_dir / "salt_down_GO_top15.csv",
    }
    tops = {k: pd.read_csv(v) for k, v in top_files.items()}

    # Load significant files
    sig_files = {
        "heat_up": export_dir / "heat_up_GO_significant.tsv",
        "heat_down": export_dir / "heat_down_GO_significant.tsv",
        "salt_up": export_dir / "salt_up_GO_significant.tsv",
        "salt_down": export_dir / "salt_down_GO_significant.tsv",
    }
    sigs = {k: pd.read_csv(v, sep="\t") for k, v in sig_files.items()}

    # Load combined top15 if present (preferred for heatmap)
    combined_path = export_dir / "all_conditions_GO_top15.csv"
    combined = pd.read_csv(combined_path) if combined_path.exists() else None

    # Plot 1-4: 4 bar charts (top 15 each)
    for cond, df in tops.items():
        plot_top15_bar(
            df,
            title=f"{cond}: Top enriched GO terms (by FDR)",
            outpath=plots_dir / f"{cond}_top15_bar.png"
        )

    # Plot 5: significant term counts per condition
    counts = plot_sig_term_counts(
        sigs,
        title="Number of Enriched GO Terms Per Condition",
        outpath=plots_dir / "sig_term_counts.png"
    )

    # Plot 6: overlap heat vs salt (for up and down)
    overlap = plot_overlap(
        sigs,
        title="Overlap of Enriched Processes (heat vs salt)",
        outpath=plots_dir / "shared_counts.png",
        outpath_overlap=plots_dir / "overlapping_genes.txt"
    )

    # Plot 7: heatmap across conditions (union of top terms)
    if combined is not None:
        plot_heatmap_from_combined(
            combined,
            title="Top GO terms across conditions (missing=0)",
            outpath=plots_dir / "top_terms_heatmap.png"
        )

    # Tiny summary CVS
    summary = pd.DataFrame([{
        "condition": k,
        "n_sig_terms": counts[k],
    } for k in ["heat_up", "heat_down", "salt_up", "salt_down"]])
    summary["shared_up_terms"] = overlap["shared_up"]
    summary["shared_down_terms"] = overlap["shared_down"]
    summary.to_csv(plots_dir / "plot_summary_numbers.csv", index=False)

    print(f"Plots saved to: {plots_dir}")

def resolve_project_paths():
    """
    Resolves paths to the export directory and plots directory.
    Assumes this script is in working-set/ and that 433_common_genes/ is a sibling or in the parent directory.
    Also assumes that top_GO_terms/ already exists.
    """
    script_dir = Path(__file__).resolve().parent  # e.g., .../working-set

    # Changes from working-set folder to 433_common_genes folder
    candidate1 = script_dir.with_stem("433_common_genes")
    # Just in case this script is in the wrong folder
    candidate2 = script_dir.parent / "433_common_genes"

    panther_root = candidate1 if candidate1.exists() else candidate2
    if not panther_root.exists():
        raise FileNotFoundError(
            "Couldn't locate the 433_common_genes folder.\n"
            f"Tried:\n- {candidate1}\n- {candidate2}\n"
            "Move this script into your working-set folder or update the path logic."
        )

    export_dir = panther_root / "top_GO_terms"
    if not export_dir.exists():
        raise FileNotFoundError(
            f"Couldn't find export folder: {export_dir}\n"
            "Make sure you already ran the export script that creates top_GO_terms/."
        )

    plots_dir = panther_root / "GO_plots"
    plots_dir.mkdir(exist_ok=True)

    return export_dir, plots_dir

def first_existing_col(df, candidates):
    """
    Returns the first column name from candidates that exists in df, or None if none found.

    """
    for c in candidates:
        if c in df.columns:
            return c
    return None


def ensure_numeric(df, col):
    """
    Ensures that the specified column in df is numeric, coercing errors to NaN.
    
    """
    df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def plot_top15_bar(df, title, outpath):
    """
    Plots a horizontal bar chart of the top 15 GO terms by FDR. 
    """
    term_col = first_existing_col(df, ["GO_term", "GO biological process complete"])
    fdr_col = first_existing_col(df, ["FDR"])
    if term_col is None or fdr_col is None:
        raise ValueError(f"Missing required columns for bar plot. Found: {list(df.columns)}")

    d = df.copy()
    d = ensure_numeric(d, fdr_col)
    d = d.dropna(subset=[term_col, fdr_col])

    # Take top 15 by smallest FDR (uses the TSV files that include all the significant terms)
    d = d.nsmallest(min(15, len(d)), fdr_col)
    d["neglog10_FDR"] = -np.log10(d[fdr_col].clip(lower=1e-300))
    d = d.sort_values("neglog10_FDR")

    plt.figure(figsize=(10, 6))
    plt.barh(d[term_col].astype(str), d["neglog10_FDR"], color="#ff509d")
    plt.xlabel("-log10(FDR)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_heatmap_from_combined(combined, title, outpath):
    """
    Plots a heatmap of -log10(FDR) for GO terms across conditions from a combined dataframe.

    """
    term_col = first_existing_col(combined, ["GO_term", "GO biological process complete"])
    fdr_col = first_existing_col(combined, ["FDR"])
    if term_col is None or fdr_col is None or "condition" not in combined.columns:
        raise ValueError(f"Combined file missing required columns. Found: {list(combined.columns)}")

    d = combined.copy()
    d = ensure_numeric(d, fdr_col)
    d = d.dropna(subset=[term_col, fdr_col, "condition"])
    d["neglog10_FDR"] = -np.log10(d[fdr_col].clip(lower=1e-300))

    pivot = d.pivot_table(
        index=term_col,
        columns="condition",
        values="neglog10_FDR",
        aggfunc="min"
    ).fillna(0)

    # Sort by max significance across conditions
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=True).index]

    plt.figure(figsize=(9, max(6, 0.28 * len(pivot))))
    plt.imshow(pivot.values, aspect="auto")
    plt.yticks(range(len(pivot.index)), pivot.index.astype(str))
    plt.xticks(range(len(pivot.columns)), pivot.columns.astype(str), rotation=45, ha="right")
    plt.colorbar(label="-log10(FDR)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_sig_term_counts(sig_tables, title, outpath):
    """
    Plots a bar chart of the number of significant GO terms per condition.

    """
    counts = {k: len(v) for k, v in sig_tables.items()}

    plt.figure(figsize=(7, 4))
    plt.bar(list(counts.keys()), list(counts.values()), color="#ff509d")
    plt.ylabel("# significant GO terms (FDR ≤ 0.05)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

    return counts

def plot_overlap(sig_tables, title, outpath, outpath_overlap):
    """
    Plots a bar chart of the overlap of significant GO terms between heat and salt conditions,
    and prints the overlapping GO term lists.
    """
    def term_set(df):
        tcol = first_existing_col(df, ["GO_term", "GO biological process complete"])
        if tcol is None:
            return set()
        # normalize: string, strip whitespace
        return set(df[tcol].astype(str).str.strip())

    hu = term_set(sig_tables["heat_up"])
    su = term_set(sig_tables["salt_up"])
    hd = term_set(sig_tables["heat_down"])
    sd = term_set(sig_tables["salt_down"])

    overlap_up = sorted(hu & su)
    overlap_down = sorted(hd & sd)

    with open(outpath_overlap, "w") as f:
        f.write(f"Overlapping UP GO terms: {len(overlap_up)}\n")
        for term in overlap_up:
            f.write(f"{term}\n")
        f.write(f"\nOverlapping DOWN GO terms: {len(overlap_down)}\n")
        for term in overlap_down:
            f.write(f"{term}\n")

    # # Print the lists
    # print(f"\nOverlapping UP GO terms (heat_up ∩ salt_up): {len(overlap_up)}")
    # for term in overlap_up:
    #     print(term)

    # print(f"\nOverlapping DOWN GO terms (heat_down ∩ salt_down): {len(overlap_down)}")
    # for term in overlap_down:
    #     print(term)

    # Plot counts
    shared_up = len(overlap_up)
    shared_down = len(overlap_down)

    plt.figure(figsize=(6, 4))
    plt.bar(["heat_up ∩ salt_up", "heat_down ∩ salt_down"], [shared_up, shared_down], color="#ff509d")
    plt.ylabel("# shared significant GO terms")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

    return {
        "shared_up": shared_up,
        "shared_down": shared_down,
        "overlap_up_terms": overlap_up,
        "overlap_down_terms": overlap_down,
    }


if __name__ == "__main__":
    main()
    print("end :)")