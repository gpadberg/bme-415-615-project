# overlap_metrics.py
# Compute normalized overlap metrics for GO-term sets (heat vs salt; up vs down)
# Reads your *_GO_significant.tsv files (exported from PANTHER) and prints + saves a summary CSV.

import pandas as pd
from pathlib import Path

script_dir = Path(__file__).resolve().parent
candidate1 = script_dir.with_stem("433_common_genes")
candidate2 = script_dir.parent / "433_common_genes"
root = candidate1 if candidate1.exists() else candidate2

export_dir = root / "top_GO_terms"
out_dir = root / "GO_plots"
out_dir.mkdir(exist_ok=True)

files = {
    "heat_up": export_dir / "heat_up_GO_significant.tsv",
    "heat_down": export_dir / "heat_down_GO_significant.tsv",
    "salt_up": export_dir / "salt_up_GO_significant.tsv",
    "salt_down": export_dir / "salt_down_GO_significant.tsv",
}

def get_term_col(df):
    for c in ["GO_term", "GO biological process complete"]:
        if c in df.columns:
            return c
    raise ValueError(f"Couldn't find GO term column. Columns: {list(df.columns)}")

def load_terms(path: Path) -> set[str]:
    df = pd.read_csv(path, sep="\t")
    tcol = get_term_col(df)
    return set(df[tcol].astype(str).str.strip())

def overlap_stats(A: set, B: set):
    inter = A & B
    union = A | B
    # use jaccard normalization method
    jaccard = len(inter) / len(union) if union else 0.0
    overlap_coeff = len(inter) / min(len(A), len(B)) if min(len(A), len(B)) else 0.0
    pct_A = len(inter) / len(A) * 100 if len(A) else 0.0
    pct_B = len(inter) / len(B) * 100 if len(B) else 0.0
    return {
        "n_A": len(A),
        "n_B": len(B),
        "n_intersection": len(inter),
        "n_union": len(union),
        "jaccard": jaccard,
        "overlap_coefficient": overlap_coeff,
        "pct_of_A_shared": pct_A,
        "pct_of_B_shared": pct_B,
    }

terms = {k: load_terms(p) for k, p in files.items()}

rows = []

up = overlap_stats(terms["heat_up"], terms["salt_up"])
up["comparison"] = "heat_up vs salt_up"
rows.append(up)

down = overlap_stats(terms["heat_down"], terms["salt_down"])
down["comparison"] = "heat_down vs salt_down"
rows.append(down)

summary = pd.DataFrame(rows)[[
    "comparison",
    "n_A", "n_B",
    "n_intersection", "n_union",
    "jaccard", "overlap_coefficient",
    "pct_of_A_shared", "pct_of_B_shared"
]]

pd.set_option("display.max_columns", None)
print(summary.to_string(index=False))

# Save CSV
out_path = out_dir / "overlap_normalization_summary.csv"
summary.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
