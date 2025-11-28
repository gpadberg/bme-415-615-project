from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TOP_N = 500
PADJ_THRESH = 0.05
LFC_THRESH = 2.0

TITLE = "Significant DEGs Across Heat/Salt and Up/Down-Regulated Sets"
OUTFILE = "DEG_heatmap_top500.png"

DEG_DIRNAME = "DESeq2_datasets"
PLOTS_DIRNAME = "plots"

def _infer_col(df, candidates, fallback_first=False):
    """
    Infers which column in DF matches one of the CANDIDATES.
    """
    for c in candidates:
        if c in df.columns:
            return c
    if fallback_first:
        return df.columns[0]
    raise ValueError(f"Could not find expected column. Tried: {candidates}. Found: {list(df.columns)}")


def load_deseq2_tabular(path: Path) -> pd.DataFrame:
    """
    Load a DESeq2-style tabular file.
    Standardizes output to columns: gene, log2FC, padj
    """
    first_data_line = None
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            first_data_line = s
            break

    if first_data_line is None:
        raise ValueError(f"{path.name} appears to be empty or only comments.")

    tokens = first_data_line.split("\t")

    header_like = any(t in {"log2FoldChange", "padj", "pvalue", "baseMean"} for t in tokens)

    if header_like:
        df = pd.read_csv(path, sep="\t", comment="#", header=0)
    else:
        df = pd.read_csv(path, sep="\t", comment="#", header=None)

        # Typical DESeq2 result tables have 7 columns:
        # gene/baseMean/log2FC/lfcSE/stat/pvalue/padj
        if df.shape[1] >= 7:
            df = df.iloc[:, :7].copy()
            df.columns = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
        elif df.shape[1] == 6:
            df = df.copy()
            df.columns = ["gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]
        else:
            raise ValueError(
                f"{path.name} has {df.shape[1]} columns; expected at least 6 or 7 for DESeq2-style output."
            )

    # Standardize column names
    gene_col = _infer_col(df, ["gene", "id", "Gene", "gene_id", "Name", "name"], fallback_first=True)
    lfc_col  = _infer_col(df, ["log2FoldChange", "logFC"])
    padj_col = _infer_col(df, ["padj", "FDR", "qvalue"])

    out = df[[gene_col, lfc_col, padj_col]].copy()
    out.columns = ["gene", "log2FC", "padj"]
    out["log2FC"] = pd.to_numeric(out["log2FC"], errors="coerce")
    out["padj"] = pd.to_numeric(out["padj"], errors="coerce")
    out = out.dropna(subset=["gene", "log2FC", "padj"])
    return out

def find_deg_files(deg_dir: Path):
    """
    Tries to locate heat and salt DEG tables inside DEG_DIRNAME using filename keywords.
    Accepts .tabular, .tsv, .txt.
    """
    candidates = []
    for ext in ("*.tabular", "*.tsv", "*.txt"):
        candidates.extend(deg_dir.glob(ext))

    def pick(keyword: str):
        hits = [p for p in candidates if keyword in p.name.lower()]
        hits.sort(key=lambda p: ("deseq" not in p.name.lower(), "result" not in p.name.lower(), len(p.name)))
        return hits[0] if hits else None

    heat = pick("heat")
    salt = pick("salt")
    return heat, salt

def plot_deg_heatmap(heat_path: Path, salt_path: Path, outpath: Path, title: str, padj_thresh=0.05, lfc_thresh=2.0, top_n=500,):
    """
    Plots a heatmap of DEGs across heat and salt conditions.
    """
    heat_df = load_deseq2_tabular(heat_path)
    salt_df = load_deseq2_tabular(salt_path)

    heat_sig = heat_df[(heat_df["padj"] <= padj_thresh) & (heat_df["log2FC"].abs() >= lfc_thresh)]
    salt_sig = salt_df[(salt_df["padj"] <= padj_thresh) & (salt_df["log2FC"].abs() >= lfc_thresh)]

    heat_up = heat_sig[heat_sig["log2FC"] > 0][["gene", "log2FC"]].copy()
    heat_down = heat_sig[heat_sig["log2FC"] < 0][["gene", "log2FC"]].copy()
    salt_up = salt_sig[salt_sig["log2FC"] > 0][["gene", "log2FC"]].copy()
    salt_down = salt_sig[salt_sig["log2FC"] < 0][["gene", "log2FC"]].copy()

    pieces = []
    for df, cond in [
        (heat_up, "Heat up"),
        (heat_down, "Heat down"),
        (salt_up, "Salt up"),
        (salt_down, "Salt down"),
    ]:
        if not df.empty:
            tmp = df.copy()
            tmp["condition"] = cond
            pieces.append(tmp)

    if not pieces:
        raise ValueError("No significant DEGs found with the current PADJ/LFC thresholds.")

    long_df = pd.concat(pieces, ignore_index=True)

    pivot = (
        long_df.pivot_table(index="gene", columns="condition", values="log2FC", aggfunc="mean")
        .fillna(0)
    )

    desired_cols = ["Heat up", "Heat down", "Salt up", "Salt down"]
    for c in desired_cols:
        if c not in pivot.columns:
            pivot[c] = 0.0
    pivot = pivot[desired_cols]

    # Pick top genes
    if top_n is not None and top_n > 0 and len(pivot) > top_n:
        strength = pivot.abs().max(axis=1)
        pivot = pivot.loc[strength.sort_values(ascending=False).head(top_n).index]

    # Sort weak to strong
    pivot = pivot.loc[pivot.abs().max(axis=1).sort_values(ascending=True).index]

    # Symmetric scaling around 0 so up/down are comparable
    vmax = float(np.nanmax(np.abs(pivot.values))) if pivot.size else 1.0
    if vmax == 0:
        vmax = 1.0

    fig = plt.figure(figsize=(7.5, 10))
    ax = fig.add_subplot(111)
    im = ax.imshow(
        pivot.values,
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )

    ax.set_title(title)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns.astype(str), rotation=45, ha="right")

    # Hide gene names
    ax.set_yticks([])

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log2 Fold Change")

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=400, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    here = Path(__file__).resolve()
    working_set = here.parent
    deg_dir = working_set / DEG_DIRNAME
    plots_dir = working_set / PLOTS_DIRNAME
    outpath = plots_dir / OUTFILE

    if not deg_dir.exists():
        raise FileNotFoundError(f"Can't find DEG folder: {deg_dir}")

    heat_path, salt_path = find_deg_files(deg_dir)
    if heat_path is None or salt_path is None:
        raise FileNotFoundError(
            f"Couldn't auto-detect heat/salt files in {deg_dir}.\n"
            f"Make sure filenames contain 'heat' and 'salt' (e.g., heat.tabular, salt.tabular)."
        )

    print(f"Using heat DEG file: {heat_path.name}")
    print(f"Using salt DEG file: {salt_path.name}")
    print(f"Output: {outpath}")

    plot_deg_heatmap(
        heat_path=heat_path,
        salt_path=salt_path,
        outpath=outpath,
        title=TITLE,
        padj_thresh=PADJ_THRESH,
        lfc_thresh=LFC_THRESH,
        top_n=TOP_N,
    )

    print("End :)")
