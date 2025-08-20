import os
import numpy as np
import pandas as pd
import scanpy as sc
import argparse

# -----------------------------
# ARGUMENTS
# -----------------------------
parser = argparse.ArgumentParser(description="DGE with Seurat-style log2FC")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

# -----------------------------
# READ DATA
# -----------------------------
adata = sc.read(myObject)

# -----------------------------
# DGE (Scanpy Wilcoxon)
# -----------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="celltype",  # directly use column
    method="wilcoxon",
    use_raw=False,
    pts=True
)

# Pick matrix for mean calculation
if adata.raw is not None:
    X = adata.raw.X
elif "counts" in adata.layers:
    X = adata.layers["counts"]
else:
    X = adata.X

if hasattr(X, "toarray"):  # convert sparse to dense
    X = X.toarray()

celltypes = adata.obs["celltype"].unique()
all_dge_list = []

for ct in celltypes:
    df = sc.get.rank_genes_groups_df(adata, group=ct)

    # Compute Seurat-style log2FC
    mask = adata.obs["celltype"] == ct
    group_mean = np.array(X[mask, :].mean(axis=0)).flatten()
    rest_mean  = np.array(X[~mask, :].mean(axis=0)).flatten()

    # avoid division by zero
    eps = 1e-9
    group_mean = np.maximum(group_mean, eps)
    rest_mean  = np.maximum(rest_mean, eps)

    logFC = np.log2(group_mean / rest_mean)
    gene_to_logFC = dict(zip(adata.var_names, logFC))

    # Assign computed logFC to both columns
    df["avg_log2FC"] = df["names"].map(gene_to_logFC)
    df["logfoldchanges"] = df["names"].map(gene_to_logFC)

    # Rename columns to Seurat-style
    df = df.rename(columns={
        "names": "gene",
        "pvals_adj": "p_val_adj"
    })

    df["cell_type"] = ct
    all_dge_list.append(df)

# -----------------------------
# SAVE RESULTS
# -----------------------------
all_dge_df = pd.concat(all_dge_list, ignore_index=True)
out_csv = f"{base_name}_all_DGE_FC.csv"
all_dge_df.to_csv(out_csv, index=False)

print(f"DGE results saved to: {out_csv}")


# -----------------------------
# HEATMAP
# -----------------------------
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# PARAMETERS
# -----------------------------
top_n = 10  # top N genes per cell type

# -----------------------------
# SELECT TOP GENES PER CELL TYPE
# -----------------------------
top_genes_all = []

for ct in celltype_order:
    # Subset DGE for this cell type
    df_ct = all_dge_df[all_dge_df["cell_type"] == ct]
    
    # Sort by logfoldchanges descending and select top_n genes
    top_genes_ct = df_ct.sort_values("logfoldchanges", ascending=False).head(top_n)["gene"].tolist()
    
    # Add to combined list
    top_genes_all.extend(top_genes_ct)

# Remove duplicates while preserving order
top_genes_all = list(dict.fromkeys(top_genes_all))

# -----------------------------
# BUILD HEATMAP MATRIX
# -----------------------------
heatmap_matrix = pd.DataFrame(index=top_genes_all, columns=celltype_order)

for ct in celltype_order:
    df_ct = all_dge_df[all_dge_df["cell_type"] == ct]
    # Map logFC values to heatmap
    heatmap_matrix[ct] = df_ct.set_index("gene")["logfoldchanges"].reindex(top_genes_all)

# Fill missing values with 0
heatmap_matrix = heatmap_matrix.fillna(0)

# -----------------------------
# PLOT HEATMAP
# -----------------------------
plt.figure(figsize=(12, max(6, len(top_genes_all)*0.3)))
sns.heatmap(
    heatmap_matrix,
    cmap="vlag",
    yticklabels=True,
    xticklabels=heatmap_matrix.columns,
    cbar_kws={"label": "log2FC"}
)
plt.title("Top marker genes per cell type (log2FC)", fontsize=16)
plt.ylabel("Genes")
plt.xlabel("Cell types")
plt.tight_layout()

heatmap_file = f"{base_name}_top_genes_heatmap_FC.png"
plt.savefig(heatmap_file, dpi=600)
plt.close()
print(f"Heatmap (logFC) saved to: {heatmap_file}")



