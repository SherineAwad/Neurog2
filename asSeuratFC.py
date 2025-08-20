import os
import numpy as np
import pandas as pd
import scanpy as sc
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
# ARGUMENTS
parser = argparse.ArgumentParser(description="DGE with Seurat-style log2FC")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
# READ DATA
adata = sc.read(myObject)
# Differential gene expression (Wilcoxon)
sc.tl.rank_genes_groups(
    adata,
    groupby="celltype",
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
# Convert sparse to dense if necessary
if hasattr(X, "toarray"):
    X = X.toarray()

celltypes = adata.obs["celltype"].unique()
all_dge_list = []

pseudocount = 1e-6  # small number to avoid divide-by-zero
# Calculate Seurat-style log2FC safely
for ct in celltypes:
    df = sc.get.rank_genes_groups_df(adata, group=ct)
    mask = adata.obs["celltype"] == ct

    # Linearize data only if it seems log-transformed
    if np.max(X) < 50:
        X_lin = np.expm1(X)  # exp(X) - 1
    else:
        X_lin = X.copy()

    # Compute mean expression per gene
    group_mean = np.array([
        np.nanmean(X_lin[mask, i][X_lin[mask, i] > 0]) + pseudocount
        for i in range(X_lin.shape[1])
    ])
    rest_mean = np.array([
        np.nanmean(X_lin[~mask, i][X_lin[~mask, i] > 0]) + pseudocount
        for i in range(X_lin.shape[1])
    ])

    # Compute log2 fold change safely
    logFC = np.log2(group_mean / rest_mean)

    # Map logFC to genes
    gene_to_logFC = dict(zip(adata.var_names, logFC))
    df["logfoldchanges"] = df["names"].map(gene_to_logFC)
    df["avg_log2FC"] = df["names"].map(gene_to_logFC)

    # Rename columns to Seurat-style
    df = df.rename(columns={"names": "gene", "pvals_adj": "p_val_adj"})
    df["cell_type"] = ct
    all_dge_list.append(df)

# Concatenate all results
all_dge_df = pd.concat(all_dge_list, ignore_index=True)
all_dge_csv = f"{base_name}_all_DGE_FC.csv"
all_dge_df.to_csv(all_dge_csv, index=False)
print(f"All DGE results saved to {all_dge_csv}")
# -----------------------------
# Heatmap of top 10 genes per cell type
# -----------------------------
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']

top_genes_all = []
for ct in celltype_order:
    df_ct = all_dge_df[all_dge_df["cell_type"] == ct]
    top_genes = df_ct.nlargest(10, "logfoldchanges")["gene"].tolist()
    top_genes_all.extend(top_genes)

# Remove duplicates while preserving order
top_genes_all = list(dict.fromkeys(top_genes_all))

# Create heatmap matrix (genes x cell types) with log2FC
heatmap_matrix = pd.DataFrame(index=top_genes_all, columns=celltype_order)
for ct in celltype_order:
    if ct in all_dge_df["cell_type"].unique():
        df_ct = all_dge_df[all_dge_df["cell_type"] == ct]
        heatmap_matrix[ct] = df_ct.set_index("gene")["logfoldchanges"].reindex(top_genes_all)

# Fill missing values with 0
heatmap_matrix = heatmap_matrix.fillna(0)

# -----------------------------
# Plot heatmap
# -----------------------------
plt.figure(figsize=(15, 15))
sns.heatmap(
    heatmap_matrix,
    cmap="vlag",         # blue-red diverging
    yticklabels=True,
    xticklabels=heatmap_matrix.columns,
    cbar_kws={"label": "log2FC"}
)
plt.title("Top marker genes per cell type (log2FC)", fontsize=16)
plt.ylabel("Genes")
plt.xlabel("Cell types")
plt.tight_layout()

heatmap_file = f"{base_name}_top_genes_logFC_heatmap.png"
plt.savefig(heatmap_file, dpi=300)
plt.close()
print(f"Heatmap (logFC) saved to: {heatmap_file}")



