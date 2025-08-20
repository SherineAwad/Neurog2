import os
import numpy as np
import pandas as pd
import scanpy as sc
import argparse
import importlib_metadata
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata
# ARGUMENTS
parser = argparse.ArgumentParser(description="Perform DGE and plot top genes heatmap")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()
myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "expression_" + base_name + ".h5ad"
# READ DATA
adata = sc.read(myObject)
# PARAMETERS
groupby_col = 'celltype'  # replace with your cell type column
# DIFFERENTIAL GENE EXPRESSION
sc.tl.rank_genes_groups(
    adata,
    groupby=groupby_col,
    method='wilcoxon',   # Seurat default
    use_raw=False,
    pts=True             # include percent expressed
)
# MERGE ALL DGE RESULTS INTO ONE CSV
celltypes = adata.obs[groupby_col].unique()
all_dge_list = []
top_genes_dict = {}
for ct in celltypes:
    df = sc.get.rank_genes_groups_df(adata, group=ct)
    
    # Rename columns to Seurat-like
    df = df.rename(columns={
        'names': 'gene',
        'logfoldchanges': 'avg_log2FC',
        'pvals_adj': 'p_val_adj'
    })
    
    # Add column for cell type
    df['cell_type'] = ct
    
    all_dge_list.append(df)
    
    # Store top 10 genes for heatmap
    top_genes_dict[ct] = df.sort_values('scores', ascending=False)['gene'].head(10).tolist()
# Concatenate all DGE results and save to current directory
all_dge_df = pd.concat(all_dge_list, ignore_index=True)
all_dge_csv = f"{base_name}_all_DGE.csv"
all_dge_df.to_csv(all_dge_csv, index=False)
print(f"All DGE results saved to {all_dge_csv}")
#Heatmaps 
import seaborn as sns

celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']

top_genes_all = []
for genes in top_genes_dict.values():
    top_genes_all.extend(genes)
top_genes_all = list(dict.fromkeys(top_genes_all))  # remove duplicates

# Subset adata to top genes
adata_top = adata[:, top_genes_all]

# Create a matrix of mean expression per cell type
mean_expr = pd.DataFrame(index=top_genes_all)
for ct in celltype_order:  # iterate over your desired order
    if ct in celltypes:     # only include existing cell types
        mean_expr[ct] = adata_top[adata_top.obs[groupby_col] == ct].X.toarray().mean(axis=0)

# Scale for visualization (z-score per gene)
mean_expr_scaled = (mean_expr - mean_expr.mean(axis=1).values[:, None]) / mean_expr.std(axis=1).values[:, None]

# Plot heatmap
plt.figure(figsize=(15, 15))
sns.heatmap(
    mean_expr_scaled, 
    cmap='vlag', 
    yticklabels=True, 
    xticklabels=mean_expr_scaled.columns,  # explicitly use ordered columns
    cbar_kws={'label': 'Z-score'}
)
plt.title("Top 10 marker genes per cell type", fontsize=16)
plt.ylabel("Genes")
plt.xlabel("Cell types")
plt.tight_layout()
heatmap_file = f"{base_name}_top10_genes_heatmap.png"
plt.savefig(heatmap_file, dpi=600)
plt.close()
print(f"Heatmap saved to {heatmap_file}")

