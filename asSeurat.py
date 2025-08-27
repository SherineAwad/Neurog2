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
groupby_col = "celltype" 
# Run Wilcoxon - FIXED: use_raw=False since we're using scaled data
sc.tl.rank_genes_groups(
    adata,
    groupby=groupby_col,
    method='wilcoxon',
    use_raw=False,  # ← CHANGED: Use scaled data for Wilcoxon
    pts=True,
)

sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.1,
    max_out_group_fraction=1.0,
    key='rank_genes_groups',
    key_added='filtered_rank_genes_groups'
)

# MERGE ALL DGE RESULTS INTO ONE CSV
celltypes = adata.obs[groupby_col].unique()
all_dge_list = []
top_genes_dict = {}

for ct in celltypes:
    # Get results - FIXED: Use correct key 'rank_genes_groups'
    df = sc.get.rank_genes_groups_df(adata, group=ct, key='rank_genes_groups')  # ← CHANGED KEY

    # Rename columns appropriately for Wilcoxon
    df = df.rename(columns={
        'names': 'gene',
        'scores': 'wilcoxon_score',  # ← RENAMED: Clear this is Wilcoxon
        'pvals_adj': 'p_val_adj'
    })

    # Add column for cell type
    df['cell_type'] = ct

    # FILTERING CRITERIA FOR WILCOXON
    df_filtered = df[
        (df['p_val_adj'] < 0.05) &              # Significant
        (df['wilcoxon_score'].abs() > 2.0)      # Effect size threshold
    ].copy()

    # Calculate mean expression differences for context
    group_mask = adata.obs[groupby_col] == ct
    other_mask = ~group_mask

    # Calculate mean expression differences from SCALED data
    group_mean = np.array(adata.X[group_mask].mean(axis=0)).flatten()
    other_mean = np.array(adata.X[other_mask].mean(axis=0)).flatten()

    # Add mean differences (NOT log2FC) - FIXED: Honest naming
    df_filtered['mean_diff'] = group_mean[df_filtered.index] - other_mean[df_filtered.index]  # FIXED
    df_filtered['direction'] = df_filtered['mean_diff'].apply(lambda x: 'up' if x > 0 else 'down')

    all_dge_list.append(df_filtered)

    # Store top 10 genes for heatmap
    top_genes_dict[ct] = df_filtered.nlargest(10, 'wilcoxon_score')['gene'].tolist()

# Concatenate all DGE results and save
all_dge_df = pd.concat(all_dge_list, ignore_index=True)

# Add method info for clarity
all_dge_df['method'] = 'wilcoxon'

# Reorder columns logically
column_order = ['gene', 'cell_type', 'wilcoxon_score', 'p_val_adj', 'mean_diff', 'direction', 'method']
all_dge_df = all_dge_df[column_order]

all_dge_csv = f"{base_name}_all_DGE_noFCwilcoxon.csv"  # ← CHANGED: Reflects method
all_dge_df.to_csv(all_dge_csv, index=False)
print(f"Wilcoxon DGE results saved to {all_dge_csv}")
print(f"Total significant genes: {len(all_dge_df)}")


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Use the same celltype order as before
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']

# Get top genes from your Wilcoxon results
top_genes_all = []
for genes in top_genes_dict.values():
    top_genes_all.extend(genes)
top_genes_all = list(dict.fromkeys(top_genes_all))  # remove duplicates

# Subset adata to top genes
adata_top = adata[:, top_genes_all]

# Create a matrix of mean expression per cell type from SCALED data
mean_expr = pd.DataFrame(index=top_genes_all)
for ct in celltype_order:
    if ct in adata.obs[groupby_col].unique():
        ct_mask = adata_top.obs[groupby_col] == ct
        if hasattr(adata_top.X, 'toarray'):
            ct_data = adata_top.X[ct_mask].toarray()
        else:
            ct_data = adata_top.X[ct_mask]
        mean_expr[ct] = ct_data.mean(axis=0)

# Scale for visualization (z-score per gene)
mean_expr_scaled = mean_expr.copy()
for gene in mean_expr.index:
    gene_values = mean_expr.loc[gene]
    mean_expr_scaled.loc[gene] = (gene_values - gene_values.mean()) / gene_values.std()

# Reorder columns to match desired order
mean_expr_scaled = mean_expr_scaled[celltype_order]

# SORT GENES TO FORM DIAGONAL PATTERN
# For each gene, find which cell type has the highest expression
gene_max_celltype = mean_expr_scaled.idxmax(axis=1)

# Create a list to store sorted genes
sorted_genes = []

# For each cell type in order, get genes that have max expression in that cell type
for ct in celltype_order:
    genes_in_celltype = gene_max_celltype[gene_max_celltype == ct].index.tolist()
    # Sort these genes by their expression in this cell type (descending)
    genes_sorted = mean_expr_scaled.loc[genes_in_celltype, ct].sort_values(ascending=False).index.tolist()
    sorted_genes.extend(genes_sorted)

# Reorder the dataframe with sorted genes
mean_expr_scaled_sorted = mean_expr_scaled.loc[sorted_genes]

# Plot heatmap with diagonal pattern
plt.figure(figsize=(12, 10))
sns.heatmap(
    mean_expr_scaled_sorted,
    cmap='vlag',
    yticklabels=True,
    xticklabels=True,
    cbar_kws={'label': 'Z-score of mean expression'},
    center=0
)
plt.title("Top marker genes per cell type - Sorted for Diagonal Pattern", fontsize=16)
plt.ylabel("Genes")
plt.xlabel("Cell types")
plt.tight_layout()
heatmap_file = f"{base_name}_NoFC_diagonal_heatmap.png"
plt.savefig(heatmap_file, dpi=600, bbox_inches='tight')
plt.close()
print(f"Diagonal heatmap saved to {heatmap_file}")

# Optional: Also save the sorted expression data
mean_expr_scaled_sorted.to_csv(f"{base_name}_NoFC_diagonal_heatmap_expression_data.csv")
print("Sorted expression data for diagonal heatmap saved")



