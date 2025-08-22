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

# REMOVED: No logfoldchanges with Wilcoxon
# logfcs = adata.uns['rank_genes_groups']['logfoldchanges']  # ← THIS WILL ERROR

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
    df_filtered['mean_diff'] = group_mean[df_filtered.index] - other_mean[df_filtered.index]
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

all_dge_csv = f"{base_name}_all_DGE_wilcoxon.csv"  # ← CHANGED: Reflects method
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

# Subset adata to top genes - FIXED: Use the main adata object (scaled data)
adata_top = adata[:, top_genes_all]

# Create a matrix of mean expression per cell type from SCALED data
mean_expr = pd.DataFrame(index=top_genes_all)
for ct in celltype_order:
    if ct in adata.obs[groupby_col].unique():  # only include existing cell types
        # Get cells of this type and calculate mean expression
        ct_mask = adata_top.obs[groupby_col] == ct
        if hasattr(adata_top.X, 'toarray'):
            ct_data = adata_top.X[ct_mask].toarray()
        else:
            ct_data = adata_top.X[ct_mask]
        mean_expr[ct] = ct_data.mean(axis=0)

# Scale for visualization (z-score per gene) - FIXED: Proper scaling
mean_expr_scaled = mean_expr.copy()
for gene in mean_expr.index:
    gene_values = mean_expr.loc[gene]
    mean_expr_scaled.loc[gene] = (gene_values - gene_values.mean()) / gene_values.std()

# Reorder columns to match desired order
mean_expr_scaled = mean_expr_scaled[celltype_order]

# Plot heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(
    mean_expr_scaled, 
    cmap='vlag', 
    yticklabels=True, 
    xticklabels=True,
    cbar_kws={'label': 'Z-score of mean expression'},
    center=0  # Center colormap at 0 for scaled data
)
plt.title("Top 10 marker genes per cell type (Wilcoxon)", fontsize=16)
plt.ylabel("Genes")
plt.xlabel("Cell types")
plt.tight_layout()
heatmap_file = f"{base_name}_top10_genes_heatmap.png"
plt.savefig(heatmap_file, dpi=600, bbox_inches='tight')
plt.close()
print(f"Heatmap saved to {heatmap_file}")

# Optional: Also save the expression data
mean_expr_scaled.to_csv(f"{base_name}_heatmap_expression_data.csv")
print("Expression data for heatmap saved")
