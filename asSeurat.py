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
parser = argparse.ArgumentParser() 
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('min_in', type=float, default=0.1) 
args = parser.parse_args()
myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "expression_" + base_name + ".h5ad"
min_in = args.min_in

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
    min_in_group_fraction=min_in,
    max_out_group_fraction=1.0,
    key='rank_genes_groups',
    key_added='filtered_rank_genes_groups'
)

# MERGE ALL DGE RESULTS INTO ONE CSV
celltypes = adata.obs[groupby_col].unique()
all_dge_list = []
top_genes_dict = {}

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------------------------------
# PARAMETERS
# -------------------------------------------------
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
groupby_col = "celltype"   # <- change if your column is named differently
base_name = "mydataset"    # <- for saving
min_in = 10                # <- placeholder, adjust as needed

all_dge_list = []
top_genes_dict = {}

# -------------------------------------------------
# DIFFERENTIAL EXPRESSION LOOP
# -------------------------------------------------
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

# -------------------------------------------------
# PARAMETERS
# -------------------------------------------------
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
groupby_col = "celltype"   # <- change if your column is named differently
base_name = "mydataset"    # <- for saving
min_in = 10                # <- placeholder, adjust as needed

all_dge_list = []
top_genes_dict = {}

# -------------------------------------------------
# DIFFERENTIAL EXPRESSION LOOP (AND logic)
# -------------------------------------------------
for ct in celltype_order:
    if ct not in adata.obs[groupby_col].unique():
        continue

    # Get DE results
    df = sc.get.rank_genes_groups_df(adata, group=ct, key='rank_genes_groups')

    # Rename columns
    df = df.rename(columns={
        'names': 'gene',
        'scores': 'wilcoxon_score',
        'pvals_adj': 'p_val_adj'
    })
    df['cell_type'] = ct

    # Calculate mean differences (fix: align with gene names)
    group_mask = adata.obs[groupby_col] == ct
    other_mask = ~group_mask

    group_mean = np.array(adata.X[group_mask].mean(axis=0)).flatten()
    other_mean = np.array(adata.X[other_mask].mean(axis=0)).flatten()
    gene_idx = {gene: i for i, gene in enumerate(adata.var_names)}
    df['mean_diff'] = df['gene'].map(lambda g: group_mean[gene_idx[g]] - other_mean[gene_idx[g]])
    df['direction'] = df['mean_diff'].apply(lambda x: 'up' if x > 0 else 'down')

    # Filter: AND logic (Wilcoxon AND mean_diff thresholds)
    df_filtered = df[
        (df['p_val_adj'] < 0.05) &
        (df['wilcoxon_score'].abs() > 2.0) &
        (df['mean_diff'].abs() > 0.25)
    ].copy()

    all_dge_list.append(df_filtered)

    # Top 10 genes by Wilcoxon score (after filtering)
    top_genes_dict[ct] = df_filtered.nlargest(10, 'wilcoxon_score')['gene'].tolist()

# Concatenate all DGE results
all_dge_df = pd.concat(all_dge_list, ignore_index=True)
all_dge_df['method'] = 'wilcoxon'
column_order = ['gene', 'cell_type', 'wilcoxon_score', 'p_val_adj',
                'mean_diff', 'direction', 'method']
all_dge_df = all_dge_df[column_order]

# Save DE results
all_dge_csv = f"{base_name}_all_DGE_noFCwilcoxon_{min_in}.csv"
all_dge_df.to_csv(all_dge_csv, index=False)
print(f"Wilcoxon DGE results saved to {all_dge_csv}")
print(f"Total significant genes: {len(all_dge_df)}")

# -------------------------------------------------
# HEATMAP OF TOP GENES
# -------------------------------------------------
# Collect unique top genes
top_genes_all = []
for genes in top_genes_dict.values():
    top_genes_all.extend(genes)
top_genes_all = list(dict.fromkeys(top_genes_all))  # remove duplicates

# Subset adata
adata_top = adata[:, top_genes_all]

# Compute mean expression per cell type
mean_expr = pd.DataFrame(index=top_genes_all)
for ct in celltype_order:
    if ct in adata.obs[groupby_col].unique():
        ct_mask = adata_top.obs[groupby_col] == ct
        if hasattr(adata_top.X, 'toarray'):
            ct_data = adata_top.X[ct_mask].toarray()
        else:
            ct_data = adata_top.X[ct_mask]
        mean_expr[ct] = ct_data.mean(axis=0)

# Scale per gene (z-score)
mean_expr_scaled = mean_expr.copy()
for gene in mean_expr.index:
    gene_values = mean_expr.loc[gene]
    mean_expr_scaled.loc[gene] = (gene_values - gene_values.mean()) / gene_values.std()
mean_expr_scaled = mean_expr_scaled[celltype_order]

# Sort genes diagonally
gene_max_celltype = mean_expr_scaled.idxmax(axis=1)
sorted_genes = []
for ct in celltype_order:
    genes_in_celltype = gene_max_celltype[gene_max_celltype == ct].index.tolist()
    genes_sorted = mean_expr_scaled.loc[genes_in_celltype, ct].sort_values(ascending=False).index.tolist()
    sorted_genes.extend(genes_sorted)

mean_expr_scaled_sorted = mean_expr_scaled.loc[sorted_genes]

# Plot heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(
    mean_expr_scaled_sorted,
    cmap='vlag',
    yticklabels=True,
    xticklabels=True,
    cbar_kws={'label': 'Z-score of mean expression'},
    center=0
)
plt.title("Top marker genes per cell type (Wilcoxon + mean_diff)", fontsize=16)
plt.ylabel("Genes")
plt.xlabel("Cell types")
plt.tight_layout()
heatmap_file = f"{base_name}_NoFC_diagonal_heatmap_{min_in}.png"
plt.savefig(heatmap_file, dpi=600, bbox_inches='tight')
plt.close()
print(f"Diagonal heatmap saved to {heatmap_file}")



