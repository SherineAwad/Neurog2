import os
import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
import argparse
import importlib_metadata
import sys
import matplotlib.pyplot as plt

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser(description="Run DE analysis without NaNs in logFC")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

# Read data
adata = sc.read(myObject)

# Ensure raw counts are available
if adata.raw is None:
    adata.raw = adata.copy()

min_frac = 0.10
cell_fraction = (adata.X > 0).sum(axis=0) / adata.n_obs
if not isinstance(cell_fraction, np.ndarray):
    cell_fraction = np.array(cell_fraction).flatten()  # for sparse matrices

genes_to_keep = cell_fraction >= min_frac
adata = adata[:, genes_to_keep]
print(f"Filtered to {adata.n_vars} genes expressed in at least {min_frac*100:.0f}% of cells.")

# Run DE using raw counts to avoid NaN logFC
sc.tl.rank_genes_groups(
    adata,
    groupby="celltype",
    use_raw=True,
    n_genes=adata.raw.shape[1]
)

# Extract DE results
df_all = sc.get.rank_genes_groups_df(adata, group=None)

#['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']
print(df_all[['group', 'names', 'scores', 'logfoldchanges', 'pvals','pvals_adj']].head())
print("Available columns in DE results:", df_all.columns.tolist())

# Show group-wise top genes
groups = df_all['group'].unique().tolist()
for group in groups:
    print(f"\nAll DE genes in group '{group}' sorted by scores (highest to lowest):")
    group_df = df_all[df_all['group'] == group].sort_values('scores', ascending=False)
    print(group_df[['names', 'scores','logfoldchanges', 'pvals_adj']])

# Save to CSV
df_all.to_csv(f"ranked_genes_{base_name}_filtered.csv", index=False)


# Get top N genes per group by absolute log fold change
top_n = 50
top_genes_combined = []

for group in groups:
    top_genes = (
        df_all[df_all['group'] == group]
        .sort_values('scores', key=abs, ascending=False)
        .head(top_n)['names']
        .tolist()
    )
    top_genes_combined.extend(top_genes)

# Remove duplicates while preserving order
top_genes_combined = list(dict.fromkeys(top_genes_combined))
print("Top genes for heatmap:", top_genes_combined)

# Set categorical ordering for cell types
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
adata.obs['celltype'] = adata.obs['celltype'].astype(str)
adata.obs['celltype'] = pd.Categorical(
    adata.obs['celltype'],
    categories=celltype_order,
    ordered=True
)
print("Confirmed category order:", adata.obs['celltype'].cat.categories)

# Filter top genes for heatmap to only those present in .raw.var_names
top_genes_filtered = [g for g in top_genes_combined if g in adata.raw.var_names]

# Take only top N genes after filtering
top_genes_filtered = top_genes_filtered[:top_n]
print(f"Top {top_n} genes for heatmap (present in raw): {top_genes_filtered}")

# Plot heatmap
sc.pl.rank_genes_groups_heatmap(
    adata,
    groups=celltype_order,
    groupby='celltype',
    var_names=top_genes_filtered,
    swap_axes=True,
    use_raw=True,
    dendrogram=False,
    cmap='bwr',
    show_gene_labels=True,
    save=f"_{base_name}_Top{top_n}Genes_fheatmap.png"
)


#adata_new.obs_names_make_unique()
#adata_new.write(newObject, compression="gzip")


