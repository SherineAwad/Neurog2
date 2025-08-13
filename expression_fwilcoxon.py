import scanpy as sc
import pandas as pd
import argparse
import os
import sys
import importlib_metadata
import numpy as np
import scipy.sparse

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "expression_" + base_name + ".h5ad"

adata = sc.read(myObject)

# Filter genes expressed in at least 10% of cells (based on log1p data)
min_frac = 0.10
cell_fraction = (adata.X > 0).sum(axis=0) / adata.n_obs
if not isinstance(cell_fraction, np.ndarray):
    cell_fraction = np.array(cell_fraction).flatten()
genes_to_keep = cell_fraction >= min_frac
adata = adata[:, genes_to_keep].copy()  # copy after slicing

# Create raw counts (unlogged) AnnData and assign to adata.raw
if scipy.sparse.issparse(adata.X):
    X_unlogged = np.expm1(adata.X.toarray())
else:
    X_unlogged = np.expm1(adata.X)

raw_adata = sc.AnnData(X=X_unlogged, obs=adata.obs.copy(), var=adata.var.copy())
adata.raw = raw_adata

# Run DE with Wilcoxon test using unlogged raw counts
sc.tl.rank_genes_groups(adata, groupby="celltype", method='wilcoxon', use_raw=True)

# Extract DE results
df_all = sc.get.rank_genes_groups_df(adata, group=None)

print(df_all[['group', 'names', 'logfoldchanges', 'pvals_adj']].head())
print("Available columns in DE results:", df_all.columns.tolist())

groups = df_all['group'].unique().tolist()

# Loop over groups and print DE genes sorted by logFC
for group in groups:
    print(f"\nAll DE genes in group '{group}' sorted by logFC (highest to lowest):")
    group_df = df_all[df_all['group'] == group].sort_values('logfoldchanges', ascending=False)
    print(group_df[['names', 'logfoldchanges', 'pvals_adj']])

df_all.to_csv(f"ranked_genes_{base_name}_fwl.csv", index=False)

# Get top N genes per group (by absolute logfc or rank)
top_n = 20
top_genes_combined = []

for group in groups:
    top_genes = (
        df_all[df_all['group'] == group]
        .sort_values('logfoldchanges', key=abs, ascending=False)  # Optional: sort by abs(logFC)
        .head(top_n)['names']
        .tolist()
    )
    top_genes_combined.extend(top_genes)

top_genes_combined = list(dict.fromkeys(top_genes_combined))
print("Top genes for heatmap:", top_genes_combined)
adata.raw = adata  # assuming `adata` holds normalized + log1p


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

# Take only top N genes after filtering (optional)
top_genes_filtered = top_genes_filtered[:top_n]

print(f"Top {top_n} genes for heatmap (present in raw): {top_genes_filtered}")


sc.pl.rank_genes_groups_heatmap(
    adata,
    groups=celltype_order,
    groupby='celltype',
    var_names=top_genes_filtered,
    swap_axes=False,            # try turning off swap_axes
    use_raw=True,
    dendrogram=False,
    cmap='bwr',
    show_gene_labels=True,      # keep True, but might need other fixes
    gene_symbols=None,          # explicitly specify if you have a gene_symbols column in adata.var
    figsize=(10,8),             # increase figure size so labels fit
    save=f"_{base_name}_Top{top_n}Genes_fwl.png"
)

adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")
