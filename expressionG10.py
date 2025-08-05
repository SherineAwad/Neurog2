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

# Prepare unlogged data with clipping to avoid overflow in expm1
adata_unlogged = adata.copy()
adata_unlogged.X = np.clip(adata_unlogged.X, a_min=None, a_max=20)
adata_unlogged.X = np.expm1(adata_unlogged.X)

# Add raw to the full data before filtering anything (important!)
adata.raw = adata_unlogged

# Assign groups
adata.obs['comparison_group'] = adata.obs['sample'].map(
    lambda x: 'control' if x == 'control_2mo' else (
        'treatment' if x in ['Neurog2_9SA_5weeks', 'Neurog2_9SA_2mo'] else 'other'
    )
)

# Filter for control/treatment cells only
adata = adata[adata.obs['comparison_group'].isin(['control', 'treatment'])].copy()

# Filter for MG only
adata = adata[adata.obs['celltype'] == 'MG'].copy()

# Calculate gene expression fraction per group
def gene_frac_expr(adata_group):
    X = adata_group.X
    if scipy.sparse.issparse(X):
        return np.asarray((X > 0).sum(axis=0)).ravel() / X.shape[0]
    else:
        return (X > 0).sum(axis=0) / X.shape[0]

min_frac = 0.1
group1 = adata[adata.obs['comparison_group'] == 'control']
group2 = adata[adata.obs['comparison_group'] == 'treatment']

frac1 = gene_frac_expr(group1)
frac2 = gene_frac_expr(group2)

gene_mask = (frac1 >= min_frac) | (frac2 >= min_frac)

# Subset adata by genes (cells already filtered)
adata = adata[:, gene_mask].copy()

# Now **re-assign** .raw to be the filtered subset of adata_unlogged (must be AnnData)
adata.raw = adata_unlogged[adata.obs_names, adata.var_names].copy()

if "log1p" not in adata.uns or "base" not in adata.uns["log1p"]:
    adata.uns["log1p"] = {"base": None}

# Run DE
sc.tl.rank_genes_groups(adata, groupby='comparison_group', method='t-test')
df_all = sc.get.rank_genes_groups_df(adata,group=None)

print(df_all[['group', 'names', 'logfoldchanges', 'pvals_adj']].head())

# Get top N genes per group (by absolute logfc or rank)
top_n = 10
top_genes_combined = []

groups = list(adata.uns['rank_genes_groups']['names'].dtype.names)

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

# Plot heatmap
sc.pl.rank_genes_groups_heatmap(
    adata,
    groupby="comparison_group",
    n_genes=top_n,
    swap_axes=True,
    use_raw=False,
    dendrogram=False,
    cmap='bwr',
    save=f"_{base_name}_Top{top_n}Genes_all_clusterttestG10.png"
)

#adata_new.obs_names_make_unique()
#adata_new.write(newObject, compression="gzip")
