import scanpy as sc
import pandas as pd
import argparse
import os
import sys
import importlib_metadata

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "expression_" + base_name + ".h5ad"


adata = sc.read(myObject)

import numpy as np 
adata_new = adata.copy()
adata_unlogged = adata.copy()
adata_unlogged.X = np.expm1(adata_unlogged.X)  # reverse log1p: exp(x) - 1
adata_new.raw = adata_unlogged

# Step 1: Define control vs treatment groups
adata_new.obs['comparison_group'] = adata_new.obs['sample'].map(
    lambda x: 'control' if x == 'control_2mo' else (
        'treatment' if x in ['Neurog2_9SA_5weeks', 'Neurog2_9SA_2mo'] else 'other'
    )
)

# Filter to only control and treatment samples
adata_new = adata_new[adata_new.obs['comparison_group'].isin(['control', 'treatment'])].copy()

sc.tl.rank_genes_groups(adata_new, groupby="comparison_group", method='t-test')
df_all = sc.get.rank_genes_groups_df(adata_new, group=None)
print(df_all[['group', 'names', 'logfoldchanges', 'pvals_adj']].head())
print("Available columns in DE results:", df_all.columns.tolist())


groups = df_all['group'].unique().tolist()
# Loop over each group and print all ranked genes sorted by logFC
for group in groups:
    print(f"\nAll DE genes in group '{group}' sorted by logFC (highest to lowest):")
    group_df = df_all[df_all['group'] == group].sort_values('logfoldchanges', ascending=False)
    print(group_df[['names', 'logfoldchanges', 'pvals_adj']])

df_all.to_csv(f"ranked_genes_{base_name}_ALLttG.csv", index=False)

# Get top N genes per group (by absolute logfc or rank)
top_n = 10
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
adata_new.raw = adata  # assuming `adata` holds normalized + log1p



# Plot heatmap
sc.pl.rank_genes_groups_heatmap(
    adata_new,
    groupby="comparison_group",
    n_genes=top_n,
    swap_axes=True,
    use_raw=True,
    dendrogram=False,
    cmap='bwr',
    save=f"_{base_name}_Top{top_n}Genes_all_clusterttestG.png"
)

#adata_new.obs_names_make_unique()
#adata_new.write(newObject, compression="gzip")
