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

sc.tl.rank_genes_groups(adata, groupby="celltype")
df_all = sc.get.rank_genes_groups_df(adata, group=None)
print(df_all[['group', 'names', 'logfoldchanges', 'pvals_adj']].head())
print("Available columns in DE results:", df_all.columns.tolist())


groups = df_all['group'].unique().tolist()
# Loop over each group and print all ranked genes sorted by logFC
for group in groups:
    print(f"\nAll DE genes in group '{group}' sorted by pvals_adj (highest to lowest):")
    group_df = df_all[df_all['group'] == group].sort_values('pvals_adj')
    print(group_df[['names', 'logfoldchanges', 'pvals_adj']])

df_all.to_csv(f"ranked_genes_{base_name}_ALL_noExp_default.csv", index=False)

# Get top N genes per group (by absolute logfc or rank)
top_n = 10
top_genes_combined = []

for group in groups:
    top_genes = (
        df_all[df_all['group'] == group]
        .sort_values('pvals_adj',ascending=True)  # Optional: sort by abs(logFC)
        .head(top_n)['names']
        .tolist()
    )
    top_genes_combined.extend(top_genes)

top_genes_combined = list(dict.fromkeys(top_genes_combined))
print("Top genes for heatmap:", top_genes_combined)

celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
adata.obs['celltype'] = adata.obs['celltype'].astype(str)
adata.obs['celltype'] = pd.Categorical(
    adata.obs['celltype'],
    categories=celltype_order,
    ordered=True
)
print("Confirmed category order:", adata.obs['celltype'].cat.categories)

# Plot heatmap
sc.pl.rank_genes_groups_heatmap(
    adata,
    groups=celltype_order,
    groupby='celltype',
    n_genes=top_n,
    swap_axes=True,
    use_raw=False,
    cmap='bwr',
    dendrogram=False,
    save=f"_{base_name}_Top{top_n}Genes_all_clusterDefault.png"
)


#adata.obs_names_make_unique()
#adata.write(newObject, compression="gzip")
