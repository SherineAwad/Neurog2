import scanpy as sc
import pandas as pd
import argparse
import os
import sys
import importlib_metadata

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

# Argument parsing
parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "expression_" + base_name + ".h5ad"

adata = sc.read(myObject)

# 1. Check your annotation column
annotation_col = 'celltype'  # change if different
if annotation_col not in adata.obs:
    raise ValueError(f"Column '{annotation_col}' not found in adata.obs")

# 2. Run differential expression analysis
sc.tl.rank_genes_groups(adata, groupby=annotation_col, method='wilcoxon')

# 3. Visualize top marker genes per cluster (optional)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# 4. Extract DE results into dict of DataFrames
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names  # list of clusters/groups

marker_genes = {}
for group in groups:
    marker_genes[group] = pd.DataFrame({
        'gene': result['names'][group],
        'logfoldchange': result['logfoldchanges'][group],
        'pvals_adj': result['pvals_adj'][group],
        'scores': result['scores'][group]
    })

# 5. Collect top N genes per cluster
top_n = 3
top_genes_combined = []
for group in groups:
    top_genes_combined.extend(marker_genes[group]['gene'].values[:top_n].tolist())
top_genes_combined = list(dict.fromkeys(top_genes_combined))  # unique list, preserve order

# 6. Plot built-in heatmap for top N genes in all clusters
sc.pl.rank_genes_groups_heatmap(
    adata, 
    groups=list(groups),  # all clusters
    n_genes=top_n,
    swap_axes=True,
    show=True,
    save=f"_heatmap_{base_name}_Top{top_n}Genes_all_clusters.png"
)

# Optional: Custom heatmap of combined unique genes across clusters
# sc.pl.heatmap(
#     adata,
#     var_names=top_genes_combined,
#     groupby=annotation_col,
#     swap_axes=True,
#     show=True,
#     save=f"_custom_heatmap_{base_name}_Top{top_n}Genes_all_clusters.png"
# )

adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")

