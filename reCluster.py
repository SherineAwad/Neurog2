import argparse
import os
import sys
import scanpy as sc
import importlib_metadata
import matplotlib.pyplot as plt


# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

# Argument parsing
parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('markers', help='Path to marker genes file (not used in this script but required for CLI compatibility)')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "reclustered_" + base_name + ".h5ad"

markers = args.markers
print(markers)

# Load data
print(f"Loading AnnData from: {myObject}")
adata = sc.read(myObject)


# Clusters to remove
clusters_to_remove = ['7', '8', '18', '27', '32', '33']
# Check for clustering column
if 'leiden' not in adata.obs:
    raise ValueError("No 'leiden' column found in .obs. Please ensure clustering was previously run.")

# Filter out unwanted clusters
initial_cells = adata.n_obs
adata = adata[~adata.obs['leiden'].isin(clusters_to_remove)].copy()

# Recompute PCA on filtered data
sc.tl.pca(adata, svd_solver='arpack')

# Recompute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Reclustering
sc.tl.leiden(adata, resolution=1.0)

sc.pl.umap(adata, color=["leiden"], save=f"_reClustered_{base_name}.png", legend_loc="on data")




marker_genes = [line.strip() for line in open(markers)]


# Save individual plots for each gene
for gene in marker_genes:
    if gene in adata.var_names:
        sc.pl.scatter(
            adata,
            color=gene,
            title=gene,
            basis='umap',
            save=f"_reClustered_{base_name}_{gene}.png"
        )

adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")

