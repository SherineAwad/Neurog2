import argparse
import os
import sys
import scanpy as sc
import importlib_metadata

# Argument parsing
parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('annotations', help='Path to annotations file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "annotated_" + base_name + ".h5ad"
sys.modules['importlib.metadata'] = importlib_metadata

annot_file = args.annotations

adata = sc.read_h5ad(myObject, backed="r")

print(adata.obs['leiden'].value_counts().sort_index())

# Load annotations
cluster_to_celltype_dict = {}
with open(annot_file, "r") as f:
    for line in f:
        cluster, celltype = line.strip().split('\t')
        cluster_to_celltype_dict[cluster] = celltype

# Ensure keys are strings
cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}

# Map cell types
adata.obs["celltype"] = adata.obs["leiden"].map(cluster_to_celltype_dict)

# Set celltype colors
celltype_colors = {
    'AC': '#e31a1c',       # Red
    'MG': '#0C727C',       # Greenish Turquoise  
    'Cones': '#026AB1',    # Blue 
    'MGPC': '#9467bd',     # Purple
    'BC': '#c2a5cf',       # Light Purple
    'Rod': '#bdbdbd'       # Grey
}

# Reorder colors to match celltype categories
celltypes_in_data = adata.obs["celltype"].astype('category').cat.categories
adata.uns["celltype_colors"] = [celltype_colors.get(ct, "#000000") for ct in celltypes_in_data]  # default to black if missing

# Make obs names unique
adata.obs_names_make_unique()

# Generate plots
figure_name = base_name + "_annotationsON.png"
sc.pl.umap(adata, color='celltype', legend_loc="on data", save=figure_name)

figure_name = base_name + "_annotations.png"
sc.pl.umap(adata, color='celltype', save=figure_name)

import pandas as pd

# Filter cells with missing celltype
unannotated = adata.obs[adata.obs['celltype'].isna()]

# Get unique cluster numbers (leiden) with missing annotations
clusters_with_na = unannotated['leiden'].unique()
print("🟡 Clusters with missing annotations:", sorted(clusters_with_na))


# Save updated object
adata.write(newObject, compression="gzip")

