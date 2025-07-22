import argparse
import os
import sys
import scanpy as sc
import importlib_metadata
import matplotlib.pyplot as plt
import numpy as np


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
clusters_to_remove = ['15','24','7', '8', '18', '27', '32', '33']
# Check for clustering column
if 'leiden' not in adata.obs:
    raise ValueError("No 'leiden' column found in .obs. Please ensure clustering was previously run.")

#Print cluster ID before removal 

original_clusters = adata.obs['leiden'].unique().tolist()
original_clusters.sort()  # Optional: sort for easier reading
print("Original cluster IDs:", original_clusters)

# Filter out unwanted clusters
initial_cells = adata.n_obs
adata = adata[~adata.obs['leiden'].isin(clusters_to_remove)].copy()


# Print remaining cluster IDs
remaining_clusters = adata.obs['leiden'].unique().tolist()
remaining_clusters.sort()  # Optional: sort for easier reading
print("Remaining cluster IDs:", remaining_clusters)


## Filter genes with zero counts in all cells
sc.pp.filter_genes(adata, min_cells=1)

## Scale
sc.pp.scale(adata)

# Recompute PCA on filtered data
sc.tl.pca(adata) #, svd_solver='arpack')

# Recompute neighborhood graph
sc.pp.neighbors(adata)  #, n_neighbors=10, n_pcs=40)

resol = 1.4 
# Reclustering
sc.tl.leiden(adata,flavor="igraph", n_iterations=2,  resolution=resol)

sc.pl.umap(adata, color=["leiden"], save=f"_reClustered_{base_name}_{resol}res.png", legend_loc="on data")

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



marker_genes  = {
    "MG": ["Rlbp1","Gfap","Apoe","Notch1","Pax6","Slc1a3","Vim"],
    "Rod": ["Rho","Nrl","Crx","Rom1"],
    "Cones": ["Opn1mw","Opn1sw","Arr3","Thrb","Gnat2"],
    "BC": ["Vsx1", "Sebox","Bhlhe23","Cabp5","Vsx1","Pcp4","Isl1"] ,
    "AC": ["Gad1","Gad2","Slc6a9","Tfap2b","Prox1","Pax6","Calb2","Pcp4","Elavl3","Isl1"],
    "HC": ["Lhx1","Cbln4","Calb1","Nefl","Nefm", "Onecut1", "Onecut2"],
    "RGC": ["Nefl","Nefm","Sncg","Thy1","Ebf3","Rbfox3","Isl1","Isl2","Pou4f1","Pou4f3","Rbpms"],
    "Microglia": ["Ptprc","Csf2rb","Sall1"],
    "Astrocytes":["Pax2","Igf2", "Gfap"]
    }
sc.pl.dotplot(adata, marker_genes, groupby="leiden", standard_scale="var", save=f"_reClustered_{base_name}_markerGenes.png")



adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")

