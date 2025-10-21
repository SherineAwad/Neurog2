import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import os

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('markers')

args = parser.parse_args()

myObject = args.myObject
newObject = "clustered_" + myObject

markers = args.markers
print(markers)


base_name = os.path.splitext(os.path.basename(newObject))[0]


combined_adata = sc.read(myObject)

figure_name = f"_{base_name}_ddClusters.png"

sc.tl.leiden(combined_adata, resolution=1.0)
sc.pl.umap(combined_adata, color=["leiden"], save=figure_name, legend_loc="on data")



# Read marker genes as a list of gene names (one per line)
marker_genes = [line.strip() for line in open(markers)]

# Use the marker_genes list directly (no .items() since it's a list)
all_genes = marker_genes

# Save individual plots for each gene
for gene in marker_genes:
    if gene in combined_adata.var_names:
        sc.pl.scatter(
            combined_adata,
            color=gene,
            title=gene,
            basis='umap',
            save=f"_{base_name}_{gene}.png"
        )

combined_adata.obs_names_make_unique()
combined_adata.write(newObject, compression="gzip")

