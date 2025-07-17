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
parser.add_argument('threshold', type=float, default=0.5, help='Doublet score threshold for filtering')

args = parser.parse_args()

markers = args.markers
threshold = args.threshold
print(f"Markers file: {markers}")
print(f"Doublet score threshold: {threshold}")

myObject = args.myObject
newObject = f"doubletsRemoved_threshold{threshold}_" + myObject

base_name = os.path.splitext(os.path.basename(newObject))[0]

# 🔹 Read in the AnnData object
adata = sc.read(myObject)


adata = adata[adata.obs['predicted_doublet'] <= threshold].copy()
print(f"✅ Kept cells with predicted_doublet <= {threshold}")



# 🔹 Save UMAP plot of new clustering
sc.pl.umap(
    adata, 
    color=["leiden"], 
    save=f"_{base_name}_doubletsRemoved_threshold{threshold}_clusters.png", 
    legend_loc="on data"
)

# Read marker genes as a list of gene names (one per line)
marker_genes = [line.strip() for line in open(markers)]

# Use the marker_genes list directly (no .items() since it's a list)
for gene in marker_genes:
    if gene in adata.var_names:
        sc.pl.scatter(
            adata,
            color=gene,
            title=gene,
            basis='umap',
            save=f"_{base_name}_threshold{threshold}_{gene}.png"
        )

adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")

