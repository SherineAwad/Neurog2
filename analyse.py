import scanpy as sc
import sys
import importlib_metadata
import argparse
import os

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
newObject = "analysed_" + myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]


combined_adata = sc.read(myObject)

# Normalize total counts per cell to 10,000 and log-transform
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(combined_adata, flavor='seurat', n_top_genes=2000)

# Scale the data (clip values to max 10)
sc.pp.scale(combined_adata, max_value=10)

# Perform PCA
sc.tl.pca(combined_adata, svd_solver='arpack')

# Compute neighborhood graph
sc.pp.neighbors(combined_adata)

# Compute UMAP embedding
sc.tl.umap(combined_adata)

# Plot UMAP colored by sample
sc.pl.umap(combined_adata, color='sample', size=2, save=f"_{base_name}.png") 

# Save the processed object
combined_adata.write(newObject)


