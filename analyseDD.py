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
newObject = "ddanalysed_" + myObject
base_name = os.path.splitext(os.path.basename(newObject))[0]


adata = sc.read(myObject)

# Normalize total counts per cell to 10,000 and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

# Scale the data (clip values to max 10)
sc.pp.scale(adata, max_value=10)

# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Compute neighborhood graph
sc.pp.neighbors(adata)

# Compute UMAP embedding
sc.tl.umap(adata)

# Plot UMAP colored by sample
sc.pl.umap(adata, color='sample', size=2, save=f"_{base_name}.png") 



cell_counts = adata.obs['sample'].value_counts()
print(cell_counts)

# Plot UMAP by sample 
samples = adata.obs['sample'].unique()
# Loop through samples and plot individual UMAPs
for sample in samples:
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color='sample',
        title=f"Sample: {sample}",
        size=20,
        save=f"_sample_{sample}.png",
        show=False
    )


sc.pl.umap(adata, color=["predicted_doublet", "doublet_score"], save=f"_{base_name}_doublets.png")


# Save the processed object
adata.write(newObject, compression="gzip")
