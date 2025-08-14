import anndata as ad

# Load original h5ad
adata = ad.read_h5ad("annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad")

# Convert categorical columns to string
for col in adata.obs.columns:
    if str(adata.obs[col].dtype) == 'category':
        adata.obs[col] = adata.obs[col].astype(str)

# Optional: drop extra metadata that may confuse SeuratDisk
adata.uns = {}
adata.layers = {}

# Ensure var index is string and unique
adata.var.index = [str(x) for x in adata.var.index]

# Save cleaned h5ad
adata.write_h5ad("neurog2_clean.h5ad")



