import os
import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
import argparse
import importlib_metadata
import sys
import matplotlib.pyplot as plt

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser(description="Run DE analysis without NaNs in logFC")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

# Read data
adata = sc.read(myObject)
# If raw is missing, create it from current adata
if adata.raw is None:
    adata.raw = adata.copy()

# Create a new AnnData object with reverted counts from raw (expm1)
X_raw = adata.raw.X
if scipy.sparse.issparse(X_raw):
    X_raw = X_raw.toarray()
counts = np.expm1(X_raw)  # revert log1p to counts

# Make a new AnnData for counts only (keep original obs/var)
adata_counts = sc.AnnData(X=counts, obs=adata.obs.copy(), var=adata.raw.var.copy())

# Filter genes expressed in >=10% cells on counts data
min_frac = 0.10
cell_fraction = (adata_counts.X > 0).sum(axis=0) / adata_counts.n_obs
genes_to_keep = cell_fraction >= min_frac
adata_counts = adata_counts[:, genes_to_keep]
print(f"Filtered to {adata_counts.n_vars} genes expressed in at least {min_frac*100:.0f}% of cells.")

# Run DE with Wilcoxon on reverted counts
sc.tl.rank_genes_groups(
    adata_counts,
    groupby="celltype",
    method="wilcoxon",
    use_raw=False,
    n_genes=adata_counts.n_vars
)

# Extract DE results
df_all = sc.get.rank_genes_groups_df(adata_counts, group=None)
print(df_all[['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']].head())
print("Available columns in DE results:", df_all.columns.tolist())

# Save to CSV
df_all.to_csv(f"ranked_genes_{base_name}_wilcoxon_counts.csv", index=False)


groups = df_all['group'].unique().tolist()


# Assuming 'df_all' and 'groups' are from DE done on adata_counts

top_n = 50
top_genes_combined = []

for group in groups:
    top_genes = (
        df_all[df_all['group'] == group]
        .sort_values('logfoldchanges', key=abs, ascending=False)
        .head(top_n)['names']
        .tolist()
    )
    top_genes_combined.extend(top_genes)

# Remove duplicates but preserve order
top_genes_combined = list(dict.fromkeys(top_genes_combined))
print("Top genes for heatmap:", top_genes_combined)

# Set categorical ordering for cell types in original adata.obs
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
adata.obs['celltype'] = adata.obs['celltype'].astype(str)
adata.obs['celltype'] = pd.Categorical(
    adata.obs['celltype'],
    categories=celltype_order,
    ordered=True
)
print("Confirmed category order:", adata.obs['celltype'].cat.categories)

# Filter top genes to those present in original adata.raw.var_names (log1p data)
top_genes_filtered = [g for g in top_genes_combined if g in adata.raw.var_names]

# Take only top_n genes after filtering
top_genes_filtered = top_genes_filtered[:top_n]
print(f"Top {top_n} genes for heatmap (present in raw): {top_genes_filtered}")

# Plot heatmap on original adata using raw data (log1p)
sc.pl.rank_genes_groups_heatmap(
    adata,
    groups=celltype_order,
    groupby='celltype',
    var_names=top_genes_filtered,
    swap_axes=True,
    use_raw=True,
    dendrogram=False,
    cmap='bwr',
    show_gene_labels=True,
    save=f"_{base_name}_Top{top_n}Genes_fwlheatmap.png"
)


#adata_new.obs_names_make_unique()
#adata_new.write(newObject, compression="gzip")


