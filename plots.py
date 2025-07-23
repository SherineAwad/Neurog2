import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse

# === 1. Parse input ===
parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()



myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "reclustered_" + base_name + ".h5ad"

adata = sc.read(myObject)
print(adata.obs.columns)
print(adata.obs.head())
list(adata.obs.columns)

'''
cell_counts = adata.obs['sample'].value_counts()
print(cell_counts)
print(adata.obs.columns)

sc.pp.pca(adata)               # if not already done
sc.pp.neighbors(adata)        # builds the neighborhood graph
sc.tl.leiden(adata)           # now Leiden will work
sc.tl.umap(adata)

figure_name = "After_"+{base_name}+".png"
sc.pl.umap(adata, color=["leiden"], save=figure_name, legend_loc="on data")


figurename1 = "figures/"+base_name+"qc_violin.png"
sc.pl.violin(
    adata,
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='leiden',     # change if you use another cluster label
    jitter=0.4,
    rotation=45,
    multi_panel=True,
    show=False
)
plt.savefig(figurename1, dpi=300, bbox_inches="tight")
plt.close()

'''


# Get unique sample values
samples = adata.obs['sample'].unique()

# Loop through samples and plot individual UMAPs
for sample in samples:
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color='sample',
        title=f"Sample: {sample}",
        size=20,
        save=f"_reclustered_{sample}.png",
        show=False
    )
# Plot UMAP colored by sample
sc.pl.umap(adata, color='sample', size=2, save=f"reclustered_{base_name}.png")

import matplotlib.pyplot as plt

# Create dataframe of counts
df = (
    adata.obs[['celltype', 'sample']]
    .value_counts()
    .reset_index(name='count')
)

# Normalize to get proportions *per sample*
df['fraction'] = df['count'] / df.groupby('sample')['count'].transform('sum')

# Pivot for plotting: index = sample, columns = celltype
pivot_df = df.pivot(index='sample', columns='celltype', values='fraction').fillna(0)

# Print available sample names to confirm
print("✅ Available sample names:", pivot_df.index.tolist())

# Manually order samples based on actual names
sample_order = ['control_2mo', 'Neurog2_9SA_5weeks', 'Neurog2_9SA_2mo']  # Adjust as needed
missing_samples = [s for s in sample_order if s not in pivot_df.index]
if missing_samples:
    raise ValueError(f"⚠️ Sample(s) not found in data: {missing_samples}")
pivot_df = pivot_df.loc[sample_order]

# Manually reorder columns to desired celltype order
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
missing_celltypes = [ct for ct in celltype_order if ct not in pivot_df.columns]
if missing_celltypes:
    raise ValueError(f"⚠️ Cell type(s) not found in data: {missing_celltypes}")
pivot_df = pivot_df[celltype_order]

# Define custom colors matching your UMAP plot
celltype_colors = {
    'AC': '#1f77b4',       # Blue
    'BC': '#ff7f0e',       # Orange
    'Cones': '#2ca02c',    # Green
    'MG': '#d62728',       # Red
    'MGPC': '#9467bd',     # Purple
    'Rod': '#8c564b'       # Brown
}

# Extract colors in order
colors = [celltype_colors[ct] for ct in celltype_order]

# Plot
ax = pivot_df.plot(kind='bar', stacked=True, figsize=(12, 6), color=colors)
plt.ylabel("Fraction of cells")
plt.xlabel("Sample")
plt.title("Cell Type Contribution per Sample")
plt.xticks(rotation=45, ha='right')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save figure
plt.savefig("figures/Reversed_stacked_bar_sample_by_celltype.png", dpi=300)
plt.close()




'''
umap_fig_name = base_name + "_umap_doublets.png"
# Option 1: Convert to string labels
adata.obs['doublet_label'] = adata.obs['predicted_doublets'].map({True: "Doublet", False: "Singlet"})

sc.pl.umap(
    adata,
    color='doublet_label',
    title='Predicted Doublets on UMAP',
    save="_doublets.png"
)


'''


