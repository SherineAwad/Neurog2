import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd


# Argument parsing
parser = argparse.ArgumentParser(description="Per-sample UMAPs and stacked bar plots")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

# Read AnnData object
adata = sc.read(myObject)
print("✅ Columns in .obs:", adata.obs.columns)
print(adata.obs.head())


# UMAP plots for each sample
samples = adata.obs['sample'].unique()
for sample in samples:
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color='sample',
        title=f"Sample: {sample}",
        size=5,
        save=f"_reclustered_{sample}.png",
        show=False
    )

# UMAP colored by sample
sc.pl.umap(adata, color='sample', size=5, save=f"reclustered_{base_name}.png")

celltype_colors = {
    'Cones': '#e31a1c',       # Red
    'MG': '#0C727C',       # Greenish Turquoise  
    'AC': '#026AB1',    # Blue 
    'MGPC': '#6E4B9E',     # Purple
    'BC': '#8A9FD1',       # Light Purple
    'Rod': '#bdbdbd'       # Grey
}

# Create dataframe of counts
df = (
    adata.obs[['celltype', 'sample']]
    .value_counts()
    .reset_index(name='count')
)

# Normalize to get proportions per sample
df['fraction'] = df['count'] / df.groupby('sample')['count'].transform('sum')

# Pivot for plotting: index = sample, columns = celltype
pivot_df = df.pivot(index='sample', columns='celltype', values='fraction').fillna(0)

# Print sample names to confirm
print("✅ Available sample names:", pivot_df.index.tolist())

# Manually define sample and celltype order
sample_order = ['control_2mo', 'Neurog2_9SA_5weeks', 'Neurog2_9SA_2mo']
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
short_sample_names = ['Control', '5weeks', '2months']

# Check and reorder samples
missing_samples = [s for s in sample_order if s not in pivot_df.index]
if missing_samples:
    raise ValueError(f"⚠️ Sample(s) not found in data: {missing_samples}")
pivot_df = pivot_df.loc[sample_order]

# Check and reorder celltypes
missing_celltypes = [ct for ct in celltype_order if ct not in pivot_df.columns]
if missing_celltypes:
    raise ValueError(f"⚠️ Cell type(s) not found in data: {missing_celltypes}")
pivot_df = pivot_df[celltype_order]

# Use the correct colors (aligned with order)
colors = [celltype_colors[ct] for ct in celltype_order]

# Plot the stacked bar chart
fig, ax = plt.subplots(figsize=(5, 4))
pivot_df.plot(kind='bar', stacked=True, color=colors, ax=ax,width=0.7)
plt.ylabel("Fraction of cells",fontsize=12)
plt.xlabel("Sample")
plt.title("Cell Type Contribution per Sample", fontsize=12)
ax.set_xticklabels(short_sample_names, rotation=0, ha='center', fontsize=10)
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save figure
plt.savefig("figures/Reversed_stacked_bar_sample_by_celltype.png", dpi=300)
plt.close()


# Make sure 'celltype' is a categorical column with this order
adata.obs["celltype"] = pd.Categorical(
    adata.obs["celltype"],
    categories=celltype_order,
    ordered=True
)

marker_genes = {
 "MG": ["Rlbp1","Slc1a3"],
  "Rod": ["Nrl", "Rho"],
  "Cones": ["Arr3","Gnat2"],
  "BC": ["Otx2", "Cabp5", "Scgn"] ,
  "AC": ["Elavl3", "Gad2", "Chat"],
  "MGPC": ["Ascl1"]
  }

import pandas as pd

# Ensure celltype column is categorical with desired order
celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'], categories=celltype_order, ordered=True)

markergenes = [
    'Rlbp1', 'Slc1a3', 'Ascl1', 'Otx2', 'Cabp5', 'Scgn',
    'Elavl3', 'Gad2', 'Chat', 'Nrl', 'Rho', 'Arr3', 'Gnat2'
]

sc.pl.dotplot(
    adata,
    markergenes,
    groupby="celltype",
    categories_order = celltype_order,
    standard_scale="var",
    figsize=(6,5),
    dot_max=1.0, show=False, dendrogram=False, save= f"_{base_name}_markerGenes.png") 





