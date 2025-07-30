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
        size=10,
        save=f"_reclustered_{sample}.png",
        show=False
    )

# UMAP colored by sample
sc.pl.umap(adata, color='sample', size=10, save=f"reclustered_{base_name}.png")

celltype_colors = {
    'AC': '#e31a1c',       # Red
    'MG': '#0C727C',       # Greenish Turquoise
    'Cones': '#026AB1',    # Blue
    'MGPC': '#9467bd',     # Purple
    'BC': '#c2a5cf',       # Light Purple
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
fig, ax = plt.subplots(figsize=(6, 6))
pivot_df.plot(kind='bar', stacked=True, color=colors, ax=ax,width=0.3)
plt.ylabel("Fraction of cells")
plt.xlabel("Sample")
plt.title("Cell Type Contribution per Sample")
plt.xticks(rotation=45, ha='right')
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

marker_genes  = {
    "MG": ["Rlbp1","Slc1a3"],
    "Rod": ["Rho","Nrl"],
    "Cones": ["Arr3","Gnat2"],
    "BC": ["Cabp5","Otx2", "Scgn"] ,
    "AC": ["Chat","Elavl3", "Gad2"],
    "MGPC": ["Ascl1"]
    }


sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="celltype",
    swap_axes=True,
    standard_scale="var",
    figsize=(8, 6),
    show=False, dendrogram=False, save= f"_{base_name}_markerGenes.png") 

plt.tight_layout()




