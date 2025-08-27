import os
import scanpy as sc
import argparse
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('myObject')

args = parser.parse_args()

myObject = args.myObject

orig_base = os.path.splitext(os.path.basename(myObject))[0]
parts = orig_base.split('_clustered_ddanalysed_doubletScores_0.8')

base_name = os.path.splitext(os.path.basename(myObject))[0]
base_name = "refined_"+ base_name 

newObject = base_name + ".h5ad"

adata = sc.read(myObject)

cell_counts = adata.obs['sample'].value_counts()
print(cell_counts)
print(adata.obs.columns)

sc.pl.violin(
         adata,
         ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
         jitter=0.4,
         multi_panel=True, save=f"_{base_name}_PreQC.png")

adata = adata[
      (adata.obs['n_genes_by_counts'] > 800) &
      (adata.obs['n_genes_by_counts'] < 8000) &
      (adata.obs['total_counts'] > 1200) &
      (adata.obs['total_counts'] < 30000) &
      (adata.obs['pct_counts_mt'] < 25), :
      ]


sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,save=f"_{base_name}_PostQC.png") 


sc.pl.violin(
    adata,
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='leiden',     # change if you use another cluster label
    jitter=0.4,
    rotation=45,
    multi_panel=True,save=f"_{base_name}_qc_violin.png") 


samples = adata.obs['sample'].unique()

# Loop through samples and plot individual UMAPs
import matplotlib.pyplot as plt

for sample in samples:
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color='sample',
        title=f"Sample: {sample}",
        size=20,
        show=False  # do not display
    )
    fig = plt.gcf()  # get the current figure
    fig.set_size_inches(10, 10)  # resize figure
    fig.savefig(f"figures/refined_{base_name}_{sample}.png", dpi=600, bbox_inches='tight')
    plt.close(fig)


adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")
