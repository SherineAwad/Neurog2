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

cleaned_name = ''.join(parts)  # 'doubletsRemoved_threshold0.8_neurog2'

base_name = "refined_" + cleaned_name

newObject = base_name + ".h5ad"

adata = sc.read(myObject)

# base_name = 'refined_doubletsRemoved_threshold0.8_neurog2'
# newObject = 'refined_doubletsRemoved_threshold0.8_neurog2.h5ad'


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
         multi_panel=True, save=f"_{base_name}_PostQC.png") 



figurename1 = "figures/"+base_name+"_qc_violin.png"
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


samples = adata.obs['sample'].unique()

# Loop through samples and plot individual UMAPs
for sample in samples:
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color='sample',
        title=f"Sample: {sample}",
        size=20,
        save=f"_refined_{base_name}_{sample}.png",
        show=False
    )



adata.obs_names_make_unique()
adata.write(newObject, compression="gzip")
