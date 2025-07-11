import scrublet as scr
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()
myObject = args.myObject

# Output file
newObject = "doublesRemoved_" + myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

# Read the full AnnData object
combined_adata = sc.read(myObject)

# Run Scrublet
counts_matrix = combined_adata.X[:]
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# To Increase expected doublet rate to 15% and simulate more doublets (3x)
## scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15, sim_doublet_ratio=3.0)


# Add prediction to .obs
combined_adata.obs['predicted_doublets'] = predicted_doublets
combined_adata.obs['doublet_status'] = ["Doublet" if x else "Singlet" for x in predicted_doublets]


import pandas as pd  # ensure this is imported if not already

combined_adata.obs['doublet_scores'] = doublet_scores
combined_adata.obs['doublet_status'] = pd.Categorical(
    combined_adata.obs['doublet_status'], categories=["Singlet", "Doublet"]
)


#  PLOT UMAP WITH DOUBLET = RED, SINGLET = LIGHT GREY
palette = {'Singlet': 'lightgrey', 'Doublet': 'red'}

sc.pl.umap(
    combined_adata,
    color='doublet_status',
    palette=palette,
    size=50,
    title='Doublets (red) vs Singlets (grey)',
    save="_doubletStatus.png"
)


scrub.plot_histogram()
plt.savefig("figures/real_simulated_histogram.png", dpi=300)
plt.close()


plt.hist(doublet_scores, bins=50)
plt.xlabel('Doublet Score')
plt.ylabel('Number of Cells')
plt.title('Scrublet Doublet Score Distribution')
plt.axvline(x=scrub.threshold_, color='r', linestyle='--', label='Default threshold')
plt.legend()
plt.tight_layout()
plt.savefig("figures/doublet_score_histogram.png", dpi=300)
plt.close()


#  NOW filter out doublets and save singlets only
combined_adata = combined_adata[~combined_adata.obs['predicted_doublets'], :]
combined_adata.write(newObject, compression="gzip")

