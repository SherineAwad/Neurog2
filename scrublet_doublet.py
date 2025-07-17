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

import scrublet as scr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

counts_matrix = combined_adata.X[:]

# Initialize Scrublet (without threshold)
scrub = scr.Scrublet(counts_matrix)

# Run doublet detection
doublet_scores, _ = scrub.scrub_doublets()

# Apply manual threshold of 0.1
threshold = 0.15
predicted_doublets = doublet_scores > threshold

# Add results to .obs
combined_adata.obs['doublet_scores'] = doublet_scores
combined_adata.obs['predicted_doublets'] = predicted_doublets
combined_adata.obs['doublet_status'] = pd.Categorical(
    ["Doublet" if x else "Singlet" for x in predicted_doublets],
    categories=["Singlet", "Doublet"]
)

# UMAP plot: Doublets in red, Singlets in light grey
palette = {'Singlet': 'lightgrey', 'Doublet': 'red'}
sc.pl.umap(
    combined_adata,
    color='doublet_status',
    palette=palette,
    size=50,
    title='Doublets (red) vs Singlets (grey)',
    save="_doubletStatus0.15.png"
)

# Plot scrublet histogram (real vs simulated)
scrub.plot_histogram()
plt.savefig("figures/real_simulated_histogram.png", dpi=300)
plt.close()

# Plot histogram of doublet scores with custom threshold
plt.hist(doublet_scores, bins=50)
plt.xlabel('Doublet Score')
plt.ylabel('Number of Cells')
plt.title('Scrublet Doublet Score Distribution')
plt.axvline(x=threshold, color='r', linestyle='--', label='Manual threshold (0.1)')
plt.legend()
plt.tight_layout()
plt.savefig("figures/doublet_score_histogram.png", dpi=300)
plt.close()


#  NOW filter out doublets and save singlets only
combined_adata = combined_adata[~combined_adata.obs['predicted_doublets'], :]
combined_adata.write(newObject, compression="gzip")

