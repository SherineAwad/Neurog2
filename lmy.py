import scanpy as sc
import numpy as np
import loompy
import scipy.sparse
import pandas as pd
import argparse
import os
import sys
import importlib_metadata

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()

myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "loom_" + base_name + ".h5ad"

# Load AnnData
adata = sc.read(myObject)

# Use adata.raw if it exists and is populated, else fallback to adata.X
X = adata.raw.X if adata.raw is not None else adata.X
var_names = adata.raw.var_names if adata.raw is not None else adata.var_names

# Convert to dense if sparse
if scipy.sparse.issparse(X):
    X = X.toarray()

# Gene filter: keep genes expressed in >=3 cells and with total counts >= 10
gene_filter = (X > 0).sum(axis=0) >= 3
gene_filter &= X.sum(axis=0) >= 10

X_filt = X[:, gene_filter]
genes_filt = np.array(var_names)[gene_filter]

# Transpose for loom (genes x cells)
M = X_filt.T

# Cell metrics
nGene = (M > 0).sum(axis=0).astype(int)
nUMI = M.sum(axis=0).astype(int)

# Replace NaNs, infs (if any)
nGene = np.nan_to_num(nGene, nan=0, posinf=0, neginf=0)
nUMI = np.nan_to_num(nUMI, nan=0, posinf=0, neginf=0)

# Create loom
row_attrs = {"Gene": genes_filt}
col_attrs = {
    "CellID": np.array(adata.obs_names),
    "nGene": nGene,
    "nUMI": nUMI
}

loomfile = base_name + ".loom"
loompy.create(loomfile, M, row_attrs, col_attrs)
print("✅ Filtered loom created successfully.")

