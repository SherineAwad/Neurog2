import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np
import os

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata

# Argument parsing
parser = argparse.ArgumentParser(description="Remove specific clusters and re-cluster")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('annotations', help='Path to annotations file')
args = parser.parse_args()
  
myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = "annotated_" + base_name + ".h5ad"
sys.modules['importlib.metadata'] = importlib_metadata

annot_file = args.annotations 

adata = sc.read_h5ad(myObject, backed="r")


cluster_to_celltype_dict = {}
with open(annot_file, "r") as f:
    for line in f:
        cluster, celltype = line.strip().split('\t')
        cluster_to_celltype_dict[cluster] = celltype

cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}

adata.obs["celltype"] = adata.obs["leiden"].map(cluster_to_celltype_dict)


figure_name = base_name+"_annotationsON.png"
adata.obs_names_make_unique()
sc.pl.umap(adata, color='celltype',legend_loc="on data", save=figure_name)

figure_name = base_name +"_annotations.png"
sc.pl.umap(adata, color='celltype',save=figure_name)

adata.write(newObject, compression="gzip")
