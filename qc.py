import os
import numpy as np
import pandas as pd
import scanpy as sc
import argparse
import importlib_metadata
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# Fix importlib.metadata for older Python environments
sys.modules['importlib.metadata'] = importlib_metadata
# ARGUMENTS
parser = argparse.ArgumentParser(description="QC")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()
myObject = args.myObject
adata = sc.read(myObject)

# Count number of cells per sample
cell_counts = adata.obs['sample'].value_counts()

print(cell_counts)


