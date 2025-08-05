import loompy
import argparse

# Create parser
parser = argparse.ArgumentParser(description="Run pySCENIC with custom input files and settings.")

# Add arguments
parser.add_argument("--loom", required=True, help="Path to input loom file for SCENIC (e.g. Rbpj_mCherry_filtered_scenic.loom)")

# Parse arguments
args = parser.parse_args()

# Assign to variables
loom_path = args.loom



with loompy.connect(loom_path, mode='r') as ds:
    print("✅ Basic Info:")
    print(f"  Shape (genes x cells): {ds.shape}")
    
    print("\n🧬 Row Attributes (Gene metadata):")
    print(f"  {list(ds.ra.keys())}")
    print(f"  Example gene names: {ds.ra['Gene'][:5]}")

    print("\n🧫 Column Attributes (Cell metadata):")
    print(f"  {list(ds.ca.keys())}")
    print(f"  Example CellIDs: {ds.ca['CellID'][:5]}")

    # Optionally look at data:
    print("\n🔢 Expression Matrix Summary:")
    expr = ds[:, :].T  # transpose to cells x genes
    print(f"  Expression matrix shape: {expr.shape}")
    print(f"  Matrix dtype: {expr.dtype}")
    print(f"  Non-zero entries: {expr.nonzero()[0].size}")

