import glob
import subprocess
import argparse
import os


# Find all feather ranking db files in current directory (or specify full path)
f_db_list = glob.glob("*feather")
if not f_db_list:
    raise FileNotFoundError("No feather files found in current directory.")

# Motif annotation file path
f_motif_path = "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# Create parser
parser = argparse.ArgumentParser(description="Run pySCENIC motif enrichment step with input files.")

# Add arguments
parser.add_argument("--loom", required=True, help="Path to input loom file for SCENIC")
parser.add_argument("--adjacency", required=True, help="Path to input adjacency CSV file")
parser.add_argument("--output", required=True, help="Path to output regulon CSV file")

# Parse arguments
args = parser.parse_args()

# Assign variables
f_loom_path_scenic = args.loom
adj_file = args.adjacency
output_file = args.output

print("Loom file:", f_loom_path_scenic)
print("Feather DB files found:", f_db_list)
print("Motif annotation file:", f_motif_path)
print("Adjacency input file:", adj_file)
print("Output regulons file:", output_file)



# Construct command as list of strings
cmd = [
    "pyscenic", "ctx", adj_file,
    *f_db_list,  # expand feather files here
    "--annotations_fname", f_motif_path,
    "--expression_mtx_fname", f_loom_path_scenic,
    "-o", output_file,
    "--mask_dropouts",
    "--num_workers", "4"
]

# Run command and capture output
result = subprocess.run(cmd, capture_output=True, text=True)

print("STDOUT:")
print(result.stdout)
print("STDERR:")
print(result.stderr)

if result.returncode != 0:
    print(f"Error: pyscenic ctx failed with code {result.returncode}")

