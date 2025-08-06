import os
import pandas as pd
import subprocess
import argparse 

# Set working directory
wdir = "/nfs/turbo/umms-thahoang/sherine/neurog2/"
os.chdir(wdir)

# Create parser
parser = argparse.ArgumentParser(description="Run pySCENIC with custom input files and settings.")

# Add arguments
parser.add_argument("--tfs", required=True, help="Path to transcription factors list (e.g. allTFs_mm.txt)")
parser.add_argument("--loom", required=True, help="Path to input loom file for SCENIC (e.g. Rbpj_mCherry_filtered_scenic.loom)")
parser.add_argument("--output", required=True, help="Path to output CSV file (e.g. Rbpj_mCherry_adj.csv)")
parser.add_argument("--workers", type=int, default=4, help="Number of workers to use (default: 4)")

# Parse arguments
args = parser.parse_args()

# Assign to variables
f_tfs = args.tfs
f_loom_path_scenic = args.loom
output_csv = args.output
num_workers = args.workers

print("Done argparsing") 

# Prepare command
cmd = [
    "pyscenic",
    "grn",
    f_loom_path_scenic,
    f_tfs,
    "-o", output_csv,
    "--num_workers", str(num_workers)
]

# Run pyscenic grn step
print("Running pyscenic GRN inference...")
result = subprocess.run(cmd, capture_output=True, text=True)

# Print outputs for debugging
print("STDOUT:\n", result.stdout)
print("STDERR:\n", result.stderr)

# Check if command ran successfully
if result.returncode != 0:
    raise RuntimeError("pyscenic failed. Check stderr for details.")

# Ensure the output file was created and is not empty
if not os.path.isfile(output_csv) or os.stat(output_csv).st_size == 0:
    raise FileNotFoundError(f"Output file '{output_csv}' not created or is empty.")

# Load adjacencies
print(f"Reading output: {output_csv}")
adjacencies = pd.read_csv(output_csv, sep='\t')  # sep='\t' if it's TSV format
print(adjacencies.head())

