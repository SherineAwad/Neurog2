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
parser = argparse.ArgumentParser(description="Perform DGE and plot top genes heatmap")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
args = parser.parse_args()
myObject = args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

# READ DATA
adata = sc.read(myObject)
groupby_col = "celltype"

print(f"Data shape: {adata.shape}")
print(f"Has raw data: {adata.raw is not None}")

# Run Wilcoxon on SCALED data (for statistical test)
sc.tl.rank_genes_groups(
    adata,
    groupby=groupby_col,
    method='wilcoxon',
    use_raw=False,  # Use scaled data for the test
    pts=True,
)

# Get unique cell types
celltypes = adata.obs[groupby_col].cat.categories

# Initialize storage
all_dge_list = []
top_genes_dict = {}

for ct in celltypes:
    print(f"\n=== Processing cell type: {ct} ===")
    
    # Get results from Wilcoxon test (using scaled data)
    df = sc.get.rank_genes_groups_df(adata, group=ct, key='rank_genes_groups')
    print(f"Initial results for {ct}: {len(df)} genes")

    # Rename columns
    df = df.rename(columns={
        'names': 'gene',
        'scores': 'wilcoxon_score',
        'pvals_adj': 'p_val_adj'
    })

    # Add column for cell type
    df['cell_type'] = ct

    # FILTERING CRITERIA
    df_filtered = df[
        (df['p_val_adj'] < 0.05) &
        (df['wilcoxon_score'].abs() > 2.0)
    ].copy()

    print(f"After filtering: {len(df_filtered)} genes")

    if len(df_filtered) > 0:
        # For log2FC calculation, use RAW counts, not scaled data
        group_mask = adata.obs[groupby_col] == ct
        other_mask = ~group_mask
        
        print(f"Group cells: {group_mask.sum()}, Other cells: {other_mask.sum()}")
        
        # Reset index to ensure proper alignment
        df_filtered = df_filtered.reset_index(drop=True)
        
        # Pre-allocate arrays for results
        mean_diffs = np.zeros(len(df_filtered))
        log2fcs = np.zeros(len(df_filtered))
        directions = [''] * len(df_filtered)
        
        valid_count = 0
        for i, row in df_filtered.iterrows():
            gene = row['gene']
            
            if gene in adata.var_names:
                gene_idx = adata.var_names.get_loc(gene)
                
                # Get RAW expression values for log2FC calculation
                if adata.raw is not None:
                    # Use raw counts if available
                    if hasattr(adata.raw.X, 'toarray'):
                        group_expr = adata.raw.X[group_mask, gene_idx].toarray().flatten()
                        other_expr = adata.raw.X[other_mask, gene_idx].toarray().flatten()
                    else:
                        group_expr = adata.raw.X[group_mask, gene_idx].flatten()
                        other_expr = adata.raw.X[other_mask, gene_idx].flatten()
                else:
                    # If no raw data, use the processed data but handle negative values
                    if hasattr(adata.X, 'toarray'):
                        group_expr = adata.X[group_mask, gene_idx].toarray().flatten()
                        other_expr = adata.X[other_mask, gene_idx].toarray().flatten()
                    else:
                        group_expr = adata.X[group_mask, gene_idx].flatten()
                        other_expr = adata.X[other_mask, gene_idx].flatten()
                
                # Convert to non-negative values if needed
                if np.min(group_expr) < 0 or np.min(other_expr) < 0:
                    # Shift all values to be positive by adding the minimum value
                    min_val = min(np.min(group_expr), np.min(other_expr))
                    group_expr = group_expr - min_val + 1e-5
                    other_expr = other_expr - min_val + 1e-5
                
                group_mean = np.mean(group_expr)
                other_mean = np.mean(other_expr)
                
                # Calculate log2FC
                pseudocount = 1e-5
                log2fc = np.log2((group_mean + pseudocount) / (other_mean + pseudocount))
                
                mean_diffs[i] = group_mean - other_mean
                log2fcs[i] = log2fc
                directions[i] = 'up' if log2fc > 0 else 'down'
                
                valid_count += 1
                
                if i < 5:  # Print first 5 for debugging
                    print(f"  Gene {i}: {gene}")
                    print(f"    Group mean: {group_mean:.4f}, Other mean: {other_mean:.4f}")
                    print(f"    log2FC: {log2fc:.4f}")
            else:
                print(f"  WARNING: Gene '{gene}' not found in var_names!")
                mean_diffs[i] = np.nan
                log2fcs[i] = np.nan
                directions[i] = 'unknown'
        
        # Add to dataframe
        df_filtered['mean_diff'] = mean_diffs
        df_filtered['log2fc'] = log2fcs
        df_filtered['direction'] = directions
        
        print(f"Successfully calculated log2FC for {valid_count} out of {len(df_filtered)} genes")
        print(f"log2FC range: {df_filtered['log2fc'].min():.3f} to {df_filtered['log2fc'].max():.3f}")
        
        # Remove rows with NaN values
        df_filtered = df_filtered.dropna(subset=['log2fc'])
        print(f"After removing NaN: {len(df_filtered)} genes")
        
        if len(df_filtered) > 0:
            all_dge_list.append(df_filtered)
            # Store top 10 upregulated genes for heatmap
            top_upregulated = df_filtered[df_filtered['log2fc'] > 0].nlargest(10, 'log2fc')['gene'].tolist()
            top_genes_dict[ct] = top_upregulated

# Concatenate all DGE results
if all_dge_list:
    all_dge_df = pd.concat(all_dge_list, ignore_index=True)
    
    # Add method info
    all_dge_df['method'] = 'wilcoxon'
    
    # Reorder columns
    column_order = ['gene', 'cell_type', 'wilcoxon_score', 'p_val_adj', 
                    'mean_diff', 'log2fc', 'direction', 'method']
    all_dge_df = all_dge_df[column_order]
    
    # Save the results
    all_dge_csv = f"{base_name}_all_DGE_wilcoxon.csv"
    all_dge_df.to_csv(all_dge_csv, index=False)
    print(f"\n=== FINAL RESULTS ===")
    print(f"Saved DGE results to {all_dge_csv}")
    print(f"Total DGE results: {len(all_dge_df)}")
    
    # Check log2FC values
    print(f"log2FC column stats:")
    print(f"  Non-null values: {all_dge_df['log2fc'].notna().sum()}")
    print(f"  Zero values: {(all_dge_df['log2fc'] == 0).sum()}")
    print(f"  Min: {all_dge_df['log2fc'].min():.6f}")
    print(f"  Max: {all_dge_df['log2fc'].max():.6f}")
    print(f"  Mean: {all_dge_df['log2fc'].mean():.6f}")
    
    # Show first few rows with log2FC
    print("\nFirst 10 rows:")
    print(all_dge_df[['gene', 'cell_type', 'log2fc']].head(10))
    
else:
    print("No differentially expressed genes found")
    all_dge_df = pd.DataFrame()

# HEATMAP GENERATION
if not all_dge_df.empty and top_genes_dict:
    print(f"\n=== GENERATING HEATMAPS ===")
    
    # Define your cell type order (update with your actual cell types)
    celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']
    
    # Get top upregulated genes from each cell type
    top_genes_all = []
    for ct in celltype_order:
        if ct in top_genes_dict:
            top_genes_all.extend(top_genes_dict[ct])
    top_genes_all = list(dict.fromkeys(top_genes_all))  # remove duplicates
    
    print(f"Number of top genes for heatmap: {len(top_genes_all)}")
    
    if top_genes_all:
        # Create log2FC matrix
        log2fc_matrix = pd.DataFrame(0.0, index=top_genes_all, columns=celltype_order)
        
        # Fill with actual log2FC values
        for ct in celltype_order:
            if ct in all_dge_df['cell_type'].unique():
                ct_data = all_dge_df[all_dge_df['cell_type'] == ct]
                for gene in top_genes_all:
                    if gene in ct_data['gene'].values:
                        log2fc_value = ct_data[ct_data['gene'] == gene]['log2fc'].values[0]
                        log2fc_matrix.loc[gene, ct] = log2fc_value
        
        print(f"log2FC matrix shape: {log2fc_matrix.shape}")
        
        # Determine appropriate color scale limits
        non_zero_values = log2fc_matrix.values[log2fc_matrix.values != 0]
        if len(non_zero_values) > 0:
            max_abs_log2fc = np.max(np.abs(non_zero_values))
            vmax = max(2, max_abs_log2fc)  # At least ±2 for good visualization
            
            # Plot log2FC heatmap
            plt.figure(figsize=(max(8, len(celltype_order)*1.5), max(6, len(top_genes_all)*0.4)))
            sns.heatmap(
                log2fc_matrix,
                cmap='RdBu_r',
                yticklabels=True,
                xticklabels=True,
                cbar_kws={'label': 'log2 Fold Change'},
                center=0,
                vmin=-vmax,
                vmax=vmax,
                square=False
            )
            plt.title("Top 10 Upregulated Marker Genes per Cell Type\n(log2 Fold Change)", fontsize=16, pad=20)
            plt.ylabel("Genes", fontsize=12)
            plt.xlabel("Cell Types", fontsize=12)
            plt.tight_layout()
            heatmap_file = f"{base_name}_top10_log2fc_heatmap.png"
            plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"log2FC heatmap saved to {heatmap_file}")
            
            # Also create a version with gene labels if there are too many genes
            if len(top_genes_all) > 30:
                plt.figure(figsize=(max(8, len(celltype_order)*1.5), 12))
                sns.heatmap(
                    log2fc_matrix,
                    cmap='RdBu_r',
                    yticklabels=False,
                    xticklabels=True,
                    cbar_kws={'label': 'log2 Fold Change'},
                    center=0,
                    vmin=-vmax,
                    vmax=vmax
                )
                plt.title("Top Upregulated Marker Genes per Cell Type\n(log2 Fold Change)", fontsize=16, pad=20)
                plt.ylabel("Genes", fontsize=12)
                plt.xlabel("Cell Types", fontsize=12)
                plt.tight_layout()
                heatmap_file_no_labels = f"{base_name}_top10_log2fc_heatmap_no_labels.png"
                plt.savefig(heatmap_file_no_labels, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"log2FC heatmap (no gene labels) saved to {heatmap_file_no_labels}")
        else:
            print("No non-zero log2FC values for heatmap")
    else:
        print("No top genes for heatmap")
else:
    print("No data for heatmap generation")

print("Analysis complete!")
