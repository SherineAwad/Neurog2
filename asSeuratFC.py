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
parser = argparse.ArgumentParser(description="Perform DGE and plot top genes heatmap with combined filtering")
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

sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.05,
    max_out_group_fraction=0.3,
    key='rank_genes_groups',
    key_added='filtered_rank_genes_groups'
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

    # For log2FC calculation, use RAW counts, not scaled data
    group_mask = adata.obs[groupby_col] == ct
    other_mask = ~group_mask

    print(f"Group cells: {group_mask.sum()}, Other cells: {other_mask.sum()}")

    # Reset index to ensure proper alignment
    df = df.reset_index(drop=True)

    # Pre-allocate arrays for log2FC results
    log2fcs = np.zeros(len(df))
    directions = [''] * len(df)

    valid_count = 0
    for i, row in df.iterrows():
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
            log2fc = np.log2((other_mean + pseudocount) / (group_mean + pseudocount))
            log2fcs[i] = log2fc
            directions[i] = 'up' if log2fc > 0 else 'down'

            valid_count += 1

        else:
            print(f"  WARNING: Gene '{gene}' not found in var_names!")
            log2fcs[i] = np.nan
            directions[i] = 'unknown'

    # Add log2FC and direction to dataframe
    df['log2fc'] = log2fcs
    df['direction'] = directions

    print(f"Calculated log2FC for {valid_count} genes")
    print(f"log2FC range: {df['log2fc'].min():.3f} to {df['log2fc'].max():.3f}")

    # Remove rows with NaN values in log2FC
    df = df.dropna(subset=['log2fc'])
    print(f"After removing NaN log2FC: {len(df)} genes")

    # COMBINED FILTERING APPROACH
    print("\nApplying combined filtering criteria...")
    
    # Option 1: Filter by both statistical significance AND effect size
    df_filtered_option1 = df[
        (df['p_val_adj'] < 0.05) &
        (df['wilcoxon_score'].abs() > 2.0) &
        (df['log2fc'].abs() > 1.0)  # |log2FC| > 2-fold
    ].copy()
    
    print(f"Option 1 (strict combined): {len(df_filtered_option1)} genes")
    
    # Option 2: Create a combined score
    significant_genes = df[df['p_val_adj'] < 0.05].copy()
    significant_genes['combined_score'] = (significant_genes['wilcoxon_score'].abs() * 0.5 + 
                                         significant_genes['log2fc'].abs() * 0.5)
    
    # Option 3: Select genes that meet either high effect OR high statistical criteria
    high_effect_genes = significant_genes[significant_genes['log2fc'].abs() > 1.0]
    high_wilcoxon_genes = significant_genes[significant_genes['wilcoxon_score'].abs() > 2.0]
    
    # Take union (genes that meet either criterion)
    df_filtered = pd.concat([high_effect_genes, high_wilcoxon_genes]).drop_duplicates()
    
    print(f"Option 3 (union): {len(df_filtered)} genes")
    print(f"  - High effect genes (|log2FC|>1): {len(high_effect_genes)}")
    print(f"  - High Wilcoxon genes (|score|>2): {len(high_wilcoxon_genes)}")
    print(f"  - Overlap between criteria: {len(high_effect_genes) + len(high_wilcoxon_genes) - len(df_filtered)}")

    if len(df_filtered) > 0:
        # Store filtered results
        all_dge_list.append(df_filtered)

        # Select top genes based on combined score (absolute log2FC + Wilcoxon)
        df_filtered_sorted = df_filtered.sort_values('combined_score', ascending=False)
        top_genes = df_filtered_sorted.head(10)['gene'].tolist()
        top_genes_dict[ct] = top_genes
        print(f"Top 10 genes for {ct} (by combined score): {top_genes}")

# Concatenate all DGE results
if all_dge_list:
    all_dge_df = pd.concat(all_dge_list, ignore_index=True)

    # Add method info
    all_dge_df['method'] = 'wilcoxon_combined'

    # Reorder columns
    column_order = ['gene', 'cell_type', 'wilcoxon_score', 'p_val_adj', 
                    'log2fc', 'direction', 'combined_score', 'method']
    all_dge_df = all_dge_df[column_order]

    # Save the results
    all_dge_csv = f"{base_name}_all_DGE_combined_filtering.csv"
    all_dge_df.to_csv(all_dge_csv, index=False)
    print(f"\n=== FINAL RESULTS ===")
    print(f"Saved DGE results to {all_dge_csv}")
    print(f"Total DGE results: {len(all_dge_df)}")

    # Summary statistics
    print(f"\nSummary statistics:")
    print(f"  Unique genes: {all_dge_df['gene'].nunique()}")
    print(f"  Up-regulated: {(all_dge_df['direction'] == 'up').sum()}")
    print(f"  Down-regulated: {(all_dge_df['direction'] == 'down').sum()}")
    
    # Check log2FC values
    print(f"\nlog2FC statistics:")
    print(f"  Min: {all_dge_df['log2fc'].min():.3f}")
    print(f"  Max: {all_dge_df['log2fc'].max():.3f}")
    print(f"  Mean: {all_dge_df['log2fc'].mean():.3f}")
    print(f"  |log2FC| > 1: {(all_dge_df['log2fc'].abs() > 1).sum()}")

    # Check Wilcoxon scores
    print(f"\nWilcoxon score statistics:")
    print(f"  Min: {all_dge_df['wilcoxon_score'].min():.3f}")
    print(f"  Max: {all_dge_df['wilcoxon_score'].max():.3f}")
    print(f"  Mean: {all_dge_df['wilcoxon_score'].mean():.3f}")
    print(f"  |score| > 2: {(all_dge_df['wilcoxon_score'].abs() > 2).sum()}")

else:
    print("No differentially expressed genes found with combined criteria")
    all_dge_df = pd.DataFrame()

# HEATMAP GENERATION
if not all_dge_df.empty and top_genes_dict:
    print(f"\n=== GENERATING HEATMAPS ===")

    # Define your cell type order
    celltype_order = ['MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones']

    # Get top genes from each cell type
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
            plt.title("Top 10 Marker Genes per Cell Type\n(Combined Wilcoxon + log2FC Filtering)", fontsize=14, pad=20)
            plt.ylabel("Genes", fontsize=12)
            plt.xlabel("Cell Types", fontsize=12)
            plt.tight_layout()
            heatmap_file = f"{base_name}_combined_filtering_heatmap.png"
            plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Combined filtering heatmap saved to {heatmap_file}")

            # Also save the heatmap data for reference
            log2fc_matrix.to_csv(f"{base_name}_heatmap_combinedfiltering.csv")
            print(f"Heatmap data saved to {base_name}_heatmap_data.csv")

print("Analysis complete!")
