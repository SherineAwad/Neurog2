# Neurog2 Project - scRNA-seq Analysis

This project focuses on the single-cell RNA sequencing (scRNA-seq) analysis of samples related to **Neurog2** expression at different stages and control. The goal is to study how MG cells develope to other cells. 

## Samples

| Sample Name           | Description          |
|-----------------------|----------------------|
| 5 weeks Neurog2_9SA   | TH1_GFP_mScarlet3    |
| 2 months control      | TH2_GFP_mScarlet3    |
| 2 months Neurog2_9SA  | TH3_GFP_mScarlet3    |

## Analysis Workflow

The analysis was performed using [Scanpy](https://scanpy.readthedocs.io/en/stable/), a scalable toolkit for analyzing single-cell gene expression data. The workflow included:


### Preprocessing 

1. **Merge Multiple Samples**
   Multiple `AnnData` objects are combined into one using their sample names as labels. This enables joint analysis while preserving sample identity.

2. **Identify Mitochondrial Genes**
   Genes that start with `"mt-"` are flagged as mitochondrial genes, which are important indicators of cell stress or damage.

3. **Calculate Quality Control (QC) Metrics**
   Standard QC metrics are computed for each cell:

   * `n_genes_by_counts`: Number of genes detected
   * `total_counts`: The total number of UMIs observed per cell
   * `pct_counts_mt`: Percent of transcripts from mitochondrial genes

4. **Visualize QC Metrics (Before Filtering)**
   Violin plots are used to visualize the distribution of these metrics to help identify low-quality cells.

5. **Filter Out Low-Quality Cells**
   Cells are removed if they have:

   * Too few or too many detected genes (e.g. <800 or >8000)
   * Extremely low or high total transcript counts
   * High mitochondrial content (e.g. >25%), indicating cell stress

6. **Further Filtering**

   * Cells with fewer than 100 genes are removed
   * Genes found in fewer than 3 cells are excluded

7. **Visualize QC Metrics (After Filtering)**
   Another set of violin plots is generated to assess the impact of filtering on the dataset.

8. **Save the Processed Data**
   The cleaned and filtered data is saved as an `.h5ad` file for downstream analysis.

---

## Figures

### 1. UMAP Visualization
![UMAP of Neurog2 samples](figures/umap_neurog2.png)  
UMAP plot colored by sample, showing clustering and distribution of single cells from different conditions.


## Per sample UMAP 

![Control 2mo](figures/umap_sample_control_2mo.png)
![Neurog2_9SA_5weeks](figures/umap_sample_Neurog2_9SA_5weeks.png)
![Neurog2_9SA_2mo](figures/umap_sample_Neurog2_9SA_2mo.png)


###  Scanpy QC Metrics — Quick Overview

#### 🔹 `n_genes_by_counts`

* **Definition**: Number of genes with **non-zero counts** in each cell.
* **Use**: Helps filter out cells with too few expressed genes (often poor quality or empty droplets).

#### 🔹 `total_counts`

* **Definition**: Total **number of counts (UMIs or reads)** in a cell.
* **Use**: Indicates cell complexity or sequencing depth. Very low values may indicate damaged cells or low capture.

#### 🔹 `pct_counts_mt`

* **Definition**: Percentage of counts from **mitochondrial genes** (e.g., genes starting with `mt-` in mouse or `MT-` in human).
* **Use**: High percentages may indicate **cell stress or apoptosis**; often used to filter out low-quality cells.


![Before Filtering QC metrics](figures/violin_QC.png)  
Violin plots displaying quality control metrics such as number of genes detected per cell, total counts, and percentage of mitochondrial gene expression.

## Filtering Criteria

Quality filtering was applied to remove low-quality cells and potential doublets. Cells were retained only if they met all the following criteria:

- Number of genes detected per cell between **800 and 8000**
- Total counts per cell between **1200 and 30000**
- Percentage of mitochondrial gene counts less than **25%**

This filtering step ensures removal of dead or dying cells and technical artifacts to improve downstream analysis quality.


### Additional Analysis Figure
![After Filtering QC metrics](figures/violin_AfterQC.png)  


## Number of cells per sample 

| Sample              | Cell Count |
|---------------------|------------|
| Neurog2_9SA_5weeks  | 27,732     |
| Neurog2_9SA_2mo     | 11,486     |
| control_2mo         | 9,701      |



## Clustering 

1. **Load the Data**
   A preprocessed `AnnData` object is loaded from disk.

2. **Normalize and Transform**

   * Normalize gene expression values 
   * Apply a logarithmic transformation to stabilize variance across genes.

3. **Feature Selection**

   * Identify the top 2,000 highly variable genes using the Seurat method. These are the most informative genes for downstream analysis.

4. **Scale the Data**

   * Standardize the expression values (mean = 0, variance = 1).
   * Clip extreme values to a maximum of 10 to reduce the impact of outliers.

5. **Dimensionality Reduction (PCA)**

   * Perform Principal Component Analysis to reduce data dimensionality and denoise the dataset.

6. **Construct the Neighborhood Graph**

   * Build a k-nearest neighbors graph based on PCA to capture the local structure of the data.

7. **UMAP Embedding**

   * Compute a 2D UMAP embedding for visualization of the dataset’s structure.

8. **Visualize UMAP by Sample**

   * Generate a UMAP plot where cells are colored by their sample origin.
   * Count how many cells belong to each sample.

9. **Per-Sample UMAP Plots**

   * Loop through each sample and generate a separate UMAP plot showing only the cells from that sample.

10. **Visualize Predicted Doublets**

* Plot a UMAP colored by predicted doublet labels and doublet scores to inspect doublet detection results.


### Marker Gene UMAP Plots
Below are the UMAP visualizations of marker gene expression across clusters. These are auto-generated from your data and saved in the figures/ directory.


### Initial Clustering 

![UMAP CLUSTERS](figures/umap_clusters.png)


### Marker Gene UMAP Plots
Below are the UMAP visualizations of marker gene expression across clusters. These are auto-generated from your data and saved in the figures/ directory.


<img src="figures/umapclustered_analysed_neurog2_Malat1.png" alt="Malat1" width="33%"><img src="figures/umapclustered_analysed_neurog2_mt-Atp6.png" alt="mt-Atp6" width="33%"><img src="figures/umapclustered_analysed_neurog2_Sox9.png" alt="Sox9" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Glul.png" alt="Glul" width="33%"><img src="figures/umapclustered_analysed_neurog2_Lhx2.png" alt="Lhx2" width="33%"><img src="figures/umapclustered_analysed_neurog2_Rlbp1.png" alt="Rlbp1" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umapclustered_analysed_neurog2_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umapclustered_analysed_neurog2_Calb2.png" alt="Calb2" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umapclustered_analysed_neurog2_Calb1.png" alt="Calb1" width="33%"><img src="figures/umapclustered_analysed_neurog2_Sebox.png" alt="Sebox" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Gad1.png" alt="Gad1" width="33%"><img src="figures/umapclustered_analysed_neurog2_Elavl3.png" alt="Elavl3" width="33%"><img src="figures/umapclustered_analysed_neurog2_Cabp5.png" alt="Cabp5" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Isl1.png" alt="Isl1" width="33%"><img src="figures/umapclustered_analysed_neurog2_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umapclustered_analysed_neurog2_Ascl1.png" alt="Ascl1" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Olig2.png" alt="Olig2" width="33%"><img src="figures/umapclustered_analysed_neurog2_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umapclustered_analysed_neurog2_Chat.png" alt="Chat" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Prdm1.png" alt="Prdm1" width="33%"><img src="figures/umapclustered_analysed_neurog2_Otx2.png" alt="Otx2" width="33%"><img src="figures/umapclustered_analysed_neurog2_Insm1.png" alt="Insm1" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Sox11.png" alt="Sox11" width="33%"><img src="figures/umapclustered_analysed_neurog2_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umapclustered_analysed_neurog2_Hes5.png" alt="Hes5" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Emx1.png" alt="Emx1" width="33%"><img src="figures/umapclustered_analysed_neurog2_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umapclustered_analysed_neurog2_GFP.png" alt="GFP" width="33%">

<img src="figures/umapclustered_analysed_neurog2_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umapclustered_analysed_neurog2_Tfap2a.png" alt="Tfap2a" width="33%"><img src="figures/umapclustered_analysed_neurog2_Bsn.png" alt="Bsn" width="33%">



### QC per Clsuter 

<img src="figures/qc_violin_by_cluster.png" width="550"/>

### Remove clusters 

| ID | Cell Type           |
|-----|--------------------|
| 7   | Bad Cells          |
| 8   | Microglia          |
| 11  | Bad Cells          |
| 20  | Microglia          |
| 28  | Monocyte           |
| 33  | RPE/Pax2           |
| 34  | SMC                |

then we reclustered and replot the marker genes as below: 


### UMAP

![UMAP RE CLUSTERS](figures/umap_reClusters.png) 


### Per sample UMAP 

![Control 2mo](figures/umap_reclustered_control_2mo.png)
![Neurog2_9SA_5weeks](figures/umap_reclustered_Neurog2_9SA_5weeks.png)
![Neurog2_9SA_2mo](figures/umap_reclustered_Neurog2_9SA_2mo.png)

<img src="figures/umap_reClustered_clustered_analysed_neurog2_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Pax2.png" alt="Pax2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Rpe65.png" alt="Rpe65" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Tie1.png" alt="Tie1" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Rho.png" alt="Rho" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Nrl.png" alt="Nrl" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Arr3.png" alt="Arr3" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Malat1.png" alt="Malat1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_mt-Atp6.png" alt="mt-Atp6" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Glul.png" alt="Glul" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Lhx2.png" alt="Lhx2" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Csf1r.png" alt="Csf1r" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Calb2.png" alt="Calb2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Calb1.png" alt="Calb1" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Gad1.png" alt="Gad1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Elavl3.png" alt="Elavl3" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Cabp5.png" alt="Cabp5" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Isl1.png" alt="Isl1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Elavl4.png" alt="Elavl4" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Olig2.png" alt="Olig2" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Chat.png" alt="Chat" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Prdm1.png" alt="Prdm1" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Otx2.png" alt="Otx2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Insm1.png" alt="Insm1" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Sox11.png" alt="Sox11" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Hes5.png" alt="Hes5" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Emx1.png" alt="Emx1" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_GFP.png" alt="GFP" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Tfap2a.png" alt="Tfap2a" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Bsn.png" alt="Bsn" width="33%">
<img src="figures/umap_reClustered_clustered_analysed_neurog2_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_reClustered_clustered_analysed_neurog2_Lhx4.png" alt="Lhx4" width="33%">


## Number of cells per sample 

| Sample              | Cell Count |
|---------------------|------------|
| Neurog2_9SA_5weeks  | 23,370     |
| Neurog2_9SA_2mo     | 10,115     |
| control_2mo         | 8,674      |

---

## Doublet Detection using `DoubletDetection`

A doublet is an artifact where two cells are captured and sequenced together, but incorrectly treated as one. Unlike `Scrublet`, which can operate effectively on clustered or preprocessed `AnnData` objects, the `DoubletDetection` tool is more sensitive to data structure and expects the **original, unclustered** `AnnData` object. Running it on a processed or subsetted object may yield suboptimal or misleading results.

In the workflow, we applied `DoubletDetection` to the original data (`adata`) to ensure it captures the full transcriptomic diversity and avoids artifacts introduced during clustering.

After running `DoubletDetection`, predicted doublets and doublet scores were stored in `adata.obs` under the keys:
- `predicted_doublet`: Boolean flag indicating whether each cell is a predicted doublet.
- `doublet_score`: Confidence score associated with doublet prediction.

The results were visualized using UMAP, colored by both prediction and score:

###  DoubletDetection Workflow

####  Data Preprocessing
- Input is a raw (or filtered) gene expression matrix.
- May optionally normalize, filter, and log-transform the data.

####  Synthetic Doublet Generation
- Creates artificial doublets by randomly pairing real cells.
- Averages their gene expression profiles to simulate doublets.

####  Clustering with Real + Synthetic Data
- Combines real and synthetic cells.
- Performs dimensionality reduction (typically **PCA**).
- Applies unsupervised clustering (usually **Phenograph**, a graph-based algorithm).

####  Voting Mechanism via Multiple Runs (Ensemble)
- Repeats the clustering multiple times (default: **50 runs**).
- Tracks how often each real cell clusters with synthetic doublets.
- Cells that frequently cluster with synthetic doublets are flagged as potential doublets.

####  Thresholding & Output
- Assigns a **doublet probability score** to each cell.
- Applies a threshold (user-defined or default) to classify each cell as a **doublet** or **singlet**.


### Doublet Scores and Conversion

   ![UMAP](figures/umap_doubletScores_neurog2_doublets.png)
   ![Doublets Thresholds](figures/threshold_test.png)
   ![Doublets Conversion](figures/convergence_test.png)


###  Understanding `doublet_score` Thresholds

The `doublet_score` typically ranges from **0 to 1**

Your filter in the code:

```python
combined_adata = combined_adata[combined_adata.obs['doublet_score'] <= threshold]
```

This means you're **keeping** cells with `doublet_score <= threshold`.


###  Interpretation of Threshold:

* **Higher threshold** (e.g., `0.9`)
  🔹 You keep **more** cells
  🔹 Less doublets are removed

* **Lower threshold** (e.g., `0.4`) 
  🔹 You keep **fewer** cells
  🔹 More potential doublets are removed 


### Doublet Removal at Threshold 0.9 


### UMAP after  clustering 

![Doublet Detection](figures/umap_clustered_ddanalysed_doubletScores_0.9_neurog2_ddClusters.png) 

### UMAP after doublet removal at threshold 0.9 
 
![Doublet Removal](figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_doubletsRemoved_threshold0.9_clusters.png)

### Marker Genes UMAP after doublet removal at threshold 0.9 

<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Cabp5.png" alt="Cabp5" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Prdx6.png" alt="Prdx6" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Vim.png" alt="Vim" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Aqp4.png" alt="Aqp4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Slc6a9.png" alt="Slc6a9" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Abca8a.png" alt="Abca8a" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Slc1a3.png" alt="Slc1a3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_GFP.png" alt="GFP" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Gfap.png" alt="Gfap" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Rbfox3.png" alt="Rbfox3" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Elavl3.png" alt="Elavl3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Otx2.png" alt="Otx2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Glul.png" alt="Glul" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Sox11.png" alt="Sox11" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Apoe.png" alt="Apoe" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Atoh7.png" alt="Atoh7" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Lhx2.png" alt="Lhx2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Prdm1.png" alt="Prdm1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Bsn.png" alt="Bsn" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Rpe65.png" alt="Rpe65" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Tfap2a.png" alt="Tfap2a" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Gad1.png" alt="Gad1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Calb2.png" alt="Calb2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Chat.png" alt="Chat" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Pax6.png" alt="Pax6" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Insm1.png" alt="Insm1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Calb1.png" alt="Calb1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Lhx4.png" alt="Lhx4" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Notch1.png" alt="Notch1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Rho.png" alt="Rho" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Ascl1.png" alt="Ascl1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Emx1.png" alt="Emx1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Arr3.png" alt="Arr3" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Pax2.png" alt="Pax2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Isl1.png" alt="Isl1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Hes1.png" alt="Hes1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Nrl.png" alt="Nrl" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Tie1.png" alt="Tie1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_mt-Atp6.png" alt="mt-Atp6" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.9_clustered_ddanalysed_doubletScores_0.9_neurog2_threshold0.9_Malat1.png" alt="Malat1" width="33%">


### Doublet removal using 0.8 threshold 

### UMAP after  clustering 

![Doublet Detection](figures/umap_clustered_ddanalysed_doubletScores_0.8_neurog2_ddClusters.png)

### UMAP after doublet removal at threshold 0.8 

![Doublet Removal](figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_doubletsRemoved_threshold0.8_clusters.png)

### Marker Genes UMAP after doublet removal at threshold 0.8  

<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Cabp5.png" alt="Cabp5" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Vim.png" alt="Vim" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Prdx6.png" alt="Prdx6" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Aqp4.png" alt="Aqp4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Abca8a.png" alt="Abca8a" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Slc1a3.png" alt="Slc1a3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Gfap.png" alt="Gfap" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_GFP.png" alt="GFP" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Neurog2.png" alt="Neurog2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Elavl3.png" alt="Elavl3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Sox11.png" alt="Sox11" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Glul.png" alt="Glul" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Lhx2.png" alt="Lhx2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Otx2.png" alt="Otx2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Apoe.png" alt="Apoe" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Tfap2a.png" alt="Tfap2a" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Prdm1.png" alt="Prdm1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Pax6.png" alt="Pax6" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Gad1.png" alt="Gad1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Rpe65.png" alt="Rpe65" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Bsn.png" alt="Bsn" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Calb2.png" alt="Calb2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Chat.png" alt="Chat" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Notch1.png" alt="Notch1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Calb1.png" alt="Calb1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Acta2.png" alt="Acta2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Insm1.png" alt="Insm1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Lhx4.png" alt="Lhx4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Emx1.png" alt="Emx1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Rho.png" alt="Rho" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Pax2.png" alt="Pax2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Arr3.png" alt="Arr3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Ccr2.png" alt="Ccr2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Isl1.png" alt="Isl1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Hes1.png" alt="Hes1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Lhx1.png" alt="Lhx1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Nrl.png" alt="Nrl" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Tie1.png" alt="Tie1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_mt-Atp6.png" alt="mt-Atp6" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2_threshold0.8_Malat1.png" alt="Malat1" width="33%">


### Pre filtering QC
![Pre QC](figures/violin_refined_doubletsRemoved_threshold0.8_neurog2_PreQC.png)


#### How we filtered 
```python
adata = adata[
    (adata.obs['n_genes_by_counts'] > 800) &
    (adata.obs['n_genes_by_counts'] < 8000) &
    (adata.obs['total_counts'] > 1200) &
    (adata.obs['total_counts'] < 30000) &
    (adata.obs['pct_counts_mt'] < 25),
    :
]

```

### Post filtering QC

![Post QC](figures/violin_refined_doubletsRemoved_threshold0.8_neurog2_PostQC.png) 


### QC per Cluster 

![Doublet QC](figures/refined_doubletsRemoved_threshold0.8_neurog2_qc_violin.png)


### UMAP per Sample 

![Control 2mo](figures/umap_refined_refined_doubletsRemoved_threshold0.8_neurog2_control_2mo.png)
![Neurog2_9SA_5weeks](figures/umap_refined_refined_doubletsRemoved_threshold0.8_neurog2_Neurog2_9SA_5weeks.png)
![Neurog2_9SA_2mo](figures/umap_refined_refined_doubletsRemoved_threshold0.8_neurog2_Neurog2_9SA_2mo.png)



### Remove clusters 

| ID | Cell Type           |
|-----|--------------------|
| 7   | Bad rod            |
| 8   | Microglia          |
| 18  | Microglia          |
| 27  | Monocyte           |
| 32  | Astrocyte          |
| 33  | Smooth muscle cells|
| 15  | Likely doublets?   |
| 24  | Likely doublets?   |

### Then reCluster 

### Different resolutions 

#### Resolution = 1.4
![UMAP resolution 1.4](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_1.4res.png)

#### Resolution = 2.0
![UMAP resoultion 2.0](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_2.0res.png)

#### Resolution = 2.5
![UMAP resoultion 2.5](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_2.5res.png)


#### Resolution = 3.0
![UMAP resoultion 3.0](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_3.0res.png)


#### Resolution = 3.5
![UMAP resoultion 3.5](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_3.5res.png)


#### Resolution = 4.0
![UMAP resoultion 4.0](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_4.0res.png)

### We will stick to resolution 2.0 

#### Resolution = 2.0
![UMAP resolution 2.0](figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_2.0res.png)

#### UMAP of marker genes  

<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Abca8a.png" alt="Abca8a" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Guca1b.png" alt="Guca1b" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pou4f1.png" alt="Pou4f1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Hes1.png" alt="Hes1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pou4f3.png" alt="Pou4f3" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Apoe.png" alt="Apoe" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Prdm1.png" alt="Prdm1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Aqp4.png" alt="Aqp4" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Igf2.png" alt="Igf2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Prdx6.png" alt="Prdx6" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Arr3.png" alt="Arr3" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Insm1.png" alt="Insm1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Prkca.png" alt="Prkca" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Isl1.png" alt="Isl1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Prox1.png" alt="Prox1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Isl2.png" alt="Isl2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Ptprc.png" alt="Ptprc" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Bhlhe23.png" alt="Bhlhe23" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Rbfox3.png" alt="Rbfox3" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Bsn.png" alt="Bsn" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Rbpms.png" alt="Rbpms" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Cabp5.png" alt="Cabp5" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Lhx2.png" alt="Lhx2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Rho.png" alt="Rho" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Calb1.png" alt="Calb1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Lhx4.png" alt="Lhx4" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Rlbp1.png" alt="Rlbp1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Calb2.png" alt="Calb2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Malat1.png" alt="Malat1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Rom1.png" alt="Rom1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Cbln4.png" alt="Cbln4" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Rpe65.png" alt="Rpe65" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_mt-Atp6.png" alt="mt-Atp6" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sall1.png" alt="Sall1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Chat.png" alt="Chat" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Nefl.png" alt="Nefl" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sebox.png" alt="Sebox" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Crx.png" alt="Crx" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Nefm.png" alt="Nefm" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sfrp2.png" alt="Sfrp2" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Slc17a7.png" alt="Slc17a7" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Cx3cr1.png" alt="Cx3cr1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Notch1.png" alt="Notch1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Slc1a3.png" alt="Slc1a3" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Ebf3.png" alt="Ebf3" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Nr2e3.png" alt="Nr2e3" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Slc6a9.png" alt="Slc6a9" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Elavl3.png" alt="Elavl3" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Nrl.png" alt="Nrl" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sncg.png" alt="Sncg" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sox11.png" alt="Sox11" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Emx1.png" alt="Emx1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Onecut1.png" alt="Onecut1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sox2.png" alt="Sox2" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Fgf15.png" alt="Fgf15" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Onecut2.png" alt="Onecut2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Sox9.png" alt="Sox9" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Opn1mw.png" alt="Opn1mw" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Tfap2a.png" alt="Tfap2a" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Gad1.png" alt="Gad1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Opn1sw.png" alt="Opn1sw" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Tfap2b.png" alt="Tfap2b" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Gad2.png" alt="Gad2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Otx2.png" alt="Otx2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Thrb.png" alt="Thrb" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Gfap.png" alt="Gfap" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pax2.png" alt="Pax2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Thy1.png" alt="Thy1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_GFP.png" alt="GFP" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pax6.png" alt="Pax6" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Tie1.png" alt="Tie1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Gli1.png" alt="Gli1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pcp4.png" alt="Pcp4" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Trpm1.png" alt="Trpm1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Glul.png" alt="Glul" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pdgfra.png" alt="Pdgfra" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Vim.png" alt="Vim" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Gnat2.png" alt="Gnat2" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Pecam1.png" alt="Pecam1" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Vsx1.png" alt="Vsx1" width="33%">
<img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Grm6.png" alt="Grm6" width="33%"><img src="figures/umap_reClustered_refined_doubletsRemoved_threshold0.8_neurog2_Vsx2.png" alt="Vsx2" width="33%">


### Dot plot for major celltypes and marker genes 

![Dot Plot](figures/dotplot__reClustered_refined_doubletsRemoved_threshold0.8_neurog2_markerGenes.png)

## Annotations 

![legend_annotations](figures/umapreclustered_refined_doubletsRemoved_threshold0.8_neurog2_annotations.png)
![ON_annotations](figures/umapreclustered_refined_doubletsRemoved_threshold0.8_neurog2_annotationsON.png)
## !!!! TO BE FIXED 
![Dot plot marker genes](figures/dotplot__annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_markerGenes.png)


## UMAP per sample 
![Control](figures/umap_reclustered_control_2mo.png)
![Neurog2_9SA_5weeks](figures/umap_reclustered_Neurog2_9SA_5weeks.png)
![Neurog2_9S1_2mo](figures/umap_reclustered_Neurog2_9SA_2mo.png)
![UMAP_All](figures/umapreclustered_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.png)

## Cell Ratio 

![CellRatio Reversed](figures/Reversed_stacked_bar_sample_by_celltype.png)


#### Gene Expression Analysis 

### 🎯🎯 Using t-test method: according to scanpy documentation, logfold change is calculated when t-test methods are used 

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top5Genes_all_clusterttest.png )

✅ ✅ [t-test Gene Expressions](https://docs.google.com/spreadsheets/d/19YY4ErDH-bcsntzXxDzc-D7wrKbziCzvH7x4rav1SbQ/edit?usp=sharing)


### 🎯🎯 Using default parameters of scanpy which doesn't calculate log foldchange

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top5Genes_all_clusterDefault.png)
 
✅ ✅ [Default Parameters Gene Expressions](https://docs.google.com/spreadsheets/d/15ME9IKEDl7INO-U6jz7JE6d_yVFePF824eoymJG3QWE/edit?usp=sharing) 



### 🎯🎯 Using Wilcoxon method 

![heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top5Genes_all_clusterwl.png)


✅ ✅ [Wilcoxon Gene Expressions](https://docs.google.com/spreadsheets/d/1Xkz7XOfQqkARZs5gnu4-vR8OwBl62RIJ0cb3jYEn1oQ/edit?usp=sharing)


### 🎯🎯 Using T-test Control vs treatments 

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top10Genes_all_clusterttestG.png) 

✅ ✅ [ttest Gene Expressions](https://docs.google.com/spreadsheets/d/1DGXnW9RvCnNbShXm3ecENHvq8_y57AAPmCFAiawWzWM/edit?usp=sharing)



### 🎯🎯 Using Scanpy default parameters Control vs treatments using MG only

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top10Genes_all_clusterDefaultG.png)


✅ ✅ [Default Gene Expressions](https://docs.google.com/spreadsheets/d/1Svw2Cc_LFwLNCNnfkRx8EegpttJcuUjn00mHoSuPZbU/edit?usp=sharing)


### 🎯🎯 Using t-test Control vs treatments 

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top10Genes_all_clusterttestG.png)

✅ ✅ [ttest Gene Expressions](https://docs.google.com/spreadsheets/d/15lF52EWJKENZb8dBTOEYmIaYLp0NJFJTW5XGaMbD4EA/edit?usp=sharing)


### 🎯🎯 Using Wilcoxon Control vs treatments 

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top10Genes_all_clusterwlG.png)

✅ ✅ [Wilcoxon Gene Expressions](https://docs.google.com/spreadsheets/d/1MV_jJk67ma8oEg4DDRkRBfIscvdIokmLa2N90W27-a8/edit?usp=sharing)


## Filter genes expressed less than 10% before proceeding with the analysis 
### Using t-test 
#### Look at `expressionG10.py` Python Script

![Heatmap Expression](figures/heatmap_annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2_Top10Genes_all_clusterttestG10.png) 

✅ ✅ [ttest Gene Expressions](https://docs.google.com/spreadsheets/d/1B6B6qDTJoKtObiHOSZRXomW-1GvsQnACPcZWqNaShXg/edit?usp=sharing)


# ✨✨✨ Running PyScenic for transcription factors network 

### ✅ Loom File Summary

**Basic Info:**  
- Shape (genes × cells): **(2354, 44755)**  

**🧬 Row Attributes (Gene metadata):**  
- Attributes: `['Gene']`  
- Example gene names: `['Rgs20', 'Adhfe1', 'Ppp1r42', 'Prex2', 'Sbspon']`  

**🧫 Column Attributes (Cell metadata):**  
- Attributes: `['CellID', 'nGene', 'nUMI']`  
- Example CellIDs: `['AAACCAAAGCCATACA-1', 'AAACCCGCAATCCGTC-1', 'AAACCCGCATCACTGC-1', 'AAACCCGCATCGTACC-1', 'AAACCCTGTTGCTGTG-1']`  

**🔢 Expression Matrix Summary:**  
- Shape: **(44755, 2354)**  
- Data type: `float32`  
- Non-zero entries: **105,353,270**


## How to run Snakemake 

For dry run to check everything before actual run:

    snakemake -j1 -p --configfile config.yaml -n

For Actual run:

    snakemake -j1 -p --configfile config.yaml


## References

- **Scanpy**  
  Wolf, F. A., Angerer, P., & Theis, F. J. (2018).  
  *Scanpy: large-scale single-cell gene expression data analysis*. Genome Biology, 19(1), 15.  
  https://doi.org/10.1186/s13059-017-1382-0

- **Scrublet**  
Wolock, S. L., Lopez, R., & Klein, A. M. (2019).  
*Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data*. Cell Systems, 8(4), 281–291.e9.  
https://doi.org/10.1016/j.cels.2018.11.005

- **DoubletDetection**  
  Gayoso, A., Shor, J., Carr, A. J., & Yosef, N. (2019).  
  *DoubletDetection: Computational doublet detection in single-cell RNA sequencing data using boosting algorithms*.  
  [GitHub Repository](https://github.com/JonathanShor/DoubletDetection)  
  *(No peer-reviewed publication; software citation based on GitHub authorship.)*



