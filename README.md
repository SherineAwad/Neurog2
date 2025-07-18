# Neurog2 Project - scRNA-seq Analysis

This project focuses on the single-cell RNA sequencing (scRNA-seq) analysis of samples related to **Neurog2** expression at different developmental stages and controls. The goal is to understand the cellular heterogeneity and transcriptional changes associated with Neurog2 in specific neuronal populations.

## Samples

| Sample Name           | Description          |
|-----------------------|----------------------|
| 5 weeks Neurog2_9SA   | TH1_GFP_mScarlet3    |
| 2 months control      | TH2_GFP_mScarlet3    |
| 2 months Neurog2_9SA  | TH3_GFP_mScarlet3    |

## Analysis Workflow

The analysis was performed using [Scanpy](https://scanpy.readthedocs.io/en/stable/), a scalable toolkit for analyzing single-cell gene expression data. The workflow included:

- Reading and normalizing raw data
- Identifying highly variable genes
- Scaling data
- Performing PCA and neighborhood graph calculation
- Computing UMAP embeddings for visualization
- Quality control and visualization

## Figures

### 1. UMAP Visualization
![UMAP of Neurog2 samples](figures/umap_neurog2.png)  
UMAP plot colored by sample, showing clustering and distribution of single cells from different conditions.


## Per sample UMAP 

![Control 2mo](figures/umap_sample_control_2mo.png)
![Neurog2_9SA_2mo](figures/umap_sample_Neurog2_9SA_2mo.png)
![Neurog2_9SA_5weeks](figures/umap_sample_Neurog2_9SA_5weeks.png)

### 2. Quality Control Violin Plot
![Before Filtering QC metrics](figures/violin_QC.png)  
Violin plots displaying quality control metrics such as number of genes detected per cell, total counts, and percentage of mitochondrial gene expression.

### 3. Additional Analysis Figure
![After Filtering QC metrics](figures/violin_AfterQC.png)  

## Filtering Criteria

Quality filtering was applied to remove low-quality cells and potential doublets. Cells were retained only if they met all the following criteria:

- Number of genes detected per cell between **800 and 8000**
- Total counts per cell between **1200 and 30000**
- Percentage of mitochondrial gene counts less than **25%**

This filtering step ensures removal of dead or dying cells and technical artifacts to improve downstream analysis quality.


## Number of cells per sample 

| Sample              | Cell Count |
|---------------------|------------|
| Neurog2_9SA_5weeks  | 27,732     |
| Neurog2_9SA_2mo     | 11,486     |
| control_2mo         | 9,701      |



### 4. Clustering 

## Marker Gene UMAP Plots
Below are the UMAP visualizations of marker gene expression across clusters. These are auto-generated from your data and saved in the figures/ directory.


![UMAP CLUSTERS](figures/umap_clusters.png)

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



## QC per Clsuter 

<img src="figures/qc_violin_by_cluster.png" width="550"/>

### 5. Removing low quality clustering and Reclustering 

We removed low quality clusters number:  ['7', '8', '11', '20', '28', '33', '34']

then we reclustered and replot the marker genes as below: 


## UMAP

![UMAP RE CLUSTERS](figures/umap_reClusters.png) 


## Per sample UMAP 

![Control 2mo](figures/umap_reclustered_control_2mo.png)
![Neurog2_9SA_2mo](figures/umap_reclustered_Neurog2_9SA_2mo.png)
![Neurog2_9SA_5weeks](figures/umap_reclustered_Neurog2_9SA_5weeks.png)

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



## QC per Cluster 

<img src="figures/qc_violin_by_reCluster.png" width="550"/>

## Number of cells per sample 

| Sample              | Cell Count |
|---------------------|------------|
| Neurog2_9SA_5weeks  | 23,370     |
| Neurog2_9SA_2mo     | 10,115     |
| control_2mo         | 8,674      |

---

## Doublet Detection with Scrublet

We are using **Scrublet**, a Python-based tool, to identify and remove potential doublets from our single-cell RNA-seq dataset.

## Understanding Doublet Scores in Scrublet

**Doublet scores in Scrublet** quantify how likely each cell is to be a **doublet**, based on how similar its gene expression profile is to simulated doublets.

---

### 🔍 In Detail

#### What is a Doublet?

A **doublet** occurs when **two cells are captured in the same droplet** during single-cell RNA sequencing. Their RNA is sequenced as if it's from one cell, producing a mixed transcriptome. This can distort downstream analyses such as clustering, dimensionality reduction, and marker gene identification.

---

### How Scrublet Works

1. **Simulates Doublets**  
   Scrublet generates **synthetic doublets** by randomly combining gene expression profiles from real cells.

2. **Embedding**  
   It runs **PCA** on both the real and synthetic cells to embed them in the same low-dimensional space.

3. **Scoring**  
   For each real cell, Scrublet calculates a **doublet score** based on its **proximity to simulated doublets** in PCA space.

---

### Interpreting the Scores

- **Doublet score range**: Typically between **0 and 1**.
- **High score (~0.5–1.0)**:  
  The cell is very similar to simulated doublets → likely a **true doublet**.
- **Low score (~0–0.2)**:  
  The cell resembles real singlets → likely a **true singlet**.

---

### Threshold for Calling Doublets

Scrublet tries to automatically find a **threshold** where the doublet score distribution separates singlets from doublets. We can:

-  Let Scrublet pick the threshold automatically (default)
- ✏️ Manually adjust the threshold based on score distribution plots


## Doublet Scores Distribution  

<img src="figures/doublet_score_histogram.png" width="550"/>


## Doublet vs Singlet UMAP using default threshold = 0.4  

<img src="figures/umap_doubletStatus.png" width="550"/>

## Doublet vs Singlet UMAP using threshold = 0.1 

<img src="figures/umap_doubletStatus0.1.png" width="550"/>


## Doublet vs Singlet UMAP using threshold = 0.18 

<img src="figures/umap_doubletStatus0.18.png" width="550"/>


## Doublet vs Singlet UMAP using threshold = 0.15 

<img src="figures/umap_doubletStatus0.15.png" width="550"/>


## Doublet Detection using `DoubletDetection`

Unlike `Scrublet`, which can operate effectively on clustered or preprocessed `AnnData` objects, the `DoubletDetection` tool is more sensitive to data structure and expects the **original, unclustered** `AnnData` object. Running it on a processed or subsetted object may yield suboptimal or misleading results.

In the workflow, we applied `DoubletDetection` to the original data (`adata`) to ensure it captures the full transcriptomic diversity and avoids artifacts introduced during clustering.

After running `DoubletDetection`, predicted doublets and doublet scores were stored in `adata.obs` under the keys:
- `predicted_doublet`: Boolean flag indicating whether each cell is a predicted doublet.
- `doublet_score`: Confidence score associated with doublet prediction.

The results were visualized using UMAP, colored by both prediction and score:

![Doublet Detection UMAP](figures/umap_doubletScores_neurog2_doublets.png)

### Doublet Scores and Conversion

   ![Doublets Thresholds](figures/threshold_test.png)
 
   ![Doublets Conversion](figures/convergence_test.png)

## Clustering after Doublet Detection 

![UMAP DOUBLET DETECTION CLUSTERS](figures/umap_ddanalysed_doubletScores_neurog2_ddClusters.png)


## Marker Genes UMAP 

### UMAPs of Gene Expression with Doublet Scores (Neurog2 Lineage)

<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Pax2.png" alt="Pax2" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Rpe65.png" alt="Rpe65" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Tie1.png" alt="Tie1" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Rho.png" alt="Rho" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Nrl.png" alt="Nrl" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Arr3.png" alt="Arr3" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Malat1.png" alt="Malat1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_mt-Atp6.png" alt="mt-Atp6" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Glul.png" alt="Glul" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Lhx2.png" alt="Lhx2" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Csf1r.png" alt="Csf1r" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Calb2.png" alt="Calb2" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Calb1.png" alt="Calb1" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Gad1.png" alt="Gad1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Elavl3.png" alt="Elavl3" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Cabp5.png" alt="Cabp5" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Isl1.png" alt="Isl1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Slc6a9.png" alt="Slc6a9" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Foxn4.png" alt="Foxn4" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Chat.png" alt="Chat" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Prdm1.png" alt="Prdm1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Otx2.png" alt="Otx2" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Insm1.png" alt="Insm1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Sox11.png" alt="Sox11" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Atoh7.png" alt="Atoh7" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Emx1.png" alt="Emx1" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_mScarlet3.png" alt="mScarlet3" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_GFP.png" alt="GFP" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Tfap2a.png" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Bsn.png" alt="Bsn" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_clustered_ddanalysed_doubletScores_neurog2_Lhx4.png" alt="Lhx4" width="33%">

Sure! Here's your explanation formatted in **Markdown**:

---

###  Understanding `doublet_score` Thresholds

The `doublet_score` typically ranges from **0 to 1**, where **higher values indicate a higher probability of a cell being a doublet**.

Your filter in the code:

```python
combined_adata = combined_adata[combined_adata.obs['doublet_score'] >= threshold]
```

This means you're **keeping** cells with `doublet_score >= threshold`.

---

###  Interpretation of Threshold:

* **Higher threshold** (e.g., `0.9`) → **Stricter filtering**
  🔹 You keep **more** cells
  🔹 Less doublets are removed

* **Lower threshold** (e.g., `0.4`) → **More relaxed filtering**
  🔹 You keep **fewer** cells
  🔹 More potential doublets are removed 

---

###  Summary:

* Relax filtering → use a **lower threshold** (e.g., `0.5 → 0.4`)
* Stricter filtering → use a **higher threshold** (e.g., `0.5 → 0.6`)



## Remove doublet cells with cutoff 0.5 

![UMAP after doublet removal](figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_doubletsRemoved_clusters.png)


## Marker Genes after doublet removal at threshold 0.5 

<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Lhx4.png" alt="Lhx4" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Pax2.png" alt="Pax2" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Rpe65.png" alt="Rpe65" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Kcnj8.png" alt="Kcnj8" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Tie1.png" alt="Tie1" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Rho.png" alt="Rho" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Nrl.png" alt="Nrl" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Arr3.png" alt="Arr3" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Malat1.png" alt="Malat1" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_mt-Atp6.png" alt="mt-Atp6" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Glul.png" alt="Glul" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Lhx2.png" alt="Lhx2" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Rbfox3.png" alt="Rbfox3" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Calb2.png" alt="Calb2" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Elavl4.png" alt="Elavl4" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Calb1.png" alt="Calb1" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Gad1.png" alt="Gad1" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Elavl3.png" alt="Elavl3" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Cabp5.png" alt="Cabp5" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Isl1.png" alt="Isl1" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Olig2.png" alt="Olig2" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Chat.png" alt="Chat" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Prdm1.png" alt="Prdm1" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Otx2.png" alt="Otx2" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Insm1.png" alt="Insm1" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Sox11.png" alt="Sox11" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Emx1.png" alt="Emx1" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_GFP.png" alt="GFP" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Neurog2.png" alt="Neurog2" width="33%">
<img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Tfap2a.png" alt="Tfap2a" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Bsn.png" alt="Bsn" width="33%"><img src="figures/umap_doubletsRemoved_clustered_doubletsDetected_neurog2_Slc17a7.png" alt="Slc17a7" width="33%">


## Remove doublet cells with cutoff 0.4

![UMAP after doublet removal](figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_doubletsRemoved_threshold0.4_clusters.png)

## Marker Genes after doublet removal at threshold 0.4 

<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Pax2.png" alt="Pax2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Rpe65.png" alt="Rpe65" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Tie1.png" alt="Tie1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Rho.png" alt="Rho" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Nrl.png" alt="Nrl" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Arr3.png" alt="Arr3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Malat1.png" alt="Malat1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_mt-Atp6.png" alt="mt-Atp6" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Glul.png" alt="Glul" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Lhx2.png" alt="Lhx2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Csf1r.png" alt="Csf1r" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Calb2.png" alt="Calb2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Calb1.png" alt="Calb1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Gad1.png" alt="Gad1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Elavl3.png" alt="Elavl3" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Cabp5.png" alt="Cabp5" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Isl1.png" alt="Isl1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Elavl4.png" alt="Elavl4" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Olig2.png" alt="Olig2" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Chat.png" alt="Chat" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Prdm1.png" alt="Prdm1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Otx2.png" alt="Otx2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Insm1.png" alt="Insm1" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Sox11.png" alt="Sox11" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Hes5.png" alt="Hes5" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Emx1.png" alt="Emx1" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_GFP.png" alt="GFP" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Tfap2a.png" alt="Tfap2a" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Bsn.png" alt="Bsn" width="33%">
<img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_doubletsRemoved_threshold0.4_clustered_doubletsDetected_neurog2_threshold0.4_Lhx4.png" alt="Lhx4" width="33%">


## Doublet Removal at Threshold 0.9 

### Checking if distribution will differ at different cutoff 


![Doublet Prediction](figures/umap_ddanalysed_doubletScores_0.9_neurog2_doublets.png)
![Doublet Conversion](figures/0.9conversion_test.png)
![Doublet Thresholds](figures/0.9threshold_test.png)


## UMAP after  clustering 

![Doublet Detection](figures/umap_clustered_ddanalysed_doubletScores_0.9_neurog2_ddClusters.png) 

## UMAP after doublet removal at threshold 0.9 
 
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


## Doublet removal using 0.8 threshold 

## UMAP after  clustering 

![Doublet Detection](figures/umap_clustered_ddanalysed_doubletScores_0.8_neurog2_ddClusters.png)

## UMAP after doublet removal at threshold 0.8 

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

## After Reclustering and QC plot  

![Doublet Reclustering](figures/umapDDreClustereddoubletScores_0.8_neurog2.png)


### Pre filtering QC
![Pre QC](figures/violin_doubletScores_0.8_neurog2_PreQC.png)

### Filter using these parameterss 

##### 📌 Filter cells based on quality thresholds

##### Keep cells with 1500–8000 detected genes
adata = adata[(adata.obs.n_genes_by_counts > 1500) & (adata.obs.n_genes_by_counts < 8000), :]

##### Keep cells with total UMI counts above 3000
adata = adata[adata.obs.total_counts > 3000, :]

#####  Keep cells with mitochondrial gene percentage below 7%
adata = adata[adata.obs.pct_counts_mt < 7, :]


### Post filtering QC

![Post QC](figures/violin_doubletScores_0.8_neurog2_PostQC.png)

## After reClustering 

![Doublet QC](figures/doubletScores_0.8_neurog2DDreClustered_qc_violin.png)
### Marker genes after reClustering 

<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Hes1.png" alt="Hes1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Otx2.png" alt="Otx2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_mScarlet3.png" alt="mScarlet3" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Sox2.png" alt="Sox2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Sox9.png" alt="Sox9" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Prkca.png" alt="Prkca" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Prox1.png" alt="Prox1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Lhx4.png" alt="Lhx4" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Thy1.png" alt="Thy1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Cabp5.png" alt="Cabp5" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Aqp4.png" alt="Aqp4" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Slc6a9.png" alt="Slc6a9" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Prdx6.png" alt="Prdx6" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Gfap.png" alt="Gfap" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Abca8a.png" alt="Abca8a" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Guca1b.png" alt="Guca1b" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Rbpms.png" alt="Rbpms" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Onecut2.png" alt="Onecut2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_GFP.png" alt="GFP" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Rho.png" alt="Rho" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Bsn.png" alt="Bsn" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Rlbp1.png" alt="Rlbp1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Sox11.png" alt="Sox11" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Opn1sw.png" alt="Opn1sw" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Neurog2.png" alt="Neurog2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Thrb.png" alt="Thrb" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Crx.png" alt="Crx" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Gnat2.png" alt="Gnat2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Cx3cr1.png" alt="Cx3cr1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Bhlhe23.png" alt="Bhlhe23" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Lhx2.png" alt="Lhx2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Opn1mw.png" alt="Opn1mw" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Igf2.png" alt="Igf2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Gad2.png" alt="Gad2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Elavl4.png" alt="Elavl4" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Notch1.png" alt="Notch1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Pecam1.png" alt="Pecam1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Calb2.png" alt="Calb2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Elavl3.png" alt="Elavl3" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Pdgfra.png" alt="Pdgfra" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Rom1.png" alt="Rom1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Sebox.png" alt="Sebox" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Vsx2.png" alt="Vsx2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Onecut1.png" alt="Onecut1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Rpe65.png" alt="Rpe65" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Prdm1.png" alt="Prdm1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Trpm1.png" alt="Trpm1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Nr2e3.png" alt="Nr2e3" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Grm6.png" alt="Grm6" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Chat.png" alt="Chat" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Tfap2a.png" alt="Tfap2a" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Acta2.png" alt="Acta2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Pou4f3.png" alt="Pou4f3" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Cldn5.png" alt="Cldn5" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Mfge8.png" alt="Mfge8" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Pvalb.png" alt="Pvalb" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Pvalb3.png" alt="Pvalb3" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Mlc1.png" alt="Mlc1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Foxp2.png" alt="Foxp2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Ttr.png" alt="Ttr" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Mpp3.png" alt="Mpp3" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Vtn.png" alt="Vtn" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Eomes.png" alt="Eomes" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Stmn2.png" alt="Stmn2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Col1a1.png" alt="Col1a1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Nckap1l.png" alt="Nckap1l" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Scgn.png" alt="Scgn" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Nefl.png" alt="Nefl" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Mmp12.png" alt="Mmp12" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Sparc.png" alt="Sparc" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Myl9.png" alt="Myl9" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Cldn4.png" alt="Cldn4" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Vtn2.png" alt="Vtn2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Gnai2.png" alt="Gnai2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Tbx21.png" alt="Tbx21" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Gjb6.png" alt="Gjb6" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_S100a10.png" alt="S100a10" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Sema3a.png" alt="Sema3a" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Atp1b2.png" alt="Atp1b2" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Pcdh15.png" alt="Pcdh15" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Slc1a2.png" alt="Slc1a2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Kcna1.png" alt="Kcna1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Slc6a1.png" alt="Slc6a1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Msn.png" alt="Msn" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Spp1.png" alt="Spp1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Acta1.png" alt="Acta1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Mgp.png" alt="Mgp" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Npy.png" alt="Npy" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Cplx2.png" alt="Cplx2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Tnnt2.png" alt="Tnnt2" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Slc32a1.png" alt="Slc32a1" width="33%">
<img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Nek7.png" alt="Nek7" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Tac1.png" alt="Tac1" width="33%"><img src="figures/umap_DDreClustered_doubletScores_0.8_neurog2_Ephx2.png" alt="Ephx2" width="33%">




### ✅ Clusters  Are More Likely to Be Low-Quality:

#### **Cluster 38**

* **n\_genes\_by\_counts**: Clear drop — low median around \~2000
* **total\_counts**: Very low, large spread
* **pct\_counts\_mt**: Slightly high (long tail to \~6%)
* **🧠 Verdict**: **Yes, should be removed**. This cluster appears consistently low-quality across all metrics.

#### **Cluster 40**

* **n\_genes\_by\_counts**: Extremely low (\~1500 median)
* **total\_counts**: Very low (\~3000–4000 median)
* **pct\_counts\_mt**: Elevated — up to \~6%
* **🧠 Verdict**: **Definitely remove** — classically low-quality.

#### **Cluster 41**

* **n\_genes\_by\_counts**: Also low
* **total\_counts**: Similar to 40, slightly higher
* **pct\_counts\_mt**: Slightly elevated
* **🧠 Verdict**: **Remove** — consistent with stressed/damaged cells.

#### **Cluster 42**

* **n\_genes\_by\_counts**: Lower end (\~2500–3000)
* **total\_counts**: Below median
* **pct\_counts\_mt**: Borderline elevated
* **🧠 Verdict**: **Optional remove** — borderline, but could be cautious and remove.



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



