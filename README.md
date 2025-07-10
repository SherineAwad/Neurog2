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


![UMAP RE CLUSTERS](figures/umap_reClusters.png)

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



## QC per Clsuter 

<img src="figures/qc_violin_by_reCluster.png" width="550"/>

## Number of cells per sample 

| Sample              | Cell Count |
|---------------------|------------|
| Neurog2_9SA_5weeks  | 23,370     |
| Neurog2_9SA_2mo     | 10,115     |
| control_2mo         | 8,674      |

---

*Generated with Scanpy for single-cell RNA-seq analysis.*

