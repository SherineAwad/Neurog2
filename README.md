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

### 2. Quality Control Violin Plot
![Before Filtering QC metrics](figures/violin_QC.png)  
Violin plots displaying quality control metrics such as number of genes detected per cell, total counts, and percentage of mitochondrial gene expression.

### 3. Additional Analysis Figure
![After Filtering QC metrics](figures/violin_AfterQC.png)  

## Filtering Criteria

Quality filtering was applied to remove low-quality cells and potential doublets. Cells were retained only if they met all the following criteria:

- Number of genes detected per cell between **800 and 6000**
- Total counts per cell between **1200 and 30000**
- Percentage of mitochondrial gene counts less than **25%**

This filtering step ensures removal of dead or dying cells and technical artifacts to improve downstream analysis quality.

📊 Marker Gene UMAP Plots
Below are the UMAP visualizations of marker gene expression across clusters. These are auto-generated from your data and saved in the figures/ directory.

<p align="left"> <img src="figures/umap_analysed_neurog2_Malat1.png" height="180"/> <img src="figures/umap_analysed_neurog2_mt-Atp6.png" height="180"/> <img src="figures/umap_analysed_neurog2_Sox9.png" height="180"/> <img src="figures/umap_analysed_neurog2_Glul.png" height="180"/> <img src="figures/umap_analysed_neurog2_Lhx2.png" height="180"/> <img src="figures/umap_analysed_neurog2_Rlbp1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Rbfox3.png" height="180"/> <img src="figures/umap_analysed_neurog2_Csf1r.png" height="180"/> <img src="figures/umap_analysed_neurog2_Calb2.png" height="180"/> <img src="figures/umap_analysed_neurog2_Elavl4.png" height="180"/> <img src="figures/umap_analysed_neurog2_Calb1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Sebox.png" height="180"/> <img src="figures/umap_analysed_neurog2_Gad1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Elavl3.png" height="180"/> <img src="figures/umap_analysed_neurog2_Cabp5.png" height="180"/> <img src="figures/umap_analysed_neurog2_Isl1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Slc6a9.png" height="180"/> <img src="figures/umap_analysed_neurog2_Ascl1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Olig2.png" height="180"/> <img src="figures/umap_analysed_neurog2_Foxn4.png" height="180"/> <img src="figures/umap_analysed_neurog2_Chat.png" height="180"/> <img src="figures/umap_analysed_neurog2_Prdm1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Otx2.png" height="180"/> <img src="figures/umap_analysed_neurog2_Insm1.png" height="180"/> <img src="figures/umap_analysed_neurog2_Sox11.png" height="180"/> <img src="figures/umap_analysed_neurog2_Atoh7.png" height="180"/> <img src="figures/umap_analysed_neurog2_Hes5.png" height="180"/> <img src="figures/umap_analysed_neurog2_Emx1.png" height="180"/> <img src="figures/umap_analysed_neurog2_mScarlet3.png" height="180"/> <img src="figures/umap_analysed_neurog2_GFP.png" height="180"/> <img src="figures/umap_analysed_neurog2_Neurog2.png" height="180"/> <img src="figures/umap_analysed_neurog2_Tfap2a.png" height="180"/> <img src="figures/umap_analysed_neurog2_Bsn.png" height="180"/> <img src="figures/umap_analysed_neurog2_Slc17a7.png" height="180"/> <img src="figures/umap_analysed_neurog2_Lhx4.png" height="180"/> </p>



---

*Generated with Scanpy for single-cell RNA-seq analysis.*

