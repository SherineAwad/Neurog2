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

- Number of genes detected per cell between **800 and 8000**
- Total counts per cell between **1200 and 30000**
- Percentage of mitochondrial gene counts less than **25%**

This filtering step ensures removal of dead or dying cells and technical artifacts to improve downstream analysis quality.

## Clustering 

#### Marker Gene UMAP Plots
Below are the UMAP visualizations of marker gene expression across clusters. These are auto-generated from your data and saved in the figures/ directory.


![UMAP CLUSTERS](figures/umap_clusters.png)

<img src="figures/umapclustered_analysed_neurog2_Malat1.png" alt="Plot 1" width="33%"> <img src="figures/umapclustered_analysed_neurog2_mt-Atp6.png" alt="Plot 2" width="33%"> <img src="figures/umapclustered_analysed_neurog2_Sox9.png" alt="Plot 3" width="33%">


<div style="text-align: left;">

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Malat1.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_mt-Atp6.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Sox9.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Glul.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Lhx2.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Rlbp1.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Rbfox3.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Csf1r.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Calb2.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Elavl4.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Calb1.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Sebox.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Gad1.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Elavl3.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Cabp5.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Isl1.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Slc6a9.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Ascl1.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Olig2.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Foxn4.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Chat.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Prdm1.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Otx2.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Insm1.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Sox11.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Atoh7.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Hes5.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Emx1.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_mScarlet3.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_GFP.png" style="width: 100%; height: auto;" />
  </div>

  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Neurog2.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-right: 1%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Tfap2a.png" style="width: 100%; height: auto;" />
  </div>
  <div style="display: inline-block; width: 31%; margin-bottom: 10px;">
    <img src="figures/umapclustered_analysed_neurog2_Bsn.png" style="width: 100%; height: auto;" />
  </div>

</div>


#### QC per Clsuter 

<img src="figures/qc_violin_by_cluster.png" width="550"/>

## Removing low quality clustering and Reclustering 



---

*Generated with Scanpy for single-cell RNA-seq analysis.*

