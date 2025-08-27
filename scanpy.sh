#!/bin/bash
#SBATCH --job-name neurog2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem=80000
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=sherinem@umich.edu
#SBATCH --account=thahoang99 


conda activate scanpy_solo_env

#python cluster.py analysed_neurog2.h5ad markers.txt
#python reCluster.py clustered_analysed_neurog2.h5ad markers.txt
#python doublet.py reclustered_clustered_analysed_neurog2.h5ad 


##New using doublet Detection 

#python detect_doublets.py  neurog2.h5ad 0.5 
python analyseDD.py doubletScores_0.5_neurog2.h5ad
python cluster.py ddanalysed_doubletScores_0.5_neurog2.h5ad markers.txt 
python removeDoublets.py clustered_ddanalysed_doubletScores_0.8_neurog2.h5ad markers.txt 0.8 

python plot_refine.py doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_neurog2.h5ad
python reCluster.py refined_doubletsRemoved_threshold0.8_neurog2.h5ad mouse_markers.txt 
python annotate.py reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad annotations.txt
#python plots.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad


#python expressionG10.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad 
#python expressionG.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python expression_defaultG.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python expression_wilcoxonG.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad 

#python expression_default.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python expression.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python expression_wilcoxon.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad

#python gexp.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python gexp2.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python gexp3.py annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#python  expression_fwilcoxon.py  annotated_reclustered_refined_doubletsRemoved_threshold0.8_neurog2.h5ad
#Rscript h5adToSeurat.R 

