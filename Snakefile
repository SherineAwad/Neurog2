with open(config['SAMPLES']) a≈s fp:
    samples = fp.read().splitlines()
SUBSET = ['Cones', 'AC'] 

rule all:
         input:
            expand("{all}.h5ad", all= config['ALL']), 
            expand("doubletScores_0.8_{all}.h5ad", all=config['ALL']),
            expand("ddanalysed_doubletScores_0.8_{all}.h5ad", all =config['ALL']),
	    expand("clustered_ddanalysed_doubletScores_0.8_{all}.h5ad", all=config['ALL']), 
            expand("doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_{all}.h5ad", all=config['ALL']), 
            expand("refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL']), 
            expand("reclustered_refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL']),
 	    expand("annotated_reclustered_refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL'])

rule preprocess: 
        input:  
            expand("{sample}_filtered_feature_bc_matrix.h5", sample = samples) 
        output: 
          expand("{all}.h5ad", all= config['ALL']), 
        params: 
          samplesFile = config['SAMPLES'],  
          name = config['ALL']
        shell: 
            """
           python preprocess.py {params.samplesFile}  {params.name}  
           """ 

rule detectDoublets: 
     input: 
         expand("{all}.h5ad", all= config['ALL'])
     params: 
         doublet_cutoff=0.8
     output: 
         expand("doubletScores_0.8_{all}.h5ad", all=config['ALL'])
     shell: 
        "python detect_doublets.py {input} {params}" 

rule analyse: 
     input:
         expand("doubletScores_0.8_{all}.h5ad", all=config['ALL']) 
     output: 
         expand("ddanalysed_doubletScores_0.8_{all}.h5ad", all =config['ALL'])
     shell: 
        """ 
        python analyseDD.py {input} 
        """ 


rule cluster: 
       input:
          expand("analysed_{all}.h5ad", all=config['ALL']) 
       params: 
          markers = config['MARKERS']  
       output:
          expand("clustered_ddanalysed_doubletScores_0.8_{all}.h5ad", all=config['ALL'])
       shell:
          """
          python cluster.py {input} {params}
          """

rule doubletRemoval: 
      input: 
          expand("reclustered_clustered_analysed_{all}.h5ad", all=config['ALL'])
      output: 
          expand("doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_{all}.h5ad", all=config['ALL'])
      params: 
         markers = config['MARKERS'],
         doublet_cutoff=0.8         
      shell: 
          """
          python removeDoublets.py clustered_ddanalysed_doubletScores_0.8_neurog2.h5ad {params.marakers} {params.doublet_cutoff}
          """

rule refine:
     input:
         expand("doubletsRemoved_threshold0.8_clustered_ddanalysed_doubletScores_0.8_{all}.h5ad", all=config['ALL'])
     output:       
        expand("refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL']),   
     shell:
       """
       python plot_refine.py {input} 
       """      


rule reCluster: 
       input:
          expand("refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL']),   
       output:    
          expand("reclustered_refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL'])
       params: 
          markers = config['MARKERS']
       shell: 
          """ 
          python reCluster.py {input} {params} 
          """

rule annotate: 
    input: 
       expand("reclustered_refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL'])
    output:
      expand("annotated_reclustered_refined_doubletsRemoved_threshold0.8_{all}.h5ad", all=config['ALL'])
    params: 
        annofile = config['ANNOFILE'] 
    shell:
         """
         python annotate.py {input} {params} 
         """ 
