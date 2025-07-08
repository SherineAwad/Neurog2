with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()
SUBSET = ['Cones', 'AC'] 

rule all:
         input:
            expand("{all}.h5ad", all= config['ALL']), 
            expand("analysed_{all}.h5ad", all=config['ALL']),
	    expand("clustered_analysed_{all}.h5ad", all=config['ALL'])
 
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

rule analyse: 
     input: 
         expand("{all}.h5ad", all=config['ALL'])
     output: 
         expand("analysed_{all}.h5ad", all=config['ALL'])
     shell: 
        """ 
        python analyse.py {input} 
        """ 


rule cluster: 
       input:
          expand("analysed_{all}.h5ad", all=config['ALL']) 
       params: 
          markers = config['MARKERS']  
       output:
          expand("clustered_analysed_{all}.h5ad", all=config['ALL'])
       shell:
          """
          python cluster.py {input} {params} 
          """

rule annotate:
       input:
          expand("clustered_{all}.h5ad", all=config['ALL'])
       params: 
          annofile = config['ANNOFILE'] 
       output:
          expand("annotated_{all}.h5ad", all=config['ALL'])
       shell:
          """
          python annotate.py {input} {params.annofile} 
          """

      
