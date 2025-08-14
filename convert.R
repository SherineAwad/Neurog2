library(SeuratDisk)
library(Seurat)


Convert("neurog2_clean.h5ad", dest = "h5seurat", overwrite = TRUE)
obj <- LoadH5Seurat("neurog2_clean.h5seurat",meta.data = FALSE, misc = FALSE)




