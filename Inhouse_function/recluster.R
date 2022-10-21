library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

recluster <- function(t,resolution) {
    t=CreateSeuratObject(counts = t@assays$RNA@counts,
                       meta.data = t@meta.data) 
    t <- NormalizeData(t, normalization.method = "LogNormalize", scale.factor = 1e4) 
    t <- FindVariableFeatures(t, selection.method = 'vst', nfeatures = 2000)
    t <- ScaleData(t)
    t <- RunPCA(t, features = VariableFeatures(object = t)) 
    
    t <- FindNeighbors(t, dims = 1:10)
    t <- FindClusters(t, resolution = resolution)
     
    set.seed(253)
    t <- RunUMAP(t, dims = 1:10,do.fast = TRUE)
    p1 = DimPlot(t,reduction = "umap",repel = T) + theme(legend.position = 'none')
    
    ggsave(paste0(outdir,"/","recluster.uamp.pdf"),width =5, height = 5)
    ggsave(paste0(outdir,"/","recluster.uamp.png"),width =5, height = 5)
    saveRDS(t,paste0(outdir,"/recluster.rds"))
}