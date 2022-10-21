library(Seurat)
library(gplots)
library(ggplot2)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
# 输入数据为不同亚群的差异基因
KEGG_anno = function(object  = NULL){
    ids=bitr(object$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
    object=merge(object,ids,by.x='gene',by.y='SYMBOL')
    gcSample=split(object$ENTREZID, object$cluster) 
    xx <- compareCluster(gcSample,
      fun = "enrichKEGG",
      organism = "hsa", pvalueCutoff = 0.05
    )
    #return(xx)
    ## xx = simplify(xx) only for GO
    options(repr.plot.height = 10,repr.plot.width = 7.5)
    dotplot(xx,showCategory=5) + theme(axis.text.x = element_text(
      angle = 45,
      vjust = 0.5, hjust = 0.5
    ))
    # return(p)
    data_GO_sim_fil <- xx@compareClusterResult
    write.table(data_GO_sim_fil,paste0(outdir,'/KEGG_result.txt',quote=F,sep='\t'))
    ggsave(paste0(outdir,"/","KEGG.pdf"),width =7.5, height = 10)
    ggsave(paste0(outdir,"/","KEGG.png"),width =7.5, height = 10)   
}

GO_anno = function(object  = NULL){
    ids=bitr(object$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
    object=merge(object,ids,by.x='gene',by.y='SYMBOL')
    gcSample=split(object$ENTREZID, object$cluster) 
   xx <- compareCluster(gcSample,
      fun = "enrichGO",
      OrgDb = "org.Hs.eg.db",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05
    )
    xx = simplify(xx)
    options(repr.plot.height = 10,repr.plot.width = 7.5)
    dotplot(xx,showCategory=5) + theme(axis.text.x = element_text(
      angle = 45,
      vjust = 0.5, hjust = 0.5
    ))
    # return(p)
    data_GO_sim_fil <- xx@compareClusterResult
    write.table(data_GO_sim_fil,paste0(outdir,'/GO_result.txt',quote=F,sep='\t'))
    ggsave(paste0(outdir,"/","GO.pdf"),width =7.5, height = 10)
    ggsave(paste0(outdir,"/","GO.png"),width =7.5, height = 10)   
    }