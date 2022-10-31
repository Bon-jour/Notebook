#A scripts to add cell anno to seurat cluster
options(bitmapType='cairo')
#setwd("~/work/test/")
print("load package...")
library(optparse)
library(Seurat)
library(clusterProfiler)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)
option_list <- list(
  make_option(c("-r", "--rds"), help="",default="All_sample_combined.rds"),
  make_option(c("-a","--cluster"),help="cluster_file",default="cluster.txt"),
  make_option(c("-o","--out"),help="out dir",default="cluster_rename_dir"),
  make_option(c("-t","--type"),help="type eg. hsa mmu ..",default="hsa")
  )
source("/Business/psn_company/sc02/1.source/enrichment2.r")
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

run_cluster<-function(immune.combined,seurat_exp_cluster_dir,type,idents){
  if(!dir.exists(seurat_exp_cluster_dir)){
    dir.create(seurat_exp_cluster_dir,recursive = T)
  }
  DefaultAssay(immune.combined) <- "RNA"
  Idents(immune.combined)<-idents
  immune.combined <- ScaleData(immune.combined, feature=rownames(immune.combined),verbose = FALSE)
  
  markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(markers,paste(seurat_exp_cluster_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
  
  #cluster summary
  marker_number<-table(markers$cluster)%>% reshape2::melt()
  colors=hue_pal()(nrow(marker_number))
  marker_number$Cluster <- as.factor(marker_number$Var1)
  ggplot(data=marker_number,mapping=aes(x=Cluster,y=value,fill=Cluster))+
    geom_bar(stat="identity",width=0.5)+theme_classic()+ylab("Number")+
    geom_text(aes(label = value),size = 3,hjust = 0.5,vjust = -0.5, position = "stack")
  
  
  ggsave(paste(seurat_exp_cluster_dir,"marker_number.pdf",sep="/"),width = 14,height = 7)
  ggsave(paste(seurat_exp_cluster_dir,"marker_number.png",sep="/"),width = 14,height = 7)
  
  
  
  all_top10_markers=markers %>%  top_n(n = 50, wt = avg_log2FC) %>%dplyr::distinct(.,gene,.keep_all = T) %>% top_n(n = 10, wt = avg_log2FC)
  VlnPlot(immune.combined, features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)
  
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7)
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7)
  
  DotPlot(immune.combined, features = all_top10_markers$gene)
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_dotplot.pdf",sep="/"),width = 20,height = 7)
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_dotplot.png",sep="/"),width = 20,height = 7)
  
  FeaturePlot(immune.combined, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)
  if(length(sample_list)>1){
    DotPlot(immune.combined, features = rev(unique(all_top10_markers$gene))) + RotatedAxis()}else{
      DotPlot(immune.combined, features = rev(unique(all_top10_markers$gene)))
    }
  
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15)
  ggsave(paste(seurat_exp_cluster_dir,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15)
  
  
  
  cluster_top10_markers=markers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)
  #immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  DoHeatmap(immune.combined, features = rev(unique(cluster_top10_markers$gene))) + NoLegend()+scale_fill_gradientn(colors = c("blue", "white", "red"))
  
  ggsave(paste(seurat_exp_cluster_dir,"all_cluster_markers_heatmap.pdf",sep="/"),width = 20,height = 14)
  ggsave(paste(seurat_exp_cluster_dir,"all_cluster_markers_heatmap.png",sep="/"),width = 20,height = 14)
  
  for( clust_num in  unique(Idents(immune.combined))){
    cluster_dir=paste(seurat_exp_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
    if(!file.exists(cluster_dir)){
      dir.create(cluster_dir)
    }
    cluster_dir_enrich=paste(cluster_dir,"enrichment",sep="/")
    if(!file.exists(cluster_dir_enrich)){dir.create(cluster_dir_enrich)}
    cluster_markers=subset(markers,cluster==clust_num)
    rownames(cluster_markers)<-cluster_markers$gene
    if(nrow(cluster_markers)>1){
      genelist=cluster_markers$gene
      enrichment(species=type,outDir=cluster_dir_enrich,geneList=genelist)
      write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
      top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)
      VlnPlot(immune.combined, features = top10_markers$gene,pt.size = 0.1 ,ncol=5)
      
      ggsave(paste(cluster_dir,"top10_vilion.pdf",sep="/"),width =20,height = 7)
      ggsave(paste(cluster_dir,"top10_vilion.png",sep="/"),width =20,height = 7)
      
      FeaturePlot(immune.combined, features = top10_markers$gene, min.cutoff = "q9",ncol=5)
      ggsave(paste(cluster_dir,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
      ggsave(paste(cluster_dir,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)
      if(length(sample_list)>1){
        DotPlot(immune.combined, features = rev(unique(top10_markers$gene)), cols = rainbow(length(sample_list))) + RotatedAxis()}else{
          DotPlot(immune.combined, features = rev(unique(top10_markers$gene)))}
      
      ggsave(paste(cluster_dir,"top10_exp_pct.pdf",sep="/"),width = 15,height =15)
      ggsave(paste(cluster_dir,"top10_exp_pct.png",sep="/"),width = 15,height =15)
      gc(TRUE)
      
    }}
}





#opt$rds="myocytes.rds"
#opt$cluster="cluster_anno.txt"

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
out_dir=opt$out
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

#load anno data
seurat_anno <- read.table(opt$cluster,sep="\t",header=F,col.names = c("cluster","anno"))

new.cluster.ids <-as.vector(seurat_anno$anno)

names(new.cluster.ids) <-seurat_anno$cluster

seurat_obj@active.ident <- seurat_obj$seurat_clusters
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$clutser_anno <- Idents(seurat_obj)
p1=DimPlot(seurat_obj,reduction = "umap")
p2=DimPlot(seurat_obj,group.by = "seurat_clusters",reduction = "umap")
plot_grid(p1,p2)
ggsave(filename = paste(out_dir,"Cluster_anno.umap.pdf",sep = "/"),width = 14,height = 7)
ggsave(filename = paste(out_dir,"Cluster_anno.umap.png",sep = "/"),width = 14,height = 7)

p1=DimPlot(seurat_obj,reduction = "tsne")
p2=DimPlot(seurat_obj,group.by = "seurat_clusters",reduction = "tsne")
plot_grid(p1,p2)
ggsave(filename = paste(out_dir,"Cluster_anno.tsne.pdf",sep = "/"),width = 14,height = 7)
ggsave(filename = paste(out_dir,"Cluster_anno.tsne.png",sep = "/"),width = 14,height = 7)

anno <- data.frame(barcode=rownames(seurat_obj@meta.data),anno=seurat_obj$clutser_anno)
write.table(anno,file=paste0(out_dir,"/anno.txt"),sep="\t",row.names = F)

sample_list=unique(seurat_obj$sample)

run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(out_dir,"/marker_diff_cluster"),type = opt$type, idents = "clutser_anno")



saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
