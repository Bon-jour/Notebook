#plot genelist 
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
  make_option(c("-a","--cluster"),help="cluster_file"),
  make_option(c("-o","--out"),help="out dir",default="show_genes_dir"),
  make_option(c("-g","--genes"),help="gene list",default = "genelist")
)

opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#1 subset cluster 
if(!dir.exists(opt$out)){
dir.create(opt$out)
}


#opt$res=0.6
#opt$rds="myocytes.rds"
#opt$cluster="cluster.txt"
#opt$workdir="workdir"

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
DefaultAssay(seurat_obj)<-"RNA"
if(is.null(opt$cluster)){
  print("use all cluster")
}else{
  
  cluster <- read.table(opt$cluster,col.names = "cluster")
  print("subset data")
  seurat_obj <-subset(seurat_obj,seurat_clusters%in%c(cluster$cluster))
  
}

#read gene list 
genes <- read.table(opt$genes,col.names = "genes")
seurat_obj <- ScaleData(seurat_obj,features = rownames(seurat_obj))
for (gene in genes$genes) {
  p <- FeaturePlot(seurat_obj,features = gene)
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".featureplot.pdf"),width = 8,height = 7)
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".featureplot.png"),width = 8,height = 7,dpi=300)
  p <- VlnPlot(seurat_obj,features = gene,group.by ="seurat_clusters")
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".VlnPlot.pdf"),width = 8,height = 7)
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".VlnPlot.png"),width = 8,height = 7,dpi=300)
  p<-VlnPlot(seurat_obj,features = gene,split.by = "group",group.by ="seurat_clusters")
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".VlnPlot.group.pdf"),width = 8,height = 7)
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".VlnPlot.group.png"),width = 8,height = 7,dpi=300)
  p <- RidgePlot(seurat_obj, features = gene)
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".RidgePlot.pdf"),width = 8,height = 7)
  ggsave(plot = p,filename = paste0(opt$out,"/",gene,".RidgePlot.png"),width = 8,height = 7,dpi=300)
  #p <- RidgePlot(seurat_obj, features = gene,split.by = "group")
  #ggsave(plot = p,filename = paste0(opt$out,"/",gene,".RidgePlot.group.pdf"),width = 8,height = 7)

}


out_width=ceiling(length(genes$genes)/20)*8
p <- DotPlot(seurat_obj,features = c(unique(genes$genes)))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(plot = p,filename = paste0(opt$out,"/","all.DotPlot.pdf"),width = out_width,height = 7)
ggsave(plot = p,filename = paste0(opt$out,"/","all.DotPlot.png"),width = out_width,height = 7,dpi=300)

p <- DoHeatmap(seurat_obj,features = c(unique(genes$genes)),size=5)+ scale_fill_gradientn(colors = c("navy","white","firebrick3"))
ggsave(plot = p,filename = paste0(opt$out,"/","all.Heatmap.pdf"),width = out_width,height = 14)
ggsave(plot = p,filename = paste0(opt$out,"/","all.Heatmap.png"),width = out_width,height = 14,dpi=300)
print("done")
