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
  make_option(c("-o","--out"),help="out dir",default="cluster_rename_dir")
  )

opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#opt$rds="myocytes.rds"
#opt$cluster="cluster_anno.txt"

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
out_dir=opt$out
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

#load anno data
seurat_anno <- read.table(opt$cluster,col.names = c("cluster","anno"))

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

saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
