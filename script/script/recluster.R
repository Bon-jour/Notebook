#run subcluster
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
source("/Business/psn_company/sc02/1.source/enrichment1.r")

option_list <- list(
  make_option(c("-r", "--rds"), help="",default="All_sample_combined.rds"),
  make_option(c("-c","--cluster"),help="cluster_file",default="cluster.txt"),
  make_option(c("-s","--res"),help="find cluster resolution",default=0.6),
  make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory",default = "Reclust"),
  make_option(c("-t", "--type"), help="species type",default = "hsa")
)

opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#opt$res=0.6
#opt$rds="myocytes.rds"
#opt$cluster="cluster.txt"
#opt$workdir="workdir"
#opt$type="hsa"

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
cluster <- read.table(opt$cluster,col.names = "cluster")
print("subset data")
seurat_obj <-subset(seurat_obj,seurat_clusters%in%c(cluster$cluster))

seurat_obj <- ScaleData(seurat_obj,feature=rownames(seurat_obj), verbose = FALSE)
seurat_obj <- FindVariableFeatures(object = seurat_obj,selection.method = 'vst', nfeatures = 2000)
seurat_obj <- RunPCA(seurat_obj,  features = VariableFeatures(object = seurat_obj) ,verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- FindClusters(seurat_obj,resolution = opt$res)
file_out=opt$workdir
#plot fig
if(dir.exists(file_out)){
  print("dir exists")
}else{
  dir.create(file_out)
}

cluster_summary=dcast(as.data.frame(table(data.frame("cluster"=Idents(seurat_obj),"sample"=seurat_obj$sample))),sample~cluster)
write.table(cluster_summary,paste(file_out,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

cluster_number<-table(seurat_obj$seurat_clusters)%>% reshape2::melt()
colors=hue_pal()(nrow(cluster_number))
cluster_number$Cluster <- as.factor(cluster_number$Var1)
p <- ggplot(data=cluster_number,mapping=aes(x=Cluster,y=value,fill=Cluster))+
  geom_bar(stat="identity",width=0.5)+theme_classic()+ylab("Number")+
  geom_text(aes(label = value),size = 3,hjust = 0.5,vjust = -0.5, position = "stack")
ggsave(plot = p,paste(file_out,"cluster_summary.pdf",sep="/"),height = 7,width = 9)
ggsave(plot = p,paste(file_out,"cluster_summary.png",sep="/"),height = 7,width = 9)

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
p2 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave(paste(file_out,"cluster_umap.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(file_out,"cluster_umap.png",sep="/"),width = 14,height = 7)

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
p2 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave(paste(file_out,"cluster_umap.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(file_out,"cluster_umap.png",sep="/"),width = 14,height = 7)

p1 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)
ggsave(paste(file_out,"cluster_tsne.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(file_out,"cluster_tsne.png",sep="/"),width = 14,height = 7)


avg_per_cluster=AverageExpression(seurat_obj,"RNA")$RNA
colnames(avg_per_cluster)=gsub("RNA.","",colnames(avg_per_cluster),perl=T)
write.table(avg_per_cluster,paste(file_out,"avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)


#find marker
print("Do find markers")
file_out_marker <- paste0(file_out,"/Marker")
#plot fig
if(dir.exists(file_out_marker)){
  print("dir exists")
}else{
  dir.create(file_out_marker)
}
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(file_out_marker,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

#cluster summary
marker_number<-table(markers$cluster)%>% reshape2::melt()
colors=hue_pal()(nrow(marker_number))
marker_number$Cluster <- as.factor(marker_number$Var1)
ggplot(data=marker_number,mapping=aes(x=Cluster,y=value,fill=Cluster))+
  geom_bar(stat="identity",width=0.5)+theme_classic()+ylab("Number")+
  geom_text(aes(label = value),size = 3,hjust = 0.5,vjust = -0.5, position = "stack")


ggsave(paste(file_out_marker,"marker_number.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(file_out_marker,"marker_number.png",sep="/"),width = 14,height = 7)



all_top10_markers=markers %>%  top_n(n = 10, wt = avg_log2FC)
VlnPlot(seurat_obj, features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)

ggsave(paste(file_out_marker,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7)
ggsave(paste(file_out_marker,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7)

DotPlot(seurat_obj, features = all_top10_markers$gene)
ggsave(paste(file_out_marker,"allmarkers_top10_dotplot.pdf",sep="/"),width = 20,height = 7)
ggsave(paste(file_out_marker,"allmarkers_top10_dotplot.png",sep="/"),width = 20,height = 7)

FeaturePlot(seurat_obj, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
ggsave(paste(file_out_marker,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
ggsave(paste(file_out_marker,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)

FeaturePlot(object = seurat_obj,features = rev(unique(all_top10_markers$gene)))

ggsave(paste(file_out_marker,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15)
ggsave(paste(file_out_marker,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15)


cluster_top10_markers=markers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = rev(unique(cluster_top10_markers$gene))) + NoLegend()

ggsave(paste(file_out_marker,"all_cluster_markers_heatmap.pdf",sep="/"),width = 20,height = 14)
ggsave(paste(file_out_marker,"all_cluster_markers_heatmap.png",sep="/"),width = 20,height = 14)

file_out_cluster=paste(file_out,"Each_cluster",sep = "/")
if(dir.exists(file_out_cluster)){
  print("dir exists")
}else{
  dir.create(file_out_cluster)
}

for( clust_num in  unique(Idents(seurat_obj))){
  cluster_dir=paste(file_out_cluster,paste("cluster",clust_num,sep="_"),sep="/")
  if(!file.exists(cluster_dir)){
    dir.create(cluster_dir)
  }
  cluster_dir_enrich=paste(cluster_dir,"enrichment",sep="/")
  if(!file.exists(cluster_dir_enrich)){dir.create(cluster_dir_enrich)}
  cluster_markers=subset(markers,cluster==clust_num)
  if(nrow(cluster_markers)>1){
    genelist=rownames(cluster_markers)
    enrichment(species=opt$type,outDir=cluster_dir_enrich,geneList=genelist)
    write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)
    VlnPlot(seurat_obj, features = top10_markers$gene,pt.size = 0.1 ,ncol=5)
    
    ggsave(paste(cluster_dir,"top10_vilion.pdf",sep="/"),width =20,height = 7)
    ggsave(paste(cluster_dir,"top10_vilion.png",sep="/"),width =20,height = 7)
    
    FeaturePlot(seurat_obj, features = top10_markers$gene, min.cutoff = "q9",ncol=5)
    ggsave(paste(cluster_dir,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
    ggsave(paste(cluster_dir,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)
    DotPlot(seurat_obj, features = rev(unique(top10_markers$gene)))
    
    ggsave(paste(cluster_dir,"top10_exp_pct.pdf",sep="/"),width = 15,height =15)
    ggsave(paste(cluster_dir,"top10_exp_pct.png",sep="/"),width = 15,height =15)
    
  }}

seurat_diff_cluster_dir=paste(file_out,"DiffAnalysis_perCluster",sep = "/")
#find different gene between sample for each cluster test
Idents(seurat_obj)<- seurat_obj$seurat_clusters
if(!file.exists(seurat_diff_cluster_dir)){
  dir.create(seurat_diff_cluster_dir)
}

seurat_obj$cluster_sample <- paste(seurat_obj$group,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
seurat_obj$celltype <- Idents(seurat_obj)
Idents(seurat_obj) <- "cluster_sample"

sample_cluster_avg=AverageExpression(seurat_obj,"RNA")$RNA
write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
cluster_sample=as.vector(unique(Idents(seurat_obj)))

sample_group <- unique(seurat_obj$group)

if(length(unique(sample_group))>1){
  sample_group_length=length(unique(sample_group))-1
  sample_group_length_1=length(unique(sample_group))
  for (idx in 1:sample_group_length) {
    idx_index=idx+1
    for(idx1 in idx_index:sample_group_length_1){
      if(idx != idx1){
        
        compare=paste(unique(sample_group)[idx],unique(sample_group)[idx1],sep="_vs_")
        compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
        if(!file.exists(compare_dir)){
          dir.create(compare_dir)
        }
        
        for(sub_cluster in unique(seurat_obj$celltype)){
          cp1=paste(unique(sample_group)[idx1],paste("cluster",sub_cluster,sep=""),sep="_")
          cp2=paste(unique(sample_group)[idx],paste("cluster",sub_cluster,sep=""),sep="_")
          if(!cp1 %in% cluster_sample || ! cp2 %in% cluster_sample){
            # print("new");
            next;
          }
          ncell1=CellsByIdentities(seurat_obj,cp1)
          ncell2=CellsByIdentities(seurat_obj,cp2)
          if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
            next;
          }
          # print(cp1)
          
          cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
          if(!file.exists(cluster_compare_dir)){
            dir.create(cluster_compare_dir)
          }
          cluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment",sep="/")
          if(!file.exists(cluster_compare_dir_enrich)){dir.create(cluster_compare_dir_enrich)}
          diff.cluster=FindMarkers(seurat_obj, ident.1 = cp1, ident.2 = cp2, verbose = FALSE)
          write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
          sub_compare_cluster1=paste(unique(sample_group)[idx],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
          sub_compare_cluster2=paste(unique(sample_group)[idx1],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
          diff.cluster$gene=row.names(diff.cluster)
          
          top3_diff_gene=diff.cluster$gene[1:3]
          diff.cluster<-subset(diff.cluster,p_val<0.05)
          diff_genes=diff.cluster$gene
          enrichment(species=opt$type,outDir=cluster_compare_dir_enrich,geneList=diff_genes)
          
          plots <- VlnPlot(seurat_obj, features = top3_diff_gene, split.by = "sample", group.by = "celltype",pt.size = 0, combine = FALSE,idents = c(sub_compare_cluster1,sub_compare_cluster2))
          CombinePlots(plots = plots, ncol = 1)
          ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 7)
          ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.png",sep="/"),width = 8,height = 7)
          
          FeaturePlot(seurat_obj, features = top3_diff_gene, split.by = "sample", max.cutoff = 3,  cols = c("grey", "red"))
          ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.pdf",sep="/"),width = 8,height = 7)
          ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.png",sep="/"),width = 8,height = 7)
          
        }
      }
    }
  }
}

saveRDS(object = seurat_obj,file = paste0(file_out,"/sub.rds"))
