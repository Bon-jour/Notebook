options(bitmapType='cairo')
library(Seurat)
library(optparse)
library(scales)
library(dplyr)
library(ggplot2)
source("/Business/psn_company/sc02/1.source/enrichment2.r")


option_list <- list(
  make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
  make_option(c("-s", "--s1"), help="s1,ident1 to compare with other" ),
  make_option(c("-c", "--s2"), help="s2,other,default ident1 compare ident2 " ),
  make_option(c("-t", "--feature"), help=" need to use",default="supp"),
  make_option(c("-v", "--verbose"), help="Shown more processing information",default="False"),
  make_option(c("-r", "--rds"), help="seurat object",default="All_sample_combined.rds"),
  make_option(c("-p", "--species"), help="oragnism",default="hsa")
  )
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

immune.combined <- readRDS(opt$rds)

if(!file.exists(opt$workdir)){
  dir.create(opt$workdir)
}
sample_name <- paste0(opt$s2,"VS",opt$s1)
results_out <- paste0(opt$workdir,"/",sample_name)

if(!file.exists(results_out)){
  dir.create(results_out)
}
DefaultAssay(immune.combined)="RNA"
markers <- Seurat::FindMarkers(object =immune.combined,ident.1 = opt$s1,ident.2=opt$s2,assay="RNA",group.by =opt$feature)
markers=subset(markers,p_val<0.05)
markers$gene <- rownames(markers)
write.table(markers,paste(results_out,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

all_top10_markers= markers %>%  top_n(n = 10, wt = avg_log2FC)
if(length(all_top10_markers)>0){
VlnPlot(immune.combined, features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)

ggsave(paste(results_out,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7)
ggsave(paste(results_out,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7)
}
DotPlot(immune.combined, features = all_top10_markers$gene)
ggsave(paste(results_out,"allmarkers_top10_dotplot.pdf",sep="/"),width = 20,height = 7)
ggsave(paste(results_out,"allmarkers_top10_dotplot.png",sep="/"),width = 20,height = 7)

FeaturePlot(immune.combined, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
ggsave(paste(results_out,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
ggsave(paste(results_out,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)


DotPlot(immune.combined, features = rev(unique(all_top10_markers$gene)))
ggsave(paste(results_out,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15)
ggsave(paste(results_out,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15)


Idents(immune.combined)<- immune.combined$seurat_clusters
seurat_diff_cluster_dir=paste(results_out,"DiffAnalysis_perCluster",sep="/")
if(!file.exists(seurat_diff_cluster_dir)){
  dir.create(seurat_diff_cluster_dir)
}

immune.combined$cluster_sample <- paste(immune.combined@meta.data[opt$feature][,1],paste("cluster",Idents(immune.combined),sep="") , sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- immune.combined$cluster_sample
cluster_sample=as.vector(unique(Idents(immune.combined)))


for(sub_cluster in unique(immune.combined$seurat_clusters)){
  cp1 <- paste0(opt$s1,"_cluster",sub_cluster)
  cp2 <- paste0(opt$s2,"_cluster",sub_cluster)
  
  if(!cp1 %in% cluster_sample || ! cp2 %in% cluster_sample){
    # print("new");
    next;
  }
  ncell1=CellsByIdentities(immune.combined,cp1)
  ncell2=CellsByIdentities(immune.combined,cp2)
  if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
    next;
  }
  cluster_compare_dir=paste(seurat_diff_cluster_dir,paste("cluster",sub_cluster,sep=""),sep="/")
  if(!file.exists(cluster_compare_dir)){
    dir.create(cluster_compare_dir)
  }
  cluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment",sep="/")
  if(!file.exists(cluster_compare_dir_enrich)){dir.create(cluster_compare_dir_enrich)}
  diff.cluster=FindMarkers(immune.combined, ident.1 = cp1, ident.2 = cp2, verbose = FALSE)
  diff.cluster_filter=subset(diff.cluster,p_val<0.05)
  write.table(diff.cluster_filter,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
  if(nrow(diff.cluster)<1){next;}
  # sub_compare_cluster1=paste(unique(sample_group)[idx],paste("cluster",unique(immune.combined$celltype),sep=""),sep="_")
  # sub_compare_cluster2=paste(unique(sample_group)[idx1],paste("cluster",unique(immune.combined$celltype),sep=""),sep="_")
  diff.cluster_filter$gene=row.names(diff.cluster_filter)
  #diff.cluster_filter=subset(diff.cluster,p_val<0.05)
  top3_diff_gene=diff.cluster_filter$gene[1:3]
  diff_genes=diff.cluster_filter$gene
  enrichment(species=opt$species,outDir=cluster_compare_dir_enrich,geneList=diff_genes)
  plots <- VlnPlot(immune.combined, features = top3_diff_gene, split.by = "sample", group.by = "celltype",pt.size = 0, combine = FALSE,idents = c(cp1,cp2))
  CombinePlots(plots = plots, ncol = 1)
  ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 7)
  ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.png",sep="/"),width = 8,height = 7)
  
  FeaturePlot(immune.combined, features = top3_diff_gene, split.by = "sample", max.cutoff = 3,  cols = c("grey", "red"))
  ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.pdf",sep="/"),width = 8,height = 7)
  ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.png",sep="/"),width = 8,height = 7)
  
  
}
