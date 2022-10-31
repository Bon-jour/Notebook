#run ST code 
#空间转录组 ST 工作流程
options(bitmapType='cairo')

library("optparse")

option_list <- list(
  make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
  make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
  #make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
 # make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
  make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT"),
  make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="10"),
  make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca  , default %default",default=30),
  make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default="False"),
  #make_option(c("-x","--matrix"),help="file matrix folder, default %default",default="False"),
  make_option(c("-t","--type"),help="organism short name,default %default",default="hsa")
)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))
#suppressPackageStartupMessages(library(SPOTlight))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
# suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(future))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(harmony))
#source("/home/Sc02/source/KEGG/enrichment.R")
source("/Business/psn_company/sc02/1.source/1.profile.R")
source("/Business/psn_company/sc02/1.source/enrichment2.r")
project_dir=getwd()
SCTransform2=function(x,y){
  SCTransform(x, assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208)
}

#load data 

# GSM4565826$percent.mt <- PercentageFeatureSet(GSM4565826, pattern = "^MT")
# 
# save.image("raw.rdata")

# 调色板
blue_green_red <- scale_color_gradient2(low = "blue", mid = "green", high = "red")
byr <- CustomPalette(low = '#0000FF', mid = "#FFFF00",high = "#FF0000", k = 7000)
bgr <- CustomPalette(low = '#0000FF', mid = "#00FF00",high = "#FF0000", k = 7000)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#getPalette(18)
#set future.globals.maxSize 
options(future.globals.maxSize=19288800000)

project_dir=getwd()

#change mt gene pattern
if(opt$type=="mmu"){
  opt$mt="^mt"
  
}

#read info
sample=read.table("sample_info.txt",sep="\t",header=T)
sample_list=c(sample$sample)
sample_group=c(sample$group)
#import 10x count matrix
print("load data...")
ifnb.list=list()
for(x in 1:length(sample_list)){
  datadir=paste0("02_cellranger/",sample_list[x]) #此处需要修改
  single.ob<-Load10X_Spatial(data.dir = datadir)
  # single.ob<-CreateSeuratObject(counts = single.data, project = "single", min.cells = 3, min.features = 200)
  single.ob$sample=sample_list[x]
  single.ob$group=sample_group[x]
  single.ob$imagerow_pos <-as.numeric(single.ob@images$slice1@coordinates$imagerow)
  single.ob$imagecol_pos <- as.numeric(single.ob@images$slice1@coordinates$imagecol)
  ifnb.list= c(ifnb.list,list(single.ob))
}


names(ifnb.list)=as.character(sample_list)

#Do SCTransform2
ifnb.list <- lapply(X =ifnb.list, FUN =SCTransform2)

#intergatedata 

print("IntegrateData...")
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000,verbose = opt$verbose)
combined <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT", anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT") 

# Choose the features to use when integrating multiple datasets. 
# This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. 
# It returns the top scoring features by this ranking.

j=0
for (i in 1:length(combined.sct@images)){
  j=j+1
  if( nrow(combined.sct@images[[j]]@coordinates)==0){
    combined.sct@images[[j]]=NULL
    j=j-1
  }
}

names(combined.sct@images)=as.character(sample_list)

seurat_qc_dir=paste(project_dir,"03_QC_SCT",sep="/")
if(!file.exists(seurat_qc_dir)){
  dir.create(seurat_qc_dir,recursive=TRUE)
}

plot2=plot4=list()

combined.sct <- PercentageFeatureSet(object = combined.sct, pattern = opt$mt,assay = "Spatial",col.name = "percent.mito")

for ( x in 1:length(sample_list) ){
  
  plot2[[x]] <- SpatialFeaturePlot(combined.sct, images=sample_list[x],features = "nCount_Spatial",crop=F) #+ theme(legend.position = "right")
  pic1=wrap_plots(plot2,ncol=2)
  
  plot4[[x]] <- SpatialFeaturePlot(combined.sct, images=sample_list[x],features = "nCount_SCT",crop=F) #+ theme(legend.position = "right")
  pic2=wrap_plots(plot4,ncol=2)
}


plot1 <- VlnPlot(combined.sct, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mito"), ncol = 3,group.by="sample", pt.size = 0.01)
plot3 <- VlnPlot(combined.sct, features = c("nFeature_SCT", "nCount_SCT", "percent.mito"), ncol = 3,group.by="sample", pt.size = 0.01)
## 图1：标准化前的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图
ggsave(plot1, width=12, height=7, filename=paste(seurat_qc_dir,"pre-QC-VlnPlot.pdf",sep="/"),limitsize = FALSE)
ggsave(plot3, width=12, height=7, filename=paste(seurat_qc_dir,"post-QC-VlnPlot.pdf",sep="/"),limitsize = FALSE)
## 图2：标准化后的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图


nrow1=ceiling(length(plot2)/2)
nrow2=ceiling(length(plot4)/2)

ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.pdf",sep="/"),limitsize = FALSE)
ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.png",sep="/"),limitsize = FALSE)
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.pdf",sep="/"),limitsize = FALSE)
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.png",sep="/"),limitsize = FALSE)


fwrite(as.data.frame(combined.sct[["Spatial"]]@counts),paste(seurat_qc_dir,"cells_GeneCounts.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["SCT"]]@data),paste(seurat_qc_dir,"cells_SCT_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["integrated"]]@data),paste(seurat_qc_dir,"cells_integrated_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

rds_dir=paste(project_dir,"seurat_Rdata",sep="/")
if(!file.exists(rds_dir)){
  dir.create(rds_dir,recursive=TRUE)
}

seurat_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
if(!file.exists(seurat_cluster_dir)){
  dir.create(seurat_cluster_dir,recursive=TRUE)
}
combined.sct <- RunPCA(combined.sct,npcs=100, verbose =  opt$verbose)
combined.sct <- RunHarmony(combined.sct, group.by.vars="orig.ident", assay.use="integrated")

## 选择后续分析使用的维度
plot=ElbowPlot(combined.sct,ndim=100,reduction="harmony")
ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir,"ElbowPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir,"ElbowPlot.png",sep="/"),dpi=300,limitsize = FALSE)
## 构建KNN图 
combined.sct <- FindNeighbors(combined.sct,reduction = "harmony", dims = 1:ndim, verbose =  opt$verbose)

save(ifnb.list,file=paste(rds_dir,"separated.sct.Rdata",sep="/"))

saveRDS(combined.sct,paste(rds_dir,"combined.sct.rds",sep="/"))

# combined.sct <- readRDS("suerat_Rdata/combined.sct.rds")

res=0.8
seurat_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
if(!file.exists(seurat_cluster_dir)){
  dir.create(seurat_cluster_dir,recursive=TRUE)
}

seurat_cluster_dir_1=paste(seurat_cluster_dir,"Integrated",sep="/")
if(!file.exists(seurat_cluster_dir_1)){
  dir.create(seurat_cluster_dir_1,recursive=TRUE)
}

## 分群
combined.sct <- FindClusters(combined.sct, resolution=res)
cluster_summary=reshape2::dcast(as.data.frame(table(data.frame("cluster"=Idents(combined.sct),"sample"=combined.sct$sample))),sample~cluster)
fwrite(cluster_summary,paste(seurat_cluster_dir_1,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

## 分群结果展示UMAP VS TSNE
plot=ElbowPlot(combined.sct,ndim=100,reduction="harmony")
ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
combined.sct <- RunUMAP(combined.sct, reduction = "harmony", dims = 1:ndim)
combined.sct <- RunTSNE(combined.sct, reduction = "harmony", dims = 1:ndim)
plot1 <- DimPlot(combined.sct, reduction = "umap",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
plot2 <- DimPlot(combined.sct, reduction = "tsne",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
plot3 <- DimPlot(combined.sct, reduction = "umap",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
plot4 <- DimPlot(combined.sct, reduction = "tsne",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
pic3=wrap_plots(plot1, plot2,plot3,plot4,ncol=2)
ggsave(pic3, width=16, height=15,  filename=paste(seurat_cluster_dir_1,"UMAP_VS_tSNE.pdf",sep="/"),dpi=1200,limitsize = FALSE)
ggsave(pic3, width=16, height=15,  filename=paste(seurat_cluster_dir_1,"UMAP_VS_tSNE.png",sep="/"), dpi=300,limitsize = FALSE)

plot5 <- DimPlot(combined.sct, reduction = "umap",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
plot6 <- DimPlot(combined.sct, reduction = "tsne",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
nrow=ceiling(length(sample_list)/2)
ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"UMAP_eachsample.pdf",sep="/"),dpi=1200,limitsize = FALSE)
ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"UMAP_eachsample.png",sep="/"),dpi=300,limitsize = FALSE)
ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"tSNE_eachsample.pdf",sep="/"), dpi=1200,limitsize = FALSE)
ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"tSNE_eachsample.png",sep="/"), dpi=300,limitsize = FALSE)

## 空间染色片上的分群结果
plot7=list()
width_add=0
for ( x in 1:length(sample_list) ){
  plot7[[x]] <- SpatialDimPlot(combined.sct,images=sample_list[x], label = TRUE, label.size = 3,repel=TRUE,crop=F)+ 
    
    labs(title=as.character(sample_list[x]))+
    theme(plot.title=element_text(color = "black", size = 20,hjust=0.5,vjust=0.5),
          legend.text=element_text(face="plain",size=12),
          legend.title = element_text(size=12,face = "bold",hjust=0.5,vjust=0.5),
          legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.8,'cm') )
  
  plot7[[x]] <- plot7[[x]] + guides(fill= guide_legend(nrow =min(15,length(unique( plot7[[x]]$data$ident))), title = "Clusters"))
  names(plot7)[x]=sample_list[x]
  width_add=max(width_add,ceiling(length(unique( plot7[[x]]$data$ident))/15)/2)
}

pic4=wrap_plots(plot7,ncol=2)
nrow=ceiling(length(sample_list)/2)
ggsave(pic4, width=13.5+width_add, height=6*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialDimPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
ggsave(pic4, width=13.5+width_add, height=6*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialDimPlot.png",sep="/"),  dpi=300,limitsize = FALSE)


##  查找cluster biomarkers特征基因，获取其表达数据
DefaultAssay(combined.sct) <- "SCT"
combined.sct <- ScaleData(combined.sct)  ##scale.data将会被替换
cluster_markers_all <- FindAllMarkers(object = combined.sct, 
                                      assay = "SCT",
                                      slot="data",
                                      only.pos = TRUE, 
                                      logfc.threshold = 0.25,
                                      min.pct = 0.1,
                                      test.use = "wilcox")   ###还有别的方法
cluster_markers_all <- cluster_markers_all %>%  filter(p_val_adj <= 0.05 )
fwrite(cluster_markers_all,paste(seurat_cluster_dir,"all_markers.xls",sep="/"),row.names=TRUE,col.names=TRUE,sep="\t")


all_top10_markers=cluster_markers_all %>% top_n(n = 10, wt = avg_log2FC)

plot8=VlnPlot(combined.sct, features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)
ggsave(plot8,filename=paste(seurat_cluster_dir,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
ggsave(plot8,filename=paste(seurat_cluster_dir,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

plot9=FeaturePlot(combined.sct, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5,reduction = "tsne")
ggsave(plot9,filename=paste(seurat_cluster_dir,"allmarkers_top10_cell_exp_distribution.tsne.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
ggsave(plot9,filename=paste(seurat_cluster_dir,"allmarkers_top10_cell_exp_distribution.tsne.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

plot9_umap=FeaturePlot(combined.sct, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5,reduction = "umap")
ggsave(plot9_umap,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_cell_exp_distribution.umap.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
ggsave(plot9_umap,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_cell_exp_distribution.umap.png",sep="/"),width = 20,height = 7,limitsize = FALSE)


plot10=DotPlot(combined.sct, features = rev(unique(all_top10_markers$gene)), cols = rainbow(length(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()
ggsave(plot10,filename=paste(seurat_cluster_dir,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15,limitsize = FALSE)
ggsave(plot10,filename=paste(seurat_cluster_dir,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15,limitsize = FALSE)

if(!file.exists(paste(seurat_cluster_dir,"SpatialFeaturePlot_all_top10",sep="/"))){
  dir.create(paste(seurat_cluster_dir,"SpatialFeaturePlot_all_top10",sep="/"),recursive=TRUE)
}

for ( y in 1:length(all_top10_markers$gene)){
  plot11 <- SpatialFeaturePlot(combined.sct,slot="data",features=all_top10_markers$gene[y],combine=F,alpha=c(0.1,1),stroke = 0.5,crop=F)
  pic5=wrap_plots(plot11,ncol=2,widths=3,heights=3)
  nrow=ceiling(length(sample_list)/2)
  ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[y],"pdf",sep="."),sep="/"),  dpi=1200,limitsize = FALSE)
  ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[y],"png",sep="."),sep="/"),  dpi=300,limitsize = FALSE)
}

DefaultAssay(combined.sct) <- "SCT"
cluster_top10_markers=cluster_markers_all %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)

plot11=DoHeatmap(combined.sct, slot="scale.data", features = rev(unique(cluster_top10_markers$gene)))

ggsave(plot11,filename=paste(seurat_cluster_dir_1,"top10_marker_echo_cluster_heatmap.pdf",sep="/"),width = 20,height = 14,limitsize = FALSE)
ggsave(plot11,filename=paste(seurat_cluster_dir_1,"top10_marker_echo_cluster_heatmap.png",sep="/"),width = 20,height = 14,limitsize = FALSE)

for( clust_num in unique(combined.sct$seurat_clusters)){
  seurat_cluster_dir_tmp=paste(seurat_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
  if(!file.exists(seurat_cluster_dir_tmp)){
    dir.create(seurat_cluster_dir_tmp,recursive=TRUE)
  }
  cluster_markers=subset(cluster_markers_all,cluster==clust_num)
  C_num=nrow(cluster_markers)
  
  seurat_cluster_dir_tmp_kegg=paste(seurat_cluster_dir_tmp,"enrichment",sep="/")
  if(!file.exists(seurat_cluster_dir_tmp_kegg)){
    dir.create(seurat_cluster_dir_tmp_kegg,recursive=TRUE)
  }
  genelist <- cluster_markers$gene
  print(length(genelist))
  enrichment(species = opt$type,geneList = genelist,outDir =seurat_cluster_dir_tmp_kegg )
  
  
  if(C_num==0){
    system(paste0("echo 'NO marker was found for this cluster!' > ",seurat_cluster_dir_tmp,paste("/cluster",clust_num,"markers.xls",sep="_")))
    next()
  }
  fwrite(cluster_markers,paste(seurat_cluster_dir_tmp,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)
  
  top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)
  
  
  plot12=VlnPlot(combined.sct, features = top10_markers$gene,pt.size = 0.1 ,ncol=5)
  
  ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.pdf",sep="/"),width =20,height = 7,limitsize = FALSE)
  ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.png",sep="/"),width =20,height = 7,limitsize = FALSE)
  
  plot13=FeaturePlot(combined.sct, features = top10_markers$gene, min.cutoff = "q9",ncol=5)
  ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
  ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7,limitsize = FALSE)
  
  plot14=DotPlot(combined.sct, features = rev(unique(top10_markers$gene)), cols = rainbow(length(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()
  
  ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.pdf",sep="/"),width = 15,height =15,limitsize = FALSE)
  ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.png",sep="/"),width = 15,height =15,limitsize = FALSE)
  
  for ( y in 1:length(top10_markers$gene)){
    plot15 <- SpatialFeaturePlot(combined.sct,slot="data",features=top10_markers$gene[y],alpha=c(0.1,1),combine=FALSE,stroke = 0.1,crop=F)
    
    pic6=wrap_plots(plot15,ncol=2,widths=3,heights=3)
    nrow=ceiling(length(sample_list)/2)
    ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[y],"pdf",sep="."),sep="/"),  dpi=1200,limitsize = FALSE)
    ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[y],"png",sep="."),sep="/"),  dpi=300,limitsize = FALSE)
  }
}

#group compare
Idents(combined.sct)<- combined.sct$seurat_clusters
seurat_diff_cluster_dir=paste(project_dir,"05_DiffAnalysis_perCluster",sep="/")
if(!file.exists(seurat_diff_cluster_dir)){
  dir.create(seurat_diff_cluster_dir)
}

combined.sct$cluster_sample <- paste(combined.sct$group,paste("cluster",Idents(combined.sct),sep="") , sep = "_")
combined.sct$celltype <- Idents(combined.sct)
Idents(combined.sct) <- "cluster_sample"

sample_cluster_avg=AverageExpression(combined.sct,"integrated")$RNA
write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
cluster_sample=as.vector(unique(Idents(combined.sct)))



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
        
        for(sub_cluster in unique(combined.sct$celltype)){
          cp1=paste(unique(sample_group)[idx1],paste("cluster",sub_cluster,sep=""),sep="_")
          cp2=paste(unique(sample_group)[idx],paste("cluster",sub_cluster,sep=""),sep="_")
          if(!cp1 %in% cluster_sample || ! cp2 %in% cluster_sample){
            # print("new");
            next;
          }
          ncell1=CellsByIdentities(combined.sct,cp1)
          ncell2=CellsByIdentities(combined.sct,cp2)
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
          diff.cluster=FindMarkers(combined.sct, ident.1 = cp1, ident.2 = cp2, verbose = FALSE)
          write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
          sub_compare_cluster1=paste(unique(sample_group)[idx],paste("cluster",unique(combined.sct$celltype),sep=""),sep="_")
          sub_compare_cluster2=paste(unique(sample_group)[idx1],paste("cluster",unique(combined.sct$celltype),sep=""),sep="_")
          diff.cluster$gene=row.names(diff.cluster)
          
          top3_diff_gene=diff.cluster$gene[1:3]
          diff_genes=diff.cluster$gene
          enrichment(species=opt$type,outDir=cluster_compare_dir_enrich,geneList=diff_genes)
          
          plots <- VlnPlot(combined.sct, features = top3_diff_gene, split.by = "sample", group.by = "celltype",pt.size = 0, combine = FALSE,idents = c(sub_compare_cluster1,sub_compare_cluster2))
          CombinePlots(plots = plots, ncol = 1)
          ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 7)
          ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.png",sep="/"),width = 8,height = 7)
          
          FeaturePlot(combined.sct, features = top3_diff_gene, split.by = "sample", max.cutoff = 3,  cols = c("grey", "red"))
          ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.pdf",sep="/"),width = 8,height = 7)
          ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.png",sep="/"),width = 8,height = 7)
          
          
          
        }
      }
    }
  }
}

