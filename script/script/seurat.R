options(bitmapType='cairo')

library("optparse")

option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
        make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
        make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
        make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT"),
        make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="10"),
        make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca, default %default",default=30),
        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default="False"),
        make_option(c("-t","--type"),help="organism short name,default %default",default="hsa"),
        make_option(c("-g","--gather"),help="remove batch effect methods",default="cca")
)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(reshape2)
library(DoubletFinder)
library(harmony)
library(dittoSeq)
library(scales)
source("/Business/psn_company/sc01/Workbase_Wangcanbiao/script/1.profile.R")
source("/Business/psn_company/sc01/Workbase_Wangcanbiao/script/enrichment2.r")
# if(! is.null(opt$workdir)){
#   print(paste("change R work directory to ",opt$workdir,sep=""))
#   setwd(opt$workdir)
#  }
project_dir=getwd()
if(!file.exists(project_dir)){
   dir.create(project_dir)
}

#change mt gene pattern
if(opt$type=="mmu"){
        opt$mt="^mt-"

}


# #read sample information matrix
# #if more than one sample, should give the list csv used in cellranger aggr
# cellranger_csv=paste(project_dir,"02_Cellranger","all.csv",sep="/")
# if(! is.null(opt$aggr_csv)){
#   print(paste("set aggr_csv to ",opt$aggr_csv,sep=""))
#   cellranger_csv=opt$aggr_csv
#  }

# sample_list=read.csv(cellranger_csv,stringsAsFactors =F)

# if(nrow(sample_list) < 2){
#   stop("this pipline need at least 2 sample\n")
# }

# #import DE contrast information
# if(! is.null(opt$contrast)){
#   print(paste("get compare contrast file ",opt$contrast,sep=""))
#   sample_contrast=read.table(opt$contrast)
#  }else{
# #sample_contrast=switch(EXPR=as.character(nrow(sample_list)),"2"=combn(sample_list[,1],2), t(combn(sample_list[,1],2)))
#  sample_contrast=t(combn(sample_list[,1],2))
#  }
sample=read.table("sample_info.txt",sep="\t",header=T)
sample_list=c(sample$sample)
sample_group=c(sample$group)
#import 10x count matrix
ifnb.list=list()
for(x in 1:length(sample_list)){
  # datadir=gsub("molecule_info.h5","filtered_feature_bc_matrix",sample_list[x,2])
   #sample=sample_list[x,1]
   #group=sample_list[x,4]
   datadir=paste0("02_cellranger/",sample_list[x],"/filtered_feature_bc_matrix")
   print(datadir)
   single.data<-Read10X(data.dir = datadir)
   single.ob<-CreateSeuratObject(counts = single.data, project = "single", min.cells = 3, min.features = 200)
   single.ob$sample=sample_list[x]
   single.ob$group=sample_group[x]
   ifnb.list= c(ifnb.list,list(single.ob))
   # saveRDS(single.ob,paste(project_dir,"/","02_Cellranger","/",sample,".rds",sep=""))
}

if(length(sample_list)>1){
single.ob=merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_list)}else{
single.ob=ifnb.list[[1]]}


cell_num_raw=table(single.ob$sample)
##QC ,filter low Q cell
seurat_qc_dir=paste(project_dir,"03_QC_Filter",sep="/")

if(!file.exists(seurat_qc_dir)){
   dir.create(seurat_qc_dir)
}
#save basic merge rds
saveRDS(single.ob,paste0(seurat_qc_dir,"/Basicinfo_inter.rds"))

#get mt percentage
single.ob[["percent.mt"]] <- PercentageFeatureSet(single.ob, pattern = opt$mt)


##output experiment raw information
out_gene_umi=as.matrix(single.ob@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file=paste(seurat_qc_dir,"Basicinfo_nGene_nUMI_mito.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=F)

# #calculate cell cycle score
# single.ob <- NormalizeData(single.ob)
# single.ob <- ScaleData(single.ob)
# single.ob <- FindVariableFeatures(object = single.ob,selection.method = 'vst', nfeatures = 2000)
# single.ob <- RunPCA(single.ob,  features = VariableFeatures(object = single.ob) ,verbose = FALSE)
# single.ob <- RunUMAP(single.ob, reduction = "pca", dims = 1:20)
# single.ob <- FindVariableFeatures(single.ob, selection.method = "vst", nfeatures = 2000)
# single.ob <- FindNeighbors(single.ob, reduction = "pca", dims = 1:20)
# single.ob <- FindClusters(single.ob)




#filter
cellqc.raw=VlnPlot(single.ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by="sample")
single.ob <- subset(single.ob, subset = nFeature_RNA > 400 & nFeature_RNA < 5000 & percent.mt <opt$mt_cutoff)
cellqc.filter=VlnPlot(single.ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by="sample")


#plot_grid(cellqc.raw,cellqc.filter,nrow=2)
cellqc.raw/cellqc.filter
ggsave(paste(seurat_qc_dir,"cells_qc_filter.pdf",sep="/"),width = 8,height = 7)
ggsave(paste(seurat_qc_dir,"cells_qc_filter.png",sep="/"),width = 8,height = 7)

cell_num_filter=table(single.ob$sample)

filter_summary=cbind(cell_num_raw,cell_num_filter)
filter_summary=data.frame("cell_num_raw"=filter_summary[,1],"cell_num_filter"=filter_summary[,2])
filter_summary[nrow(filter_summary)+1,]=apply(filter_summary,2,sum)
rownames(filter_summary)[nrow(filter_summary)]="total"
filter_summary$percent=apply(filter_summary,1,function(x){round(x[2]*100/x[1],2)})

write.table(filter_summary,paste(seurat_qc_dir,"cells_filter_stat.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

#for more than one sample
#split data by sample
#names(ifnb.list)=sample_list
seurat_exp_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
if(!file.exists(seurat_exp_cluster_dir)){dir.create(seurat_exp_cluster_dir)}
if(!file.exists(paste(seurat_exp_cluster_dir,"All_sample_combined.rds",sep="/"))){
if(length(sample_list)>1){
ifnb.list <- SplitObject(single.ob, split.by = "sample")


ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- ScaleData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
        x <- RunPCA(x)
        x <- RunUMAP(x, dims = 1:10)
        x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
        x <- FindClusters(x)
        ##doublet finder need test
        annotations <- x@meta.data$seurat_clusters
        sweep.res.list <- paramSweep_v3(x, PCs = 1:10)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
        homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
        double_rate <- 0.0005272+0.000007589*length(x$sample)
        nExp_poi <- round(double_rate*length(x$sample))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        x <- my_doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
        x <- my_doubletFinder_v3_1(x, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = "DF.classifications_raw", sct = FALSE)
        x@meta.data[,"Double_status"] <- x@meta.data$DF.classifications_raw
#       meta_data<-x@meta.data
        #meta_data$DF_hi.lo[which(meta_data$DF_hi.lo == "Singlet" & meta_data$DF.classifications_raw== "Doublet")] <- "Doublet_lo"
#       meta_data$DF_hi.lo[which(meta_data$Double_status == "Doublet")] <- "Doublet"
        #x@meta.data$DF_hi.lo[which(x@meta.data$DF_hi.lo == "Doublet" & x@meta.data$DF.classifications_raw == "Singlet")] <- "Doublet_lo"
        #x@meta.data$DF_hi.lo[which(x@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
        #x@meta.data$DF_hi.lo<-meta_data$DF_hi.lo
        sample_name=unique(x$sample) #get sample name
        out_gene_umi=as.matrix(x@meta.data)
        row_name=c("Cell",rownames(out_gene_umi))
        out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
        out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
        write.table(out_gene_umi,file=paste0(seurat_qc_dir,"/",sample_name,".cells_doublet_info.tsv"),sep="\t")
        p1<-DimPlot(x, group.by="seurat_clusters")
        p2<-DimPlot(x, group.by="Double_status",order=c("Singlet","Doublet"),cols =c("gray","red"))
        p1+p2
        ggsave(paste0(seurat_qc_dir,"/",sample_name,".cells_qc_filter.pdf"),width =16,height = 7)
        #x<- subset(x,DF_hi.lo=="Singlet")
        # return(x)
        x<-subset(x,Double_status=="Singlet")
        x<-x


})


#integrat sample object
if(opt$gather=="harmony"){
immune.combined<-merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_list)%>%NormalizeData() %>% FindVariableFeatures()
#immune.combined<-subset(immune.combined,Double_status=="Singlet") 
immune.combined<-ScaleData(immune.combined,feature=rownames(immune.combined)) %>% RunPCA(verbose = FALSE) %>% RunHarmony( group.by.vars = "sample")
immune.combined<-RunUMAP(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined<-RunTSNE(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:20) %>% FindClusters()
}else{
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"
#immune.combined<-subset(immune.combined,Double_status=="Singlet")
immune.combined <- ScaleData(immune.combined,feature=rownames(immune.combined), verbose = FALSE)
immune.combined <- RunPCA(immune.combined,  features = VariableFeatures(object = immune.combined) ,verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined)
}}else{
immune.combined<-single.ob
#immune.combined<-subset(immune.combined,Double_status=="Singlet")
immune.combined <-NormalizeData(immune.combined)
immune.combined <- ScaleData(immune.combined,feature=rownames(immune.combined), verbose = FALSE)
immune.combined <- FindVariableFeatures(object = immune.combined,selection.method = 'vst', nfeatures = 2000)
immune.combined <- RunPCA(immune.combined,  features = VariableFeatures(object = immune.combined) ,verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined)


}



seurat_exp_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
if(!file.exists(seurat_exp_cluster_dir)){
   dir.create(seurat_exp_cluster_dir)
}

#write.table(immune.combined[["RNA"]]@data,paste(seurat_exp_cluster_dir,"cells_logNormalize_expression.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

write.table(immune.combined[["RNA"]]@counts,paste(seurat_exp_cluster_dir,"cells_GeneCounts.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

#run scale pca cluster and umap

cluster_summary=dcast(as.data.frame(table(data.frame("cluster"=Idents(immune.combined),"sample"=immune.combined$sample))),sample~cluster)
write.table(cluster_summary,paste(seurat_exp_cluster_dir,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

#cell cycling
immune.combined <- CellCycleScoring(immune.combined, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
DimPlot(immune.combined,group.by="Phase")
ggsave(paste(seurat_qc_dir,"cells_cell_cycle.pdf",sep="/"),width = 8,height = 7)
ggsave(paste(seurat_qc_dir,"cells_cell_cycle.png",sep="/"),width = 8,height = 7)
#cell cycle bar plot
p<-dittoBarPlot(immune.combined, "Phase", group.by = "sample")
ggsave(p,filename=paste(seurat_qc_dir,"cell_cycle_per_sample.pdf",sep="/"))

saveRDS(object=immune.combined,paste(seurat_exp_cluster_dir,"All_sample_combined.rds",sep="/"))}else{
immune.combined<-readRDS(paste(seurat_exp_cluster_dir,"All_sample_combined.rds",sep="/"))
}


#ElbowPlot
ElbowPlot(immune.combined)
ggsave(paste(seurat_exp_cluster_dir,"pca.ElbowPlot.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(seurat_exp_cluster_dir,"pca.ElbowPlot.png",sep="/"),width = 14,height = 7)

#cluster summary
cluster_number<-table(immune.combined$seurat_clusters)%>% reshape2::melt()
colors=hue_pal()(nrow(cluster_number))
cluster_number$Cluster <- as.factor(cluster_number$Var1)
ggplot(data=cluster_number,mapping=aes(x=Cluster,y=value,fill=Cluster))+
  geom_bar(stat="identity",width=0.5)+theme_classic()+ylab("Number")+
  geom_text(aes(label = value),size = 3,hjust = 0.5,vjust = -0.5, position = "stack")


ggsave(paste(seurat_exp_cluster_dir,"cluster_number.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(seurat_exp_cluster_dir,"cluster_number.png",sep="/"),width = 14,height = 7)


p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.png",sep="/"),width = 14,height = 7)

out_width=ceiling(sqrt(length(sample_list)))+1
height_for_save <- length(sample_list)/out_width*7
length_for_save <-ceiling(out_width/2)*7

DimPlot(immune.combined, reduction = "umap", split.by = "sample",ncol=out_width)
ggsave(paste(seurat_exp_cluster_dir,"cluster_umap_eachsample.pdf",sep="/"),width = length_for_save,height = height_for_save)
ggsave(paste(seurat_exp_cluster_dir,"cluster_umap_eachsample.png",sep="/"),width = length_for_save,height = height_for_save)

p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)
ggsave(paste(seurat_exp_cluster_dir,"cluster_tsne.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(seurat_exp_cluster_dir,"cluster_tsne.png",sep="/"),width = 14,height = 7)

DimPlot(immune.combined, reduction = "tsne", split.by = "sample",ncol=out_width)
ggsave(paste(seurat_exp_cluster_dir,"cluster_tsne_eachsample.pdf",sep="/"),width = length_for_save,height = height_for_save)
ggsave(paste(seurat_exp_cluster_dir,"cluster_tsne_eachsample.png",sep="/"),width = length_for_save,height = height_for_save)



avg_per_cluster=AverageExpression(immune.combined,"RNA")$RNA
colnames(avg_per_cluster)=gsub("RNA.","",colnames(avg_per_cluster),perl=T)
write.table(avg_per_cluster,paste(seurat_exp_cluster_dir,"avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)


#find marker
DefaultAssay(immune.combined) <- "RNA"
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
DoHeatmap(immune.combined, features = rev(unique(cluster_top10_markers$gene))) + NoLegend()

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
        enrichment(species=opt$type,outDir=cluster_dir_enrich,geneList=genelist)
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

}}




# #find different gene between sample for each cluster
# seurat_diff_cluster_dir=paste(project_dir,"05_DiffAnalysis_perCluster",sep="/")
# if(!file.exists(seurat_diff_cluster_dir)){
#    dir.create(seurat_diff_cluster_dir)
# }


# immune.combined$cluster_sample <- paste(immune.combined$sample,paste("cluster",Idents(immune.combined),sep="") , sep = "_")
# immune.combined$celltype <- Idents(immune.combined)
# Idents(immune.combined) <- "cluster_sample"

# sample_cluster_avg=AverageExpression(immune.combined,"RNA")$RNA
# write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)


# cluster_sample=as.vector(unique(Idents(immune.combined)))
# for(x in 1:nrow(sample_contrast)){

#       compare=paste(sample_contrast[x,1],sample_contrast[x,2],sep="_vs_")
#       compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
#       if(!file.exists(compare_dir)){
#          dir.create(compare_dir)
#       }
#       for(sub_cluster in unique(immune.combined$celltype)){
#               cp1=paste(sample_contrast[x,2],paste("cluster",sub_cluster,sep=""),sep="_")
#               cp2=paste(sample_contrast[x,1],paste("cluster",sub_cluster,sep=""),sep="_")
#               if(!cp1 %in% cluster_sample || ! cp2 %in% cluster_sample){
#                       next;
#               }
#               ncell1=CellsByIdentities(immune.combined,cp1)
#               ncell2=CellsByIdentities(immune.combined,cp2)
#               if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
#                       next;
#               }
#               cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
#               if(!file.exists(cluster_compare_dir)){
#                  dir.create(cluster_compare_dir)
#               }


#               diff.cluster=FindMarkers(immune.combined, ident.1 = cp1, ident.2 = cp2, verbose = FALSE)
#               write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
#               sub_compare_cluster1=paste(sample_contrast[x,1],paste("cluster",unique(immune.combined$celltype),sep=""),sep="_")
#               sub_compare_cluster2=paste(sample_contrast[x,2],paste("cluster",unique(immune.combined$celltype),sep=""),sep="_")
#               diff.cluster$gene=row.names(diff.cluster)

#               top3_diff_gene=diff.cluster$gene[1:3]
#               plots <- VlnPlot(immune.combined, features = top3_diff_gene, split.by = "sample", group.by = "celltype",pt.size = 0, combine = FALSE,idents = c(sub_compare_cluster1,sub_compare_cluster2))
#               CombinePlots(plots = plots, ncol = 1)
#               ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 7)
#               ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.png",sep="/"),width = 8,height = 7)

#               FeaturePlot(immune.combined, features = top3_diff_gene, split.by = "sample", max.cutoff = 3,  cols = c("grey", "red"))
#               ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.pdf",sep="/"),width = 8,height = 7)
#               ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.png",sep="/"),width = 8,height = 7)

#       }


# }


#find different gene between sample for each cluster test
Idents(immune.combined)<- immune.combined$seurat_clusters
seurat_diff_cluster_dir=paste(project_dir,"05_DiffAnalysis_perCluster",sep="/")
if(!file.exists(seurat_diff_cluster_dir)){
   dir.create(seurat_diff_cluster_dir)
}

immune.combined$cluster_sample <- paste(immune.combined$group,paste("cluster",Idents(immune.combined),sep="") , sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "cluster_sample"

sample_cluster_avg=AverageExpression(immune.combined,"RNA")$RNA
write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
cluster_sample=as.vector(unique(Idents(immune.combined)))



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

                        for(sub_cluster in unique(immune.combined$celltype)){
                        cp1=paste(unique(sample_group)[idx1],paste("cluster",sub_cluster,sep=""),sep="_")
                        cp2=paste(unique(sample_group)[idx],paste("cluster",sub_cluster,sep=""),sep="_")
                        if(!cp1 %in% cluster_sample || ! cp2 %in% cluster_sample){
                                # print("new");
                                next;
                        }
                        ncell1=CellsByIdentities(immune.combined,cp1)
                        ncell2=CellsByIdentities(immune.combined,cp2)
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
                        diff.cluster=FindMarkers(immune.combined, ident.1 = cp1, ident.2 = cp2, verbose = FALSE)
                        write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
                        sub_compare_cluster1=paste(unique(sample_group)[idx],paste("cluster",unique(immune.combined$celltype),sep=""),sep="_")
                        sub_compare_cluster2=paste(unique(sample_group)[idx1],paste("cluster",unique(immune.combined$celltype),sep=""),sep="_")
                        diff.cluster$gene=row.names(diff.cluster)
			diff.cluster=subset(diff.cluster,p_val<0.05)
                        top3_diff_gene=diff.cluster$gene[1:3]
                        diff_genes=diff.cluster$gene
                        enrichment(species=opt$type,outDir=cluster_compare_dir_enrich,geneList=diff_genes)

                        plots <- VlnPlot(immune.combined, features = top3_diff_gene, split.by = "sample", group.by = "celltype",pt.size = 0, combine = FALSE,idents = c(sub_compare_cluster1,sub_compare_cluster2))
                        CombinePlots(plots = plots, ncol = 1)
                        ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 7)
                        ggsave(paste(cluster_compare_dir,"top3_diffgene_exp_vilion.png",sep="/"),width = 8,height = 7)

                        FeaturePlot(immune.combined, features = top3_diff_gene, split.by = "sample", max.cutoff = 3,  cols = c("grey", "red"))
                        ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.pdf",sep="/"),width = 8,height = 7)
                        ggsave(paste(cluster_compare_dir,"top3_diffgene_umap.png",sep="/"),width = 8,height = 7)



                        }
                        }
                }
        }
}





# getFeaturePct<-function(object,features,celltype,assay="RNA",slot="data"){
#       data <- GetAssayData(object = new.ob[[assay]], slot = slot)
#       cells=WhichCells(object = object, idents = celltype)
#       round(x = Matrix::rowSums(x = data[features,cells, drop = FALSE] > 0) /length(cells),digits = 3)
# }


#library(SingleR,lib.loc="/Business/psn_company/t01/public/scRNA/R_library/")
# load("/Business/psn_company/t01/public/scRNA/script/HumanPrimaryCellAtlasData.rds")
# sce_for_SingleR <- GetAssayData(immune.combined, slot="counts")
# pred.mouseImmu <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.main,method="cluster",clusters=immune.combined@meta.data$seurat_clusters)
# write.table(pred.mouseImmu,"singleR.result.xls",sep="\t",quote=F,row.names=T,col.names=NA)




# makeScRNASampleCountMatrix<-function(object){
#    samples=unique(object$sample)
#    countMatrix=data.frame()
#    for(x in samples){
#       sub.ob=subset(object,subset= sample== x)
#       countMatix=cbind(countMatrix,Matrix::rowSums(sub.ob[["RNA"]]@counts))
#    }
#    colnames(countMatrix)=samples
#    return(countMatrix)
# }


###monocle 2 拟时序分析
Idents(immune.combined)<- immune.combined$seurat_clusters
trajectory_cluster_dir=paste(project_dir,"06_trajectory_analysis",sep="/")
if(!file.exists(trajectory_cluster_dir)){
   dir.create(trajectory_cluster_dir)
}

sample_name<-unique(immune.combined$sample)[1]
subset_data_for_monocle<-subset(immune.combined,sample==sample_name)
immune.combined_monocle <- run_trajectory(subset_data_for_monocle)
saveRDS(object= immune.combined_monocle,paste0(trajectory_cluster_dir,"/subset_data_for_monocle.rds"))


##singleR
singler_results_dir=paste(project_dir,"08_singleR_analysis",sep="/")
if(!file.exists(singler_results_dir)){
   dir.create(singler_results_dir)
}

cell_anno_obj<-run_singleR_anno(immune.combined,type=opt$type)

##cell chat
cell_inter_dir=paste(project_dir,"07_cellchat_analysis",sep="/")
if(!file.exists(cell_inter_dir)){
   dir.create(cell_inter_dir)
}
cell_chat_obj<-run_cellchat(immune.combined,type=opt$type)
