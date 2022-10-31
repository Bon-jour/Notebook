
suppressPackageStartupMessages(library(optparse))
options(bitmapType='cairo')
option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
	    make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
        make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
        make_option(c("-b", "--RMbE"), help="remove batch effect by function IntegrateData , default %default",default=TRUE),
        make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT"),
        make_option(c("-r", "--mt_regress"), help="regree the mt gene in the SCTransform, default %default",default=FALSE),
        make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="10"),

        make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca  , default %default",default=30),

        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default=FALSE)

)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

################################


################################

if(! is.null(opt$workdir)){
  print(paste("change R work directory to ",opt$workdir,sep=""))
  setwd(opt$workdir)
 }
project_dir=getwd()


#read sample information matrix
#if more than one sample, should give the list csv used in cellranger aggr
cellranger_csv=paste(project_dir,"02_Cellranger","all.csv",sep="/")
if(! is.null(opt$aggr_csv)){
  print(paste("set aggr_csv to ",opt$aggr_csv,sep=""))
  cellranger_csv=opt$aggr_csv
 }

if(exists("cellranger_csv")){
    sample_list=read.csv(cellranger_csv,stringsAsFactors =F)
}else{
    stop("no cellranger_csv file was found\n")
}

if(nrow(sample_list) < 2){
  stop("this pipline need at least 2 sample\n")
}


#import DE contrast information
if(! is.null(opt$contrast)){
  print(paste("get compare contrast file ",opt$contrast,sep=""))
  sample_contrast=read.table(opt$contrast)
 }else{
#sample_contrast=switch(EXPR=as.character(nrow(sample_list)),"2"=combn(sample_list[,1],2), t(combn(sample_list[,1],2)))
 sample_contrast=t(combn(sample_list[,1],2))
 }

##################################
suppressPackageStartupMessages(library(SPOTlight))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(harmony))
options(future.globals.maxSize = 10000 * 1024^2)
###############################

if(opt$mt_regress){
    SCTransform2=function(x,y){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208,vars.to.regress ="percent.mt", verbose = opt$verbose)
    }
}else{
    SCTransform2=function(x,y){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208, verbose = opt$verbose)
    }
}
# 调色板
blue_green_red <- scale_color_gradient2(low = "blue", mid = "green", high = "red")
byr <- CustomPalette(low = '#0000FF', mid = "#FFFF00",high = "#FF0000", k = 7000)
bgr <- CustomPalette(low = '#0000FF', mid = "#00FF00",high = "#FF0000", k = 7000)
#br
#rblack
#whitered
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#getPalette(18)

#########基础空转分析#############

# 读取空转数据

spatial.ob=list()
#VariableFeatures=c()
for(x in 1:nrow(sample_list)){
    datadir=gsub("filtered_feature_bc_matrix.h5","",sample_list[x,2])
    spatial.ob[[x]]<-suppressWarnings(Load10X_Spatial(data.dir = datadir))
    #spatial.ob <- Seurat::SCTransform(spatial.ob, assay ='Spatial',verbose = opt$verbose)
    spatial.ob[[x]]@meta.data$orig.ident<-as.character(sample_list[x,1])
    spatial.ob[[x]]$sample=as.character(sample_list[x,1])     #添加在meta.data slot里面
    spatial.ob[[x]]$orig.ident=as.character(sample_list[x,1]) 
    spatial.ob[[x]]$group=as.character(sample_list[x,4])
    names(spatial.ob[[x]]@images)=as.character(sample_list[x,1])

    spatial.ob[[x]][["percent.mt"]] <- PercentageFeatureSet(spatial.ob[[x]], pattern = opt$mt)
   # saveRDS(spatial.ob[[x]],paste(project_dir,"/","02_Cellranger","/",sample,".rds",sep=""))
}

#spatial.ob=merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_list[,1])
#names(spatial.ob@images)=as.character(sample_list[,1])
#spatial.ob$sample=as.character(sample_list[,1])
#VariableFeatures(spatial.ob) <-VariableFeatures
#cell_num_raw=table(spatial.ob$sample)
#spatial.ob[["percent.mt"]] <- PercentageFeatureSet(spatial.ob, pattern = opt$mt)
names(spatial.ob)=as.character(sample_list[,1])

# 3 QC与SCTransform去除批次效应
spatial.ob <- lapply(X =spatial.ob, FUN =SCTransform2)
if( opt$RMbE ){
    #combined <- SplitObject(spatial.ob, split.by = "sample")

    features <- SelectIntegrationFeatures(object.list = spatial.ob, nfeatures = 3000,verbose = opt$verbose)
    combined <- PrepSCTIntegration(object.list = spatial.ob, anchor.features = features,verbose = opt$verbose)
    anchors <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT", anchor.features = features,verbose = opt$verbose)
    combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",verbose = opt$verbose) 
    # Choose the features to use when integrating multiple datasets. 
    # This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. 
    # It returns the top scoring features by this ranking.
    VariableFeatures(combined.sct)=rownames(combined.sct@assays$integrated@scale.data)
    DefaultAssay(combined.sct) <- "integrated"
}else{

    combined.sct=merge(spatial.ob[[1]],spatial.ob[2:length(spatial.ob)],add.cell.ids=sample_list[,1])
    VariableFeatures(combined.sct)=rownames(combined.sct@assays$SCT@data)
    DefaultAssay(combined.sct) <- "SCT"
}

j=0
for (i in 1:length(combined.sct@images)){
    j=j+1
   if( nrow(combined.sct@images[[j]]@coordinates)==0){
       combined.sct@images[[j]]=NULL
       j=j-1
   }
}
names(combined.sct@images)=as.character(sample_list[,1])

seurat_qc_dir=paste(project_dir,"03_QC_SCT",sep="/")
if(!file.exists(seurat_qc_dir)){
    dir.create(seurat_qc_dir,recursive=TRUE)
}

plot2=plot4=list()

for ( x in 1:nrow(sample_list) ){

    plot2[[x]] <- SpatialFeaturePlot(combined.sct, images=sample_list[x,1],features = "nCount_Spatial") #+ theme(legend.position = "right")
    pic1=wrap_plots(plot2,ncol=2)

    plot4[[x]] <- SpatialFeaturePlot(combined.sct, images=sample_list[x,1],features = "nCount_SCT") #+ theme(legend.position = "right")
    pic2=wrap_plots(plot4,ncol=2)
}

plot1 <- VlnPlot(combined.sct, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3,group.by="sample", pt.size = 0.01)
plot3 <- VlnPlot(combined.sct, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), ncol = 3,group.by="sample", pt.size = 0.01)
## 图1：标准化前的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图
ggsave(plot1, width=12, height=7, filename=paste(seurat_qc_dir,"pre-QC-VlnPlot.pdf",sep="/"),limitsize = FALSE)
ggsave(plot3, width=12, height=7, filename=paste(seurat_qc_dir,"post-QC-VlnPlot.pdf",sep="/"),limitsize = FALSE)
## 图2：标准化后的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图

nrow1=ceiling(length(plot2)/2);
nrow2=ceiling(length(plot4)/2)

ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.pdf",sep="/"),limitsize = FALSE)
ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.png",sep="/"),limitsize = FALSE)
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.pdf",sep="/"),limitsize = FALSE)
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.png",sep="/"),limitsize = FALSE)

fwrite(as.data.frame(combined.sct[["Spatial"]]@counts),paste(seurat_qc_dir,"cells_GeneCounts.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["SCT"]]@data),paste(seurat_qc_dir,"cells_SCT_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["integrated"]]@data),paste(seurat_qc_dir,"cells_integrated_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

rds_dir=paste(project_dir,"suerat_Rdata",sep="/")
if(!file.exists(rds_dir)){
    dir.create(rds_dir,recursive=TRUE)
}
seurat_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
if(!file.exists(seurat_cluster_dir)){
    dir.create(seurat_cluster_dir,recursive=TRUE)
}

combined.sct <- RunPCA(combined.sct,npcs=100, verbose =  opt$verbose)
combined.sct <- RunHarmony(combined.sct, group.by.vars="orig.ident", assay.use="integrated", verbose =  opt$verbose)
## 选择后续分析使用的维度
plot=ElbowPlot(combined.sct,ndim=100,reduction="harmony")
ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir,"ElbowPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir,"ElbowPlot.png",sep="/"),dpi=300,limitsize = FALSE)
## 构建KNN图 
combined.sct <- FindNeighbors(combined.sct,reduction = "harmony", dims = 1:ndim, verbose =  opt$verbose)

save(spatial.ob,file=paste(rds_dir,"separated.sct.Rdata",sep="/"))

saveRDS(combined.sct,paste(rds_dir,"combined.sct.rds",sep="/"))
