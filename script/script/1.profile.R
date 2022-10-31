



my_doubletFinder_v3<-function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE,sct = FALSE)
{
    require(Seurat)
    require(fields)
    require(KernSmooth)
    if (reuse.pANN != FALSE) {
        pANN.old <- seu@meta.data[, reuse.pANN]
        classifications <- rep("Singlet", length(pANN.old))
        classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("DF.classifications", pN, pK, nExp,
            sep = "_")] <- classifications
        return(seu)
    }
    if (reuse.pANN == FALSE) {
        real.cells <- rownames(seu@meta.data)
        data <- seu@assays$RNA@counts[, real.cells]
        n_real.cells <- length(real.cells)
        n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
        print(paste("Creating", n_doublets, "artificial doublets...",
            sep = " "))
        real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
        real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
        doublets <- (data[, real.cells1] + data[, real.cells2])/2
        colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
        data_wdoublets <- cbind(data, doublets)
        orig.commands <- seu@commands
        if (sct == FALSE) {
            print("Creating Seurat object...")
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
            print("Normalizing Seurat object...")
            seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                margin = orig.commands$NormalizeData.RNA@params$margin)
            print("Finding variable genes...")
            seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
            print("Scaling data...")
            seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features,
                model.use = orig.commands$ScaleData.RNA$model.use,
                do.scale = orig.commands$ScaleData.RNA$do.scale,
                do.center = orig.commands$ScaleData.RNA$do.center,
                scale.max = orig.commands$ScaleData.RNA$scale.max,
                block.size = orig.commands$ScaleData.RNA$block.size,
                min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
            print("Running PCA...")
            seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features,
                npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca,
                weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                verbose = FALSE)
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,
                PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets)
            gc()
        }
        if (sct == TRUE) {
            require(sctransform)
            print("Creating Seurat object...")
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
            print("Running SCTransform...")
            seu_wdoublets <- SCTransform(seu_wdoublets)
            print("Running PCA...")
            seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,
                PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets)
            gc()
        }
        print("Calculating PC distance matrix...")
        dist.mat <- fields::rdist(pca.coord)
        print("Computing pANN...")
        pANN <- as.data.frame(matrix(0L, nrow = n_real.cells,
            ncol = 1))
        rownames(pANN) <- real.cells
        colnames(pANN) <- "pANN"
        k <- round(nCells * pK)
        for (i in 1:n_real.cells) {
            neighbors <- order(dist.mat[, i])
            neighbors <- neighbors[2:(k + 1)]
            neighbor.names <- rownames(dist.mat)[neighbors]
            pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
        }
        print("Classifying doublets..")
        classifications <- rep("Singlet", n_real.cells)
        classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data),1]
        seu@meta.data[, paste("DF.classifications", "raw",
            sep = "_")] <- classifications
        return(seu)
    }
}





my_doubletFinder_v3_1<- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE,sct = FALSE)
{
    require(Seurat)
    require(fields)
    require(KernSmooth)
    if (reuse.pANN != FALSE) {
        pANN.old <- seu@meta.data[, reuse.pANN]
        classifications <- rep("Singlet", length(pANN.old))
        classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("DF.classifications", "filter",
            sep = "_")] <- classifications
        return(seu)
    }
    if (reuse.pANN == FALSE) {
        real.cells <- rownames(seu@meta.data)
        data <- seu@assays$RNA@counts[, real.cells]
        n_real.cells <- length(real.cells)
        n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
        print(paste("Creating", n_doublets, "artificial doublets...",
            sep = " "))
        real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
        real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
        doublets <- (data[, real.cells1] + data[, real.cells2])/2
        colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
        data_wdoublets <- cbind(data, doublets)
        orig.commands <- seu@commands
        if (sct == FALSE) {
            print("Creating Seurat object...")
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
            print("Normalizing Seurat object...")
            seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                margin = orig.commands$NormalizeData.RNA@params$margin)
            print("Finding variable genes...")
            seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
            print("Scaling data...")
            seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features,
                model.use = orig.commands$ScaleData.RNA$model.use,
                do.scale = orig.commands$ScaleData.RNA$do.scale,
                do.center = orig.commands$ScaleData.RNA$do.center,
                scale.max = orig.commands$ScaleData.RNA$scale.max,
                block.size = orig.commands$ScaleData.RNA$block.size,
                min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
            print("Running PCA...")
            seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features,
                npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca,
                weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                verbose = FALSE)
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,
                PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets)
            gc()
        }
        if (sct == TRUE) {
            require(sctransform)
            print("Creating Seurat object...")
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
            print("Running SCTransform...")
            seu_wdoublets <- SCTransform(seu_wdoublets)
            print("Running PCA...")
            seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,
                PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets)
            gc()
        }
        print("Calculating PC distance matrix...")
        dist.mat <- fields::rdist(pca.coord)
        print("Computing pANN...")
        pANN <- as.data.frame(matrix(0L, nrow = n_real.cells,
            ncol = 1))
        rownames(pANN) <- real.cells
        colnames(pANN) <- "pANN"
        k <- round(nCells * pK)
        for (i in 1:n_real.cells) {
            neighbors <- order(dist.mat[, i])
            neighbors <- neighbors[2:(k + 1)]
            neighbor.names <- rownames(dist.mat)[neighbors]
            pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
        }
        print("Classifying doublets..")
        classifications <- rep("Singlet", n_real.cells)
        classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data),1]
        seu@meta.data[, paste("DF.classifications",pN,pK,"filter",
            sep = "_")] <- classifications
        return(seu)
    }
}

run_trajectory <- function(ser.obj){
  #library needed packages 
  suppressWarnings(library(Seurat))
  suppressWarnings(library(monocle))
  #format data
  print("format data...")
  data <- as(as.matrix(ser.obj@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = ser.obj@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  #Construct monocle cds
  clustered_spleen_monocle <- newCellDataSet(data,
                                             phenoData = pd,
                                             featureData = fd,
                                             lowerDetectionLimit = 0.5,
                                             expressionFamily = negbinomial.size())
  
  #准备分组
  print("format data done , filter select genes ")
  pData(clustered_spleen_monocle)$Cluster<-as.factor(pData(clustered_spleen_monocle)$seurat_clusters) 
  pData(clustered_spleen_monocle)['Cluster_anno']=ser.obj@active.ident
  clustered_spleen_monocle <- estimateSizeFactors(clustered_spleen_monocle)
  clustered_spleen_monocle <- estimateDispersions(clustered_spleen_monocle)
  clustered_spleen_monocle <- detectGenes(clustered_spleen_monocle, min_expr = 0.1)
  expressed_genes2 <- row.names(subset(fData(clustered_spleen_monocle), num_cells_expressed >= 3))
  print("do differentialGeneTest accroding Cluster")
  diff_test_res3 <- differentialGeneTest(clustered_spleen_monocle[expressed_genes2,], fullModelFormulaStr = "~Cluster")
  diff_test_order<-diff_test_res3[order(diff_test_res3$qval),]
  diff_test_order<-diff_test_order[,c("gene_short_name", "pval", "qval")]
  write.table(diff_test_order,paste(trajectory_cluster_dir,"diff_gene_over_cluster.txt",sep="/"),sep="\t",quote=F,row.names=F)
  
  ordering_genes2 <- row.names(subset(diff_test_res3, qval < 0.01))
  clustered_spleen_monocle <- setOrderingFilter(clustered_spleen_monocle, ordering_genes2)
  # clustered_spleen_monocle <- reduceDimension(clustered_spleen_monocle, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
  # clustered_spleen_monocle <- clusterCells(clustered_spleen_monocle, num_clusters = 2)
  clustered_spleen_monocle <- reduceDimension(clustered_spleen_monocle, method = 'DDRTree')
  clustered_spleen_monocle <- orderCells(clustered_spleen_monocle)
  
  p1<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cluster")
  ggsave(p1,file=paste(trajectory_cluster_dir,"monocle_Cluster.pdf",sep="/"),width=12)
  ggsave(p1,file=paste(trajectory_cluster_dir,"monocle_Cluster.png",sep="/"),width=12)

  p2<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "State")
  ggsave(p2,file=paste(trajectory_cluster_dir,"monocle_State.pdf",sep="/"))
  ggsave(p2,file=paste(trajectory_cluster_dir,"monocle_State.png",sep="/"))

  p3<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "Pseudotime")
  ggsave(p3,file=paste(trajectory_cluster_dir,"monocle_Pseudotime.pdf",sep="/"))
  ggsave(p3,file=paste(trajectory_cluster_dir,"monocle_Pseudotime.png",sep="/"))

  p<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "sample")
  ggsave(p,file=paste(trajectory_cluster_dir,"monocle_sample.pdf",sep="/"),width=10)
  ggsave(p,file=paste(trajectory_cluster_dir,"monocle_sample.png",sep="/"),width=10)
  
  pn1<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "State")+facet_wrap(~State, nrow = 1)
  ggsave(pn1,file=paste(trajectory_cluster_dir,"monocle_State_split.pdf",sep="/"))
  ggsave(pn1,file=paste(trajectory_cluster_dir,"monocle_State_split.png",sep="/"))

  cluster_num=length(unique(ser.obj$seurat_clusters))
  if(cluster_num <=10){
    width=10
    height=8
    nrow=1
  }else if(cluster_num >10 && cluster_num <=20){
    nrow=3
    width=10
    height=9
  }else if(cluster_num >20 && cluster_num <=30){
    nrow=4
    width=10
    height=12
  }else{
    nrow=5
    width=10
    height=14
  }
  pn2<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cluster")+facet_wrap(~Cluster, nrow = nrow)
  ggsave(pn2,file=paste(trajectory_cluster_dir,"monocle_Cluster_split.pdf",sep="/"),width=width,height=height)
  ggsave(pn2,file=paste(trajectory_cluster_dir,"monocle_Cluster_split.png",sep="/"),width=width,height=height)
  pn3<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "orig.ident")+facet_wrap(~orig.ident, nrow = nrow)
  ggsave(pn3,file=paste(trajectory_cluster_dir,"monocle_sample_split.pdf",sep="/"),width=10)
  ggsave(pn3,file=paste(trajectory_cluster_dir,"monocle_sample_split.png",sep="/"),width=10)
  #
  print("do differentialGeneTest accroding Pseudotime")
  diff_test_Pseudotime <- differentialGeneTest(clustered_spleen_monocle[expressed_genes2,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_order <- diff_test_Pseudotime[order(diff_test_Pseudotime$qval),]  ###按照q值进行排序
  diff_order <- diff_order[,c("gene_short_name", "pval", "qval")]
  write.table(diff_order,paste(trajectory_cluster_dir,"diff_test_Pseudotime.txt",sep="/"),sep="\t",row.names=F,quote=F)
  
  
  
  select.gene=rownames(diff_order)[1:6]
  cds_subset <- clustered_spleen_monocle[select.gene,]
  p4<-plot_genes_in_pseudotime(cds_subset, color_by = "State")
  ggsave(p4,file=paste(trajectory_cluster_dir,"monocle_gene_State_top6.pdf",sep="/"))
  ggsave(p4,file=paste(trajectory_cluster_dir,"monocle_gene_State_top6.png",sep="/"))
  p7<-plot_genes_in_pseudotime(cds_subset, color_by = "Cluster")
  ggsave(p7,file=paste(trajectory_cluster_dir,"monocle_gene_Cluster_top6.pdf",sep="/"))
  ggsave(p7,file=paste(trajectory_cluster_dir,"monocle_gene_Cluster_top6.png",sep="/"))
  select.gene.top50=rownames(diff_order)[1:50]
  cds_subset_top50=clustered_spleen_monocle[select.gene.top50,]
  p5<-plot_pseudotime_heatmap(cds_subset_top50,num_clusters = 3,cores = 1,show_rownames = T,return_heatmap=TRUE)
  ggsave(p5,file=paste(trajectory_cluster_dir,"monocle_gene_heatmap_top50.pdf",sep="/"))
   ggsave(p5,file=paste(trajectory_cluster_dir,"monocle_gene_heatmap_top50.png",sep="/"))

  p8<-plot_complex_cell_trajectory(clustered_spleen_monocle,color_by = "Cluster")
  ggsave(file=paste(trajectory_cluster_dir,"cell_trajectory_cluster.pdf",sep="/"),p8)
  ggsave(file=paste(trajectory_cluster_dir,"cell_trajectory_cluster.png",sep="/"),p8)
  
  p9<-plot_complex_cell_trajectory(clustered_spleen_monocle,color_by = "State")
  ggsave(file=paste(trajectory_cluster_dir,"cell_trajectory_state.pdf",sep="/"),p9)
  ggsave(file=paste(trajectory_cluster_dir,"cell_trajectory_state.png",sep="/"),p9)
  
  p10<-plot_complex_cell_trajectory(clustered_spleen_monocle,color_by = "Pseudotime")
  ggsave(file=paste(trajectory_cluster_dir,"cell_trajectory_Pseudotime.pdf",sep="/"),p10)
  ggsave(file=paste(trajectory_cluster_dir,"cell_trajectory_Pseudotime.png",sep="/"),p10)  
  
  BEAM_test <- BEAM(clustered_spleen_monocle[expressed_genes2,], branch_point = 1, cores = 1)
  BEAM_order <- BEAM_test[order(BEAM_test$qval),]
  BEAM_order <- BEAM_order[,c("gene_short_name", "pval", "qval")]
  write.table(BEAM_order,paste(trajectory_cluster_dir,"BEAM_test.txt",sep="/"),sep="\t",row.names=F,quote=F)
  
  select.genes=rownames(BEAM_order)[1:100]
  p<-plot_genes_branched_heatmap(clustered_spleen_monocle[select.genes,],branch_point = 1,cores = 4,num_clusters = 3,show_rownames = T,return_heatmap = T)
  pdf(paste(trajectory_cluster_dir,"genes_branched_heatmap_top100.pdf",sep="/"),height=14,width=12)
  p$ph_res
  dev.off()
  
  save(clustered_spleen_monocle,file =paste(trajectory_cluster_dir,"monocle.rda",sep="/"))
  
  return(clustered_spleen_monocle)
}


run_cellchat <- function(sc_object,type="hsa"){
  library(CellChat)
  library(ComplexHeatmap)
  if(type=="hsa"){
    CellChatDB <- CellChatDB.human # set CellChatDB <- CellChatDB.human if working on the human dataset
  }else if(type=="mmu"){
    CellChatDB <- CellChatDB.mouse
  }else{next;
  }
  library(ggalluvial)
  suppressWarnings(library(CellChat))
  suppressWarnings(library(ggalluvial))
  suppressWarnings(library(svglite))
  #prepare info txt
  sc_object$supp<-Idents(sc_object)
  data.input  <- sc_object@assays$RNA@data
  meta<-sc_object@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "supp")

  #add meta data
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "supp") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  cellchat@DB <- CellChatDB
  #calculate
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)

  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  #Calculate the aggregated cell-cell communication network
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  #plot
  #net view
  groupSize <- as.numeric(table(cellchat@idents))
  pdf(paste0(cell_inter_dir,"/","ggregated-cell-cell-communication-network.pdf"))
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  png(paste0(cell_inter_dir,"/","ggregated-cell-cell-communication-network.png"))
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()


  #net view per idents
  mat <- cellchat@net$weight
  circos_out=paste0(cell_inter_dir,"/","circos_plot/")
  if(!dir.exists(circos_out)){
    dir.create(circos_out)
  }
  for (i in 1:nrow(mat)) {

    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    pdf(paste0(circos_out,"/",rownames(mat)[i],".cell-cell-communication-network.pdf"))
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()

    png(paste0(circos_out,"/",rownames(mat)[i],".cell-cell-communication-network.png"))
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
  }


  # Access all the signaling pathways showing significant communications
  pathways.show.all <- cellchat@netP$pathways
  # check the order of cell identity to set suitable vertex.receiver
  # levels(cellchat@idents)
  number_choose=floor(length(unique(cellchat@idents))/2)
  vertex.receiver = seq(1,number_choose)
  netAnalysis_contribution_out=paste0(cell_inter_dir,"/","netAnalysis_contribution_out/")
  if(!dir.exists(netAnalysis_contribution_out)){
    dir.create(netAnalysis_contribution_out)
  }

  hierarchy_out=paste0(cell_inter_dir,"/","hierarchy_plot/")
  if(!dir.exists(hierarchy_out)){
    dir.create(hierarchy_out)
  }

  Chord_out=paste0(cell_inter_dir,"/","Chord_plot/")
  if(!dir.exists(Chord_out)){
    dir.create(Chord_out)
  }

  Circle_out=paste0(cell_inter_dir,"/","Circle_plot/")
  if(!dir.exists(Circle_out)){
    dir.create(Circle_out)
  }


df.net <- subsetCommunication(cellchat)

write.table(df.net,file = paste0(cell_inter_dir,"/cellchat.Communication.net.tsv"),sep = "\t",row.names = F)


for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(netAnalysis_contribution_out,"/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 7, units = 'in', dpi = 300)
  ggsave(filename=paste0(netAnalysis_contribution_out,"/",pathways.show.all[i], "_L-R_contribution.png"), plot=gg, width = 8, height = 7, units = 'in', dpi = 300)

  #netVisual_aggregate plot
  pdf(paste0(hierarchy_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.pdf"))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i],layout = "hierarchy" ,vertex.receiver = vertex.receiver,)
  dev.off()

  png(paste0(hierarchy_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.png"))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i],layout = "hierarchy" ,vertex.receiver = vertex.receiver)
  dev.off()


  #netVisual_aggregate circle plot
  pdf(paste0(Circle_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.circle.pdf"))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i],layout = "circle")
  dev.off()

  png(paste0(Circle_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.circle.png"))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i],layout = "circle")
  dev.off()

  #Chord diagram
  pdf(paste0(Chord_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.Chord.pdf"))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i],layout = "chord")
  dev.off()

  png(paste0(Chord_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.Chord.png"))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i],layout = "chord")
  dev.off()
}

heatmap_out=paste0(cell_inter_dir,"/","heatmap_plot/")
if(!dir.exists(heatmap_out)){
  dir.create(heatmap_out)
}



#plot alldotplot dir
#pdf(paste0(cell_inter_dir,"/","significant-interactions.all.heatmap.pdf"),width = 24, height = 14)
netVisual_bubble(cellchat, remove.isolate = FALSE)
ggsave(paste0(cell_inter_dir,"/","significant-interactions.all.heatmap.pdf"),width = 30, height = 10)

#png(paste0(cell_inter_dir,"/","significant-interactions.all.heatmap.png"))
netVisual_bubble(cellchat, remove.isolate = FALSE)
#dev.off()
ggsave(paste0(cell_inter_dir,"/","significant-interactions.all.heatmap.png"),width = 30, height = 10)

exp_out=paste0(cell_inter_dir,"/","expression_plot/")
if(!dir.exists(exp_out)){
  dir.create(exp_out)
}

for(id in cellchat@netP$pathways){
  # pdf(paste0(exp_out,"/","Pathway.",id,".GeneExpression.pdf"),width =8, height = 12)
  gg <- plotGeneExpression(cellchat, signaling =id)
  ggsave(filename=paste0(exp_out,"/","Pathway.",id,".GeneExpression.pdf"), plot=gg, width = 8, height = 7, units = 'in', dpi = 300)
  ggsave(filename=paste0(exp_out,"/","Pathway.",id,".GeneExpression.png"), plot=gg, width = 8, height = 7, units = 'in', dpi = 300)

}

#Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
mcs_out=paste0(cell_inter_dir,"/","major_contributing_signalingt/")
if(!dir.exists(mcs_out)){
  dir.create(mcs_out)
}

for(id in cellchat@netP$pathways){
  pdf(paste0(mcs_out,"/","Pathway.",id,".netAnalysis_signalingRole_network.pdf"),width =8, height = 12)
  netAnalysis_signalingRole_network(cellchat, signaling = id, width = 10, height = 5, font.size = 10)
  dev.off()

  png(paste0(mcs_out,"/","Pathway.",id,".netAnalysis_signalingRole_network.png"))
  netAnalysis_signalingRole_network(cellchat, signaling = id, width = 10, height = 5, font.size = 10)
  dev.off()

}

  print("basic plot done ")
  # cell_chat_files=list.files(pattern="png$|svg$")
  # dir.create("07_cellchat_analysis/hierarchy")
  # file.copy(from=cell_chat_files,to="07_cellchat_analysis/hierarchy/")
  # file.remove(cell_chat_files)
  #Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest
  ggsave(filename=paste0(cell_inter_dir,"/","dominant_senders_receivers.pdf"), plot=gg1, width = 14, height = 7, units = 'in', dpi = 300)
  ggsave(filename=paste0(cell_inter_dir,"/","dominant_senders_receivers.png"), plot=gg1, width = 14, height = 7, units = 'in', dpi = 300)
  print("netAnalysis_signalingRole_scatter done")

  #Identify signals contributing most to outgoing or incoming signaling of certain cell groups

   # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  pdf(paste0(cell_inter_dir,"/","signals_incoming_outcoming_groups.pdf"), width = 20, height = 14)
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  draw(ht1 + ht2)
  dev.off()

  png(paste0(cell_inter_dir,"/","signals_incoming_outcoming_groups.png"))
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  draw(ht1 + ht2)
  dev.off()
  #Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
  print("Cell Chat NMF")
  suppressWarnings(library(NMF))
  suppressWarnings(library(ggalluvial))
  #此处看图需要定义 nPatterns
  selectK(cellchat, pattern = "outgoing")
  ggsave(paste0(cell_inter_dir,"/","Outgoing.NMF.pattern_select.pdf"),width =18, height = 12)
  ggsave(paste0(cell_inter_dir,"/","Outgoing.NMF.pattern_select.pdf"),width =18, height = 12)
 

  nPatterns = 3
  pdf(paste0(cell_inter_dir,"/","Outgoing.pattern_cluster.pdf"),width =8, height = 12)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  dev.off()
  png(paste0(cell_inter_dir,"/","Outgoing.pattern_cluster.png"))
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,width = 6,height =10)
  dev.off()
  # river plot
  netAnalysis_river(cellchat, pattern = "outgoing")
  ggsave(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_river.pdf"),width =18, height = 12)
  ggsave(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_river.pdf"),width =18, height = 12)
  
  #dot plot
  netAnalysis_dot(cellchat, pattern = "outgoing")
  ggsave(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_dot.pdf"),width =8, height = 12)
  ggsave(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_dot.png"),width =8, height = 12)
  
  #此处看图需要定义 nPatterns
  selectK(cellchat, pattern = "incoming")
  ggsave(paste0(cell_inter_dir,"/","Incoming.NMF.pattern_select.pdf"),width =18, height = 12)
  ggsave(paste0(cell_inter_dir,"/","Incoming.NMF.pattern_select.png"),width =18, height = 12)

  nPatterns = 3
  pdf(paste0(cell_inter_dir,"/","Incoming.pattern_cluster.pdf"),width =8, height = 12)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  dev.off()
  png(paste0(cell_inter_dir,"/","Incoming.pattern_cluster.png"))
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width = 6,height =10)
  dev.off()
  # river plot
  netAnalysis_river(cellchat, pattern = "incoming")
  ggsave(paste0(cell_inter_dir,"/","Incoming.netAnalysis_river.pdf"),width =18, height = 12)
  ggsave(paste0(cell_inter_dir,"/","Incoming.netAnalysis_river.png"),width =18, height = 12)
  

  #dot plot
  netAnalysis_dot(cellchat, pattern = "incoming")
  ggsave(paste0(cell_inter_dir,"/","Incoming.netAnalysis_dot.pdf"),width =8, height = 12)
  ggsave(paste0(cell_inter_dir,"/","Incoming.netAnalysis_dot.png"),width =8, height = 12)

  for (i in 1:length(pathways.show.all)) {

	pdf(paste0(heatmap_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.heatmap.pdf"),width = 8, height = 7)
	# par(mfrow=c(1,1))
	gg <- netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds")
	draw(gg)
	dev.off()

	png(paste0(heatmap_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.heatmap.png"))
	# par(mfrow=c(1,1))
	gg <- netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds")
	draw(gg)
	dev.off()
	}
  df.net <- subsetCommunication(cellchat)
  write.table(df.net,file =paste0(cell_inter_dir,"/", "cellchat.Communication.net.tsv"),sep = "\t",row.names = F) 
  #Identify signaling groups based on their functional similarity
  # cellchat <- computeNetSimilarity(cellchat, type = "functional")
  # cellchat <- netEmbedding(cellchat, type = "functional")
  #> Manifold learning of the signaling networks for a single dataset
  #  cellchat <- netClustering(cellchat, type = "functional")
  #> Classification learning of the signaling networks for a single dataset
  # Visualization in 2D-space
  # pdf(paste0(cell_inter_dir,"/","functional_similarity.netVisual_embedding.pdf"),width =8, height = 12)
  # netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  # dev.off()

  #cellchat <- computeNetSimilarity(cellchat, type = "structural")
  #cellchat <- netEmbedding(cellchat, type = "structural")
  #> Manifold learning of the signaling networks for a single dataset
  # cellchat <- netClustering(cellchat, type = "structural")
  #> Classification learning of the signaling networks for a single dataset
  # Visualization in 2D-space
  #pdf(paste0(cell_inter_dir,"/","structure_similarity.netVisual_embedding.pdf"),width =8, height = 12)
  #netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  #dev.off()
  saveRDS(cellchat, file = paste0(cell_inter_dir,"/","cellchat.rds"))

  return(cellchat)
}





run_singleR_anno <- function(sc_object,type="hsa"){
  if(type=="hsa"){
    ref.se=readRDS("/Business/psn_company/sc02/1.source/hpca.se.rds")
  }else if(type=="mmu"){
    ref.se=readRDS("/Business/psn_company/sc02/1.source/mmu.se.rds")
  }else{next;
	}
  suppressWarnings(library(SingleR))
  #counts <- GetAssayData(sc_object, slot="counts")
  counts<-as.matrix(sc_object[["RNA"]]@counts)
  cluster=sc_object$seurat_clusters
  singler_res <- SingleR(test = counts, ref = ref.se, labels = ref.se$label.main,method="cluster",clusters=cluster)
  print("break1")  
  write.table(singler_res,paste0(singler_results_dir,"/","singleR.result.xls"),sep="\t",quote=F,row.names=T,col.names=NA)
  #plot
  
  cluster_len=length(unique(sc_object$seurat_clusters))-1
  cluster_anno=data.frame(cluster=c(0:cluster_len),label=singler_res$labels)
  
  write.table(cluster_anno,paste0(singler_results_dir,"/","singleR.cell.anno.xls"),sep="\t",quote=F,row.names=F)
  
  new.cluster.ids <-as.vector(singler_res$labels)
  names(new.cluster.ids) <-c(0:cluster_len)
  sc_object@active.ident <- sc_object$seurat_clusters
  sc_object <- RenameIdents(sc_object, new.cluster.ids)
  sc_object$clutser_anno <- Idents(sc_object)
  DimPlot(sc_object,reduction = "umap")
  ggsave(paste0(singler_results_dir,"/","singleR_cell_anno.umap.pdf"),width =8, height = 7)
  ggsave(paste0(singler_results_dir,"/","singleR_cell_anno.umap.png"),width =8, height = 7)
  p1=DimPlot(sc_object,reduction = "umap")
  p2=DimPlot(sc_object,reduction = "umap",group.by = "seurat_clusters")
  p2+p1
  ggsave(paste0(singler_results_dir,"/","singleR_cell_cluster_anno.umap.pdf"),width =15, height = 7)
  ggsave(paste0(singler_results_dir,"/","singleR_cell_cluster_anno.umap.png"),width =15, height = 7)
  p1=DimPlot(sc_object,reduction = "tsne")
  p2=DimPlot(sc_object,reduction = "tsne",group.by = "seurat_clusters")
  p2+p1
  ggsave(paste0(singler_results_dir,"/","singleR_cell_cluster_anno.tsne.pdf"),width =15, height = 7)
  ggsave(paste0(singler_results_dir,"/","singleR_cell_cluster_anno.tsne.png"),width =15, height = 7)
  p1 <- dittoBarPlot(object = 
                 sc_object,var = "clutser_anno",
               group.by = "sample")
  p1
  ggsave(paste0(singler_results_dir,"/","cell_anno_per_sample.pdf"),width =8, height = 7)
  ggsave(paste0(singler_results_dir,"/","cell_anno_per_sample.png"),width =8, height = 7)    
                
}


