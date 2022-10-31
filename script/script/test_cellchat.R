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
  
  

for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(netAnalysis_contribution_out,"/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 7, units = 'in', dpi = 300)
  ggsave(filename=paste0(netAnalysis_contribution_out,"/",pathways.show.all[i], "_L-R_contribution.png"), plot=gg, width = 8, height = 7, units = 'in', dpi = 300)
  
  #netVisual_aggregate plot
  pdf(paste0(hierarchy_out,"/",pathways.show.all[i], "_L-R_netVisual_aggregate.pdf"),width=18,height=12)
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
  
#plot alldotplot dir 
pdf(paste0(cell_inter_dir,"/","significant-interactions.all.heatmap.pdf"),width = 24, height = 14)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()

png(paste0(cell_inter_dir,"/","significant-interactions.all.heatmap.png"))
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()

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
  pdf(paste0(cell_inter_dir,"/","Outgoing.NMF.pattern_select.pdf"),width =8, height = 12)
  selectK(cellchat, pattern = "outgoing")
  dev.off()
  png(paste0(cell_inter_dir,"/","Outgoing.NMF.pattern_select.png"))
  selectK(cellchat, pattern = "outgoing")
  dev.off()
  
  nPatterns = 3
  pdf(paste0(cell_inter_dir,"/","Outgoing.pattern_cluster.pdf"),width =8, height = 12)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  dev.off()
  png(paste0(cell_inter_dir,"/","Outgoing.pattern_cluster.png"))
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,width = 6,height =10)
  dev.off()
  # river plot
  pdf(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_river.pdf"),width =8, height = 12)
  netAnalysis_river(cellchat, pattern = "outgoing")
  dev.off()
  png(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_river.png"))
  netAnalysis_river(cellchat, pattern = "outgoing")
  dev.off()
  #dot plot
  pdf(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_dot.pdf"),width =8, height = 12)
  netAnalysis_dot(cellchat, pattern = "outgoing")
  dev.off()
  png(paste0(cell_inter_dir,"/","Outgoing.netAnalysis_dot.png"))
  netAnalysis_dot(cellchat, pattern = "outgoing")
  dev.off()
  #此处看图需要定义 nPatterns
  pdf(paste0(cell_inter_dir,"/","Incoming.NMF.pattern_select.pdf"),width =8, height = 12)
  selectK(cellchat, pattern = "incoming")
  dev.off()
  png(paste0(cell_inter_dir,"/","Incoming.NMF.pattern_select.png"))
  selectK(cellchat, pattern = "incoming")
  dev.off()
  nPatterns = 3
  pdf(paste0(cell_inter_dir,"/","Incoming.pattern_cluster.pdf"),width =8, height = 12)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  dev.off()
  png(paste0(cell_inter_dir,"/","Incoming.pattern_cluster.png"))
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width = 6,height =10)
  dev.off()
  # river plot
  pdf(paste0(cell_inter_dir,"/","Incoming.netAnalysis_river.pdf"),width =8, height = 12)
  netAnalysis_river(cellchat, pattern = "incoming")
  dev.off()
  png(paste0(cell_inter_dir,"/","Incoming.netAnalysis_river.png"))
  netAnalysis_river(cellchat, pattern = "incoming")
  dev.off()
  
  #dot plot
  pdf(paste0(cell_inter_dir,"/","Incoming.netAnalysis_dot.pdf"),width =8, height = 12)
  netAnalysis_dot(cellchat, pattern = "incoming")
  dev.off()
  png(paste0(cell_inter_dir,"/","Incoming.netAnalysis_dot.png"))
  netAnalysis_dot(cellchat, pattern = "incoming")
  dev.off()
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
