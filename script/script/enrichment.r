enrichment<-function(species,outDir,geneList){
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(org.Rn.eg.db)
    library(ggplot2)
    library(dplyr)
    species.org=switch(species,
                        "hsa"=org.Hs.eg.db,
                        "mmu"=org.Mm.eg.db,
                        "rno"=org.Rn.eg.db,
                    )
    species.kegg=switch(species,
                        "hsa"="/Business/psn_company/sc01/local/KEGG.db/hsa/",
                        "mmu"="/Business/psn_company/sc01/local/KEGG.db/mmu/",
			"rno"="/Business/psn_company/sc01/local/KEGG.db/rno/",
                    )
    library(KEGG.db,lib.loc=species.kegg)
    if(species=="hsa"){
    EG2symbol=toTable(org.Hs.egSYMBOL)
}else if(species=="mmu"){
    EG2symbol=toTable(org.Mm.egSYMBOL)
}else if(species=="rno"){
    EG2symbol=toTable(org.Rn.egSYMBOL)
}
    #EG2symbol=toTable(paste0(species.org,"SYMBOL"))
    if("TRUE" %in% unique(geneList%in%EG2symbol$symbol)){
    gene.ncbi=bitr(geneList, fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = species.org)
    
    processGO(geneList,species.org,outDir)
    processKEGG(gene.ncbi[,2],species,outDir)
    processREACTOME(gene.ncbi[,2],species,outDir)
}}


processGO<- function(geneList,species.org,outDir){
    en <- enrichGO(gene       = geneList,
                OrgDb         = species.org,
                keyType       = 'SYMBOL',
                ont           = "ALL",
                minGSSize     = 1 ,
                maxGSSize     = 100000 ,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
    if(!is.null(en)){
    res=en@result
    res$Total=apply(res,1,function(x){getBGnumber(x[5])})
    df=data.frame( Category=res$ONTOLOGY,
                GO=res$ID, 
				Term=res$Description, 
				List=res$Count, 
				Total=res$Total,
				Pvalue=res$pvalue,
				adjustPvalue=res$p.adjust,
				Gene=res$geneID
				)
    df<-df[order(df$Pvalue),]
    write.table(df,file=paste(outDir,"/","GO_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')

    if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]}
    out_df$Term=factor(out_df$Term,levels=unique(out_df$Term))
    out_df$Term=substring(out_df$Term,1,50)
    out_df<-out_df[!duplicated(out_df[,"Term"]),]
    out_df$rich=out_df$List / out_df$Total
    out_df$Number=out_df$List
    q<-qplot(rich,Term,data=out_df,colour=Pvalue,size=Number,main="GO Enrichment")+
	scale_colour_gradient(low="red",high="green",limits=c(0,1))+
	theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey80",linetype="dotted"))
    ggsave(paste(outDir,"GO.richfactor.pdf",sep="/"),width = 10, height = 10)
    ggsave(paste(outDir,"GO.richfactor.png",sep="/"),width = 10, height = 10)
   
    out_df=out_df[order(out_df$Category),]

    out_df$Category <- factor(out_df$Category, levels = unique(out_df$Category))
    out_df$Term <- factor(out_df$Term, levels = out_df$Term)
    p<- ggplot(data = out_df) +
    geom_col(aes(x = Term, y = -log10(as.numeric(Pvalue)), fill = Category),width =0.8) +
    scale_color_brewer(type="seq",palette="Dark2")+
    theme( strip.text.y = element_text(angle = 0),axis.text.x=element_text(size=10,angle=80,hjust=1)) +
    ylab("-log10(P-value)")
    ggsave(paste(outDir,"GO_enrichment_pvalue_barplot.pdf",sep="/"), width = 10, height = 10)
    ggsave(paste(outDir,"GO_enrichment_pvalue_barplot.png",sep="/"), width = 10, height = 10)

    library(topGO)
    for(catergory in c("CC","BP",'MF')){
        ego <- enrichGO(gene       = geneList,
                OrgDb         = species.org,
                keyType       = 'SYMBOL',
                ont           = catergory,
                minGSSize     = 1 ,
                maxGSSize     = 100000 ,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
        ego@result$p.adjust[ego@result$p.adjust==0]=1e-300
        pdf(paste(outDir,"/","GODAG_",catergory,"_top10",".pdf",sep=""))
        plotGOgraph(ego)
        dev.off()
        png(paste(outDir,"/","GODAG_",catergory,"_top10",".png",sep=""))
        plotGOgraph(ego)
        dev.off()
    
    }
}
}


processKEGG<- function(geneList,species,outDir){
    species.kegg=switch(species,
                        "hsa"="/Business/psn_company/sc01/local/KEGG.db/hsa/",
                        "mmu"="/Business/psn_company/sc01/local/KEGG.db/mmu/",
			"rno"="/Business/psn_company/sc01/local/KEGG.db/rno/",
                    )
    organism=strsplit(species.kegg,split="/")[[1]][length(strsplit(species.kegg,split="/")[[1]])]
    library(KEGG.db,lib.loc=species.kegg)
	en <- enrichKEGG(gene     = geneList,
                organism      = organism,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                minGSSize     = 1 ,
                maxGSSize     = 100000 ,
                qvalueCutoff  = 1,
                use_internal_data =T)

    kegg_level=read.table("/Business/psn_company/sc01/local/KEGG.db/pathway_level",header=F,sep="\t",stringsAsFactors=F,colClasses="character",row.names=1)
    #if(length(unique(en@result$p.adjust<0.05))>1)
    if(!is.null(en)){
    res=en@result
    res$Total=apply(res,1,function(x){getBGnumber(x[4])})
    df=data.frame( Category=kegg_level[gsub(species,"",res$ID),1] ,
                PathwayID=res$ID, 
				Pathway=res$Description, 
				List=res$Count, 
				Total=res$Total,
				Pvalue=res$pvalue,
				adjustPvalue=res$p.adjust,
				Gene=res$geneID
				)
    write.table(df,file=paste(outDir,"/","KEGG_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
    if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]}
    out_df$Pathway=factor(out_df$Pathway,levels=unique(out_df$Pathway))
    out_df$Pathway=substring(out_df$Pathway,1,50)
    out_df<-out_df[!duplicated(out_df[,"Pathway"]),]
    out_df$rich=out_df$List / out_df$Total
    out_df$Number=out_df$List
    q<-qplot(rich,Pathway,data=out_df,colour=Pvalue,size=Number,main="KEGG Enrichment")+
	scale_colour_gradient(low="red",high="green",limits=c(0,1))+
	theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey80",linetype="dotted"))
    ggsave(paste(outDir,"KEGG.richfactor.pdf",sep="/"),width = 10, height = 10)
    ggsave(paste(outDir,"KEGG.richfactor.png",sep="/"),width = 10, height = 10)
   
    out_df=out_df[order(out_df$Category),]
    out_df$Category <- factor(out_df$Category, levels = unique(out_df$Category))
    out_df$Pathway <- factor(out_df$Pathway, levels = out_df$Pathway)
    p<- ggplot(data = out_df) +
    geom_col(aes(x = Pathway, y = -log10(as.numeric(Pvalue)), fill = Category),width =0.8) +
    scale_color_brewer(type="seq",palette="Dark2")+
    theme( strip.text.y = element_text(angle = 0),axis.text.x=element_text(size=10,angle=80,hjust=1)) +
    ylab("-log10(P-value)")
    ggsave(paste(outDir,"KEGG_enrichment_pvalue_barplot.pdf",sep="/"), width = 10, height = 10)
    ggsave(paste(outDir,"KEGG_enrichment_pvalue_barplot.png",sep="/"), width = 10, height = 10)

}}

processREACTOME <- function(geneList,species,outDir){
    library(ReactomePA)
    species.reactome=switch(species,
                            'hsa'='human', 
                            'mmu'='mouse', 
                            'rno'='rat'
                            )
    en <- enrichPathway(gene       = geneList,
                organism      = species.reactome,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                minGSSize     = 1 ,
                maxGSSize     = 100000 ,
                qvalueCutoff  = 1)
 #   if(length(unique(en@result$p.adjust<0.05))>1){  
    if(!is.null(en)){
    res=en@result
    res$Total=apply(res,1,function(x){getBGnumber(x[4])})
    df=data.frame( PathwayID=res$ID, 
				Pathway=res$Description, 
				List=res$Count, 
				Total=res$Total,
				Pvalue=res$pvalue,
				adjustPvalue=res$p.adjust,
				Gene=res$geneID
				)
    write.table(df,file=paste(outDir,"/","Reactome_enrichment.xls",sep=""),col.names=T,row.names=F,quote=F,sep='\t')

    #Plot
    if(nrow(df) < 20){out_df <- df[seq(1,nrow(df),1),]}else{out_df <- df[seq(1,20,1),]}
    out_df$Pathway=factor(out_df$Pathway,levels=unique(out_df$Pathway))
    out_df$Pathway=substring(out_df$Pathway,1,50)
    out_df<-out_df[!duplicated(out_df[,"Pathway"]),]
    out_df$rich=out_df$List / out_df$Total
    out_df$Number=out_df$List
    q<-qplot(rich,Pathway,data=out_df,colour=Pvalue,size=Number,main="Reactome Enrichment")+
	scale_colour_gradient(low="red",high="green",limits=c(0,1))+
	theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey80",linetype="dotted"))
    ggsave(paste(outDir,"Reactome.richfactor.pdf",sep="/"),width = 10, height = 10)
    ggsave(paste(outDir,"Reactome.richfactor.png",sep="/"),width = 10, height = 10)


    p<- ggplot(data = out_df) +
    geom_col(aes(y = Pathway, x = -log10(as.numeric(Pvalue))),width =0.8) +
    theme( strip.text.y = element_text(angle = 0),axis.text.x=element_text(size=10,angle=80,hjust=1)) +
    xlab("-log10(P-value)")
    ggsave(paste(outDir,"Reactome_enrichment_pvalue_barplot.pdf",sep="/"), width = 10, height = 10)
    ggsave(paste(outDir,"Reactome_enrichment_pvalue_barplot.png",sep="/"), width = 10, height = 10)

}}

richfactor<-function(en){
    total <- apply(en,1,function(x){getBGnumber(x[4])})
    return(en$Count/total)
}
getBGnumber <- function(ratio,split="/"){
	list<-strsplit(ratio, split = split)[[1]]
	return(as.numeric(list[1]))
}

