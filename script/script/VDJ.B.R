library(immunarch)
library(limma)
#setwd("E:\\1.work\\vdj\\BCR/")
#LOAD DATA
immdata_10x <- repLoad("1.raw_data/")
immdata_10x$meta$Sample<-strsplit2(immdata_10x$meta$Sample,split = "_annotations")[,1]
names(immdata_10x$data)<-immdata_10x$meta$Sample
exp_vol <- repExplore(immdata_10x$data, .method = "volume")

p1 <- vis(exp_vol)
#gene usage
if(!dir.exists("Geneusage")){
  dir.create("Geneusage")
}

imm_gu <- geneUsage(immdata_10x$data, "macmul.IGHV",.ambig = "exc",.norm = T)
write.table(imm_gu,"Geneusage/gene_usage.xls",sep="\t",row.names = F)
vis(imm_gu[c(1, 2)])+theme(axis.text.x = element_text(size = 7,color="black"))
ggsave("Geneusage/Gene_usage.barplot.pdf",width =14,height = 7)
ggsave("Geneusage/Gene_usage.barplot.png",width = 14,height = 7,dpi = 300)

vis(imm_gu[c(1, 2)])+theme(axis.text.x = element_text(size = 6,color="black"))+coord_polar()
ggsave("Geneusage/Gene_usage.circleplot.png",width = 14,height = 7,dpi = 300)
ggsave("Geneusage/Gene_usage.circleplot.png",width = 14,height = 7,dpi = 300)

imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)

p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 4)+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 4)+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))

p1 + p2

ggsave("Geneusage/Gene_usage.corrplot.pdf",width =14,height = 7)
ggsave("Geneusage/Gene_usage.corrplot.png",width = 14,height = 7,dpi = 300)



#CDR3 length
if(!dir.exists("CDR3_feature")){
  dir.create("CDR3_feature")
}
exp_len <- repExplore(immdata_10x$data, .method = "len", .col = "aa")
write.table("CDR3_feature/CDR3_length_distribution.xls",row.names = F)
exp_cnt <- repExplore(immdata_10x$data, .method = "count")
write.table("CDR3_feature/CDR3_abundances_distribution.xls",row.names = F)
exp_vol <- repExplore(immdata_10x$data, .method = "volume")
write.table("CDR3_feature/clonotypes_summary.xls",row.names = F)

p1 <- vis(exp_len)+ scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))
p2 <- vis(exp_cnt)+ scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))
p3 <- vis(exp_vol)+ scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p1
p2
p3
ggsave(p1,filename = "CDR3_feature/CDR3_length_distribution.pdf",width = 8,height = 6)
ggsave(p1,filename = "CDR3_feature/CDR3_length_distribution.png",width = 8,height = 6)

ggsave(p2,filename = "CDR3_feature/CDR3_abundances_distribution.pdf",width = 8,height = 6)
ggsave(p2,filename = "CDR3_feature/CDR3_abundances_distribution.png",width = 8,height = 6)

ggsave(p3,filename = "CDR3_feature/clonotypes_summary.pdf",width = 8,height = 6)
ggsave(p3,filename = "CDR3_feature/clonotypes_summary.png",width = 8,height = 6)



names(immdata_10x$data)
for (ix in c(1:length(immdata_10x$data))) {
  p1 <- vis(spectratype(immdata_10x$data[[ix]], .quant = "id", .col = "nt"))+theme_bw()
  p2 <- vis(spectratype(immdata_10x$data[[ix]], .quant = "count", .col = "aa+v"))+theme_bw()

  p1 + p2
  outnames <- paste0("CDR3_feature/",names(immdata_10x$data)[ix],".Spectratyping.pdf")
  ggsave(filename = outnames,width = 14,height =7 )
  outnames <- paste0("CDR3_feature/",names(immdata_10x$data)[ix],".Spectratyping.png")
  ggsave(filename = outnames,width = 14,height =7 )
}


#Clonality
if(!dir.exists("Clonality")){
  dir.create("Clonality")
}

imm_pr <- repClonality(immdata_10x$data, .method = "clonal.prop")
imm_pr_out <- imm_pr%>%data.frame()%>%tibble::rownames_to_column("sample")
write.table(x = imm_pr_out,file = "Clonality/clonal.prop.xls",sep="\t",row.names = F)

imm_rare <- repClonality(immdata_10x$data, .method = "rare")
imm_rare_out <- imm_rare%>%data.frame()%>%tibble::rownames_to_column("sample")
write.table(x = imm_rare_out,file = "Clonality/clonal.rare_proportion.xls",sep="\t",row.names = F)


imm_hom <- repClonality(immdata_10x$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

if(length(immdata_10x$data)>1){
imm_top <- repClonality(immdata_10x$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top_out <- imm_top%>%data.frame()%>%tibble::rownames_to_column("sample")
write.table(x = imm_top_out,file = "Clonality/clonal.top_proportion.xls",sep="\t",row.names = F)

vis(imm_top)+scale_fill_brewer(palette = "Set3")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggsave(filename = "Clonality/clonal.top_proportion.pdf",width = 8,height = 5)
ggsave(filename = "Clonality/clonal.top_proportion.png",width = 8,height = 5)

}


vis(imm_pr)+scale_fill_brewer(palette = "Set3")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggsave(filename = "Clonality/clonal.prop.pdf",width = 8,height = 5)
ggsave(filename = "Clonality/clonal.prop.png",width = 8,height = 5)



vis(imm_rare)+scale_fill_brewer(palette = "Set3")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))

ggsave(filename = "Clonality/clonal.rare_proportion.pdf",width = 8,height = 5)
ggsave(filename = "Clonality/clonal.rare_proportion.png",width = 8,height = 5)



vis(imm_hom) +scale_fill_brewer(palette = "Set3")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggsave(filename = "Clonality/clonal.hom_proportion.pdf",width = 8,height = 5)
ggsave(filename = "Clonality/clonal.hom_proportion.png",width = 8,height = 5)


# Compute statistics and visualise them
# Chao1 diversity measure
if(!dir.exists("Diversity")){
  dir.create("Diversity")
}

div_chao <- repDiversity(immdata_10x$data, "chao1")
div_chao_out <- div_chao%>%data.frame()%>%tibble::rownames_to_column("sample")
write.table(x = div_chao_out,file = "Diversity/chao.xls",sep="\t",row.names = F)

# Hill numbers
div_hill <- repDiversity(immdata_10x$data, "hill")
div_hill_out <- div_hill%>%data.frame()%>%tibble::rownames_to_column("sample")
write.table(x = div_hill_out,file = "Diversity/hill.xls",sep="\t",row.names = F)

# D50
div_d50 <- repDiversity(immdata_10x$data, "d50")
div_d50_out <- div_hill%>%data.frame()%>%tibble::rownames_to_column("sample")
write.table(x = div_d50_out,file = "Diversity/d50.xls",sep="\t",row.names = F)



# Ecological diversity measure
div_div <- repDiversity(immdata_10x$data, "div")



p1 <- vis(div_chao)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggsave(plot = p1,filename = "Diversity/chao.pdf",width = 8,height = 5)
ggsave(plot = p1,filename = "Diversity/chao.png",width = 8,height = 5,dpi = 300)

p2 <- vis(div_hill)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))
ggsave(plot = p2,filename = "Diversity/hill.pdf",width = 8,height = 5)
ggsave(plot = p2,filename = "Diversity/hill.png",width = 8,height = 5,dpi = 300)


p3 <- vis(div_d50)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))

ggsave(plot = p3,filename = "Diversity/d50.pdf",width = 8,height = 5)
ggsave(plot = p3,filename = "Diversity/d50.png",width = 8,height = 5,dpi = 300)



p4<-vis(div_div)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.text.x= element_text(angle = 45, hjust = 0.5, vjust = 0.5))


p1 +p2+p3+p4
ggsave(filename = "Diversity/diversity_all.pdf",width = 10,height =8)
ggsave(filename = "Diversity/diversity_all.png",width = 10,height = 8,dpi = 300)


##
imm_raref <- repDiversity(immdata_10x$data, "raref", .verbose = F)

p1 <- vis(imm_raref)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))

p1
ggsave(plot = p1,filename = "rarefaction.pdf",width = 10,height = 7)
ggsave(plot = p1,filename = "rarefaction.png",width = 10,height = 7,dpi = 300)

##

if(length(immdata_10x$data)>1){
  
if(!dir.exists("Kmer")){
  dir.create("Kmer")
}


kmers <- getKmers(immdata_10x$data, 3)
kmers <- kmers[- grep(";", kmers$Kmer),]

vis(kmers)+
  theme(text=element_text(family="serif"))+
  theme(axis.text.x= element_text(size=6,angle = 45, hjust = 0.5, vjust = 0.5))
ggsave("Kmer/Kmer_distribution.pdf",width = 12,height = 6)
ggsave("Kmer/Kmer_distribution.png",width = 12,height = 6,dpi = 300)


p1 <- vis(kmers, .head = 5)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))
p2 <- vis(kmers, .head = 10)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))
p3 <- vis(kmers, .head = 30)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=16,  family="serif"))

(p1 + p2) / p3
ggsave("Kmer/Kmer_top.pdf",width = 10,height = 6)
ggsave("Kmer/Kmer_top.png",width = 10,height = 6,dpi = 300)


##Tracking of clonotypes
if(!dir.exists("Clonotypes_Tracking")){
  dir.create("Clonotypes_Tracking")
}

#for (ix in c(1:length(immdata_10x$data))) {}
tc1 <- trackClonotypes(immdata_10x$data, list(1, 5), .col = "nt")
tc1$CDR3.nt<-paste0("clonotype",c(1:5))

p1 <- vis(tc1)+scale_fill_brewer(palette = "Set2")+
  theme(text=element_text(size=8,  family="serif"),legend.position = "bottom")
ggsave(plot = p1,"Clonotypes_Tracking/top5_clonotypes_tracking.pdf",width = 10,height = 7)
ggsave(plot = p1,"Clonotypes_Tracking/top5_clonotypes_tracking.png",width = 10,height = 7,dpi = 300)

#p2 <- vis(tc2)+scale_fill_brewer(palette = "Set2")+
 # theme(text=element_text(size=8,  family="serif"))
#ggsave(plot = p1,"Clonotypes_Tracking/top5_clonotypes_tracking.pdf",width = 10,height = 7)
#ggsave(plot = p1,"Clonotypes_Tracking/top5_clonotypes_tracking.png",width = 10,height = 7,dpi = 300)

for (ix in c(1:length(immdata_10x$data))) {
  tc <- trackClonotypes(immdata_10x$data, list(names(immdata_10x$data)[ix], 5), .col = "aa")
  out_pdf <- paste0("Clonotypes_Tracking/",names(immdata_10x$data)[ix],".top5_clonotypes_tracking.pdf")
  ggsave(out_pdf,width = 10,height = 7)
  out_png <- paste0("Clonotypes_Tracking/",names(immdata_10x$data)[ix],".top5_clonotypes_tracking.png")
  ggsave(out_png,width = 10,height = 7,dpi=300)
}




#db anno


#vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb")
vdjdb=readRDS("/Business/psn_company/sc02/1.source/vdjdb.rds")
#saveRDS(vdjdb,file = "vdjdb.rds")
anno_data <- dbAnnotate(immdata_10x$data, vdjdb, "CDR3.aa", "cdr3")
anno_data_res<-vdjdb[vdjdb$cdr3  %in%  anno_data$CDR3.aa,]
anno_data_res<-as.data.frame(anno_data_res)
anno_data_merge <- dplyr::left_join(anno_data_res,anno_data,c("cdr3" = "CDR3.aa"))
write.table(x = anno_data_merge,"CDR3_annoation.xls",row.names = F,sep="\t")
}

