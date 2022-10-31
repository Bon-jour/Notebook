#mkdir -p summary/02_cellranger
#cp -r SPL202103794/Hru_musc/outs/ summary/02_cellranger/Hru_musc
#cp -r SPL202103794/Hru_skin/outs/ summary/02_cellranger/Hru_skin/
#sed '1d' summary/sample_info.txt | cut -f 1 | while read i
#do
#	mkdir -p summary/02_cellranger/${i}
#	cp -r   ${1}/${i}/outs/cloupe.cloupe  ${1}/${i}/outs/web_summary.html  ${1}/${i}/outs/filtered_feature_bc_matrix  ${1}/${i}/outs/filtered_feature_bc_matrix.h5 summary/02_cellranger/${i}
#done

sample=$(sed -n '2p' summary/sample_info.txt |cut -f 1)
sample_num=$(sed '1d' summary/sample_info.txt |wc -l)
cp -r /Business/psn_company/sc02/1.source/Report .
perl /Business/psn_company/sc02/1.source/10x_webSummary.pl -i CellrangerOut/ -o summary/02_cellranger/
if [ $sample_num -eq 1 ];then cp CellrangerOut/${sample}/outs/web_summary.html Report/table/
else cp aggr/outs/web_summary.html Report/table/
fi
cp -r  summary/02_cellranger/table*txt Report/table/
head -n 100  summary/04_Cluster/cells_GeneCounts.xls |cut -f 1-50 | sed '1s/^/Gene/'> Report/table/table4.txt
sed '1s/^/Sample/' summary/03_QC_Filter/cells_filter_stat.xls >Report/table/table5.txt
cp summary/03_QC_Filter/cells_qc_filter.png Report/pictures/seurat/
cp summary/04_Cluster/pca.ElbowPlot.png Report/pictures/seurat/pca.png
cp summary/04_Cluster/cluster_number.png Report/pictures/seurat/
cp summary/04_Cluster/cluster_umap.png Report/pictures/seurat/umap1.png
cp summary/04_Cluster/cluster_umap_eachsample.png Report/pictures/seurat/umap2.png
cp summary/04_Cluster/cluster_tsne.png Report/pictures/seurat/tsne1.png
cp summary/04_Cluster/cluster_tsne_eachsample.png Report/pictures/seurat/tsne2.png
cp summary/04_Cluster/marker_number.png Report/pictures/seurat/marker_number.png
head -n 20 summary/04_Cluster/all_markers.xls | sed '1s/^/Gene/'> Report/table/table6.txt
cp summary/04_Cluster/allmarkers_top10_vilion.png Report/pictures/seurat/top10_gene.png
cp summary/04_Cluster/allmarkers_top10_dotplot.png Report/pictures/seurat/top10_gene_dotplot.png

cp summary/04_Cluster/allmarkers_top10_cell_exp_distribution.png  Report/pictures/seurat/top10_gene_umap.png

#cp summary/04_Cluster/all_cluster_markers_heatmap.png Report/pictures/seurat/top10_gene_heatmap.png
convert summary/04_Cluster/all_cluster_markers_heatmap.pdf Report/pictures/seurat/top10_gene_heatmap.png
#chose sample data 
for i in summary/04_Cluster/cluster_*/enrichment/; do  num=$(ls ${i}|wc -l); if [ ${num} -gt 9 ];then  cluster_cho=${i}; break;fi; done
head -n 20  ${cluster_cho}/GO_enrichment.xls|cut -f 1-7 > Report/table/table7.txt
cp ${cluster_cho}/GO_enrichment_pvalue_barplot.png Report/pictures/seurat/top_gene_goenrich.png
cp ${cluster_cho}/GO.richfactor.png Report/pictures/seurat/top_gene_goenrich_dot.png

cp ${cluster_cho}/GODAG_MF_top10.png  Report/pictures/seurat/top_gene_goenrich_net.png

head -n 20  ${cluster_cho}/KEGG_enrichment.xls|cut -f 1-7 > Report/table/table8.txt

cp ${cluster_cho}/KEGG_enrichment_pvalue_barplot.png Report/pictures/seurat/top_gene_keggenrich.png

cp ${cluster_cho}/KEGG.richfactor.png Report/pictures/seurat/top_gene_keggenrich_dot.png

head -n 20  ${cluster_cho}/Reactome_enrichment.xls|cut -f 1-6 > Report/table/table9.txt

cp ${cluster_cho}/Reactome_enrichment_pvalue_barplot.png Report/pictures/seurat/top_gene_reactome.png

sample_number=$(sed '1d' summary/sample_info.txt |wc -l)
if [ ${sample_number} -gt 1 ];then
contrast_dir=$(ls summary/05_DiffAnalysis_perCluster/| grep -v xls| cut -f 1|head -n 1 )

head -n 20 summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/diff_gene.xls |sed '1s/^/Gene/'>Report/table/table10.txt

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/top3_diffgene_exp_vilion.png Report/pictures/seurat/diff_gene_expression.png 

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/top3_diffgene_umap.png Report/pictures/seurat/diff_gene_umap.png

head -n 20 summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GO_enrichment.xls|cut -f 1-7 >Report/table/table11.txt

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GO_enrichment_pvalue_barplot.png Report/pictures/seurat/diff_gene_gobarplot.png

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GO.richfactor.png Report/pictures/seurat/diff_gene_godotplot.png

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GODAG_CC_top10.png Report/pictures/seurat/diff_gene_gonet.png

head -n 20 summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/KEGG_enrichment.xls |cut -f 1-7> Report/table/table12.txt

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/KEGG_enrichment_pvalue_barplot.png Report/pictures/seurat/diff_gene_keggbarplot.png

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/KEGG.richfactor.png  Report/pictures/seurat/diff_gene_keggdotplot.png

head -n 20 summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/Reactome_enrichment.xls |cut -f 1-6 > Report/table/table13.txt

cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/Reactome_enrichment_pvalue_barplot.png Report/pictures/seurat/diff_gene_reactome.png
fi

cp summary/08_singleR_analysis/singleR_cell_cluster_anno.umap.png Report/pictures/step_ana/singler.png

cp summary/06_trajectory_analysis/monocle_Pseudotime.png Report/pictures/step_ana/pseudotime1.png

cp summary/06_trajectory_analysis/monocle_Cluster.png Report/pictures/step_ana/pseudotime2.png

cp summary/06_trajectory_analysis/monocle_gene_State_top6.png Report/pictures/step_ana/pseudotime3.png

cp summary/06_trajectory_analysis/monocle_gene_heatmap_top50.png Report/pictures/step_ana/pseudotime4.png

cp summary/07_cellchat_analysis/ggregated-cell-cell-communication-network.png Report/pictures/step_ana/cellchat1-0.png

ls summary/07_cellchat_analysis/hierarchy_plot/*netVisual_aggregate.pdf|sed -n '1p' | while read i ; do id=${i%%.pdf}; convert ${i} Report/pictures/step_ana/cellchat2.png ; done
