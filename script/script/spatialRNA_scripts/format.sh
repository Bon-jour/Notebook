#!/bin/bash
cp -r /Business/psn_company/sc02/1.source/10xST_Report Report/
mkdir temp
cat sample_info.txt |cut -f 1| sed '1d'|while read i
do
	~/1.source/1.soft/csvtk transpose  St_out/${i}/outs/metrics_summary.csv >temp/${i}.temp1
done

~/1.source/1.soft/csvtk join -f 1 temp/*.temp1 >temp/all.temp2

sed -n '1p;2p;3p;9p;11p;12p;13p;14p' temp/all.temp2 | sed 's/,/\t/g' > Report/table/table1.txt

sed -n '1p;22p;4p;6p;23p;7p' temp/all.temp2 | sed 's/,/\t/g' > Report/table/table2.txt

sed -n '1p;16,21p' temp/all.temp2 | sed 's/,/\t/g' > Report/table/table3.txt

cp summary/04_Cluster/Integrated/cluster_summary.xls Report/table/table4.txt

cp summary/04_Cluster/Integrated/UMAP_VS_tSNE.png Report/pictures/analysis/UMAP_VS_tSNE.png

cp summary/04_Cluster/Integrated/tSNE_eachsample.png Report/pictures/analysis/tSNE_eachsample.png

cp summary/04_Cluster/Integrated/SpatialDimPlot.png Report/pictures/analysis/SpatialDimPlot.png

cp summary/04_Cluster/all_markers.xls Report/table/table5.txt

cp summary/04_Cluster/allmarkers_top10_exp_pct.png Report/pictures/analysis/allmarkers_top10_exp_pct.png

cp summary/04_Cluster/Integrated/allmarkers_top10_cell_exp_distribution.umap.png Report/pictures/analysis/allmarkers_top10_cell_exp_distribution.umap.png

ls  summary/04_Cluster/SpatialFeaturePlot_all_top10/*png|sed -n '1p'|while read  i; do cp ${i} Report/pictures/analysis/all_top10.ST.png; done

cp summary/04_Cluster/Integrated/allmarkers_top10_cell_exp_distribution.umap.png Report/pictures/analysis/allmarkers_top10_cell_exp_distribution.umap.png

cp summary/04_Cluster/Integrated/top10_marker_echo_cluster_heatmap.png Report/pictures/analysis/top10_marker_echo_cluster_heatmap.png

cp summary/04_Cluster/cluster_0/enrichment/GO_enrichment.xls Report/table/table6.txt

cp summary/04_Cluster/cluster_0/enrichment/GO_enrichment_pvalue_barplot.png Report/pictures/analysis/cluster_GO_enrichment_pvalue_barplot.png

cp summary/04_Cluster/cluster_0/enrichment/GO.richfactor.png  Report/pictures/analysis/cluster_GO.richfactor.png

cp summary/04_Cluster/cluster_0/enrichment/GODAG_BP_top10.png Report/pictures/analysis/cluster_GODAG_BP_top10.png

cp summary/04_Cluster/cluster_0/enrichment/KEGG_enrichment.xls Report/table/table7.txt

cp summary/04_Cluster/cluster_0/enrichment/KEGG_enrichment_pvalue_barplot.png Report/pictures/analysis/cluster_KEGG_enrichment_pvalue_barplot.png

cp summary/04_Cluster/cluster_0/enrichment/KEGG.richfactor.png Report/pictures/analysis/cluster_KEGG.richfactor.png

sample_number=$(sed '1d' summary/sample_info.txt |wc -l)
if [ ${sample_number} -gt 1 ];then
	contrast_dir=$(ls summary/05_DiffAnalysis_perCluster/| grep -v xls| cut -f 1|head -n 1 )
	cat  summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/diff_gene.xls| sed '1s/^/Gene/' > Report/table/table8.txt
	cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/top3_diffgene_umap.png Report/pictures/analysis/sample_top3_diffgene_umap.png
	cp  summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GO_enrichment.xls Report/table/table9.txt
	cp  summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GO_enrichment_pvalue_barplot.png Report/pictures/analysis/GO_enrichment_pvalue_barplot_group.png
	cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GO.richfactor.png Report/pictures/analysis/GO.richfactor_group.png
	cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/GODAG_BP_top10.png Report/pictures/analysis/GODAG_BP_top10_group.png
	cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/KEGG_enrichment.xls Report/table/table10.txt
	cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/KEGG_enrichment_pvalue_barplot.png Report/pictures/analysis/KEGG_enrichment_pvalue_barplot_group.png
	cp summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster0/enrichment/KEGG.richfactor.png Report/pictures/analysis/KEGG.richfactor.png

fi

