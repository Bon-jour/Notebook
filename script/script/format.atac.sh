#!/bin/bash
cp -r ~/1.source/report_atac/  report
mkdir temp
sed '1d' sample_info.txt |cut -f 1 |while read i
do
        csvtk transpose ATAC/${i}/outs/summary.csv > temp/${i}.raw.csv
done

sample_num=$(sed '1d' sample_info.txt |wc -l)
if [ ${sample_num} -gt 1 ];then
        csvtk join -f 1 temp/*raw.csv -T >temp/all.table
else
        sed 's/,/\t/g' temp/*raw.csv > temp/all.table
fi
sed -n '1p;27p;19,23p' temp/all.table >report/table/table1.txt

sed -n '1p;4p;8p;11p;14p;15p' temp/all.table >report/table/table2.txt

sed -n '1p;5p;12p;13p;16p;26p' temp/all.table >report/table/table3.txt

sed '1s/^/Sample/' summary/03_QC_Filter/cells_filter_stat.xls >report/table/table4.txt

cp summary/03_QC_Filter/cells_qc_filter.png report/pictures/signac/cells_qc_filter.png

cp summary/03_QC_Filter/TSS_enrichment_score.png report/pictures/signac/peaks_enrich_score.png

cp summary/03_QC_Filter/nuleosomeSignal_fragmentLengthHistogram.png report/pictures/signac/nucl_signal.png

cp summary/04_Cluster/cluster_number.png report/pictures/signac/cluster_number.png

head -n 20 summary/04_Cluster/cluster_special_peaks.xls > report/table/table5.txt

cp summary/04_Cluster/allpeaks_top4_vilion.png report/pictures/signac/top4_peaks.png

cp summary/04_Cluster/allpeaks_top4_cell_exp_distribution.png report/pictures/signac/top4_peaks_umap.png

for i in summary/04_Cluster/cluster_*
do
        if [ -f ${i}/enrichment/GO_enrichment.xls ];then
                choose_dir=${i}
                break
        fi
done

cp ${choose_dir}/topcoveragePlot.png report/pictures/signac/top4_peaks_region.png

head -n 20  ${choose_dir}/*peaks_closetGene.xls > report/table/table6.txt

head -n 20  ${choose_dir}/enrichment/GO_enrichment.xls |cut -f  1-7 > report/table/table7.txt

cp ${choose_dir}/enrichment/GO_enrichment_pvalue_barplot.png report/pictures/signac/top_gene_goenrich.png

cp ${choose_dir}/enrichment/GO.richfactor.png report/pictures/signac/top_gene_goenrich_dot.png

cp ${choose_dir}/enrichment/GODAG_BP_top10.png report/pictures/signac/top_gene_goenrich_net.png

head -n 20  ${choose_dir}/enrichment/KEGG_enrichment.xls |cut -f  1-7 > report/table/table8.txt

cp ${choose_dir}/enrichment/KEGG_enrichment_pvalue_barplot.png report/pictures/signac/top_gene_keggenrich.png

cp ${choose_dir}/enrichment/KEGG.richfactor.png report/pictures/signac/top_gene_keggenrich_dot.png

head -n 20  ${choose_dir}/enrichment/Reactome_enrichment.xls|cut -f 1-6 > report/table/table9.txt

cp ${choose_dir}/enrichment/Reactome_enrichment_pvalue_barplot.png report/pictures/signac/top_gene_reactome.png

head -n 20 ${choose_dir}/*enriched.motifs.xls.xls > report/table/table10.txt

cp ${choose_dir}/top4_enrichedMotif.png report/pictures/signac/peaks_motif.png

num=$(sed '1d' sample_info.txt|cut -f 2 |sort |uniq |wc -l)

if [ $num -gt 1 ];then
        contrast_dir=$(ls summary/05_DiffAnalysis_perCluster/| grep -v xls| cut -f 1|head -n 1 )
        for i in summary/05_DiffAnalysis_perCluster/${contrast_dir}/cluster*
        do
                if [ -f ${i}/enrichment/GO_enrichment.xls ];then
                        contrast_dir=${i}
                        break
                fi
        done

        head -n 20 ${contrast_dir}/diff_peaks.xls > report/table/table11.txt

        cp ${contrast_dir}/top4_vilion.png report/pictures/signac/diff_peaks_vilion.png

        cp ${contrast_dir}/top4_cell_exp_distribution.png report/pictures/signac/diff_peaks_umap.png

        head -n 20  ${contrast_dir}/diff_peaks_closetGene.xls > report/table/table12.txt

        head -n 20 ${contrast_dir}/enrichment/GO_enrichment.xls|cut -f 1-7 >report/table/table13.txt

        cp ${contrast_dir}/enrichment/GO_enrichment_pvalue_barplot.png report/pictures/signac/diff_gene_gobarplot.png

        cp ${contrast_dir}/enrichment/GO.richfactor.png report/pictures/signac/diff_gene_godotplot.png

        cp ${contrast_dir}/enrichment/GODAG_CC_top10.png report/pictures/signac/diff_gene_gonet.png

        head -n 20 ${contrast_dir}/enrichment/KEGG_enrichment.xls |cut -f 1-7> report/table/table14.txt

        cp ${contrast_dir}/enrichment/KEGG_enrichment_pvalue_barplot.png report/pictures/signac/diff_gene_keggbarplot.png

        cp ${contrast_dir}/enrichment/KEGG.richfactor.png  report/pictures/signac/diff_gene_keggdotplot.png

        head -n 20 ${contrast_dir}/enrichment/Reactome_enrichment.xls |cut -f 1-6 > report/table/table15.txt

        cp ${contrast_dir}/enrichment/Reactome_enrichment_pvalue_barplot.png report/pictures/signac/diff_gene_reactome.png

        head -n 20 ${contrast_dir}/diff_enriched.motifs.xls > report/table/table16.txt

        cp ${contrast_dir}//top4_enrichedMotif.png report/pictures/signac/diff_peaks_motif.png

fi

