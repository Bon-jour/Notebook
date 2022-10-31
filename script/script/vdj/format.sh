type=$1
cp -r  /Business/psn_company/sc01/Workbase_Wangcanbiao/script/vdj/report/ .
num=$(sed '1d' sample_info.txt |wc -l)
mkdir temp
if [ ${num} -eq 1 ];then
        sample=$(cut -f 1 sample_info.txt | sed '1d' )
        /Business/psn_company/sc01/Workbase_Wangcanbiao/script/vdj/csvtk -T transpose ${type}/${sample}/outs/metrics_summary.csv | sed "1iCategory\t${sample}">temp/all.table
else
        cut -f 1  sample_info.txt | sed '1d' | while read sample
        do
                /Business/psn_company/sc01/Workbase_Wangcanbiao/script/vdj/csvtk transpose ${type}/${sample}/outs/metrics_summary.csv| sed "1iCategory,${sample}" >temp/${sample}.table.temp
                /Business/psn_company/sc01/Workbase_Wangcanbiao/script/vdj/csvtk join  -T -f 1 temp/*.table.temp  > temp/all.table
        done
fi

sed -n '1p;6,10p' temp/all.table  > report/table/table1.txt
sed -n '1,4p;11,28p' temp/all.table >report/table/table2.txt

sample=$(cut -f 1 sample_info.txt| sed -n '2p')
if [ ${num} -gt 1 ];then
        cp aggr/outs/web_summary.html report/table/web_summary.html
else
        cp ${type}/${sample}/outs/web_summary.html report/table/web_summary.html
fi
num=$(cat summary/02_cellranger/${sample}/filtered_contig_annotations.csv|wc -l)

if [ ${num} -gt 500 ];then
        head -n 500 summary/02_cellranger/${sample}/filtered_contig_annotations.csv|sed 's/,/\t/g' >  report/table/table3.txt
else
        cat summary/02_cellranger/${sample}/filtered_contig_annotations.csv|sed 's/,/\t/g' >  report/table/table3.txt
fi

num=$(cat summary/02_cellranger/${sample}/consensus_annotations.csv|wc -l)

if [ ${num} -gt 500 ];then
        head -n 500 summary/02_cellranger/${sample}/consensus_annotations.csv|sed 's/,/\t/g' >  report/table/table4.txt
else
        cat summary/02_cellranger/${sample}/consensus_annotations.csv|sed 's/,/\t/g' >  report/table/table4.txt
fi


num=$(cat summary/02_cellranger/${sample}/clonotypes.csv|wc -l)

if [ ${num} -gt 500 ];then
        head -n 500 summary/02_cellranger/${sample}/clonotypes.csv|sed 's/,/\t/g' >  report/table/table5.txt
else
        cat summary/02_cellranger/${sample}/clonotypes.csv|sed 's/,/\t/g' >  report/table/table5.txt
fi

num=$(sed '1d' sample_info.txt |wc -l)
if [ ${num} -eq 1 ];then
        cp summary/03_analysis/vdjtools/rareplot.rarefaction.strict.png report/pictures/analysis/rarefaction.png

        cp summary/03_analysis/Clonality/clonal.top_proportion.png report/pictures/analysis/

        cp summary/03_analysis/Clonality/clonal.rare_proportion.png report/pictures/analysis/

        cp summary/03_analysis/Geneusage/Gene_usage.barplot.png report/pictures/analysis/

        convert  summary/03_analysis/vdjtools/${sample}/V_J_gene.pairs.heatmap.pdf report/pictures/analysis/gene.pairs.heatmap.png

        convert  summary/03_analysis/vdjtools/${sample}/V_J_gene.pairs.circos.pdf report/pictures/analysis/gene.pairs.circos.png  

        cp summary/03_analysis/Diversity/diversity_all.png report/pictures/analysis/diversity_all.png

        cp summary/03_analysis/CDR3_feature/CDR3_length_distribution.png report/pictures/analysis/CDR3_length_distribution.png

        cp summary/03_analysis/CDR3_feature/CDR3_abundances_distribution.png report/pictures/analysis/CDR3_abundances_distribution.png
        
        cp summary/03_analysis/CDR3_feature/${sample}.Spectratyping.png report/pictures/analysis/Spectratyping.png
else
	cp summary/03_analysis/Geneusage/Gene_usage.corrplot.png report/pictures/analysis/

        cp summary/03_analysis/rarefaction.png report/pictures/analysis/

        cp summary/03_analysis/Clonality/clonal.prop.png report/pictures/analysis/

        cp summary/03_analysis/Clonality/clonal.top_proportion.png report/pictures/analysis/

        cp summary/03_analysis/Clonality/clonal.rare_proportion.png report/pictures/analysis/

        cp summary/03_analysis/Geneusage/Gene_usage.barplot.png report/pictures/analysis/

        convert  summary/03_analysis/vdjtools/${sample}/V_J_gene.pairs.heatmap.pdf report/pictures/analysis/gene.pairs.heatmap.png

        convert  summary/03_analysis/vdjtools/${sample}/V_J_gene.pairs.circos.pdf report/pictures/analysis/gene.pairs.circos.png  

        cp summary/03_analysis/Diversity/diversity_all.png report/pictures/analysis/diversity_all.png

        cp summary/03_analysis/CDR3_feature/CDR3_length_distribution.png report/pictures/analysis/CDR3_length_distribution.png

        cp summary/03_analysis/CDR3_feature/${sample}.Spectratyping.png report/pictures/analysis/Spectratyping.png

        cp summary/03_analysis/CDR3_feature/CDR3_abundances_distribution.png report/pictures/analysis/CDR3_abundances_distribution.png

        cp summary/03_analysis/Kmer/Kmer_distribution.png report/pictures/analysis/Kmer_distribution.png

        cp summary/03_analysis/Clonotypes_Tracking/top5_clonotypes_tracking.png report/pictures/analysis/top5_clonotypes_tracking.png

        sed 's/\"//g' summary/03_analysis/CDR3_annoation.xls>report/table/table6.txt
fi

