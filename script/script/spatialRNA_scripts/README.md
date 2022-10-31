
# 空转项目

## 标准分析流程

```bash
# 先执行 注意参数需要设置
Rscript spatialRNA.R 

# 完成后执行
Rscript ../spRNA_clustering_intergrated.R -l 0.1 &
Rscript ../spRNA_clustering_intergrated.R -l 0.2 &
......
Rscript ../spRNA_clustering_intergrated.R -l 1 &

Rscript ../spRNA_clustering_separated.R -l 0.1 &
Rscript ../spRNA_clustering_separated.R -l 0.2 &
......
Rscript ../spRNA_clustering_separated.R -l 1 &


```

## spotlight分析

```bash
#  注意参数需要设置
Rscript spotlight.R 

```

## spotlight subcluster分析

```bash
#  注意参数需要设置
Rscript spotlight_new2.R -d suerat_Rdata/Y18_sc_celltype.rds &

Rscript spotlight_new2.R -d suerat_Rdata/Y18_sc_celltype.rds -t Granulosa  &

```
