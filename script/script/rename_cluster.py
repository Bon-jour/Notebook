#coding:utf-8
# A scripts to run sc-rna recluster 
# Author: Personalbio idsz.Dai
import argparse
import datetime
import os
def execCmd(cmd):
    try:
        print("命令{}开始运行{}".format(cmd, datetime.datetime.now()))
        os.system(cmd)
        print("命令{}结束运行{}" .format(cmd, datetime.datetime.now()))
    except:
        print("{}\t运行失败".format(cmd))

def run_reclust(rds_file,cluster_anno,out_dir):
    rds_file=os.path.abspath(rds_file)
    cluster_anno=os.path.abspath(cluster_anno)
    out_dir=os.path.abspath(out_dir)
    os.system("df1")
    os.system("cp /Business/psn_company/sc02/1.source/rename_cluster.R ./")
    cmd="bsub -q psnpublic -I -e run.err -o run.log  \
        Rscript rename_cluster.R -r {} -a {} -o {}"\
        .format(rds_file,cluster_anno,out_dir)
    execCmd(cmd)

def main():
    parser = argparse.ArgumentParser(description='Recluster seurat object ')

    parser.add_argument('--rds', type=str,default="A rds data contain processed seurat obj \
        always can be found in 04_Cluster/All_sample_combined.rds",required=True,
                        help='A rds data contain processed seurat obj')
    parser.add_argument('--anno', required=True,default='A filelist contain choose cluster  \
        to do recluster,one cluster per line',type=str,
                        help='raw data dir')
    parser.add_argument('--out',required=True,type=str,
                        default='out_file',
                        help='Out dir')      
    args = parser.parse_args()

    run_reclust(rds_file=args.rds,cluster_anno=args.anno,out_dir=args.out)

if __name__ == '__main__':
    main()
