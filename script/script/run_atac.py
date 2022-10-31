#coding:utf-8
# A scripts to run multi sample vdj
# Author: Personalbio idsz.Dai
import argparse
import os
from pickle import FALSE
import threading
import shutil
import datetime
import json
import time
run_cmds=open("all.commands",'w')

def execCmd(cmd):
    try:
        print("命令{}开始运行{}".format(cmd, datetime.datetime.now()))
        time.sleep(10)
        os.system(cmd)
        print("命令{}结束运行{}" .format(cmd, datetime.datetime.now()))
    except:
        print("{}\t运行失败".format(cmd))
   # os.system(cmd)

def runQC(data_dir):
    data_dir=os.path.abspath(data_dir)
    cmd="mkdir QC && bsub -q psnpublic -e qc.log -e qc.err /Business/psn_company/t01/local/bin/fastqc -j /Business/psn_company/t01/local/jdk1.8.0_172/bin/java -t 10 -o QC {}/*/*gz".format(data_dir)
    execCmd(cmd)


def run_cellranger(sample_info,data_dir,species):
    #database choose
    # if species=="hsa":
    #     db_choose="/Business/psn_company/t01/public/Database/SingleCell_Reference/refdata-gex-GRCh38-2020-A"
    # elif species=="mmu":
    #     db_choose="/Business/psn_company/t01/public/Database/SingleCell_Reference/refdata-gex-mm10-2020-A"
    genome_config="/Business/psn_company/sc02/1.source/ref.json"
    genome_dic=json.load(open(genome_config))
    db_choose=genome_dic[species]['dir']
    #parse sample info
    sample_file=os.path.abspath(sample_info)
    data_dir=os.path.abspath(data_dir)
    sample_info_dic={}
    cmds=[]
    for line in open(sample_file):
        if not line.startswith("sample"):
            line=line.strip().split()
            sample_info_dic[line[0]]=line[1]
            cmd="bsub -q psnpublic  -e {}.err -o {}.log \
                /Business/psn_company/t01/public/software/cellranger-6.0.0/bin/cellranger count \
                    --id={} --fastqs={} --sample={} \
                        --transcriptome={}"\
                            .format(line[0],line[0],line[0],data_dir,line[0],db_choose)
            cmds.append(cmd)
            run_cmds.write("{}\n".format(cmd))
    threads = []
    for cmd in cmds:
        th = threading.Thread(target=execCmd, args=(cmd,))
        th.start()
        threads.append(th)

    # waiting threading down
    for th in threads:
        th.join()
    return(sample_info_dic)

def run_format_file(sample_info_dic):
    if os.path.exists("CellrangerOut"):
        pass
    else:
        os.mkdir("CellrangerOut")
    for sample in sample_info_dic.keys():
        if os.path.exists("{}/{}".format("CellrangerOut",sample)):
            pass
        else:
            shutil.move(sample,"CellrangerOut")

def check_dir(sample_info_dic):
    check_flag=True
    for sample in sample_info_dic.keys():
        if os.path.exists("{}/outs".format(sample)): ##check cellranger outs results 
            pass
        else:
            check_flag=False
    return(check_flag)                

#library_id,molecule_h5,batch
#LV123,/opt/runs/LV123/outs/molecule_info.h5,v2_lib

def run_cellranger_aggr(sample_info_dic):
    #generate aggr csv
    aggr_out=open("aggr.csv",'w')
    aggr_out.write("sample_id,molecule_h5,batch\n")
    for sample in sample_info_dic.keys():
        h5_dir=os.path.abspath("./CellrangerOut/{}/outs/molecule_info.h5".format(sample))
        aggr_out.write("{},{},{}\n".format(sample,h5_dir,sample_info_dic[sample]))
    aggr_out.close()
    cmd="bsub -q psnpublic -I  -e aggr.cellranger.err -o aggr.cellranger.log \
        /Business/psn_company/t01/public/software/cellranger-6.0.0/bin/cellranger \
            aggr --id=aggr --csv=aggr.csv"
    run_cmds.write("{}\n".format(cmd))
    execCmd(cmd)


def run_step_analysis(sample_info_dic,sample_type,mt_type,mt_cutoff):
    if os.path.exists("summary/02_cellranger"):
        pass
    else:
        os.makedirs("summary/02_cellranger")
    for sample in sample_info_dic.keys():
        if os.path.exists("summary/02_cellranger/{}".format(sample)):
            pass
        else:
            os.makedirs("summary/02_cellranger/{}".format(sample))
        cmd="cp -r  \
            CellrangerOut/{}/outs/cloupe.cloupe  CellrangerOut/{}/outs/web_summary.html  \
                CellrangerOut/{}/outs/filtered_feature_bc_matrix  \
                    CellrangerOut/{}/outs/filtered_feature_bc_matrix.h5 summary/02_cellranger/{}".format(\
                        sample,sample,sample,sample,sample)
        execCmd(cmd)
        #shutil.copyfile("{}/{}/outs/filtered_contig_annotations.csv".format(immune_type,sample),"summary/03_analysis/{}_annotations.csv".format(sample))
    shutil.copyfile("sample_info.txt","summary/sample_info.txt")
    shutil.copyfile("/Business/psn_company/sc02/1.source/seurat.R","summary/seurat.R")
    cmd="cd summary/;bsub -q psnpublic -I  -e seurat.err -o seurat.log Rscript seurat.R -t {} -m {} -f {}".format(sample_type,mt_type,mt_cutoff)
    run_cmds.write("{}\n".format(cmd))
    execCmd(cmd)

def format_file():
    cmd="cp /Business/psn_company/sc02/1.source/format.sh . && bash format.sh"
    execCmd(cmd)


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--sample_info', type=str,default="sample_info.txt",required=True,
                        help='sample info')
    parser.add_argument('--data', required=True,default='1.raw_data',type=str,
                        help='raw data dir')
    parser.add_argument('--species',required=True,type=str,
                        default='hsa',
                        help='choose species,related to genome and database to use ')
    parser.add_argument('--mt',required=True,type=str,
                        default='^MT',
                        help='Single cell mt gene pattern')
    parser.add_argument('--mt_cut',required=True,type=str,
                        default='10',
                        help='Single cell mt gene cut off value')




    # parser.add_argument('--species', dest='accumulate',required=True,
    #                     default='hsa',
    #                     help='choose species,related to genome and database to use , hsa or mmu will be accept')

    args = parser.parse_args()
    runQC(args.data)
    sample_info_dic=run_cellranger(args.sample_info,args.data,args.species)
    time.sleep(10)
    check_flag=check_dir(sample_info_dic)
    while check_flag == False:
        print("Cellranger results have not been parpred...")
        time.sleep(50*60)
        check_flag=check_dir(sample_info_dic)
    run_format_file(sample_info_dic)
    run_cellranger_aggr(sample_info_dic)
    run_step_analysis(sample_info_dic,args.species,args.mt,args.mt_cut)
    format_file()
if __name__ == '__main__':
    main()
