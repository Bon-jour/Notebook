#coding:utf-8
# A scripts to run multi sample vdj
# Author: Personalbio idsz.Dai
import argparse
import os
import threading
import shutil
import datetime
run_cmds=open("all.commands",'w')
def execCmd(cmd):
    try:
        print("命令{}\n开始运行{}\n".format(cmd, datetime.datetime.now()))
        os.system(cmd)
        print("命令{}\n结束运行{}\n" .format(cmd, datetime.datetime.now()))
    except:
        print("{}\t运行失败".format(cmd))
    #os.system(cmd)

def run_cellranger(sample_info,data_dir,species,immune_type):
    #database choose
    if species=="hsa":
        db_choose="/Business/psn_company/sc02/1.source/vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
    elif species=="mmu":
        db_choose="/Business/psn_company/sc02/1.source/vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0/"
    if immune_type=="TCR":
        chain_type="TR"
    elif immune_type=="BCR":
        chain_type="IG"
    #parse sample info
    sample_file=os.path.abspath(sample_info)
    data_dir=os.path.abspath(data_dir)
    sample_info_dic={}
    cmds=[]
    for line in open(sample_file):
        if not line.startswith("sample"):
            line=line.strip().split()
            sample_info_dic[line[0]]=line[1]
            cmd="bsub -q psnpublic -I  -e {}.cellranger.err -o {}.cellranger.log /Business/psn_company/t01/public/software/cellranger-6.0.0/bin/cellranger vdj\
            --id={} --fastqs={}/{} --sample={}  --reference={}\
                --localcores=6 --chain={}".format(line[0],line[0],line[0],data_dir,line[0],line[0],db_choose,chain_type)
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

def run_format_file(sample_info_dic,immune_type):
    if os.path.exists(immune_type):
        pass
    else:
        os.mkdir(immune_type)
    if os.path.exists("summary/02_cellranger"):
        pass
    else:
        os.makedirs("summary/02_cellranger")
    for sample in sample_info_dic:
        if os.path.exists("{}/{}".format(immune_type,sample)):
            pass
        else:
            shutil.move(sample,immune_type)
        if os.path.exists("{}/{}".format(immune_type,sample)):
            pass
        else:
            os.makedirs("{}/{}".format(immune_type,sample))

        if os.path.exists("summary/02_cellranger/{}".format(sample)):
            pass
        else:
            os.makedirs("summary/02_cellranger/{}".format(sample))
        cmd="cp -r  \
                {}/{}/outs/consensus_annotations.csv  {}/{}/outs/clonotypes.csv  \
                    {}/{}/outs/filtered_contig_annotations.csv  \
                        summary/02_cellranger/{}".format(\
                            immune_type,sample,immune_type,sample,immune_type,sample,sample)
        execCmd(cmd)
def run_vdj_aggr(sample_info_dic,immune_type):
    #generate aggr csv
    aggr_out=open("aggr.csv",'w')
    aggr_out.write("sample_id,vdj_contig_info,donor,origin\n")
    for sample in sample_info_dic.keys():
        pb_dir=os.path.abspath("./{}/{}/outs/vdj_contig_info.pb".format(immune_type,sample))
        aggr_out.write("{},{},{},{}\n".format(sample,pb_dir,sample_info_dic[sample],sample))
    aggr_out.close()
    cmd="bsub -q psnpublic -I  -e aggr.cellranger.err -o aggr.cellranger.log /Business/psn_company/t01/public/software/cellranger-6.0.0/bin/cellranger aggr --id=aggr --csv=aggr.csv"
    run_cmds.write("{}\n".format(cmd))
    execCmd(cmd)

def run_step_analysis(sample_info_dic,immune_type):
    if os.path.exists("summary/03_analysis"):
        pass
    else:
        os.makedirs("summary/03_analysis")
    if os.path.exists("summary/03_analysis/1.raw_data"):
        pass
    else:
         os.makedirs("summary/03_analysis/1.raw_data")
    for sample in sample_info_dic.keys():
        shutil.copyfile("{}/{}/outs/filtered_contig_annotations.csv".format(immune_type,sample),"summary/03_analysis/1.raw_data/{}_annotations.csv".format(sample))
    if immune_type=="TCR":
        shutil.copyfile("/Business/psn_company/sc02/1.source/VDJ.T.R","summary/03_analysis/VDJ.R")
    elif immune_type=="BCR":
        shutil.copyfile("/Business/psn_company/sc02/1.source/VDJ.B.R","summary/03_analysis/VDJ.R")
    cmd="cd summary/03_analysis;/Business/psn_company/sc01/software/R-4.0.5/bin/Rscript VDJ.R"
    run_cmds.write("{}\n".format(cmd))
    execCmd(cmd)

def run_vdjtools_step(sample_info_dic,immune_type):
    if os.path.exists("summary/03_analysis/vdjtools"):
        pass
    else:
        os.makedirs("summary/03_analysis/vdjtools")
    for sample in sample_info_dic:
        cmd="cp -r {}/{}/outs/consensus_annotations.csv summary/03_analysis/vdjtools/{}.consensus_annotations.csv && \
            cp -r {}/{}/outs/clonotypes.csv summary/03_analysis/vdjtools/{}.clonotypes.csv && \
            cp -r {}/{}/outs/filtered_contig_annotations.csv summary/03_analysis/vdjtools/{}.filtered_contig_annotations.csv && cp -r {}/{}/outs/web_summary.html  summary/02_cellranger/{}/ && cp -r sample_info.txt summary/03_analysis/vdjtools/"\
                .format(immune_type,sample,sample,immune_type,sample,sample,immune_type,sample,sample,immune_type,sample,sample)
        execCmd(cmd)
    if immune_type=="BCR":
        run_src="run_vdjtools.B.sh"
    elif immune_type=="TCR":
        run_src="run_vdjtools.T.sh"
    cmd="cd summary/03_analysis/vdjtools &&\
        cp /Business/psn_company/sc02/1.source/1.soft/{} ./run_vdjtools.sh &&\
            bash  run_vdjtools.sh".format(run_src)
    execCmd(cmd)

def format(immune_type):
    cmd="cp /Business/psn_company/sc01/Workbase_Wangcanbiao/script/vdj/format.sh . && bash format.sh {} ".format(immune_type)
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
    parser.add_argument('--immune',required=True,type=str,
                        default='TCR',
                        help='choose immune rep databases,related to immune rep databases to use,\
                            default database TCR. Other option is BCR')


    # parser.add_argument('--species', dest='accumulate',required=True,
    #                     default='hsa',
    #                     help='choose species,related to genome and database to use , hsa or mmu will be accept')

    args = parser.parse_args()

    sample_info_dic=run_cellranger(args.sample_info,args.data,args.species,args.immune)
    run_format_file(sample_info_dic,args.immune)
    run_vdj_aggr(sample_info_dic,args.immune)
    run_step_analysis(sample_info_dic,args.immune)
    run_vdjtools_step(sample_info_dic,args.immune)
    format(args.immune)
if __name__ == '__main__':
    main()

