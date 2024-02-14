#!/usr/bin/env python

import pandas as pd
import os
import numpy as np
import subprocess
import sys

os.chdir("/data/BB_Bioinformatics/Kevin/BCAC/code")
os.getcwd()

#work on sumdat
#EURO
# sumdatfile="../result/PRS1/euro_training_sumstats.txt"
# sumdat=pd.read_csv(sumdatfile,delimiter="\t")
# sumdat1=pd.DataFrame({"SNP":sumdat['rsid'],"A1":sumdat['a1'],"A2":sumdat['a0'],"BETA":sumdat["beta"],"P":sumdat['p']})
# sumdat1.to_csv("../result/PRS1/prscsx/euro.sumdat.prscsx",index=False,sep="\t")


#EAS
# sumdatfile="../result/PRS1/asian_training_sumstats.txt"
# sumdat=pd.read_csv(sumdatfile,delimiter="\t")
# sumdat1=pd.DataFrame({"SNP":sumdat['rsid'],"A1":sumdat['a1'],"A2":sumdat['a0'],"BETA":sumdat["beta"],"P":sumdat['p']})
# sumdat1.to_csv("../result/PRS1/prscsx/asian.sumdat.prscsx",index=False,sep="\t")

print("start")
print(os.environ['OMP_NUM_THREADS'])
print(os.environ['HOSTNAME'])
print(os.environ['PWD'])
opt=int(sys.argv[1])
print(opt)
#to pick n_gwas: PRSCSx uses total sample size (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8652106/)
#euro 33505+33961+26699+35150=129315
#asian 5940+4882+6859+7545=25226
prscrsdir="/data/BB_Bioinformatics/Kevin/tools/PRScsx/"
def run_prscsx (ref_dir=prscrsdir+"1kg",bim_prefix="../result/PRS1/euro_onco_tuning",
sst_file="../result/PRS1/prscsx/euro.sumdat.prscsx,../result/PRS1/prscsx/asian.sumdat.prscsx",
n_gwas="129315,25226", pop="EUR,EAS",OUTPUT_DIR="../result/PRS1/prscsx",
OUTPUT_FILE_PREFIX="euro_onco_tuning",
PARAM_PHI="1e-6",PHI="1e-06",SEED="1000"):
    #run on each chrom
    for chr in range(22):
        file1out=OUTPUT_DIR+OUTPUT_FILE_PREFIX+"_EUR_pst_eff_a1_b0.5_phi"+PHI+"_chr"+str(chr+1)+".txt"
        file2out=OUTPUT_DIR+OUTPUT_FILE_PREFIX+"_EAS_pst_eff_a1_b0.5_phi"+PHI+"_chr"+str(chr+1)+".txt"
        print(file1out)
        print(file2out)
        if (not os.path.exists(file1out)) or (not os.path.exists(file2out)):
            cmd=[prscrsdir+"PRScsx.py", "--ref_dir="+ref_dir, "--bim_prefix="+bim_prefix,"--chrom="+str(chr+1), "--sst_file="+sst_file, "--n_gwas="+n_gwas, "--pop="+pop, "--out_dir="+OUTPUT_DIR, "--out_name="+OUTPUT_FILE_PREFIX, "--phi="+PARAM_PHI, "--seed="+SEED]
            #cmd=prscrsdir+"PRScsx.py --ref_dir="+ref_dir+" --bim_prefix="+bim_prefix+"--chrom="+str(chr+1)+" --sst_file="+sst_file+" --n_gwas="+n_gwas+" --pop="+pop+" --out_dir="+OUTPUT_DIR+" --out_name="+OUTPUT_FILE_PREFIX +" --phi="+PARAM_PHI+" --seed="+SEED
            print(cmd)
            subprocess.call(cmd,shell=False)
    print("done")

if (opt==1):
    run_prscsx()

if (opt==2):
    run_prscsx(PARAM_PHI="1e-4",PHI="1e-04")

if (opt==3):
    run_prscsx(PARAM_PHI="1e-2",PHI="1e-02")

if (opt==4):
    run_prscsx(PARAM_PHI="1",PHI="1e+00")

if (opt==5):
    run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning")

if (opt==6):
    run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning",PARAM_PHI="1e-4",PHI="1e-04")

if (opt==7):
    run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning",PARAM_PHI="1e-2",PHI="1e-02")

if (opt==8):
    run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning",PARAM_PHI="1",PHI="1e+00")

# if (opt==9):
#     run_prscsx(bim_prefix="../result/PRS1/african_onco_validation",OUTPUT_FILE_PREFIX="african_onco_validation")
# 
# if (opt==10):
#     run_prscsx(bim_prefix="../result/PRS1/african_onco_validation",OUTPUT_FILE_PREFIX="african_onco_validation",PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==11):
#     run_prscsx(bim_prefix="../result/PRS1/african_onco_validation",OUTPUT_FILE_PREFIX="african_onco_validation",PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==12):
#     run_prscsx(bim_prefix="../result/PRS1/african_onco_validation",OUTPUT_FILE_PREFIX="african_onco_validation",PARAM_PHI="1",PHI="1e+00")
# 
# if (opt==13):
#     run_prscsx(bim_prefix="../result/PRS1/african_icogs_validation",OUTPUT_FILE_PREFIX="african_icogs_validation")
# 
# if (opt==14):
#     run_prscsx(bim_prefix="../result/PRS1/african_icogs_validation",OUTPUT_FILE_PREFIX="african_icogs_validation",PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==15):
#     run_prscsx(bim_prefix="../result/PRS1/african_icogs_validation",OUTPUT_FILE_PREFIX="african_icogs_validation",PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==16):
#     run_prscsx(bim_prefix="../result/PRS1/african_icogs_validation",OUTPUT_FILE_PREFIX="african_icogs_validation",PARAM_PHI="1",PHI="1e+00")
# 
# if (opt==17):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_validation",OUTPUT_FILE_PREFIX="asian_onco_validation")
# 
# if (opt==18):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_validation",OUTPUT_FILE_PREFIX="asian_onco_validation",PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==19):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_validation",OUTPUT_FILE_PREFIX="asian_onco_validation",PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==20):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_validation",OUTPUT_FILE_PREFIX="asian_onco_validation",PARAM_PHI="1",PHI="1e+00")
# 
# if (opt==21):
#     run_prscsx(bim_prefix="../result/PRS1/euro_onco_validation",OUTPUT_FILE_PREFIX="euro_onco_validation")
# 
# if (opt==22):
#     run_prscsx(bim_prefix="../result/PRS1/euro_onco_validation",OUTPUT_FILE_PREFIX="euro_onco_validation",PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==23):
#     run_prscsx(bim_prefix="../result/PRS1/euro_onco_validation",OUTPUT_FILE_PREFIX="euro_onco_validation",PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==24):
#     run_prscsx(bim_prefix="../result/PRS1/euro_onco_validation",OUTPUT_FILE_PREFIX="euro_onco_validation",PARAM_PHI="1",PHI="1e+00")
# 
# if (opt==25):
#     run_prscsx(bim_prefix="../result/PRS1/hispanic_onco_validation",OUTPUT_FILE_PREFIX="hispanic_onco_validation")
# 
# if (opt==26):
#     run_prscsx(bim_prefix="../result/PRS1/hispanic_onco_validation",OUTPUT_FILE_PREFIX="hispanic_onco_validation",PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==27):
#     run_prscsx(bim_prefix="../result/PRS1/hispanic_onco_validation",OUTPUT_FILE_PREFIX="hispanic_onco_validation",PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==28):
#     run_prscsx(bim_prefix="../result/PRS1/hispanic_onco_validation",OUTPUT_FILE_PREFIX="hispanic_onco_validation",PARAM_PHI="1",PHI="1e+00")

#swarm --sbatch '--export=MKL_NUM_THREADS=1,NUMEXPR_NUM_THREADS=1,OMP_NUM_THREADS=1' -f /data/BB_Bioinformatics/Kevin/BCAC/code/prscsx.swarm --module python/3.9 -g 64 --time=3-00:00:00 --gres=lscratch:64 
