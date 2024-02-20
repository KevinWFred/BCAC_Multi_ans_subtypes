#!/usr/bin/env python

import pandas as pd
import os
import numpy as np
import subprocess
import sys

os.chdir("/data/BB_Bioinformatics/Kevin/BCAC/code")
os.getcwd()
#sumstaf files created in QC_subtype_PRS.R

subtypes=["LumA","LumB","LumB_HN","Her2E","TripN"]
PARAM_PHIs=["1e-6","1e-4","1e-2",1]
PHIs=["1e-06","1e-04","1e-02","1e+00"]
print("start")
print(os.environ['OMP_NUM_THREADS'])
print(os.environ['HOSTNAME'])
print(os.environ['PWD'])
pop=str(sys.argv[1])
print(pop)
subtypeidx=int(sys.argv[2]) #0-4
subtype=subtypes[subtypeidx]
print(subtype)
PHIidx=int(sys.argv[3]) #0-3
PARAM_PHI=PARAM_PHIs[PHIidx]
print(PARAM_PHI)
PHI=PHIs[PHIidx]
print(PHI)
bim_prefix=sys.argv[4] #"../result/PRS1/euro_onco_tuning"
print(bim_prefix)
OUTPUT_FILE_PREFIX=sys.argv[5] #euro_onco_tuning
print(OUTPUT_FILE_PREFIX)

# d1=pd.DataFrame()
# for i in range(5):
#   for j in range(4):
#     temp = pd.DataFrame(
#       {
#         'code':["/data/BB_Bioinformatics/Kevin/BCAC/code/run_subtype_prscsx.py"],
#         'pop':["euro"],
#         'subtypeidx': [i],
#         'PHIidx': [j],
#         'bim_prefix': ["../result/PRS1/euro_onco_tuning"],
#         'OUTPUT_FILE_PREFIX': ["euro_onco_tuning"]
#       }
# )
#     d1=pd.concat([d1, temp])
# 
# d2=pd.DataFrame()
# for i in range(5):
#   for j in range(4):
#     temp = pd.DataFrame(
#       {
#         'code':["/data/BB_Bioinformatics/Kevin/BCAC/code/run_subtype_prscsx.py"],
#         'pop':["asian"],
#         'subtypeidx': [i],
#         'PHIidx': [j],
#         'bim_prefix': ["../result/PRS1/asian_onco_tuning"],
#         'OUTPUT_FILE_PREFIX': ["asian_onco_tuning"]
#       }
# )
#     d2=pd.concat([d2, temp])      
# 
# d=pd.concat([d1, d2])
# d.to_csv("run_subtype_prscsx.swarm",sep="\t",index=False,header=False)
#to pick n_gwas: PRSCSx uses total sample size (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8652106/)
#euro
# LumA    LumB LumB_HN   Her2E   TripN 
# 77252   64117   64515   61926   65149
#asian
# LumA    LumB LumB_HN   Her2E   TripN 
# 15218   14468   13662   13833   13754
samplesize=["77252,15218","64117,14468","64515,13662","61926,13833","65149,13754"]
sst_files=["../result/PRS_subtype/euro/euro_LumA_PRScsx_sumstats.txt,../result/PRS_subtype/asian/asian_LumA_PRScsx_sumstats.txt",
"../result/PRS_subtype/euro/euro_LumB_PRScsx_sumstats.txt,../result/PRS_subtype/asian/asian_LumB_PRScsx_sumstats.txt",
"../result/PRS_subtype/euro/euro_LumB_HN_PRScsx_sumstats.txt,../result/PRS_subtype/asian/asian_LumB_HN_PRScsx_sumstats.txt",
"../result/PRS_subtype/euro/euro_Her2E_PRScsx_sumstats.txt,../result/PRS_subtype/asian/asian_Her2E_PRScsx_sumstats.txt",
"../result/PRS_subtype/euro/euro_TripN_PRScsx_sumstats.txt,../result/PRS_subtype/asian/asian_TripN_PRScsx_sumstats.txt"]
prscrsdir="/data/BB_Bioinformatics/Kevin/tools/PRScsx/"
def run_prscsx (ref_dir=prscrsdir+"1kg",bim_prefix="../result/PRS1/euro_onco_tuning",
sst_file="../result/PRS_subtype/euro/euro_LumA_PRScsx_sumstats.txt,../result/PRS_subtype/asian/asian_LumA_PRScsx_sumstats.txt",
n_gwas="77252,15218", pop="EUR,EAS",OUTPUT_DIR="../result/PRS_subtype/prscsx/euro/LumA",
OUTPUT_FILE_PREFIX="euro_onco_tuning",
PARAM_PHI="1e-6",PHI="1e-06",SEED="1000"):
    #run on each chrom
    for chr in range(22):
        file1out=OUTPUT_DIR+OUTPUT_FILE_PREFIX+"_EUR_pst_eff_a1_b0.5_phi"+PHI+"_chr"+str(chr+1)+".txt"
        file2out=OUTPUT_DIR+OUTPUT_FILE_PREFIX+"_EAS_pst_eff_a1_b0.5_phi"+PHI+"_chr"+str(chr+1)+".txt"
        print(file1out)
        print(file2out)
        if (not os.path.exists(file1out)) or (not os.path.exists(file2out)):
            cmd=[prscrsdir+"PRScsx.py", "--ref_dir="+ref_dir, "--bim_prefix="+bim_prefix,"--chrom="+str(chr+1), "--sst_file="+sst_file, "--n_gwas="+n_gwas, "--pop="+pop, "--out_dir="+OUTPUT_DIR, "--out_name="+OUTPUT_FILE_PREFIX, "--phi="+str(PARAM_PHI), "--seed="+SEED]
            #cmd=prscrsdir+"PRScsx.py --ref_dir="+ref_dir+" --bim_prefix="+bim_prefix+"--chrom="+str(chr+1)+" --sst_file="+sst_file+" --n_gwas="+n_gwas+" --pop="+pop+" --out_dir="+OUTPUT_DIR+" --out_name="+OUTPUT_FILE_PREFIX +" --phi="+PARAM_PHI+" --seed="+SEED
            print(cmd)
            subprocess.call(cmd,shell=False)
    print("done")

outdir="../result/PRS_subtype/prscsx/"+pop+"/"+subtype+"/"
print(outdir)
print(prscrsdir+"1kg")
print(bim_prefix)
print(sst_files[subtypeidx])
print(samplesize[subtypeidx])
print(OUTPUT_FILE_PREFIX)
os.makedirs(outdir, exist_ok=True)
run_prscsx(ref_dir=prscrsdir+"1kg",bim_prefix=bim_prefix,
sst_file=sst_files[subtypeidx],
n_gwas=samplesize[subtypeidx], pop="EUR,EAS",OUTPUT_DIR=outdir,
OUTPUT_FILE_PREFIX=OUTPUT_FILE_PREFIX,
PARAM_PHI=PARAM_PHI,PHI=PHI,SEED="1000")


# if (opt==1):
#     run_prscsx()
# 
# if (opt==2):
#     run_prscsx(PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==3):
#     run_prscsx(PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==4):
#     run_prscsx(PARAM_PHI="1",PHI="1e+00")
# 
# if (opt==5):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning")
# 
# if (opt==6):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning",PARAM_PHI="1e-4",PHI="1e-04")
# 
# if (opt==7):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning",PARAM_PHI="1e-2",PHI="1e-02")
# 
# if (opt==8):
#     run_prscsx(bim_prefix="../result/PRS1/asian_onco_tuning",OUTPUT_FILE_PREFIX="asian_onco_tuning",PARAM_PHI="1",PHI="1e+00")


#swarm --sbatch '--export=MKL_NUM_THREADS=1,NUMEXPR_NUM_THREADS=1,OMP_NUM_THREADS=1' -f /data/BB_Bioinformatics/Kevin/BCAC/code/run_subtype_prscsx.swarm --module python/3.9 -g 64 --time=3-00:00:00 --gres=lscratch:64 
