#!/usr/bin/env Rscript
#apply HWE/info score cutoff on BCAC data

library(data.table)
library(tidyverse)
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
tmp=read.table("../result/imp_onco/asian/asian.psam")
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

#MAF>0.005
#QCed SNP lists:
african_onco=as.data.frame(fread("../result/imp_onco/african/african.pvar"))
print(dim(african_onco)) #18816456
asian_onco=as.data.frame(fread("../result/imp_onco/asian/asian.pvar"))
print(dim(asian_onco)) #10083629
euro_onco=as.data.frame(fread("../result/imp_onco/euro/euro.pvar"))
print(dim(euro_onco)) #10347052
hispanic_onco=as.data.frame(fread("../result/imp_onco/hispanic/hispanic.pvar"))
print(dim(hispanic_onco)) #11133499
african_icogs=as.data.frame(fread("../result/imp_icogs/african/african.pvar"))
print(dim(african_icogs)) #18629686
asian_icogs=as.data.frame(fread("../result/imp_icogs/asian/asian.pvar"))
print(dim(asian_icogs)) #10589245
euro_icogs=as.data.frame(fread("../result/imp_icogs/euro/euro.pvar"))
print(dim(euro_icogs)) #10339918

#--hard-call-threshold 0.4999
filter_hwe=function(prefix="../result/imp_onco/asian/asian",outprefix="../result/imp_QC/onco/asian_hwe")
{
  #plink2 --pfile ... --hwe --keep-if PHENO1 == control --write-snplist --out hwe_pass
  cmd=paste0(plink2," --pfile ",prefix," --hwe 1e-5 --keep-if pheno == 1 --write-snplist --out ",outprefix," --memory 200000 --threads 14")
  system(cmd)
  cmd=paste0(plink2," --pfile ",prefix," --hwe 1e-6 --keep-if pheno == 1 --write-snplist --out ",outprefix,"1e6"," --memory 200000 --threads 14")
  system(cmd)
}
filter_hwe(prefix="../result/imp_onco/african/african",outprefix="../result/imp_QC/onco/african_hwe")
filter_hwe(prefix="../result/imp_onco/asian/asian",outprefix="../result/imp_QC/onco/asian_hwe")
filter_hwe(prefix="../result/imp_onco/euro/euro",outprefix="../result/imp_QC/onco/euro_hwe")
filter_hwe(prefix="../result/imp_onco/hispanic/hispanic",outprefix="../result/imp_QC/onco/hispanic_hwe")

filter_hwe(prefix="../result/imp_icogs/african/african",outprefix="../result/imp_QC/icogs/african_hwe")
filter_hwe(prefix="../result/imp_icogs/asian/asian",outprefix="../result/imp_QC/icogs/asian_hwe")
filter_hwe(prefix="../result/imp_icogs/euro/euro",outprefix="../result/imp_QC/icogs/euro_hwe")

filter_info=function(prefix="../result/imp_QC/onco/asian_hwe1e6",infoprefix="../data/imp_infoscore/data/asian_onco/",
                     outprefix="../result/imp_QC/onco/asian_info")
{
  print(outprefix)
  info=as.data.frame(fread(paste0(infoprefix,"infoscore.txt")))
  hwe=as.data.frame(fread(paste0(prefix,".snplist"),header=F))
  print(paste0("hwe: ",nrow(hwe)))
  info1=info %>% filter(R2>0.3) %>% select(ID)
  info2=intersect(hwe[,1],info1$ID)
  print(paste0("info: ",length(info2)))
  fwrite(as.matrix(info2),file=paste0(outprefix,".snplist"),row.names = F,col.names = F,quote=F)
}

filter_info(prefix="../result/imp_QC/onco/african_hwe1e6",infoprefix="../data/imp_infoscore/data/african_onco/",
            outprefix="../result/imp_QC/onco/african_info")
filter_info(prefix="../result/imp_QC/onco/asian_hwe1e6",infoprefix="../data/imp_infoscore/data/asian_onco/",
                     outprefix="../result/imp_QC/onco/asian_info")
filter_info(prefix="../result/imp_QC/onco/euro_hwe1e6",infoprefix="../data/imp_infoscore/data/euro_onco/",
            outprefix="../result/imp_QC/onco/euro_info")
filter_info(prefix="../result/imp_QC/onco/hispanic_hwe1e6",infoprefix="../data/imp_infoscore/data/hispanic_onco/",
            outprefix="../result/imp_QC/onco/hispanic_info")

filter_info(prefix="../result/imp_QC/icogs/african_hwe1e6",infoprefix="../data/imp_infoscore/data/african_icogs/",
            outprefix="../result/imp_QC/icogs/african_info")
filter_info(prefix="../result/imp_QC/icogs/asian_hwe1e6",infoprefix="../data/imp_infoscore/data/asian_icogs/",
            outprefix="../result/imp_QC/icogs/asian_info")
filter_info(prefix="../result/imp_QC/icogs/euro_hwe1e6",infoprefix="../data/imp_infoscore/data/euro_icogs/",
            outprefix="../result/imp_QC/icogs/euro_info")

# [1] "../result/imp_QC/onco/african_info"
# [1] "hwe: 18801351"
# [1] "info: 18755874"
# 
# [1] "../result/imp_QC/onco/asian_info"
# [1] "hwe: 9035988"
# [1] "info: 8405662"
# 
# [1] "../result/imp_QC/onco/euro_info"
# [1] "hwe: 10134904"
# [1] "info: 9928681"
# 
# [1] "../result/imp_QC/onco/hispanic_info"
# [1] "hwe: 11116135"
# [1] "info: 10958799"
# 
# [1] "../result/imp_QC/icogs/african_info"
# [1] "hwe: 18594971"
# [1] "info: 18456758"
# 
# [1] "../result/imp_QC/icogs/asian_info"
# [1] "hwe: 8118922"
# [1] "info: 6221720"
# 
# [1] "../result/imp_QC/icogs/euro_info"
# [1] "hwe: 8319711"
# [1] "info: 7922354"

#prefix: only based on MAF of 0.005
get_qc_plinkdata=function(prefix="../result/imp_onco/euro/euro",outprefix="../result/imp_QC/onco/euro/euro",
                          snpfile="../result/imp_QC/onco/euro_info.snplist")
{
  cmd=paste0(plink2," --pfile ",prefix," --extract ",snpfile," --make-pgen --out ",outprefix," --memory 128000 --threads 8")
  system(cmd)
}

get_qc_plinkdata(prefix="../result/imp_onco/euro/euro",outprefix="../result/imp_QC/onco/euro/euro",
                          snpfile="../result/imp_QC/onco/euro_info.snplist")

get_qc_plinkdata(prefix="../result/imp_onco/asian/asian",outprefix="../result/imp_QC/onco/asian/asian",
                 snpfile="../result/imp_QC/onco/asian_info.snplist")

get_qc_plinkdata(prefix="../result/imp_onco/african/african",outprefix="../result/imp_QC/onco/african/african",
                 snpfile="../result/imp_QC/onco/african_info.snplist")

get_qc_plinkdata(prefix="../result/imp_icogs/euro/euro",outprefix="../result/imp_QC/icogs/euro/euro",
                 snpfile="../result/imp_QC/icogs/euro_info.snplist")

get_qc_plinkdata(prefix="../result/imp_icogs/asian/asian",outprefix="../result/imp_QC/icogs/asian/asian",
                 snpfile="../result/imp_QC/icogs/asian_info.snplist")

