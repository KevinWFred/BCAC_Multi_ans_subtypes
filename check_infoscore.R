#!/usr/bin/env Rscript

library(data.table)

# The columns in the European files are:
#   CHROM POS END REF ALT ID AF MAF R2
# 
# Columns in the African and Asian files are:
#   ID CHROM POS REF ALT AF R2
# 
# The European and Asian R2 are from just one imputation batch not an average across all batches.
get_infoscore=function(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_icogs/",opt=0)
{
  print(folder)
  print(Sys.time())
  allres=NULL
  for (chr in 1:23)
  {
    if(opt==1)
    {
      tmp=as.data.frame(fread(paste0(folder,"topmed_icogs_index_chr",chr,".txt")))
    }else
    {
      tmp=as.data.frame(fread(paste0(folder,"topmed_index_chr",chr,".txt")))
    }
    
    if (ncol(tmp)==9)
    {
      colnames(tmp)=c("CHR","POS","END","REF","ALT","ID","AF","MAF","R2")
    }else
    {
      colnames(tmp)=c("ID","CHR","POS","REF","ALT","AF","R2")
    }
    allres=rbind(allres,tmp)
  }
  print(dim(allres))
  fwrite(allres,file=paste0(folder,"infoscore.txt"),row.names = F,sep="\t",quote=F)
  #return(allres)
  print(Sys.time())
}
#euro_icogs=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_icogs/",opt=1)
# asian_icogs=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_icogs/")
# african_icogs=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/african_icogs/")
# euro_onco=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_onco/")
# asian_onco=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_onco/")
# african_onco=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/african_onco/")
# hispanic_onco=get_infoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/hispanic_onco/")

#get info score only on snps in plink files
get_myinfoscore=function(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_icogs/",prefix="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian")
{
  tmp=as.data.frame(fread(paste0(folder,"infoscore.txt")))
  pvar=as.data.frame(fread(paste0(prefix,".pvar")))
  idx=match(pvar$ID,tmp$ID)
  if (sum(is.na(idx)>0)) warning("SNPs not found!")
  tmp1=tmp[idx,]
  fwrite(tmp1,file=paste0(folder,"myinfoscore.txt"),row.names = F,sep="\t",quote=F)
}
folders=c("/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_icogs/",
          "/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/african_icogs/",
          "/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_icogs/",
          "/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_onco/",
          "/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/african_onco/",
          "/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/euro_onco/",
          "/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/hispanic_onco/")
prefixes=c("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian",
           "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african",
           "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro",
           "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian",
           "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african",
           "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro",
           "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/hispanic/hispanic")
#get_myinfoscore(folder="/data/BB_Bioinformatics/Kevin/BCAC/data/imp_infoscore/data/asian_icogs/",prefix="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian")
args = commandArgs(trailingOnly=TRUE)
i1=as.numeric(args[1])
print(i1)
folder=folders[i1]
print(folder)
prefix=prefixes[i1]
print(prefix)
get_myinfoscore(folder=folder,prefix=prefix)
