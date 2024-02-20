#!/usr/bin/env Rscript
# to select SNPs for EUR/ASN separately. The SNPs are in icogs and onco; having MAF>0.01 in icogs/onco.
#create genotype data for tuning/validation, snpid use rsid.
#use snps from onco

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
library(data.table)
library(dplyr)
subtypes=c("LumA","LumB","LumB_HN","Her2E","TripN")

#compute effective sample sizes for each pop/subtype, can be used for LDpred2,not used though
#training samples
trainingsamples_icogs=read.table("../result/PRS_icogs_trainingpheno.txt",sep="\t",header=T)
trainingsamples_onco=read.table("../result/PRS_onco_trainingpheno.txt",sep="\t",header=T)
phenoicogs=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
phenoicogs$ID=phenoicogs$SG_ID
idx=which(colnames(phenoicogs) %in% paste0("pc",1:10))
phenoicogs[,idx]=NULL
tmp=read.table("../result/imp_icogs/merged1.eigenvec")
colnames(tmp)=c("ID",paste0("pc",1:20))
tmp1=unlist(strsplit(tmp$ID,"_"))
tmp$ID=tmp1[seq(1,length(tmp1),2)]
idx=match(phenoicogs$ID,tmp$ID)
phenoicogs=cbind(phenoicogs,tmp[idx,2:11])
phenoicogs$y=NA
phenoicogs$y[which(phenoicogs$Behaviour1==1)]=1
phenoicogs$y[which(is.na(phenoicogs$Behaviour1))]=0

phenoonco=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
phenoonco$ID=phenoonco$Onc_ID
colnames(phenoonco)=gsub("PC_","pc",colnames(phenoonco))
phenoonco$y=NA
phenoonco$y[which(phenoonco$Behaviour1==1)]=1
phenoonco$y[which(is.na(phenoonco$Behaviour1))]=0

call_subtype=function(pheno=phenoicogs)
{
  pheno$LumA=pheno$LumB=pheno$LumB_HN=pheno$Her2E=pheno$TripN=NA
  idx=which(is.na(pheno$Behaviour1))
  pheno$LumA[idx]=pheno$LumB[idx]=pheno$LumB_HN[idx]=pheno$Her2E[idx]=pheno$TripN[idx]=0
  y.pheno.mis1 <- pheno[,c("ER_status1","PR_status1","HER2_status1","Grade1")]
  
  idx.1 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &(y.pheno.mis1[,4]==1|y.pheno.mis1[,4]==2))
  pheno$LumA[idx.1]=1
  #define Luminal-B like
  idx.2 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==1)
  pheno$LumB[idx.2]=1
  #for Luminal B HER2 negative-like
  idx.3 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &y.pheno.mis1[,4]==3)
  pheno$LumB_HN[idx.3]=1
  #for HER2 enriched-like
  idx.4 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==1)
  pheno$Her2E[idx.4]=1
  #for Triple negative
  idx.5 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==0)
  pheno$TripN[idx.5]=1
  return(pheno)
  
}
phenoicogs=call_subtype(pheno=phenoicogs)
phenoonco=call_subtype(pheno=phenoonco)
compute_train_size_subtype=function(pop="euro")
{
  phenoicogs1=phenoicogs[phenoicogs$ID %in% trainingsamples_icogs$SG_ID,]
  phenoonco1=phenoonco[phenoonco$ID %in% trainingsamples_onco$Onc_ID,]
  if (pop=="euro")
  {
    phenoicogs1=phenoicogs1[phenoicogs1$EthnicityGeno %in% "European",]
    phenoonco1=phenoonco1[phenoonco1$EthnicityGeno %in% "European",]
  }
  if (pop=="asian")
  {
    phenoicogs1=phenoicogs1[phenoicogs1$EthnicityGeno %in% "Asian",]
    phenoonco1=phenoonco1[phenoonco1$EthnicityGeno %in% "Asian",]
  }
  effsize=rep(NA,length(subtypes))
  names(effsize)=subtypes
  totalsize=effsize
  for (i in 1:length(subtypes))
  {
    subtype=subtypes[i]
    effsize[i]=1/(1/sum(phenoicogs1[,subtype]==0,na.rm = T)+1/(sum(phenoicogs1[,subtype]==1,na.rm = T)))+
      1/(1/sum(phenoonco1[,subtype]==0,na.rm = T)+1/(sum(phenoonco1[,subtype]==1,na.rm = T)))
    totalsize[i]=(sum(phenoicogs1[,subtype]==0,na.rm = T))+(sum(phenoicogs1[,subtype]==1,na.rm = T))+
      (sum(phenoonco1[,subtype]==0,na.rm = T))+(sum(phenoonco1[,subtype]==1,na.rm = T))
  }
  return(list(effsize,totalsize))
}
euro_size=compute_train_size_subtype(pop="euro")
# [[1]]
# LumA      LumB   LumB_HN     Her2E     TripN 
# 13095.591  3642.110  3997.814  1670.667  4562.831 
# 
# [[2]]
# LumA    LumB LumB_HN   Her2E   TripN 
# 77252   64117   64515   61926   65149
asian_size=compute_train_size_subtype(pop="asian")
# [[1]]
# LumA      LumB   LumB_HN     Her2E     TripN 
# 1982.8726 1435.8404  796.8670  944.8316  882.3290 
# 
# [[2]]
# LumA    LumB LumB_HN   Her2E   TripN 
# 15218   14468   13662   13833   13754
#training SNPs after QC (logistic)
eurallsumstat=as.data.frame(fread("../result/PRS1/euro_training_sumstats.txt"))
eurallctslebsumstat=as.data.frame(fread("../result/PRS1/ctsleb/EUR/EUR_CTSLEB.sumdat"))

pop="euro"
load(paste0("../result/PRS_subtype/",pop,"/subtype_sumstats.RData"))
table(allres$SNP %in% eurallsumstat$SNP)
idx=match(eurallsumstat$SNP,allres$SNP)
eursumstat=allres[idx,]
tmp=unlist(strsplit(eursumstat$SNP,":"))
alt=tmp[seq(3,length(tmp),4)]
ref=tmp[seq(4,length(tmp),4)]
all(eursumstat$a1==alt) #T
table(eursumstat$a1==eurallsumstat$a1)
table(eurallctslebsumstat$A1==alt)

#sumstat for PRScsx
#sumdat1=pd.DataFrame({"SNP":sumdat['rsid'],"A1":sumdat['a1'],"A2":sumdat['a0'],"BETA":sumdat["beta"],"P":sumdat['p']})
idx=match(eursumstat$SNP,eurallsumstat$SNP)
eursumstat$rsid=eurallsumstat$rsid[idx]

form_PRScsxsumstat=function(pop="euro",sumstat=eursumstat)
{
  sumstat0=data.frame(SNP=sumstat$rsid,A1=sumstat$a1,A2=sumstat$a0)
  for (i in 1:length(subtypes))
  {
    sumstat1=data.frame(BETA=sumstat[,paste0("beta_",subtypes[i])],P=sumstat[,paste0("p_",subtypes[i])])
    fwrite(cbind.data.frame(sumstat0,sumstat1),file=paste0("../result/PRS_subtype/",pop,"/",pop,"_",subtypes[i],"_PRScsx_sumstats.txt"),
           row.names = F,sep="\t",quote=F)
  }
}
form_PRScsxsumstat(pop="euro",sumstat=eursumstat)

form_CTSLEB_sumdat=function(sumdat=eursumstat,pop="euro",bimprefix="../result/PRS1/euro_onco_tuning")
{
  bim=as.data.frame(fread(paste0(bimprefix,".bim")))
  idx=match(sumdat$rsid,bim$V2)
  print(sum(is.na(idx)))
  bim=bim[idx,]
  bim$name=paste0("chr",bim$V1,":",bim$V4,":",bim$V5,":",bim$V6)
  print(all(sumdat$SNP==bim$name)) #alleles are consistent
  
  idx=which(colnames(sumdat)=="chr")
  if (length(idx)>0) colnames(sumdat)[idx]="CHR"
  idx=which(colnames(sumdat) %in% c("pos","position"))
  if (length(idx)>0) colnames(sumdat)[idx]="BP"
  idx=which(colnames(sumdat) %in% c("effect_allele","a1"))
  if (length(idx)>0) colnames(sumdat)[idx]="A1"
  idx=which(colnames(sumdat) %in% c("reference_allele","a0"))
  if (length(idx)>0) colnames(sumdat)[idx]="A2"
  sumdat$SNP=sumdat$rsid
  sumdat0=sumdat[,c("CHR","BP","A1","A2","SNP")]
  for (i in 1:length(subtypes))
  {
    subtype=subtypes[i]
    sumdat1=data.frame(P=sumdat[,paste0("p_",subtype)],BETA=sumdat[,paste0("beta_",subtype)],
                       SE=sumdat[,paste0("beta_se_",subtype)])
    fwrite(cbind.data.frame(sumdat0,sumdat1),file=paste0("../result/PRS_subtype/",pop,"/",pop,"_",subtypes[i],"_CTSLEB_sumstats.txt"),
           row.names = F,sep="\t",quote=F)
  }
}
form_CTSLEB_sumdat(sumdat=eursumstat,pop="euro",bimprefix="../result/PRS1/euro_onco_tuning")


pop="asian"
easallsumstat=as.data.frame(fread("../result/PRS1/asian_training_sumstats.txt"))
easallctslebsumstat=as.data.frame(fread("../result/PRS1/ctsleb/EAS/EAS_CTSLEB.sumdat"))

load(paste0("../result/PRS_subtype/",pop,"/subtype_sumstats.RData"))
table(allres$SNP %in% easallsumstat$SNP)
table(easallsumstat$SNP %in% allres$SNP)
idx=match(easallsumstat$SNP,allres$SNP)
eassumstat=allres[idx,]
tmp=unlist(strsplit(eassumstat$SNP,":"))
alt=tmp[seq(3,length(tmp),4)]
ref=tmp[seq(4,length(tmp),4)]
all(eassumstat$a1==alt) #T
table(eassumstat$a1==easallsumstat$a1)
table(easallctslebsumstat$A1==alt)

idx=match(eassumstat$SNP,easallsumstat$SNP)
eassumstat$rsid=easallsumstat$rsid[idx]  
form_PRScsxsumstat(pop="asian",sumstat=eassumstat)
form_CTSLEB_sumdat(sumdat=eassumstat,pop="asian",bimprefix="../result/PRS1/asian_onco_tuning")
  