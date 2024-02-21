#!/usr/bin/env Rscript

setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")

library(data.table)
library(pROC)
library(boot)

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
phenoicogs=phenoicogs[,c("ID",paste0("pc",1:10),c("Behaviour1","ER_status1","PR_status1","HER2_status1","Grade1"))]

phenoonco=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
phenoonco$ID=phenoonco$Onc_ID
colnames(phenoonco)=gsub("PC_","pc",colnames(phenoonco))
phenoonco=phenoonco[,c("ID",paste0("pc",1:10),c("Behaviour1","ER_status1","PR_status1","HER2_status1","Grade1"))]

addsubtype=function(pheno=phenoicogs)
{
  idx=which(pheno==888,arr.ind = T)
  pheno[idx]=NA
  pheno$lumA=NA
  pheno$lumA[which(is.na(pheno$Behaviour1))]=0
  pheno1=pheno[,c("ER_status1","PR_status1","HER2_status1","Grade1")]

  idx.1 <- which((pheno1[,1]==1|pheno1[,2]==1)
                 &pheno1[,3]==0
                 &(pheno1[,4]==1|pheno1[,4]==2))
  pheno$lumA[idx.1]=1
  
  pheno$lumB=NA
  pheno$lumB[which(is.na(pheno$Behaviour1))]=0
  idx.2 <- which((pheno1[,1]==1|pheno1[,2]==1)
                 &pheno1[,3]==1)
  pheno$lumB[idx.2]=1
  
  #for Luminal B HER2 negative-like
  pheno$LumB_HN=NA
  pheno$LumB_HN[which(is.na(pheno$Behaviour1))]=0
  idx.3 <- which((pheno1[,1]==1|pheno1[,2]==1)
                 &pheno1[,3]==0
                 &pheno1[,4]==3)
  pheno$LumB_HN[idx.3]=1
  
  #for HER2 enriched-like
  pheno$HER2E=NA
  pheno$HER2E[which(is.na(pheno$Behaviour1))]=0
  idx.4 <- which(pheno1[,1]==0&pheno1[,2]==0
                 &pheno1[,3]==1)
  pheno$HER2E[idx.4]=1
  #for Triple negative
  pheno$TripN=NA
  pheno$TripN[which(is.na(pheno$Behaviour1))]=0
  idx.5 <- which(pheno1[,1]==0&pheno1[,2]==0
                 &pheno1[,3]==0)
  pheno$TripN[idx.5]=1
  return(pheno)
}
phenoicogs=addsubtype(pheno=phenoicogs)
phenoonco=addsubtype(pheno=phenoonco)

AUCBoot = function(data,indices){
  if (max(data$y,na.rm=T)!=1) data$y=data$y-1 #y:1,2-->0,1
  boot_data = data[indices, ]
  model4 <- glm(y~prs, data=boot_data,family = "binomial")
  predicted4 <- predict(model4,boot_data, type="response")
  auc=as.numeric(pROC::auc(boot_data$y,predicted4,quiet=T))
  return(c(auc))
}
#load data
outprefix="../result/PRS1/euro/euro"
methodprefix="CT"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
allres1=list()
allres1$allprs=list()
allres1$allfamval=list()
#combine 2 AF
allres1$allprs[[1]]=rbind(allres$allprs[[1]],allres$allprs[[2]])
allres1$allfamval[[1]]=rbind(allres$allfamval[[1]],allres$allfamval[[2]])
allres1$allprs[[2]]=allres$allprs[[3]]
allres1$allfamval[[2]]=allres$allfamval[[3]]
allres1$allprs[[3]]=allres$allprs[[4]]
allres1$allfamval[[3]]=allres$allfamval[[4]]
allres1$allprs[[4]]=allres$allprs[[5]]
allres1$allfamval[[4]]=allres$allfamval[[5]]

aucres=data.frame(matrix(NA,nrow=4,ncol=5))
colnames(aucres)=colnames(phenoicogs)[17:21]
rownames(aucres)=c("AFR","ASA","EUR","HIS")
aucreslow=aucreshigh=aucres
for (i in 1:4)
{
  famval=allres1$allfamval[[i]]
  prs=allres1$allprs[[i]]
  if (!"numeric" %in% class(prs)) #matrix
  {
    prs=prs[,ncol(prs)]
  }
  
  pheno=phenoonco
  if (i==1) pheno=rbind(phenoicogs,phenoonco)
  idx=match(famval$V1,pheno$ID)
  for (j in 1:5)
  {
    pheno.prs=data.frame(y=pheno[idx,16+j],prs=prs)
    model1 <- glm(I(y==1)~prs, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    aucres[i,j]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    boot_auc = boot(data =pheno.prs, statistic = AUCBoot, R = 1000)
    tmp=boot.ci(boot_auc,type="perc")
    aucreslow[i,j]=tmp$percent[4] #low
    aucreshigh[i,j]=tmp$percent[5] #high
  }
}
allaucres=list(aucres=aucres,aucreslow=aucreslow,aucreshigh=aucreshigh)

get_subtypeauc=function(outprefix="../result/PRS1/Weighted",methodprefix="CT")
{
  load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
  allres1=list()
  allres1$allprs=list()
  allres1$allfamval=list()
  #combine 2 AF, change prs to array if need
  if ("numeric" %in% class(allres$allprs[[1]])) #array
  {
    allres1$allprs[[1]]=c(allres$allprs[[1]],allres$allprs[[2]])
    allres1$allprs[[2]]=allres$allprs[[3]]
    allres1$allprs[[3]]=allres$allprs[[4]]
    allres1$allprs[[4]]=allres$allprs[[5]]
    
  }else
  {
    allres1$allprs[[1]]=rbind(allres$allprs[[1]],allres$allprs[[2]])
    allres1$allprs[[1]]=allres1$allprs[[1]][,ncol(allres$allprs[[1]])]
    allres1$allprs[[2]]=allres$allprs[[3]]
    allres1$allprs[[2]]=allres1$allprs[[2]][,ncol(allres$allprs[[1]])]
    allres1$allprs[[3]]=allres$allprs[[4]]
    allres1$allprs[[3]]=allres1$allprs[[3]][,ncol(allres$allprs[[1]])]
    allres1$allprs[[4]]=allres$allprs[[5]]
    allres1$allprs[[4]]=allres1$allprs[[4]][,ncol(allres$allprs[[1]])]
  }
  allres1$allfamval[[1]]=rbind(allres$allfamval[[1]],allres$allfamval[[2]])
  allres1$allfamval[[2]]=allres$allfamval[[3]]
  allres1$allfamval[[3]]=allres$allfamval[[4]]
  allres1$allfamval[[4]]=allres$allfamval[[5]]
  
  aucres=data.frame(matrix(NA,nrow=4,ncol=5))
  colnames(aucres)=colnames(phenoicogs)[17:21]
  rownames(aucres)=c("AFR","ASN","EUR","HIS")
  aucreslow=aucreshigh=aucres
  for (i in 1:4)
  {
    famval=allres1$allfamval[[i]]
    prs=allres1$allprs[[i]]
    if (!"numeric" %in% class(prs)) #matrix
    {
      prs=prs[,ncol(prs)]
    }
    pheno=phenoonco
    if (i==1) pheno=rbind(phenoicogs,phenoonco)
    idx=match(famval$V1,pheno$ID)
    for (j in 1:5)
    {
      pheno.prs=data.frame(y=pheno[idx,16+j],prs=prs)
      model1 <- glm(I(y==1)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      aucres[i,j]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
      boot_auc = boot(data =pheno.prs, statistic = AUCBoot, R = 1000)
      tmp=boot.ci(boot_auc,type="perc")
      aucreslow[i,j]=tmp$percent[4] #low
      aucreshigh[i,j]=tmp$percent[5] #high
    }
  }
  allaucres=list(aucres=aucres,aucreslow=aucreslow,aucreshigh=aucreshigh)
  save(allaucres,file=paste0(outprefix,"_",methodprefix,"_valauc_onsubtypes.RData"))
}

get_subtypeauc(outprefix="../result/PRS1/Weighted",methodprefix="LDpred")
get_subtypeauc(outprefix="../result/PRS1/TargetEUR",methodprefix="PRScsx")
get_subtypeauc(outprefix="../result/PRS1/TargetEAS",methodprefix="PRScsx")

get_auctable=function(outprefix="../result/PRS1/Weighted",methodprefix="CT",optwrite=T)
{
  load(paste0(outprefix,"_",methodprefix,"_valauc_onsubtypes.RData"))
  allaucres$aucres=round(allaucres$aucres,3)
  allaucres$aucreslow=round(allaucres$aucreslow,3)
  allaucres$aucreshigh=round(allaucres$aucreshigh,3)
  alltable=allaucres$aucres
  for (i in 1:nrow(alltable))
  {
    for (j in 1:ncol(alltable))
    {
      alltable[i,j]=paste0(alltable[i,j],"(",allaucres$aucreslow[i,j],",",allaucres$aucreshigh[i,j],")")
    }
  }
  if (optwrite==T)
  {
    write.csv(alltable,file=paste0(outprefix,"_",methodprefix,"_valauc_onsubtypes.csv"))
  }else
  {
    return(alltable)
  }
  
}
get_auctable(outprefix="../result/PRS1/Weighted",methodprefix="LDpred")
get_auctable(outprefix="../result/PRS1/TargetEUR",methodprefix="PRScsx")
get_auctable(outprefix="../result/PRS1/TargetEAS",methodprefix="PRScsx")

#use EUR and EAS tables to get only one (most populations use EUR results; EAS use EAS )
get_2auctable=function(outprefix1="../result/PRS1/TargetEUR",
                       outprefix2="../result/PRS1/TargetEAS",
                       outprefix="../result/PRS1/Target2", #output
                       methodprefix="PRScsx")
{
  alltable1=get_auctable(outprefix=outprefix1,methodprefix=methodprefix,optwrite=F)
  alltable2=get_auctable(outprefix=outprefix2,methodprefix=methodprefix,optwrite=F)
  idx=which(rownames(alltable1)=="ASN")
  alltable=alltable1
  alltable[idx,]=alltable2[idx,]
  
  if (optwrite==T)
  {
    write.csv(alltable,file=paste0(outprefix,"_",methodprefix,"_valauc_onsubtypes.csv"))
  }else
  {
    return(alltable)
  }
}
