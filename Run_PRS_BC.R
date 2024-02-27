#!/usr/bin/env Rscript
#subtype=NULL: canceroutcome; nonNULL:subtype
.libPaths(c("/data/wangx53",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")

set.seed(1000)
library(data.table)
library(pROC)
#library(SuperLearner)
library(RISCA)
library(boot)
library(tidyverse)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
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

subtypes=c("LumA","LumB","LumB_HN","Her2E","TripN")
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

#validation plink files
prefix_val1="../result/PRS1/african_onco_validation"
prefix_val2="../result/PRS1/african_icogs_validation"
prefix_val3="../result/PRS1/asian_onco_validation"
prefix_val4="../result/PRS1/euro_onco_validation"
prefix_val5="../result/PRS1/hispanic_onco_validation"
validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)

get_pheno.prs=function(prs=prs$V6,famval,subtype=NULL)
{
  phenotype=phenoicogs
  if (!all(famval$V2 %in% phenoicogs$ID))
  {
    phenotype=phenoonco
  }
  phenotype=phenotype[match(famval$V2,phenotype$ID),]
  if (!is.null(subtype))
  {
    phenotype$y=NA
    phenotype$y[which(phenotype[,subtype]==1)]=1
    phenotype$y[which(is.na(phenotype$Behaviour1))]=0
  }else
  {
    phenotype$y=famval$V6-1 #cancer
  }
  pheno.prs=phenotype
  pheno.prs$prs=prs
  pheno.prs=pheno.prs[,c("ID",paste0("pc",1:10),"age","Behaviour1","y","prs")]
  return(pheno.prs)
}
myscale=function(pheno.prs)
{
  idx=which(is.na(pheno.prs$Behaviour1)) #control
  pheno.prs$prs=(pheno.prs$prs-mean(pheno.prs$prs[idx]))/sd(pheno.prs$prs[idx])
  return(pheno.prs)
}

RISCA_AUC=function(pheno.prs)
{
  set.seed(1000)
  pheno.prs=pheno.prs[!is.na(pheno.prs$y),]
  pheno.prs=myscale(pheno.prs)
  # roc_obj_pv = roc.binary(status = "y", estimator = "pv",
  #                         variable = "prs",
  #                         confounders = ~EV1+EV2+EV3+EV4+EV5+EV6+
  #                           EV7+EV8+EV9+EV10+AGE,
  #                         data = pheno.prs,
  #                         precision=seq(0.05,0.95, by=0.05))
  confounders = c(paste0("pc",1:10),"age")
  idx=match(confounders,colnames(pheno.prs))
  idx=complete.cases(pheno.prs[,idx])
  pheno.prs=pheno.prs[idx,]
  roc_obj_ipw = roc.binary(status = "y", estimator = "ipw",
                           variable = "prs",
                           #confounders = ~1,
                           confounders = ~pc1+pc2+pc3+pc4+pc5+pc6+
                             pc7+pc8+pc9+pc10+age,
                           data = pheno.prs,
                           precision=seq(0.05,0.95, by=0.05))
  #print(roc_obj_pv$auc)
  #print(roc_obj_ipw$auc)
  return(roc_obj_ipw$auc)
}

AUCBoot = function(data,indices){
  if (max(data$y,na.rm = T)!=1) data$y=data$y-1 #y:1,2-->0,1
  boot_data = data[indices, ]
  model4 <- glm(y~prs, data=boot_data,family = "binomial")
  predicted4 <- predict(model4,boot_data, type="response")
  auc=as.numeric(pROC::auc(boot_data$y,predicted4,quiet=T))
  return(c(auc))
}

AUCadjBoot = function(data,indices){
  data=data[!is.na(data$y),]
  boot_data = data[indices, ]
  confounders = c(paste0("pc",1:10),"age")
  idx=match(confounders,colnames(data))
  idx=complete.cases(boot_data[,idx])
  boot_data=boot_data[idx,]
  roc_obj_ipw = roc.binary(status = "y", estimator = "ipw",
                           variable = "prs",
                           confounders = ~pc1+pc2+pc3+pc4+pc5+pc6+
                             pc7+pc8+pc9+pc10+age,
                           data = boot_data,
                           precision=seq(0.05,0.95, by=0.05))
  auc=roc_obj_ipw$auc
  return(c(auc))
}

#allprs is a list of 5 val PRS dataframe
#allfamval is a list of 5 val fam dataframe
# get_valauc=function(allprs,allfamval,outprefix,methodprefix="CT",subtype=NULL)
# {
#   #validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
#   #validationscores=c(score_val1,score_val2,score_val3,score_val4,score_val5)
#   
#   CTauc_val=data.frame(matrix(NA,nrow=6,ncol=6))
#   colnames(CTauc_val)=c("onco_african","icogs_african","onco_asian","onco_euro","onco_hispanic","african")
#   rownames(CTauc_val)=c("AUC","AUClow","AUChigh","AUCadj","AUCadjlow","AUCadjhigh")
#   afrpheno.prs=NULL
#   for (i in 1:5)
#   {
#     prs=allprs[[i]]
#     famval=allfamval[[i]]
#     if ("numeric" %in% class(prs))
#     {
#       pheno.prs=get_pheno.prs(prs=prs,famval,subtype = subtype)  #prs should align with famval
#     }else
#     {
#       all(prs[,1]==famval[,1])
#       pheno.prs=get_pheno.prs(prs=prs[,ncol(prs)],famval,subtype = subtype) #prs should align with famval
#     }
#     
#     if (i %in% c(1,2))
#     {
#       afrpheno.prs=rbind(afrpheno.prs,pheno.prs)
#     }
#     if (length(unique(pheno.prs$y[!is.na(pheno.prs$y)]))>1) #for Her2E,icogsafrican, there is no cases
#     {
#       model1 <- glm(I(y==1)~prs, data=pheno.prs,family = "binomial")
#       predicted1 <- predict(model1,pheno.prs, type="response")
#       CTauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
#       boot_auc = boot(data =pheno.prs, statistic = AUCBoot, R = 1000)
#       tmp=boot.ci(boot_auc,type="perc")
#       CTauc_val[2,i]=tmp$percent[4] #low
#       CTauc_val[3,i]=tmp$percent[5] #high
#       CTauc_val[4,i]=RISCA_AUC(pheno.prs)
#       boot_aucadj = boot(data =pheno.prs, statistic = AUCadjBoot, R = 1000)
#       tmp=boot.ci(boot_aucadj,type="perc")
#       CTauc_val[5,i]=tmp$percent[4] #adjlow
#       CTauc_val[6,i]=tmp$percent[5] #adjhigh
#     }
#     
#   }
#   
#   model1 <- glm(I(y==1)~prs, data=afrpheno.prs,family = "binomial")
#   predicted1 <- predict(model1,afrpheno.prs, type="response")
#   CTauc_val[1,6]=as.numeric(pROC::auc(afrpheno.prs$y,predicted1,quiet=T))
#   boot_auc = boot(data =afrpheno.prs, statistic = AUCBoot, R = 1000)
#   tmp=boot.ci(boot_auc,type="perc")
#   CTauc_val[2,6]=tmp$percent[4] #low
#   CTauc_val[3,6]=tmp$percent[5] #high
#   CTauc_val[4,6]=RISCA_AUC(afrpheno.prs)
#   boot_aucadj = boot(data =afrpheno.prs, statistic = AUCadjBoot, R = 1000)
#   tmp=boot.ci(boot_aucadj,type="perc")
#   CTauc_val[5,6]=tmp$percent[4] #adjlow
#   CTauc_val[6,6]=tmp$percent[5] #adjhigh
#   write.table(CTauc_val,file=paste0(outprefix,"_",methodprefix,"_valauc.txt"),row.names=F,sep="\t",quote=F)
#   allres=list(auc=CTauc_val,allprs=allprs,allfamval=allfamval)
#   save(allres,file=paste0(outprefix,"_",methodprefix,"_valauc.RData"))
#   return(allres)
# }

#combine icogs/onco african
#allprs is a list of 4 val PRS dataframe
#allfamval is a list of 4 val fam dataframe
get_valauc=function(allprs,allfamval,outprefix,methodprefix="CT",subtype=NULL)
{
  #validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  #validationscores=c(score_val1,score_val2,score_val3,score_val4,score_val5)
  if (length(allprs)==5)
  {
    allprs1=list()
    if (class(allprs[[1]])!="numeric")
    {
      allprs1[[1]]=rbind(allprs[[1]],allprs[[2]])
    }else
    {
      allprs1[[1]]=c(allprs[[1]],allprs[[2]])
    }
    allprs1[[2]]=allprs[[3]]
    allprs1[[3]]=allprs[[4]]
    allprs1[[4]]=allprs[[5]]
    allfamval1=list()
    allfamval1[[1]]=rbind(allfamval[[1]],allfamval[[2]])
    allfamval1[[2]]=allfamval[[3]]
    allfamval1[[3]]=allfamval[[4]]
    allfamval1[[4]]=allfamval[[5]]
    allprs=allprs1
    allfamval=allfamval1
  }
  CTauc_val=data.frame(matrix(NA,nrow=6,ncol=4))
  colnames(CTauc_val)=c("african","onco_asian","onco_euro","onco_hispanic")
  rownames(CTauc_val)=c("AUC","AUClow","AUChigh","AUCadj","AUCadjlow","AUCadjhigh")

  for (i in 1:4)
  {
    prs=allprs[[i]]
    famval=allfamval[[i]]
    if ("numeric" %in% class(prs))
    {
      pheno.prs=get_pheno.prs(prs=prs,famval,subtype = subtype)  #prs should align with famval
    }else
    {
      all(prs[,1]==famval[,1])
      pheno.prs=get_pheno.prs(prs=prs[,ncol(prs)],famval,subtype = subtype) #prs should align with famval
    }
    
    if (length(unique(pheno.prs$y[!is.na(pheno.prs$y)]))>1)
    {
      model1 <- glm(I(y==1)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      CTauc_val[1,i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
      boot_auc = boot(data =pheno.prs, statistic = AUCBoot, R = 1000)
      tmp=boot.ci(boot_auc,type="perc")
      CTauc_val[2,i]=tmp$percent[4] #low
      CTauc_val[3,i]=tmp$percent[5] #high
      CTauc_val[4,i]=RISCA_AUC(pheno.prs)
      boot_aucadj = boot(data =pheno.prs, statistic = AUCadjBoot, R = 1000)
      tmp=boot.ci(boot_aucadj,type="perc")
      CTauc_val[5,i]=tmp$percent[4] #adjlow
      CTauc_val[6,i]=tmp$percent[5] #adjhigh
    }
  }
  
  write.table(CTauc_val,file=paste0(outprefix,"_",methodprefix,"_valauc.txt"),row.names=F,sep="\t",quote=F)
  allres=list(auc=CTauc_val,allprs=allprs,allfamval=allfamval)
  save(allres,file=paste0(outprefix,"_",methodprefix,"_valauc.RData"))
  return(allres)
}

# #for EUR
# metafile="../result/PRS1/euro_training_sumstats.txt"
# metadat=as.data.frame(fread(metafile))

#CT
#euro_training_sumstats.txt came from QC_logistic_PRS.R
runCT=function(sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro")
{
  #CT method
  #parameters for clumping
  pthr=1
  r2thr=0.1
  kbpthr=500
  
  # cmd=paste0(plink," --bfile ",prefix_tun," --clump ",sumstatfile," --clump-p1 ",
  #            pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field rsid --clump-field p --out ",outprefix," --memory 128000")
  # system(cmd)
  clumpsnp=read.table(paste0(outprefix,".clumped"),header=T)
  # write.table(clumpsnp$SNP,file=paste0(outprefix,".clumpedsnp"),row.names=F,col.names = F,quote=F)
  # sumstat=as.data.frame(fread(sumstatfile))
  # tmp=data.frame(SNP=sumstat$rsid,A1=sumstat$a1,beta=sumstat$beta)
  # write.table(tmp,file=paste0(outprefix,".score"),row.names=F,col.names=T,sep=" ",quote=F)
  # tmp=data.frame(SNP=sumstat$rsid,P=sumstat$p)
  # write.table(tmp,file=paste0(outprefix,".pvalue"),row.names=F,col.names=T,sep=" ",quote=F)
  # cmd=paste0(plink2," --bfile ",prefix_tun," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",outprefix,"_tun"," --memory 128000")
  # system(cmd)
  prefix_val1="../result/PRS1/african_onco_validation"
  score_val1=paste0(outprefix,"_validation_onco_african")
  # cmd=paste0(plink2," --bfile ",prefix_val1," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",score_val1," --memory 128000")
  # system(cmd)
  prefix_val2="../result/PRS1/african_icogs_validation"
  score_val2=paste0(outprefix,"_validation_icogs_african")
  # cmd=paste0(plink2," --bfile ",prefix_val2," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",score_val2," --memory 128000")
  # system(cmd)
  prefix_val3="../result/PRS1/asian_onco_validation"
  score_val3=paste0(outprefix,"_validation_onco_asian")
  # cmd=paste0(plink2," --bfile ",prefix_val3," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",score_val3," --memory 128000")
  # system(cmd)
  prefix_val4="../result/PRS1/euro_onco_validation"
  score_val4=paste0(outprefix,"_validation_onco_euro")
  # cmd=paste0(plink2," --bfile ",prefix_val4," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",score_val4," --memory 128000")
  # system(cmd)
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  score_val5=paste0(outprefix,"_validation_onco_hispanic")
  # cmd=paste0(plink2," --bfile ",prefix_val5," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",score_val5," --memory 128000")
  # system(cmd)
  
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-05,5E-03,5E-02,5E-01,1) 
  famtun=read.table(paste0(prefix_tun,".fam"))
  auc_tun=rep(0,length(pthres))
  for (i in 1:length(pthres))
  {
    prs=read.table(paste0(outprefix,"_tun.p_value_",i,".sscore"))
    if (any(!is.na(prs$V6)))
    {
      #all(prs$V1==famtun$V1)
      pheno.prs=data.frame(y=famtun$V6,prs=prs$V6)
      model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
      predicted1 <- predict(model1,pheno.prs, type="response")
      auc_tun[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    }
  }
  
  idx_optimal=which.max(auc_tun) #euro:3,asian 5
  pvalue_opitmal=pthres[idx_optimal]
  print(paste0("pvalue_optimal: ",pvalue_opitmal))
  print(paste0("number of snps: ",sum(clumpsnp$P<pvalue_opitmal)))
  
  # #add tuning PRS on the other tuning data(for example tuning eur PRS on tuning asian data),used for weighted PRS
  # other_prefix_tun="../result/PRS1/euro_onco_tuning"
  # other_score_tun="../result/PRS1/asian/asian_tun_euro_tun" #asian PRS on euro onco
  # if (prefix_tun=="../result/PRS1/euro_onco_tuning")
  # {
  #   other_prefix_tun="../result/PRS1/asian_onco_tuning"
  #   other_score_tun="../result/PRS1/euro/euro_tun_asian_tun"
  # }
  # cmd=paste0(plink2," --bfile ",other_prefix_tun," --score ",outprefix,".score ",
  #            "--q-score-range /data/BB_Bioinformatics/Kevin/PRS_EASLC/result/range_list ",outprefix,".pvalue --extract ",outprefix,".clumpedsnp ",
  #            "--out ",other_score_tun," --memory 128000")
  # system(cmd)
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  validationscores=c(score_val1,score_val2,score_val3,score_val4,score_val5)
  allprs=allfamval=list()
  for (i in 1:5)
  {
    prs=read.table(paste0(validationscores[i],".p_value_",idx_optimal,".sscore"))
    allprs[[i]]=prs
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
  }
  CTauc_val=get_valauc(allprs,allfamval,outprefix,methodprefix="CT")
  return(CTauc_val)
}
print(Sys.time())
euroCT=runCT()
print(Sys.time())
#euro
#5e-6, 387 SNPs
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic african
# AUC         0.5597353     0.5522557    0.5973346 6.274929e-01     0.5676048 0.5574493
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000

print(Sys.time())
asianCT=runCT(sumstatfile="../result/PRS1/asian_training_sumstats.txt",prefix_tun="../result/PRS1/asian_onco_tuning",outprefix="../result/PRS1/asian/asian")
print(Sys.time())
#asian
#5e-5, 127 SNPs
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic afircan
# AUC         0.5233717     0.5236688    0.5600331 5.339368e-01     0.5219497 0.524
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000

Weighted_CTprs=function()
{
  famtun1=as.data.frame(fread("../result/PRS1/euro_onco_tuning.fam"))
  famtun2=as.data.frame(fread("../result/PRS1/asian_onco_tuning.fam"))
  famtun=rbind(famtun1,famtun2)
  idx_optimal1=3
  idx_optimal2=5
  prs1_tun1=read.table(paste0("../result/PRS1/euro/euro_tun.p_value_",idx_optimal1,".sscore"))
  prs1_tun2=read.table(paste0("../result/PRS1/euro/euro_tun_asian_tun.p_value_",idx_optimal1,".sscore"))
  prs1_tun=rbind(prs1_tun1,prs1_tun2) #euro first,asian second
  prs2_tun1=read.table(paste0("../result/PRS1/asian/asian_tun.p_value_",idx_optimal2,".sscore"))
  prs2_tun2=read.table(paste0("../result/PRS1/asian/asian_tun_euro_tun.p_value_",idx_optimal2,".sscore"))
  prs2_tun=rbind(prs2_tun2,prs2_tun1) #euro first,asian second
  all(prs1_tun$V2==famtun$V2)
  all(prs2_tun$V2==famtun$V2)
  dat_tun=cbind.data.frame(y=famtun$V6-1,prs1=prs1_tun$V6,prs2=prs2_tun$V6)
  model1 <- glm(y~prs1+prs2, data=dat_tun,family = "binomial")
  
  #print(summary(model1))
  w1=summary(model1)$coefficients[2,1]
  w2=summary(model1)$coefficients[3,1]
  
  #apply on validations
  outprefix1="../result/PRS1/euro/euro"
  score_val1_1=paste0(outprefix1,"_validation_onco_african")
  score_val2_1=paste0(outprefix1,"_validation_icogs_african")
  score_val3_1=paste0(outprefix1,"_validation_onco_asian")
  score_val4_1=paste0(outprefix1,"_validation_onco_euro")
  score_val5_1=paste0(outprefix1,"_validation_onco_hispanic")
  outprefix2="../result/PRS1/asian/asian"
  score_val1_2=paste0(outprefix2,"_validation_onco_african")
  score_val2_2=paste0(outprefix2,"_validation_icogs_african")
  score_val3_2=paste0(outprefix2,"_validation_onco_asian")
  score_val4_2=paste0(outprefix2,"_validation_onco_euro")
  score_val5_2=paste0(outprefix2,"_validation_onco_hispanic")
  
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  validationscores1=c(score_val1_1,score_val2_1,score_val3_1,score_val4_1,score_val5_1) #euro prs
  validationscores2=c(score_val1_2,score_val2_2,score_val3_2,score_val4_2,score_val5_2) #asian prs
  
  allprs=allfamval=list()
  for (i in 1:5)
  {
    prs1=read.table(paste0(validationscores1[i],".p_value_",idx_optimal1,".sscore"))
    prs2=read.table(paste0(validationscores2[i],".p_value_",idx_optimal2,".sscore"))
    prs=prs1$V6*w1+prs2$V6*w2
    allprs[[i]]=prs
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
    all(prs1$V1==famval$V1)
    all(prs2$V1==famval$V1)
  }
  Weighted_CTauc_val=get_valauc(allprs,allfamval,outprefix="../result/PRS1/Weighted",methodprefix="CT")
  
  return(Weighted_CTauc_val)
}
WeightedCT=Weighted_CTprs()
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic   african
# AUC         0.5590266     0.5518596    0.5967572 6.269937e-01     0.5671613 0.5568423


#LDpred2
get_PRS_genotype=function(betafile="../result/PRS1/euro/euro_LDpred_weights.txt",betaprefix="../result/PRS1/euro_ldpredweight",prefix="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro")
{
  ldpredprefix=paste0(outprefix,"_",basename(prefix),"_ldpred")
  scorefile=paste0(ldpredprefix,".score")
  if (!file.exists(scorefile)) #if file exist, no need to generate again
  {
    beta=as.data.frame(fread(betafile))
    bim=as.data.frame(fread(paste0(prefix,".bim")))
    comsnps=intersect(bim$V2,beta$snp)
    idx1=match(comsnps,bim$V2)
    idx2=match(comsnps,beta$snp)
    bim=bim[idx1,]
    beta=beta[idx2,]
    idx1=which(beta$effect==bim$V5)
    idx2=which(beta$effect==bim$V6)
    a2=bim$V6
    a2[idx2]=bim$V5[idx2]
    tmp=data.frame(snp=beta$snp,a2=a2)
    a2file=paste0(betaprefix,".a2")
    write.table(tmp,file=a2file,col.names = F,row.names = F,sep="\t",quote=F)
    snpfile=paste0(betaprefix,".snp")
    write.table(tmp$snp,file=snpfile,col.names = F,row.names = F,sep="\t",quote=F)
    
    cmd=paste0(plink2," --bfile ",prefix, " --ref-allele ",a2file," --extract ",snpfile," --make-bed --out ",ldpredprefix)
    system(cmd)
    cmd=paste0(plink," --bfile ",ldpredprefix," --allow-no-sex --score ",betafile," header sum --out ",ldpredprefix," --memory 64000 --threads 8")
    system(cmd)
    PRS=read.table(paste0(ldpredprefix,".profile"),header=T)
    system(paste0("rm ",ldpredprefix,".bed"))
    system(paste0("rm ",ldpredprefix,".bim"))
    system(paste0("rm ",ldpredprefix,".fam"))
    
    write.table(PRS,file=scorefile,row.names = F,sep="\t",quote=F)
  }else
  {
    PRS=read.table(scorefile,header=T)
  }
  
  return(PRS)
}
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# 1. Read in the phenotype and covariate files
# read in phenotype and covariates
library(data.table)
library(magrittr)
run_ldpred2=function(tmpdir0="euro",sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro",opfiltertuning=T)
{
  print(sumstatfile)
  print(prefix_tun)
  #print(prefix_val)
  #2. Obtain HapMap3 SNPs
  #LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs
  #load HapMap3 SNPs
  # info =readRDS("/data/BB_Bioinformatics/Kevin/tools/ldpred2/map_hm3_ldpred2.rds")
  #write.table(info$rsid,file="../result/PRS1/ldpred2_hw3snp.txt",row.names = F,col.names = F,quote=F)
  #3. Load and transform the summary statistic file
  #Load summary statistic file
  # Read in the summary statistic file
  # sumdat=as.data.frame(fread(sumstatfile)) #8242741,6753434
  # print("sumdat dim:")
  # print(dim(sumdat))
  # # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  # sumdat1=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=sumdat$a0,a1=sumdat$a1,n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  # sumdat1=sumdat1[sumdat1$beta!=0,]
  # fwrite(sumdat1,file=paste0(sumstatfile,".ldpred"),row.names=F,sep="\t")
  # sumstats <- bigreadr::fread2(paste0(sumstatfile,".ldpred")) 
  # 
  # # Filter out hapmap SNPs
  # sumstats <- sumstats[sumstats$rsid%in% info$rsid,] #1012780,811530
  # write.table(sumstats$rsid,file=paste0(outprefix,"_sumstats_hw3snp.txt"),row.names = F,col.names = F,quote=F)
  # print("sumstat in HM3:")
  # print(nrow(sumstats))
  # #3. Calculate the LD matrix
  # # Get maximum amount of cores
  # NCORES <- nb_cores()
  # # Open a temporary file
  # tmpdir <- tempfile(tmpdir = paste0("tmp_data_",tmpdir0))
  # #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  # 
  # #to remove mono snps in tuning data (which result NAs in computing corr)
  # if (opfiltertuning)
  # {
  #   cmd=paste0(plink," --bfile ",prefix_tun," --maf 0.01 --extract ../result/PRS1/ldpred2_hw3snp.txt"," --make-bed --out ",prefix_tun,"_ldpred_maf01 --memory 128000 --threads 8")
  #   system(cmd)
  #   cmd=paste0(plink," --bfile ",prefix_tun,"_ldpred_maf01 --extract ",outprefix,"_sumstats_hw3snp.txt --make-bed --out ",prefix_tun,"_ldpred_maf01 --memory 128000 --threads 8")
  #   system(cmd)
  # }
  # 
  # # preprocess the bed file (only need to do once for each data set)
  # if (opfiltertuning)
  # {
  #   file.remove(paste0(prefix_tun,"_ldpred_maf01.bk"))
  #   #if (!file.exists(paste0(prefix_tun,"_ldpred_maf01.rds")))
  #   snp_readBed(paste0(prefix_tun,"_ldpred_maf01.bed"))
  #   # now attach the genotype object
  #   obj.bigSNP <- snp_attach(paste0(prefix_tun,"_ldpred_maf01.rds"))
  # }else
  # {
  #   #if (!file.exists(paste0(prefix_tun,".rds")))
  #   snp_readBed(paste0(prefix_tun,".bed"))
  #   obj.bigSNP <- snp_attach(paste0(prefix_tun,".rds"))
  # }
  # 
  # # extract the SNP information from the genotype
  # map <- obj.bigSNP$map[-3]
  # names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # table(map$rsid %in% sumstats$rsid)
  # # perform SNP matching
  # info_snp <- snp_match(sumstats, map)
  # table(info_snp$rsid==map$rsid)
  # write.table(info_snp,file=paste0(outprefix,"_LDpred_","info_snp.txt"),row.names=F,sep="\t",quote=F)
  # # Assign the genotype to a variable for easier downstream analysis
  # genotype <- obj.bigSNP$genotypes
  # genotype1 = snp_fastImputeSimple(genotype)
  # # Rename the data structures
  # CHR <- map$chr
  # POS <- map$pos
  # # get the CM information from 1000 Genome
  # # will download the 1000G file to the current directory (".")
  # POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  # 
  # # Initialize variables for storing the LD score and LD matrix
  # corr <- NULL
  # ld <- NULL
  # 
  # # calculate LD
  # for (chr in 1:22) {
  #   print(chr)
  #   # Extract SNPs that are included in the chromosome
  #   ind.chr <- which(info_snp$chr == chr)
  #   ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  #   # Calculate the LD
  #   corr0 <- snp_cor(
  #     genotype,
  #     ind.col = ind.chr2,
  #     ncores = NCORES,
  #     infos.pos = POS2[ind.chr2],
  #     size = 3 / 1000
  #   )
  #   if (chr == 1) {
  #     ld <- Matrix::colSums(corr0^2)
  #     corr <- as_SFBM(corr0, tmpdir)
  #   } else {
  #     ld0=Matrix::colSums(corr0^2)
  #     ld <- c(ld, Matrix::colSums(corr0^2))
  #     corr$add_columns(corr0, nrow(corr))
  #   }
  #   if (sum(is.na(ld))>0) stop(chr)
  # }
  # save(ld,corr,file=paste0(outprefix,"_LDpred_","ld.RData"))
  # #4. Perform LD score regression
  # df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  # rownames(df_beta)=info_snp$rsid.ss
  # ldsc <- snp_ldsc(   ld,
  #                     length(ld),
  #                     chi2 = (df_beta$beta / df_beta$beta_se)^2,
  #                     sample_size = df_beta$n_eff,
  #                     blocks = NULL)
  # # ldsc <- snp_ldsc2(corr,df_beta)
  # h2_est <- ldsc[["h2"]] #asn:0.058
  # print(paste0("h2_est:",h2_est))
  # if (ldsc[['h2']] < 0) print('h2 negative')
  # 
  # 
  # #6 grid model
  # # Prepare data for grid model
  # p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  # h2_seq <- round(h2_est * c(0.3,0.7, 1, 1.4), 4)
  # grid.param <-
  #   expand.grid(p = p_seq,
  #               h2 = h2_seq,
  #               sparse = c(FALSE, TRUE))
  # # Get adjusted beta from grid model
  # set.seed(1000) # to get the same result every time
  # beta_grid <-
  #   snp_ldpred2_grid(corr, df_beta, grid.param, ncores = 1)
  # 
  # 
  # famtun=read.table(paste0(prefix_tun,".fam"))
  # pred_grid <- big_prodMat( genotype1, 
  #                           beta_grid, 
  #                           ind.col = info_snp$`_NUM_ID_`)
  # rownames(pred_grid)=famtun$V2
  # 
  # #to get AUC on tunning
  # auc_tun=rep(0,ncol(pred_grid))
  # 
  # for (i in 1:ncol(pred_grid))
  # {
  #   if (any(!is.na(pred_grid[,i])))
  #   {
  #     pheno.prs=cbind.data.frame(y=famtun$V6,prs=pred_grid[,i])
  #     
  #     model1 <- glm(I(y==2)~prs, data=pheno.prs,family = "binomial")
  #     predicted1 <- predict(model1,pheno.prs, type="response")
  #     auc_tun[i]= as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
  #   }
  # }
  # idx_optimal=which.max(auc_tun)
  # param_optimal=grid.param[idx_optimal,]
  # ldpred_beta=data.frame(snp=info_snp$rsid.ss,effect=info_snp$a1,beta=beta_grid[,idx_optimal])
  betafile=paste0(outprefix,"_LDpred_weights.txt")
  # write.table(ldpred_beta,file=betafile,row.names = F,sep="\t",quote=F)
  
  # #save PRS on tuning data
  # #euro:../result/PRS1/euro/euro_tuning_ldpredprs.txt
  # #asian:../result/PRS1/asian/asian_tuning_ldpredprs.txt
  # tun.prs=cbind.data.frame(famtun,prs=pred_grid[,idx_optimal])
  # write.table(tun.prs,file=paste0(outprefix,"_tuning_ldpredprs.txt"),row.names = F,sep="\t",quote=F)
  # save(ld,corr,beta_grid,pred_grid,ldsc,df_beta,grid.param,auc_tun,idx_optimal,param_optimal,file=paste0(outprefix,"LDpred_pred.RData"))
  # #add tuning PRS on the other tuning data(for example tuning eur PRS on tuning asian data),used for weighted PRS
  # other_prefix_tun="../result/PRS1/euro_onco_tuning" #asian PRS on euro onco
  # if (prefix_tun=="../result/PRS1/euro_onco_tuning")
  # {
  #   other_prefix_tun="../result/PRS1/asian_onco_tuning"
  # }
  # #euro prs on onco asian tuning: ../result/PRS1/euro/euro_asian_onco_tuning_ldpred.score
  # #asian prs on onco euro tuning: ../result/PRS1/asian/asian_euro_onco_tuning_ldpred.score
  # other_prs_tun=get_PRS_genotype(betafile=betafile,betaprefix=paste0(outprefix,"_ldpredweight"),prefix=other_prefix_tun,outprefix=outprefix)
  # 
  #validation
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  
  allprs=allfamval=list()
  for (i in 1:5)
  {
    prs=get_PRS_genotype(betafile=betafile,betaprefix=paste0(outprefix,"_ldpredweight"),prefix=validationprefix[i],outprefix=outprefix)
    allprs[[i]]=prs
    #prs=read.table(paste0(validationscores[i],".p_value_",idx_optimal,".sscore"))
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
  }
  LDpredauc_val=get_valauc(allprs,allfamval,outprefix=outprefix,methodprefix="LDpred")
  
  return(LDpredauc_val)
}
euroLDpred=run_ldpred2(tmpdir0="euro",sumstatfile="../result/PRS1/euro_training_sumstats.txt",prefix_tun="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/euro/euro",opfiltertuning=T)

#euro:
#         onco_african icogs_african   onco_asian    onco_euro onco_hispanic   african
# AUC         0.5680855     0.5579379    0.6117766 6.582231e-01     0.6067224 0.5628277
# totaln   5569.0000000  1758.0000000 5406.0000000 2.771000e+04  2413.0000000        NA
# casen    3481.0000000   941.0000000 2663.0000000 1.480900e+04  1195.0000000        NA
# controln 2088.0000000   817.0000000 2743.0000000 1.290100e+04  1218.0000000        NA

asianLDpred=run_ldpred2(tmpdir0="asian",sumstatfile="../result/PRS1/asian_training_sumstats.txt",prefix_tun="../result/PRS1/asian_onco_tuning",outprefix="../result/PRS1/asian/asian",opfiltertuning=T)
#asian
#          onco_african icogs_african   onco_asian   onco_euro onco_hispanic   african
# AUC         0.5505593     0.5346977    0.5874594 5.75309e-01     0.5638608 0.5518476
# totaln   5569.0000000  1758.0000000 5406.0000000 2.77100e+04  2413.0000000        NA
# casen    3481.0000000   941.0000000 2663.0000000 1.48090e+04  1195.0000000        NA
# controln 2088.0000000   817.0000000 2743.0000000 1.29010e+04  1218.0000000        NA


Weighted_ldpredprs=function()
{
  famtun1=as.data.frame(fread("../result/PRS1/euro_onco_tuning.fam"))
  famtun2=as.data.frame(fread("../result/PRS1/asian_onco_tuning.fam"))
  famtun=rbind(famtun1,famtun2)
  prs1_tun1=read.table("../result/PRS1/euro/euro_tuning_ldpredprs.txt",header=T)
  prs1_tun2=read.table("../result/PRS1/euro/euro_asian_onco_tuning_ldpred.score",header = T)
  prs1_tun=rbind(data.frame(IID=prs1_tun1$V2,prs=prs1_tun1$prs),data.frame(IID=prs1_tun2$IID,prs=prs1_tun2$SCORESUM)) #euro first,asian second
  prs2_tun1=read.table("../result/PRS1/asian/asian_tuning_ldpredprs.txt",header=T)
  prs2_tun2=read.table("../result/PRS1/asian/asian_euro_onco_tuning_ldpred.score",header = T)
  prs2_tun=rbind(data.frame(IID=prs2_tun2$IID,prs=prs2_tun2$SCORESUM),data.frame(IID=prs2_tun1$V2,prs=prs2_tun1$prs)) #euro first,asian second
  all(prs1_tun$IID==famtun$V2)
  all(prs2_tun$IID==famtun$V2)
  dat_tun=cbind.data.frame(y=famtun$V6-1,prs1=prs1_tun$prs,prs2=prs2_tun$prs)
  model1 <- glm(y~prs1+prs2, data=dat_tun,family = "binomial")
  
  #print(summary(model1))
  w1=summary(model1)$coefficients[2,1]
  w2=summary(model1)$coefficients[3,1]
  
  #apply on validations
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  #euro prs
  validationprefix1=paste0("../result/PRS1/euro/euro_",basename(validationprefix),"_ldpred")
  #asian prs
  validationprefix2=paste0("../result/PRS1/asian/asian_",basename(validationprefix),"_ldpred")
  
  allprs=allfamval=list()
  for (i in 1:5)
  {
    prs1=read.table(paste0(validationprefix1[i],".score"),header=T)
    prs2=read.table(paste0(validationprefix2[i],".score"),header=T)
    prs=prs1$SCORESUM*w1+prs2$SCORESUM*w2
    allprs[[i]]=prs
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
  }
  Weighted_LDpredauc_val=get_valauc(allprs,allfamval,outprefix="../result/PRS1/Weighted",methodprefix="LDpred")
  
  return(Weighted_LDpredauc_val)
}
WeightedLDpred=Weighted_ldpredprs()
#          onco_african icogs_african   onco_asian    onco_euro onco_hispanic   african
# AUC         0.5755086     0.5609777    0.6187005 6.530507e-01     0.6060151 0.5717066

#collect the above val auc results
outprefix="../result/PRS1/euro/euro"
methodprefix="CT"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
euroCT=allres
outprefix="../result/PRS1/asian/asian"
methodprefix="CT"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
asianCT=allres
outprefix="../result/PRS1/Weighted"
methodprefix="CT"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
WeightedCT=allres

outprefix="../result/PRS1/euro/euro"
methodprefix="LDpred"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
euroLDpred=allres
outprefix="../result/PRS1/asian/asian"
methodprefix="LDpred"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
asianLDpred=allres
outprefix="../result/PRS1/Weighted"
methodprefix="LDpred"
load(paste0(outprefix,"_",methodprefix,"_valauc.RData"))
WeightedLDpred=allres
save(euroCT,asianCT,WeightedCT,euroLDpred,asianLDpred,WeightedLDpred,file="../result/Run_PRS_BC.RData")

get_auctable=function(aucres=euroCT$auc)
{
  aucres=round(aucres,digits = 3)
  res=data.frame(matrix(NA,1,4))
  colnames(res)=c("AFR","ASN","EUR","HIS")
  res$AFR=paste0(aucres$african[1],"(",aucres$african[2],",",aucres$african[3],")")
  res$ASN=paste0(aucres$onco_asian[1],"(",aucres$onco_asian[2],",",aucres$onco_asian[3],")")
  res$EUR=paste0(aucres$onco_euro[1],"(",aucres$onco_euro[2],",",aucres$onco_euro[3],")")
  res$HIS=paste0(aucres$onco_hispanic[1],"(",aucres$onco_hispanic[2],",",aucres$onco_hispanic[3],")")
  return(res)
}
tmp=list(euroCT$auc,asianCT$auc,WeightedCT$auc,euroLDpred$auc,asianLDpred$auc,WeightedLDpred$auc)
allauctable=NULL
for (i in 1:length(tmp))
{
  tmp1=get_auctable(aucres=tmp[[i]])
  allauctable=rbind(allauctable,tmp1)
}
rownames(allauctable)=c("EURCT","ASNCT","WeightedCT","EURLDpred2","ASNLDpred2","WeightedLDpred2")
write.csv(allauctable,file="../result/validationAUC.csv")  

#PRSCSx
phi = c("phi1e+00","phi1e-02","phi1e-04","phi1e-06")

#get PRS based on tuning data
#targetEAS_EURprs_asian_onco_tuning, targetEAS,EUR PRS on asian_onco_tuning
# Asian target PRS: ../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_tuningPRS.sscore, ../result/PRS1/prscsx/targetEAS_EASprs_asian_onco_tuningPRS.sscore
# beta: ../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_.prs_coeff, ../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_.prs_coeff
get_PRScsxprs_tuning=function(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix="../result/PRS1/asian_onco_tuning",outprefix="../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_tuning")
{
  if (!file.exists(paste0(prefix,"PRS.sscore")))
  {
    #read posterior SNP effect size
    betadat=NULL
    for (chr in 1:22)
    {
      datlist=list()
      for (i in 1:4)
      {
        myphi=phi[i]
        myfile=paste0(prefix,myphi,"_chr",chr,".txt")
        if (file.exists(myfile))
        {
          datlist[[i]]=fread(myfile)
          
        }else
        {
          warning(paste0(myfile," is not existed"))
        }
      }
      #all(datlist[[1]]$V4==datlist[[4]]$V4)
      betadat_chr=cbind.data.frame(SNP=datlist[[1]]$V2,A1=datlist[[1]]$V4,phi1=datlist[[1]]$V6,phi2=datlist[[2]]$V6,phi4=datlist[[3]]$V6,phi6=datlist[[4]]$V6)
      betadat=rbind(betadat,betadat_chr)
    }
    betaoutfile=paste0(prefix,".prs_coeff")
    write.table(betadat,file=betaoutfile,row.names = F,sep="\t",quote=F)
    
    #use plink to get PRS
    cmd=paste0(plink2," --pfile ",bimprefix," --score-col-nums 3-",ncol(betadat),
               " --score ",betaoutfile," cols=+scoresums,-scoreavgs header no-mean-imputation ",
               " --out ",outprefix,"PRS"," --memory 64000 --threads 8")
    
    system(cmd)
  }
  
  prsdat=read.table(paste0(outprefix,"PRS.sscore"))
  prsdat=prsdat %>% select(V1,V6,V7,V8,V9) %>%
    rename(ID=V1,phi1=V6,phi2=V7,phi4=V8,phi6=V9)
  write.table(prsdat,file=paste0(outprefix,"PRS.sscore"),row.names = F,sep="\t",quote=F)
  #return(prsdat)
}
#Target EAS
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix="../result/PRS1/asian_onco_tuning",outprefix="../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_tuning")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix="../result/PRS1/asian_onco_tuning",outprefix="../result/PRS1/prscsx/targetEAS_EASprs_asian_onco_tuning")
#get prs on validation
prefix_val1="../result/PRS1/african_onco_validation"
prefix_val2="../result/PRS1/african_icogs_validation"
prefix_val3="../result/PRS1/asian_onco_validation"
prefix_val4="../result/PRS1/euro_onco_validation"
prefix_val5="../result/PRS1/hispanic_onco_validation"
validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
#../result/PRS1/prscsx/targetEAS_EURprs_african_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val1,outprefix="../result/PRS1/prscsx/targetEAS_EURprs_african_onco_validation")
#../result/PRS1/prscsx/targetEAS_EASprs_african_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val1,outprefix="../result/PRS1/prscsx/targetEAS_EASprs_african_onco_validation")

#../result/PRS1/prscsx/targetEAS_EURprs_african_icogs_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val2,outprefix="../result/PRS1/prscsx/targetEAS_EURprs_african_icogs_validation")
#../result/PRS1/prscsx/targetEAS_EASprs_african_icogs_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val2,outprefix="../result/PRS1/prscsx/targetEAS_EASprs_african_icogs_validation")

#../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val3,outprefix="../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_validation")
#../result/PRS1/prscsx/targetEAS_EASprs_asian_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val3,outprefix="../result/PRS1/prscsx/targetEAS_EASprs_asian_onco_validation")

#../result/PRS1/prscsx/targetEAS_EURprs_euro_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val4,outprefix="../result/PRS1/prscsx/targetEAS_EURprs_euro_onco_validation")
#../result/PRS1/prscsx/targetEAS_EASprs_euro_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val4,outprefix="../result/PRS1/prscsx/targetEAS_EASprs_euro_onco_validation")

#../result/PRS1/prscsx/targetEAS_EURprs_hispanic_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val5,outprefix="../result/PRS1/prscsx/targetEAS_EURprs_hispanic_onco_validation")
#../result/PRS1/prscsx/targetEAS_EASprs_hispanic_onco_validationPRS.sscore
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/asian_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val5,outprefix="../result/PRS1/prscsx/targetEAS_EASprs_hispanic_onco_validation")

#targetEUR
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/prscsx/targetEUR_EURprs_asian_onco_tuning")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix="../result/PRS1/euro_onco_tuning",outprefix="../result/PRS1/prscsx/targetEUR_EASprs_asian_onco_tuning")

get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val1,outprefix="../result/PRS1/prscsx/targetEUR_EURprs_african_onco_validation")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val1,outprefix="../result/PRS1/prscsx/targetEUR_EASprs_african_onco_validation")

get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val2,outprefix="../result/PRS1/prscsx/targetEUR_EURprs_african_icogs_validation")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val2,outprefix="../result/PRS1/prscsx/targetEUR_EASprs_african_icogs_validation")

get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val3,outprefix="../result/PRS1/prscsx/targetEUR_EURprs_asian_onco_validation")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val3,outprefix="../result/PRS1/prscsx/targetEUR_EASprs_asian_onco_validation")

get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val4,outprefix="../result/PRS1/prscsx/targetEUR_EURprs_euro_onco_validation")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val4,outprefix="../result/PRS1/prscsx/targetEUR_EASprs_euro_onco_validation")

get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EUR_pst_eff_a1_b0.5_",bimprefix=prefix_val5,outprefix="../result/PRS1/prscsx/targetEUR_EURprs_hispanic_onco_validation")
get_PRScsxprs_tuning(prefix="../result/PRS1/prscsx/euro_onco_tuning_EAS_pst_eff_a1_b0.5_",bimprefix=prefix_val5,outprefix="../result/PRS1/prscsx/targetEUR_EASprs_hispanic_onco_validation")

#get AUC on validation
targetEASEURprsfiles=paste0("../result/PRS1/prscsx/targetEAS_EURprs_",basename(validationprefix),"PRS.sscore")
targetEASEASprsfiles=paste0("../result/PRS1/prscsx/targetEAS_EASprs_",basename(validationprefix),"PRS.sscore")
PRScsxprs=function(PRS_EUR_tunfile="../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_tuningPRS.sscore",
                   PRS_EAS_tunfile="../result/PRS1/prscsx/targetEAS_EASprs_asian_onco_tuningPRS.sscore",
                   EURprsfiles=targetEASEURprsfiles,EASprsfiles=targetEASEASprsfiles,
                   target="EAS",subtype=NULL)
{
  #get optimal hyperparameter phi
  PRS_EUR_tun=read.table(PRS_EUR_tunfile,header=T)
  PRS_EAS_tun=read.table(PRS_EAS_tunfile,header=T)
  auc_tun=data.frame(EUR=rep(0,4),EAS=rep(0,4),weighted=rep(0,4))
  w1=w2=rep(0,4)
  phenotype1=phenoonco[match(PRS_EUR_tun$ID,phenoonco$ID),]
  if (!is.null(subtype))
  {
    print(subtype)
    phenotype1$y=NA
    phenotype1$y[which(phenotype1[,subtype]==1)]=1
    phenotype1$y[which(is.na(phenotype1$Behaviour1))]=0
  }
  for(i in 1:4)
  {
    pheno.prs=cbind.data.frame(phenotype1,prseur=PRS_EUR_tun[,i+1],prseas=PRS_EAS_tun[,i+1])
    model1 <- glm(y~prseur, data=pheno.prs,family = "binomial")
    predicted1 <- predict(model1,pheno.prs, type="response")
    auc_tun$EUR[i]=as.numeric(pROC::auc(pheno.prs$y,predicted1,quiet=T))
    model2 <- glm(y~prseas, data=pheno.prs,family = "binomial")
    predicted2 <- predict(model2,pheno.prs, type="response")
    auc_tun$EAS[i]=as.numeric(pROC::auc(pheno.prs$y,predicted2,quiet=T))
    model3 <- glm(y~prseur+prseas, data=pheno.prs,family = "binomial")
    
    w1[i]=summary(model3)$coefficients[2,1]
    w2[i]=summary(model3)$coefficients[3,1]
    pheno.prs=cbind.data.frame(phenotype1,weightedprs=w1[i]*PRS_EUR_tun[,i+1]+w2[i]*PRS_EAS_tun[,i+1])
    model4 <- glm(y~weightedprs, data=pheno.prs,family = "binomial")
    predicted4 <- predict(model4,pheno.prs, type="response")
    auc_tun$weighted[i]=as.numeric(pROC::auc(pheno.prs$y,predicted4,quiet=T))
  }
  
  idx_optimal=which.max(auc_tun$weighted)
  print(idx_optimal)
  allprs=allfamval=list()
  for (i in 1:5)
  {
    prs1=read.table(EURprsfiles[i],header=T)
    prs1=prs1[,idx_optimal+1]
    prs2=read.table(EASprsfiles[i],header=T)
    prs2=prs2[,idx_optimal+1]
    prs=prs1*w1[idx_optimal]+prs2*w2[idx_optimal]
    allprs[[i]]=prs
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
  }
  if (is.null(subtype))
  {
    Target_PRScsx_val=get_valauc(allprs,allfamval,outprefix=paste0("../result/PRS1/Target",target),methodprefix="PRScsx",subtype=subtype)
  }else
  {
    Target_PRScsx_val=get_valauc(allprs,allfamval,outprefix=paste0("../result/PRS_subtype/prscsx/",subtype,"/Target",target),methodprefix="PRScsx",subtype=subtype)
  }
  
  
  return(Target_PRScsx_val)
}
EASPRScsxprs=PRScsxprs(PRS_EUR_tunfile="../result/PRS1/prscsx/targetEAS_EURprs_asian_onco_tuningPRS.sscore",
                   PRS_EAS_tunfile="../result/PRS1/prscsx/targetEAS_EASprs_asian_onco_tuningPRS.sscore",
                   EURprsfiles=targetEASEURprsfiles,EASprsfiles=targetEASEASprsfiles,
                   target="EAS")
load("../result/PRS1/TargetEAS_PRScsx_valauc.RData")
asianPRScsx=allres
load("../result/PRS1/TargetEUR_PRScsx_valauc.RData")
euroPRScsx=allres
save(euroCT,asianCT,WeightedCT,euroLDpred,asianLDpred,WeightedLDpred,euroPRScsx,asianPRScsx,file="../result/Run_PRS_BC.RData")

targetEUREURprsfiles=paste0("../result/PRS1/prscsx/targetEUR_EURprs_",basename(validationprefix),"PRS.sscore")
targetEUREASprsfiles=paste0("../result/PRS1/prscsx/targetEUR_EASprs_",basename(validationprefix),"PRS.sscore")
EURPRScsxprs=PRScsxprs(PRS_EUR_tunfile="../result/PRS1/prscsx/targetEUR_EURprs_asian_onco_tuningPRS.sscore",
          PRS_EAS_tunfile="../result/PRS1/prscsx/targetEUR_EASprs_asian_onco_tuningPRS.sscore",
          EURprsfiles=targetEUREURprsfiles,EASprsfiles=targetEUREASprsfiles,
          target="EUR")

EASPRScsxprs=read.table("../result/PRS1/TargetEAS_PRScsx_valauc.txt",header=T)
EURPRScsxprs=read.table("../result/PRS1/TargetEUR_PRScsx_valauc.txt",header=T)
EASPRScsx_table=get_auctable(aucres=EASPRScsxprs)
EURPRScsx_table=get_auctable(aucres=EURPRScsxprs)
PRScsx_table=rbind(EASPRScsx_table,EURPRScsx_table)
rownames(PRScsx_table)=c("TargetEAS","TargetEUR")
write.csv(PRScsx_table,file="../result/PRS1/prscsx/PRScsx_table.csv")

#CTSLEB
#results from CTSLEB_super
CTSLEBprs=function(superresfile="/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EAS/EAS_CTSLEB_runsuper_linear_comb.RData",
       target="EAS",outdir="../result/PRS1/",subtype=NULL)
{
  load(superresfile)
  #find the best PRS from learners
  auctun=rep(0,6)
  for (i in 1:6)
  {
    auctun[i]=res[[8+i]]$auc[1]
  }
  idx_optimal=which.max(auctun)
  print(paste0("idx_optimal: ",idx_optimal))
  #including all the prs
  prsmat=res[[idx_optimal+8]]$allprs
  allprs=allfamval=list()
  for (i in 1:5)
  {
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
    idx=match(famval[,1],prsmat[,1])
    if (any(is.na(idx))) warning("some valiation samples are not in PRS matrix!")
    prs=prsmat[idx,2]
    allprs[[i]]=prs
  }
  CTsleb=get_valauc(allprs,allfamval,outprefix=paste0(outdir,target),methodprefix="CTSLEB",subtype=subtype)
}
EASCTslebprs=read.table("../result/PRS1/EAS_CTSLEB_valauc.txt",header=T)
EASCTsleb_table=get_auctable(aucres=EASCTslebprs)
write.csv(EASCTsleb_table,file="../result/PRS1/ctsleb/EASCTsleb_table.csv",row.names = F)
CTSLEBprs(superresfile="/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/ctsleb/EUR/EUR_CTSLEB_runsuper_linear_comb.RData",
                   target="EUR",outdir="../result/PRS1/",subtype=NULL)
EURCTslebprs=read.table("../result/PRS1/EUR_CTSLEB_valauc.txt",header=T)
EURCTsleb_table=get_auctable(aucres=EURCTslebprs)
write.csv(EURCTsleb_table,file="../result/PRS1/ctsleb/EURCTsleb_table.csv",row.names = F)

#PRSCSx on subtypes
PRScsx_subtype=function(subtype="LumA",outfolder="../result/PRS_subtype/prscsx/LumA/")
{
  if (!dir.exists(outfolder)) dir.create(outfolder)
  
  #Target EAS
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EURprs_asian_onco_tuningPRS.sscore
  prefix1=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/asian_onco_tuning_EUR_pst_eff_a1_b0.5_")
  get_PRScsxprs_tuning(prefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/asian_onco_tuning_EUR_pst_eff_a1_b0.5_"),bimprefix="../result/PRS1/asian_onco_tuning",outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_asian_onco_tuning"))
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EASprs_asian_onco_tuningPRS.sscore
  prefix2=prefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/asian_onco_tuning_EAS_pst_eff_a1_b0.5_")
  get_PRScsxprs_tuning(prefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/asian_onco_tuning_EAS_pst_eff_a1_b0.5_"),bimprefix="../result/PRS1/asian_onco_tuning",outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_asian_onco_tuning"))
  
  #get prs on validation
  prefix_val1="../result/PRS1/african_onco_validation"
  prefix_val2="../result/PRS1/african_icogs_validation"
  prefix_val3="../result/PRS1/asian_onco_validation"
  prefix_val4="../result/PRS1/euro_onco_validation"
  prefix_val5="../result/PRS1/hispanic_onco_validation"
  validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EURprs_african_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val1,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_african_onco_validation"))
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EASprs_african_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val1,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_african_onco_validation"))
  
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EURprs_african_icogs_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val2,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_african_icogs_validation"))
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EASprs_african_icogs_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val2,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_african_icogs_validation"))
  
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EURprs_asian_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val3,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_asian_onco_validation"))
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EASprs_asian_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val3,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_asian_onco_validation"))
  
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EURprs_euro_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val4,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_euro_onco_validation"))
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EASprs_euro_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val4,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_euro_onco_validation"))
  
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EURprs_hispanic_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val5,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_hispanic_onco_validation"))
  #../result/PRS_subtype/prscsx/asian/LumA/targetEAS_EASprs_hispanic_onco_validationPRS.sscore
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val5,outprefix=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_hispanic_onco_validation"))
  
  #targetEUR
  prefix1=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/euro_onco_tuning_EUR_pst_eff_a1_b0.5_")
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix="../result/PRS1/euro_onco_tuning",outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_asian_onco_tuning"))
  prefix2=prefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/euro_onco_tuning_EAS_pst_eff_a1_b0.5_")
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix="../result/PRS1/euro_onco_tuning",outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_asian_onco_tuning"))
  
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val1,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_african_onco_validation"))
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val1,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_african_onco_validation"))
  
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val2,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_african_icogs_validation"))
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val2,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_african_icogs_validation"))
  
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val3,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_asian_onco_validation"))
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val3,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_asian_onco_validation"))
  
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val4,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_euro_onco_validation"))
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val4,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_euro_onco_validation"))
  
  get_PRScsxprs_tuning(prefix=prefix1,bimprefix=prefix_val5,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_hispanic_onco_validation"))
  get_PRScsxprs_tuning(prefix=prefix2,bimprefix=prefix_val5,outprefix=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_hispanic_onco_validation"))
  
  #get AUC on validation
  targetEASEURprsfiles=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_",basename(validationprefix),"PRS.sscore")
  targetEASEASprsfiles=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_",basename(validationprefix),"PRS.sscore")
  EASPRScsxprs=PRScsxprs(PRS_EUR_tunfile=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EURprs_asian_onco_tuningPRS.sscore"),
                         PRS_EAS_tunfile=paste0("../result/PRS_subtype/prscsx/asian/",subtype,"/targetEAS_EASprs_asian_onco_tuningPRS.sscore"),
                         EURprsfiles=targetEASEURprsfiles,EASprsfiles=targetEASEASprsfiles,
                         target="EAS",subtype=subtype)
  targetEUREURprsfiles=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_",basename(validationprefix),"PRS.sscore")
  targetEUREASprsfiles=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_",basename(validationprefix),"PRS.sscore")
  EURPRScsxprs=PRScsxprs(PRS_EUR_tunfile=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EURprs_asian_onco_tuningPRS.sscore"),
                         PRS_EAS_tunfile=paste0("../result/PRS_subtype/prscsx/euro/",subtype,"/targetEUR_EASprs_asian_onco_tuningPRS.sscore"),
                         EURprsfiles=targetEUREURprsfiles,EASprsfiles=targetEUREASprsfiles,
                         target="EUR",subtype=subtype)
}
for (i in 5:length(subtypes))
{
  subtype=subtypes[i]
  print(Sys.time())
  PRScsx_subtype(subtype=subtype,outfolder=paste0("../result/PRS_subtype/prscsx/",subtype,"/"))
  print(Sys.time())
}

#PROSPER result
plink2score=function(bimprefix,scorefile,outprefix)
{
  scoredat=as.data.frame(fread(scorefile,header=T))
  cmd=paste0(plink2," --pfile ",bimprefix," --score-col-nums ",ncol(scoredat),
             " --score ",scorefile," cols=+scoresums,-scoreavgs header no-mean-imputation ",
             " --out ",outprefix," --memory 64000 --threads 8")
  system(cmd)
}

PROSPERprs=function(scorefile="/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/prosper/PROSPER/after_ensemble_EAS/PROSPER_prs_file.txt",
                   target="EAS",outdir="../result/PRS1/prosper/",subtype=NULL)
{
  score=as.data.frame(fread(scorefile,header=T))
  scorefile1=paste0(outdir,target,"_prosper_score.txt")
  fwrite(score[,c("rsid","a1","weight")],file=scorefile1,row.names = F,sep="\t",quote=F)
  score_val1=paste0(outdir,basename(prefix_val1))
  plink2score(bimprefix=prefix_val1,scorefile=scorefile1,outprefix=score_val1)
  score_val2=paste0(outdir,basename(prefix_val2))
  plink2score(bimprefix=prefix_val2,scorefile=scorefile1,outprefix=score_val2)
  score_val3=paste0(outdir,basename(prefix_val3))
  plink2score(bimprefix=prefix_val3,scorefile=scorefile1,outprefix=score_val3)
  score_val4=paste0(outdir,basename(prefix_val4))
  plink2score(bimprefix=prefix_val4,scorefile=scorefile1,outprefix=score_val4)
  score_val5=paste0(outdir,basename(prefix_val5))
  plink2score(bimprefix=prefix_val5,scorefile=scorefile1,outprefix=score_val5)
  validationscoreprefix=c(score_val1,score_val2,score_val3,score_val4,score_val5)
  allprs=allfamval=list()
  for (i in 1:5)
  {
    famval=read.table(paste0(validationprefix[i],".fam"))
    allfamval[[i]]=famval
    prsmat=as.data.frame(fread(paste0(validationscoreprefix[i],".sscore"),header=T))
    idx=match(famval[,1],prsmat[,1])
    if (any(is.na(idx))) warning("some valiation samples are not in PRS matrix!")
    prs=prsmat[idx,ncol(prsmat)]
    allprs[[i]]=prs
  }
  prospers=get_valauc(allprs,allfamval,outprefix=paste0(outdir,target),methodprefix="PROSPERS",subtype=subtype)
}
PROSPERprs(scorefile="/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/prosper/PROSPER/after_ensemble_EUR/PROSPER_prs_file.txt",
                    target="EUR",outdir="../result/PRS1/prosper/",subtype=NULL)
EASPROSPERprs=read.table("../result/PRS1/prosper/EAS_PROSPERS_valauc.txt",header=T)
EURPROSPERprs=read.table("../result/PRS1/prosper/EUR_PROSPERS_valauc.txt",header=T)
EASPROSPER_table=get_auctable(aucres=EASPROSPERprs)
EURPROSPER_table=get_auctable(aucres=EURPROSPERprs)
PROSPER_table=rbind(EASPROSPER_table,EURPROSPER_table)
rownames(PROSPER_table)=c("TargetEAS","TargetEUR")
write.csv(PROSPER_table,file="../result/PRS1/prosper/PROSPER_table.csv")  
