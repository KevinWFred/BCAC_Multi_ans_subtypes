#!/usr/bin/env Rscript

library(data.table)
library(readxl)
pheno_icogs_749=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/iCOGS_pheno_v10_zhang749.txt",header=T,sep="\t")
table(pheno_icogs_749$Behaviour1,pheno_icogs_749$status,useNA="ifany")
dim(pheno_icogs_749)
# [1] 99316    51
pheno_onco_749=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/OncoArray_pheno_v10_zhang749.txt",header=T,sep="\t")
dim(pheno_onco_749)
# [1] 138971     50

#only include controls and invasives Behaviour1(1,888,NA):NAs are controls
pheno_icogs_749=pheno_icogs_749[which(is.na(pheno_icogs_749$Behaviour1)| pheno_icogs_749$Behaviour1==1),]
pheno_onco_749=pheno_onco_749[which(is.na(pheno_onco_749$Behaviour1)| pheno_onco_749$Behaviour1==1),]
dim(pheno_icogs_749) #92938
dim(pheno_onco_749) #130791

#form age variable
pheno_icogs_749$age=NA
idx=which(is.na(pheno_icogs_749$Behaviour1))
pheno_icogs_749$age[idx]=pheno_icogs_749$ageInt[idx]
idx=which(pheno_icogs_749$Behaviour1==1)
pheno_icogs_749$age[idx]=pheno_icogs_749$AgeDiagIndex[idx]
pheno_icogs_749$age[which(pheno_icogs_749$age==888)]=NA
table(pheno_icogs_749$Behaviour1,is.na(pheno_icogs_749$age),useNA="ifany")
#      FALSE  TRUE
# 1    48218   262
# <NA> 41437  3021
pheno_onco_749$age=NA
idx=which(is.na(pheno_onco_749$Behaviour1))
pheno_onco_749$age[idx]=pheno_onco_749$ageInt[idx]
idx=which(pheno_onco_749$Behaviour1==1)
pheno_onco_749$age[idx]=pheno_onco_749$AgeDiagIndex[idx]
pheno_onco_749$age[which(pheno_onco_749$age==888)]=NA
table(pheno_onco_749$Behaviour1,is.na(pheno_onco_749$age),useNA = "ifany")
#      FALSE  TRUE
# 1    72674    91
# <NA> 56671  1355
#keep samples only have age info
pheno_icogs_749=pheno_icogs_749[!is.na(pheno_icogs_749$age),]
pheno_onco_749=pheno_onco_749[!is.na(pheno_onco_749$age),]
dim(pheno_icogs_749) #89655
dim(pheno_onco_749) #129345

#remove overlap samples in icogs
olapsamples=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/iCOGS_Onco_overlap_v10_list_for_investigators.txt",header=T,sep="\t")
length(unique(olapsamples$ICOGS_ID)) #19018
length(unique(olapsamples$Onco_ID)) #19018
idx=which(olapsamples$ICOGS_ID %in% pheno_icogs_749$SG_ID)
icogssamples_749_torm=olapsamples$ICOGS_ID[idx]
length(icogssamples_749_torm) #16715
sum(duplicated(icogssamples_749_torm)) #0
icogssamples_749_torm=unique(icogssamples_749_torm)
length(icogssamples_749_torm) #16715
pheno_icogs_749=pheno_icogs_749[!pheno_icogs_749$SG_ID %in% icogssamples_749_torm,]
dim(pheno_icogs_749) #72940


oldsup1=as.data.frame(read_excel("/data/BB_Bioinformatics/ProjectData/BCAC/doc/41588_2020_609_MOESM3_ESM.xlsx",sheet = 1,skip=2))

allstudy_749=unique(c(pheno_icogs_749$study,pheno_onco_749$study))
allcountry_749=rep(NA,length(allstudy_749))
for (i in 1:length(allstudy_749))
{
  idx=match(allstudy_749[i],pheno_icogs_749$study)
  if (sum(!is.na(idx))>0)
  {
    allcountry_749[i]=pheno_icogs_749$StudyCountry[idx[1]]
    
  }else
  {
    idx=match(allstudy_749[i],pheno_onco_749$study)
    allcountry_749[i]=pheno_onco_749$StudyCountry[idx[1]]
    
  }
}

#Sample size tables for 749

#only remove Norway data in onco
study2rm=c("ICICLE","GLACIER","LMBC") #remove in icogs
sum(study2rm %in% allstudy_749) #3

sup1=data.frame(Acronym=allstudy_749,Country=allcountry_749,icogs_ctrl="",icogs_invasive="",onco_ctrl="",onco_invasive="")
for (i in 1:nrow(sup1))
{
  idx=which(pheno_icogs_749$study==sup1$Acronym[i] & is.na(pheno_icogs_749$Behaviour1) )
  if (length(idx)>0) sup1$icogs_ctrl[i]=length(idx)
  idx=which(pheno_icogs_749$study==sup1$Acronym[i] & pheno_icogs_749$Behaviour1==1 )
  if (length(idx)>0) sup1$icogs_invasive[i]=length(idx)
  idx=which(pheno_onco_749$study==sup1$Acronym[i] & is.na(pheno_onco_749$Behaviour1) ) 
  if (length(idx)>0) sup1$onco_ctrl[i]=length(idx)
  idx=which(pheno_onco_749$study==sup1$Acronym[i] & pheno_onco_749$Behaviour1==1 )
  if (length(idx)>0) sup1$onco_invasive[i]=length(idx)
  if (sup1$Country[i]=="Norway")
  {
    sup1$onco_ctrl[i]=""
    sup1$onco_invasive[i]=""
  }
  if (sup1$Acronym[i] %in% study2rm)
  {
    sup1$icogs_ctrl[i]=""
    sup1$icogs_invasive[i]=""
  }
}
#remove blank studies
tmp=sup1
idx=which(sup1=="",arr.ind = T)
tmp[idx]=0
tmp1=as.numeric(tmp[,3])+as.numeric(tmp[,4])+as.numeric(tmp[,5])+as.numeric(tmp[,6])
idx=which(tmp1==0)
sup1=sup1[-idx,]

sup1$Acronym[!sup1$Acronym %in% oldsup1$Acronym]
tmp=oldsup1$Acronym[! oldsup1$Acronym %in% sup1$Acronym][1:5]
#studies dropped out
# [1] "BCFR"                                                                                                                                 
# [2] "BPC3"                                                                                                                                 
# [3] "HEBCS"                                                                                                                                
# [4] "UK2"                                                                                                                                  
# [5] "VUMC"
idx=which(oldsup1$Acronym %in% tmp)
sup1=sup1[order(sup1$Acronym),]
write.csv(sup1,file="../result/supplement749_v10_1.csv",row.names = F,quote=F)

sup2=data.frame(Acronym=allstudy_749,Country=allcountry_749,icogs_ctrl="",icogs_ERP="",icogs_ERN="",icogs_ERU="",icogs_PRP="",icogs_PRN="",icogs_PRU="",
                icogs_HER2P="",icogs_HER2N="",icogs_HER2U="",icogs_grade1="",icogs_grade2="",icogs_grade3="",icogs_gradeU="",
                onco_ctrl="",onco_ERP="",onco_ERN="",onco_ERU="",onco_PRP="",onco_PRN="",onco_PRU="",
                onco_HER2P="",onco_HER2N="",onco_HER2U="",onco_grade1="",onco_grade2="",onco_grade3="",onco_gradeU="")
for (i in 1:nrow(sup2))
{
  if (!sup2$Acronym[i] %in% study2rm)
  {
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & is.na(pheno_icogs_749$Behaviour1) )
    if (length(idx)>0) sup2$icogs_ctrl[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$ER_status1==1 )
    if (length(idx)>0) sup2$icogs_ERP[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$ER_status1==0 )
    if (length(idx)>0) sup2$icogs_ERN[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$ER_status1==888 )
    if (length(idx)>0) sup2$icogs_ERU[i]=length(idx)
    
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$PR_status1==1 )
    if (length(idx)>0) sup2$icogs_PRP[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$PR_status1==0 )
    if (length(idx)>0) sup2$icogs_PRN[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$PR_status1==888 )
    if (length(idx)>0) sup2$icogs_PRU[i]=length(idx)
    
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$HER2_status1==1 )
    if (length(idx)>0) sup2$icogs_HER2P[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$HER2_status1==0 )
    if (length(idx)>0) sup2$icogs_HER2N[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$HER2_status1==888 )
    if (length(idx)>0) sup2$icogs_HER2U[i]=length(idx)
    
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$Grade1==1 )
    if (length(idx)>0) sup2$icogs_grade1[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$Grade1==2 )
    if (length(idx)>0) sup2$icogs_grade2[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$Grade1==3 )
    if (length(idx)>0) sup2$icogs_grade3[i]=length(idx)
    idx=which(pheno_icogs_749$study==sup2$Acronym[i] & pheno_icogs_749$Behaviour1==1 & pheno_icogs_749$Grade1==888 )
    if (length(idx)>0) sup2$icogs_gradeU[i]=length(idx)
  }
  
  
  if (!sup2$Country[i] %in% "Norway")
  {
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & is.na(pheno_onco_749$Behaviour1) ) 
    if (length(idx)>0) sup2$onco_ctrl[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$ER_status1==1 )
    if (length(idx)>0) sup2$onco_ERP[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$ER_status1==0 )
    if (length(idx)>0) sup2$onco_ERN[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$ER_status1==888 )
    if (length(idx)>0) sup2$onco_ERU[i]=length(idx)
    
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$PR_status1==1 )
    if (length(idx)>0) sup2$onco_PRP[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$PR_status1==0 )
    if (length(idx)>0) sup2$onco_PRN[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$PR_status1==888 )
    if (length(idx)>0) sup2$onco_PRU[i]=length(idx)
    
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$HER2_status1==1 )
    if (length(idx)>0) sup2$onco_HER2P[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$HER2_status1==0 )
    if (length(idx)>0) sup2$onco_HER2N[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$HER2_status1==888 )
    if (length(idx)>0) sup2$onco_HER2U[i]=length(idx)
    
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$Grade1==1 )
    if (length(idx)>0) sup2$onco_grade1[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$Grade1==2 )
    if (length(idx)>0) sup2$onco_grade2[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$Grade1==3 )
    if (length(idx)>0) sup2$onco_grade3[i]=length(idx)
    idx=which(pheno_onco_749$study==sup2$Acronym[i] & pheno_onco_749$Behaviour1==1 & pheno_onco_749$Grade1==888 )
    if (length(idx)>0) sup2$onco_gradeU[i]=length(idx)
  }
}

tmp=rep(NA,nrow(sup2))
for (i in 1:nrow(sup2))
{
  tmp[i]=sum(as.numeric(unlist(sup2[i,3:ncol(sup2)])),na.rm=T)
}
idx=which(tmp==0)
sup2=sup2[-idx,]
#sup2=sup2[sup2$Acronym %in% oldsup1$Acronym,]
sup2=sup2[order(sup2$Acronym),]
#We excluded the iCOGS  and OncoArray data from Belarus (the Hannover-Minsk Breast Cancer Study) because it lacked PR and HER2 data
# idx=which(sup2$Country=="Belarus")
# sup2=sup2[-idx,]
write.csv(sup2,file="../result/supplement749_v10_2.csv",row.names = F,quote=F)

#sample used in analysis
icogs_samples=pheno_icogs_749$SG_ID[pheno_icogs_749$study %in% sup1$Acronym & !pheno_icogs_749$study %in% study2rm] #67776
onco_samples=pheno_onco_749$Onc_ID[pheno_onco_749$study %in% sup1$Acronym & pheno_onco_749$StudyCountry!="Norway"] #128140
pheno_icogs_749=pheno_icogs_749[pheno_icogs_749$SG_ID %in% icogs_samples,]
write.table(pheno_icogs_749,file="../data/concept_749_zhang_icogs_pheno_v10_age.txt",row.names = F,sep="\t",quote=F)
pheno_onco_749=pheno_onco_749[pheno_onco_749$Onc_ID %in% onco_samples,]
write.table(pheno_onco_749,file="../data/concept_750_zhang_onco_pheno_v10_02_age.txt",row.names = F,sep="\t",quote=F)

