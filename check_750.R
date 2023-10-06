#!/usr/bin/env Rscript

# Concept #750
# Use all n=172,338 samples in the Oncoarray dataset. iCOGS overlap =19,311 so iCOGS only samples = 94,236 (113,547-19,311)
# Therefore total samples that can be used in this analyses = 172,338 + 94,236 = 266,574.
library(data.table)
library(readxl)

#to find out samples with genotype data
basedir="/data/BB_Bioinformatics/ProjectData/BCAC/"
#onco
tmp=list.files(paste0(basedir,"onco/"),"*.sample")
onco_asian=tmp[grepl("asian",tmp)]
onco_african=tmp[grepl("african",tmp)]
onco_hispanic=tmp[grepl("hispanic",tmp)]
onco_european=tmp[!tmp %in% c(onco_asian,onco_african,onco_hispanic)]

get_sampelnames=function(myfile=paste0(basedir,"onco/",onco_asian[1]))
{
  tmp=read.table(myfile,header = T)
  tmp=tmp$ID_1
  tmp=tmp[tmp!=0]
  tmp=unlist(strsplit(tmp,"_"))
  return(tmp[seq(1,length(tmp),2)])
}
onco_asiansamples=get_sampelnames(paste0(basedir,"onco/",onco_asian[2])) #27515
onco_africansamples=get_sampelnames(paste0(basedir,"onco/",onco_african[2])) #5804
onco_hispanicsamples=get_sampelnames(paste0(basedir,"onco/",onco_hispanic[2])) #2560
onco_europeansamples=get_sampelnames(paste0(basedir,"onco/",onco_european[2])) #133236
onco_gensamples=c(onco_asiansamples,onco_africansamples,onco_hispanicsamples,onco_europeansamples)
onco_gensamples=onco_gensamples[onco_gensamples!=0]
length(onco_gensamples) #169115

#icogs
tmp=list.files(paste0(basedir,"icogs/"),"*.sample")
icogs_asian=tmp[grepl("asian",tmp)]
icogs_african=tmp[grepl("african",tmp)]
icogs_european=tmp[grepl("euro",tmp)]
length(icogs_asian)+length(icogs_african)+length(icogs_european)

icogs_asiansamples=get_sampelnames(paste0(basedir,"icogs/",icogs_asian[2])) #12500, first row is 0 0
icogs_africansamples=get_sampelnames(paste0(basedir,"icogs/",icogs_african[2])) #2047
icogs_europeansamples=get_sampelnames(paste0(basedir,"icogs/",icogs_european[2])) #99000

length(intersect(onco_europeansamples,icogs_europeansamples)) #34157
length(intersect(onco_asiansamples,icogs_asiansamples)) #1393
length(intersect(onco_africansamples,icogs_africansamples)) #0
icogs_gensamples=c(icogs_asiansamples,icogs_africansamples,icogs_europeansamples)
icogs_gensamples=icogs_gensamples[icogs_gensamples!=0]
length(icogs_gensamples) #113547


#check sample sizes based on phenotype data
pheno_icogs=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
dim(pheno_icogs)
# [1] 113547     47
pheno_onco=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_onco_pheno_v15_02_corrected.txt",header=T,sep="\t")
dim(pheno_onco)
# [1] 172338     46
table(pheno_onco$Onc_ID %in% onco_gensamples)

#only include controls and invasives: Behaviour1(1,888,NA),NAs are controls
pheno_icogs=pheno_icogs[which(is.na(pheno_icogs$Behaviour1)| pheno_icogs$Behaviour1==1),]
dim(pheno_icogs) #106604
pheno_onco=pheno_onco[which(is.na(pheno_onco$Behaviour1)| pheno_onco$Behaviour1==1),]
dim(pheno_onco) #161229

#form age variable
pheno_icogs$age=NA
idx=which(is.na(pheno_icogs$Behaviour1))
pheno_icogs$age[idx]=pheno_icogs$ageInt[idx]
idx=which(pheno_icogs$Behaviour1==1)
pheno_icogs$age[idx]=pheno_icogs$AgeDiagIndex[idx]
pheno_icogs$age[which(pheno_icogs$age==888)]=NA
table(pheno_icogs$Behaviour1,is.na(pheno_icogs$age),useNA="ifany")
# FALSE  TRUE
# 1    54855   224
# <NA> 48501  3024
pheno_onco$age=NA
idx=which(is.na(pheno_onco$Behaviour1))
pheno_onco$age[idx]=pheno_onco$ageInt[idx]
idx=which(pheno_onco$Behaviour1==1)
pheno_onco$age[idx]=pheno_onco$AgeDiagIndex[idx]
pheno_onco$age[which(pheno_onco$age==888)]=NA
table(pheno_onco$Behaviour1,is.na(pheno_onco$age),useNA = "ifany")
# FALSE  TRUE
# 1    83489   131
# <NA> 76013  1596
#keep samples only have age info
pheno_icogs=pheno_icogs[!is.na(pheno_icogs$age),]
pheno_onco=pheno_onco[!is.na(pheno_onco$age),]
dim(pheno_icogs) #103356
dim(pheno_onco) #159502

#remove overlap samples in icogs
olapsamples=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/iCOGS_Onco_overlap_v15_list_for_investigators.txt",header=T,sep="\t")
length(unique(olapsamples$ICOGS_ID)) #19608
length(unique(olapsamples$Onco_ID)) #19942
sum(duplicated(olapsamples$ICOGS_ID)) #402
sum(duplicated(olapsamples$Onco_ID)) #68

idx=which(olapsamples$ICOGS_ID %in% pheno_icogs$SG_ID)
icogssamples_torm=olapsamples$ICOGS_ID[idx]
length(icogssamples_torm) #18376
sum(duplicated(icogssamples_torm)) #253
icogssamples_torm=unique(icogssamples_torm)
length(icogssamples_torm) #18123
pheno_icogs=pheno_icogs[!pheno_icogs$SG_ID %in% icogssamples_torm,]
dim(pheno_icogs)
#85233    48

table(pheno_icogs$SG_ID %in% icogs_gensamples) #all samples with pheno have genotype data
# TRUE 
# 85233
table(pheno_onco$Onc_ID %in% onco_gensamples) #7887 pheno not have genotype data
# FALSE   TRUE 
# 7887 151615 
#keep samples with genotype data
pheno_onco=pheno_onco[pheno_onco$Onc_ID %in% onco_gensamples,]
dim(pheno_onco)
#151615     47

#write.table(pheno_icogs,file="../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",row.names = F,sep="\t",quote=F)
#write.table(pheno_onco,file="../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",row.names = F,sep="\t",quote=F)

#remove icogs African studies
table(pheno_icogs$EthnicityGeno)
# African    Asian European 
# 1758    10822    72653 
unique(pheno_icogs$study[pheno_icogs$EthnicityGeno=="African"])
#"NBHS" "SCCS"
pheno_icogs=pheno_icogs[pheno_icogs$EthnicityGeno!="African",]
dim(pheno_icogs)
#83475    48
#remove onco other studies
table(pheno_onco$EthnicityGeno)
# African    Asian European    other 
# 5569    25214   118419     2413
unique(pheno_onco$study[pheno_onco$EthnicityGeno=="other"])
#"CAMA"    "COLBCCC"
pheno_onco=pheno_onco[pheno_onco$EthnicityGeno!="other",]
dim(pheno_onco)
#149202     47

#To generate supplementary table1
#only remove Norway data in onco
study2rm=c("ICICLE","GLACIER","LMBC") #remove in icogs
# LMBC excluded because inflation was not reduced by the PCs. 
# ICICLE mainly Ductal in situ
# GLACIER mainly lobular invasive

#stratify by ancestries
allstudy=unique(c(pheno_icogs$study,pheno_onco$study))
allcountry=rep(NA,length(allstudy))
for (i in 1:length(allstudy))
{
  idx=match(allstudy[i],pheno_icogs$study)
  if (sum(!is.na(idx))>0)
  {
    allcountry[i]=pheno_icogs$StudyCountry[idx[1]]
    
  }else
  {
    idx=match(allstudy[i],pheno_onco$study)
    allcountry[i]=pheno_onco$StudyCountry[idx[1]]
    
  }
}

sup11=data.frame(Acronym=allstudy,Country=allcountry,icogs_ancestry="",icogs_ctrl="",icogs_invasive="",onco_ancestry="",onco_ctrl="",onco_invasive="")
for (i in 1:nrow(sup11))
{
  {
    idx=which(pheno_icogs$study==sup11$Acronym[i])
    n=length(idx)
    if(length(idx)>0) sup11$icogs_ancestry[i]=paste0(unique(pheno_icogs$EthnicityGeno[idx]),collapse = "_")
    idx=which(pheno_icogs$study==sup11$Acronym[i] & is.na(pheno_icogs$Behaviour1))
    if (length(idx)>0) sup11$icogs_ctrl[i]=length(idx)
    n1=length(idx)
    idx=which(pheno_icogs$study==sup11$Acronym[i] & pheno_icogs$Behaviour1==1)
    if (length(idx)>0) sup11$icogs_invasive[i]=length(idx)
    n2=length(idx)
    if(n!=n1+n2) stop(i)
    idx=which(pheno_onco$study==sup11$Acronym[i])
    if(length(idx)>0) sup11$onco_ancestry[i]=paste0(unique(pheno_onco$EthnicityGeno[idx]),collapse = "_")
    idx=which(pheno_onco$study==sup11$Acronym[i] & is.na(pheno_onco$Behaviour1)) 
    if (length(idx)>0) sup11$onco_ctrl[i]=length(idx)
    idx=which(pheno_onco$study==sup11$Acronym[i] & pheno_onco$Behaviour1==1)
    if (length(idx)>0) sup11$onco_invasive[i]=length(idx)
    if (sup11$Country[i]=="Norway")
    {
      sup11$onco_ctrl[i]=""
      sup11$onco_invasive[i]=""
    }
    if (sup11$Acronym[i] %in% study2rm)
    {
      sup11$icogs_ctrl[i]=""
      sup11$icogs_invasive[i]=""
    }
  } 
}

sup11=sup11[order(sup11$Acronym),]
#remove empty rows
idx=which(sup11$icogs_ctrl!=""|sup11$icogs_invasive!=""|sup11$onco_ctrl!=""|sup11$onco_invasive!="")
sup11=sup11[idx,]
#write.csv(sup11,file="../result/supplement750_1.csv",row.names = F,quote=F)
sum(as.numeric(sup11$icogs_ctrl),na.rm=T)+sum(as.numeric(sup11$icogs_invasive),na.rm=T)
#[1] 78288 icogs
sum(as.numeric(sup11$icogs_ctrl),na.rm=T) #39445 icogs control
sum(as.numeric(sup11$icogs_invasive),na.rm=T) #38843 icogs case

#sample used in analysis
# icogs_samples=pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm]
# write.table(icogs_samples,file="../result/icogs_samples_750.txt",row.names = F,col.names = F,quote=F)
#this is the final sample table after removing studies
pheno_icogs_final=pheno_icogs[!pheno_icogs$study %in% study2rm,]
dim(pheno_icogs_final) #78288
idx=match(icogs_samples,pheno_icogs$SG_ID)
tmp=pheno_icogs$EthnicityGeno[idx]

sum(as.numeric(sup11$onco_ctrl),na.rm=T)+sum(as.numeric(sup11$onco_invasive),na.rm=T) #148051 onco samples
sum(as.numeric(sup11$onco_ctrl),na.rm=T) #66932 onco control
sum(as.numeric(sup11$onco_invasive),na.rm=T) #81119 onco case
# onco_samples=pheno_onco$Onc_ID[pheno_onco$study %in% sup11$Acronym & pheno_onco$StudyCountry!="Norway"]
# write.table(onco_samples,file="../result/onco_samples_750.txt",row.names = F,col.names = F,quote=F)
#final onco sample table 
pheno_onco_final=pheno_onco[pheno_onco$StudyCountry!="Norway",]
dim(pheno_onco_final) #148051

#write.table(pheno_icogs_final,file="../data/concept_750_zhang_icogs_pheno_v15_02_age_used.txt",row.names = F,sep="\t",quote=F)
#write.table(pheno_onco_final,file="../data/concept_750_zhang_onco_pheno_v15_02_corrected_age_used.txt",row.names = F,sep="\t",quote=F)
#check correlation among ancestries


#Supplemental table 2
sup22=data.frame(Acronym=allstudy,Country=allcountry,icogs_ctrl="",icogs_ERP="",icogs_ERN="",icogs_ERU="",icogs_PRP="",icogs_PRN="",icogs_PRU="",
                icogs_HER2P="",icogs_HER2N="",icogs_HER2U="",icogs_grade1="",icogs_grade2="",icogs_grade3="",icogs_gradeU="",
                onco_ctrl="",onco_ERP="",onco_ERN="",onco_ERU="",onco_PRP="",onco_PRN="",onco_PRU="",
                onco_HER2P="",onco_HER2N="",onco_HER2U="",onco_grade1="",onco_grade2="",onco_grade3="",onco_gradeU="")
for (i in 1:nrow(sup22))
{
  if (!sup22$Acronym[i] %in% study2rm)
  {
    idx=which(pheno_icogs$study==sup22$Acronym[i] & is.na(pheno_icogs$Behaviour1) )
    if (length(idx)>0) sup22$icogs_ctrl[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$ER_status1==1 )
    if (length(idx)>0) sup22$icogs_ERP[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$ER_status1==0 )
    if (length(idx)>0) sup22$icogs_ERN[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$ER_status1==888 )
    if (length(idx)>0) sup22$icogs_ERU[i]=length(idx)
    
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$PR_status1==1 )
    if (length(idx)>0) sup22$icogs_PRP[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$PR_status1==0 )
    if (length(idx)>0) sup22$icogs_PRN[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$PR_status1==888 )
    if (length(idx)>0) sup22$icogs_PRU[i]=length(idx)
    
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$HER2_status1==1 )
    if (length(idx)>0) sup22$icogs_HER2P[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$HER2_status1==0 )
    if (length(idx)>0) sup22$icogs_HER2N[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$HER2_status1==888 )
    if (length(idx)>0) sup22$icogs_HER2U[i]=length(idx)
    
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$Grade1==1 )
    if (length(idx)>0) sup22$icogs_grade1[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$Grade1==2 )
    if (length(idx)>0) sup22$icogs_grade2[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$Grade1==3 )
    if (length(idx)>0) sup22$icogs_grade3[i]=length(idx)
    idx=which(pheno_icogs$study==sup22$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$Grade1==888 )
    if (length(idx)>0) sup22$icogs_gradeU[i]=length(idx)
  }
  
  if (!sup22$Country[i] %in% "Norway")
  {
    idx=which(pheno_onco$study==sup22$Acronym[i] & is.na(pheno_onco$Behaviour1) ) 
    if (length(idx)>0) sup22$onco_ctrl[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$ER_status1==1 )
    if (length(idx)>0) sup22$onco_ERP[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$ER_status1==0 )
    if (length(idx)>0) sup22$onco_ERN[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$ER_status1==888 )
    if (length(idx)>0) sup22$onco_ERU[i]=length(idx)
    
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$PR_status1==1 )
    if (length(idx)>0) sup22$onco_PRP[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$PR_status1==0 )
    if (length(idx)>0) sup22$onco_PRN[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$PR_status1==888 )
    if (length(idx)>0) sup22$onco_PRU[i]=length(idx)
    
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$HER2_status1==1 )
    if (length(idx)>0) sup22$onco_HER2P[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$HER2_status1==0 )
    if (length(idx)>0) sup22$onco_HER2N[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$HER2_status1==888 )
    if (length(idx)>0) sup22$onco_HER2U[i]=length(idx)
    
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Grade1==1 )
    if (length(idx)>0) sup22$onco_grade1[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Grade1==2 )
    if (length(idx)>0) sup22$onco_grade2[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Grade1==3 )
    if (length(idx)>0) sup22$onco_grade3[i]=length(idx)
    idx=which(pheno_onco$study==sup22$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Grade1==888 )
    if (length(idx)>0) sup22$onco_gradeU[i]=length(idx)
  }
}
#remove empty studies 
sup22=sup22[sup22$Acronym %in% sup11$Acronym,]
sup22=sup22[order(sup22$Acronym),]
#write.csv(sup22,file="../result/supplement750_2.csv",row.names = F,quote=F)


#Supplementary table 3: other is Hispanic
sup33=data.frame(ancestry=c("European","Asian","African","other"),icogs_ctrl="",icogs_invasive="",onco_ctrl="",onco_invasive="")
for (i in 1:4)
{
  idx=which(pheno_icogs_final$EthnicityGeno==sup33$ancestry[i] & is.na(pheno_icogs_final$Behaviour1))
  if (length(idx)>0)
  sup33$icogs_ctrl[i]=length(idx)
  idx=which(pheno_icogs_final$EthnicityGeno==sup33$ancestry[i] & pheno_icogs_final$Behaviour1==1)
  if (length(idx)>0)
    sup33$icogs_invasive[i]=length(idx)
  idx=which(pheno_onco_final$EthnicityGeno==sup33$ancestry[i] & is.na(pheno_onco_final$Behaviour1))
  if (length(idx)>0)
  sup33$onco_ctrl[i]=length(idx)
  idx=which(pheno_onco_final$EthnicityGeno==sup33$ancestry[i] & pheno_onco_final$Behaviour1==1)
  if (length(idx)>0)
  sup33$onco_invasive[i]=length(idx)
}
write.csv(sup33,file="../result/supplement750_3.csv",row.names = F,quote=F)

#subtype sample size
# we could define five intrinsic breast cancer subtypes based on the four tumor characteristics: ER, PR, HER2,
# grade. The five intrinsic subtypes are: 1. (ER or PR)+, HER2-, grade 1 or 2; 
# 2. (ER or PR)+, HER2+; # 3. (ER or PR)+, HER2-, grade 3; 
# 4. (ER & PR)-, HER2-; 5. ER-PR-HER2-.

sup44=function(myphenoicogs=pheno_icogs_final,myphenoonco=pheno_onco_final)
{
  sup44=data.frame(matrix(NA,nrow=2,ncol=5))
  rownames(sup44)=c("icogs","onco")
  colnames(sup44)=c("HRP_HER2N_lowgrade","HRP_HER2P","HRP_HER2N_highgrade","HRN_HER2P","HRN_HER2N")
  idx=which((myphenoicogs$ER_status1==1|myphenoicogs$PR_status1==1) & myphenoicogs$HER2_status1==0 & myphenoicogs$Grade1 %in% c(1,2))
  sup44[1,1]=length(idx)
  idx=which((myphenoicogs$ER_status1==1|myphenoicogs$PR_status1==1) & myphenoicogs$HER2_status1==1)
  sup44[1,2]=length(idx)
  idx=which((myphenoicogs$ER_status1==1|myphenoicogs$PR_status1==1) & myphenoicogs$HER2_status1==0 & myphenoicogs$Grade1 %in% c(3))
  sup44[1,3]=length(idx)
  idx=which((myphenoicogs$ER_status1==0 & myphenoicogs$PR_status1==0) & myphenoicogs$HER2_status1==1)
  sup44[1,4]=length(idx)
  idx=which((myphenoicogs$ER_status1==0 & myphenoicogs$PR_status1==0) & myphenoicogs$HER2_status1==0)
  sup44[1,5]=length(idx)
  
  idx=which((myphenoonco$ER_status1==1|myphenoonco$PR_status1==1) & myphenoonco$HER2_status1==0 & myphenoonco$Grade1 %in% c(1,2))
  sup44[2,1]=length(idx)
  idx=which((myphenoonco$ER_status1==1|myphenoonco$PR_status1==1) & myphenoonco$HER2_status1==1)
  sup44[2,2]=length(idx)
  idx=which((myphenoonco$ER_status1==1|myphenoonco$PR_status1==1) & myphenoonco$HER2_status1==0 & myphenoonco$Grade1 %in% c(3))
  sup44[2,3]=length(idx)
  idx=which((myphenoonco$ER_status1==0 & myphenoonco$PR_status1==0) & myphenoonco$HER2_status1==1)
  sup44[2,4]=length(idx)
  idx=which((myphenoonco$ER_status1==0 & myphenoonco$PR_status1==0) & myphenoonco$HER2_status1==0)
  sup44[2,5]=length(idx)
  return(sup44)
}
tmp=sup44()
#write.csv(tmp,file="../result/supplement750_4.csv",row.names = F,quote=F)

sup44_euro=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="European",],
                 myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="European",])
#write.csv(sup44_euro,file="../result/supplement750_4_euro.csv",row.names = F,quote=F)
sup44_asian=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="Asian",],
                 myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="Asian",])
write.csv(sup44_asian,file="../result/supplement750_4_asian.csv",row.names = F,quote=F)
sup44_african=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="African",],
                 myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="African",])
#write.csv(sup44_african,file="../result/supplement750_4_african.csv",row.names = F,quote=F)
sup44_other=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="other",],
                  myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="other",])
# write.csv(sup44_other,file="../result/supplement750_4_other.csv",row.names = F,quote=F)

