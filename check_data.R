#!/usr/bin/env Rscript


# For concept #749, iCOGS  (n=99,000), Onco (n=136,459). You can only use the data as per the pheno IDs in these files for analyses under this concept.
# For concept #750, iCOGS  (n=113,547), Onco (n=172,338). You can only use the data as per the pheno IDs in these files for analyses under this concept.
# 
# With regards to the overlapping samples for these analyses using the corrected files:
#   
#   Concept #749
# Use all n=136,459 samples in the Oncoarray dataset. iCOGS overlap =17,824 so iCOGS only samples = 81,176 (99,000-17,824)
# Therefore total samples that can be used in this analyses = 136,459 + 81,176 = 217,635.
# 
# Concept #750
# Use all n=172,338 samples in the Oncoarray dataset. iCOGS overlap =19,311 so iCOGS only samples = 94,236 (113,547-19,311)
# Therefore total samples that can be used in this analyses = 172,338 + 94,236 = 266,574.
library(data.table)
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

#plink="/usr/local/apps/plink/1.9/plink"
#plink2="/usr/local/apps/plink/2.3-alpha/plink2"

# prefix=gsub(".sample","",paste0(basedir,"icogs/",asian[1]),fixed = T)
# cmd=paste0(plink2," --bgen ",prefix,".bgen ref-unknown --sample ",prefix,".sample --maf 0.01 --geno 0.1 --max-alleles 2  -make-pgen --out ../result/",basename(prefix))
# system(cmd)
# tmp=read.table("../result/zhang_750_asian_icogs_topmed_1.psam")
# cmd=paste0(plink2," --bgen ",prefix,".bgen ref-unknown --sample ",prefix,".sample --maf 0.01 --geno 0.1 --max-alleles 2 -make-bed --out ../result/",basename(prefix))
# system(cmd)
# tmp=read.table("../result/zhang_750_asian_icogs_topmed_1.fam")
# cmd=paste0(plink2," --pfile ../result/",basename(prefix)," --recode A-transpose --out ../result/",basename(prefix))
# system(cmd)
# tmp=as.data.frame(fread("../result/zhang_750_asian_icogs_topmed_1.traw"))


#check sample size
library(readxl)
pheno_icogs_749=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_749_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
table(pheno_icogs_749$Behaviour1,pheno_icogs_749$status,useNA="ifany")
dim(pheno_icogs_749)
# [1] 99000    50
tmp=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_749_zhang_onco_pheno_v15_02.txt",header=T,sep="\t")
sum(tmp$study %in% c("GC-HBOC")) #5227
pheno_onco_749=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_749_zhang_onco_pheno_v15_02_corrected.txt",header=T,sep="\t")
dim(pheno_onco_749)
# [1] 136459     49
pheno_icogs=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
dim(pheno_icogs)
# [1] 113547     47
pheno_onco=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_onco_pheno_v15_02_corrected.txt",header=T,sep="\t")
dim(pheno_onco)
# [1] 172338     46
table(pheno_onco$Onc_ID %in% onco_gensamples)

#only include controls and invasives Behaviour1(1,888,NA):NAs are controls
pheno_icogs_749=pheno_icogs_749[which(is.na(pheno_icogs_749$Behaviour1)| pheno_icogs_749$Behaviour1==1),]
pheno_onco_749=pheno_onco_749[which(is.na(pheno_onco_749$Behaviour1)| pheno_onco_749$Behaviour1==1),]
pheno_icogs=pheno_icogs[which(is.na(pheno_icogs$Behaviour1)| pheno_icogs$Behaviour1==1),]
dim(pheno_icogs) #106604
pheno_onco=pheno_onco[which(is.na(pheno_onco$Behaviour1)| pheno_onco$Behaviour1==1),]
dim(pheno_onco) #161229

#form age variable
pheno_icogs_749$age=NA
idx=which(is.na(pheno_icogs_749$Behaviour1))
pheno_icogs_749$age[idx]=pheno_icogs_749$ageInt[idx]
idx=which(pheno_icogs_749$Behaviour1==1)
pheno_icogs_749$age[idx]=pheno_icogs_749$AgeDiagIndex[idx]
pheno_icogs_749$age[which(pheno_icogs_749$age==888)]=NA
table(pheno_icogs_749$Behaviour1,is.na(pheno_icogs_749$age),useNA="ifany")
# FALSE  TRUE
# 1    48191   223
# <NA> 41203  3020
pheno_onco_749$age=NA
idx=which(is.na(pheno_onco_749$Behaviour1))
pheno_onco_749$age[idx]=pheno_onco_749$ageInt[idx]
idx=which(pheno_onco_749$Behaviour1==1)
pheno_onco_749$age[idx]=pheno_onco_749$AgeDiagIndex[idx]
pheno_onco_749$age[which(pheno_onco_749$age==888)]=NA
table(pheno_onco_749$Behaviour1,is.na(pheno_onco_749$age),useNA = "ifany")
# FALSE  TRUE
# 1    65943    93
# <NA> 60363  1362

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
pheno_icogs_749=pheno_icogs_749[!is.na(pheno_icogs_749$age),]
pheno_onco_749=pheno_onco_749[!is.na(pheno_onco_749$age),]
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
idx=which(olapsamples$ICOGS_ID %in% pheno_icogs_749$SG_ID)
icogssamples_749_torm=olapsamples$ICOGS_ID[idx]
length(icogssamples_749_torm) #16993
sum(duplicated(icogssamples_749_torm)) #252
icogssamples_749_torm=unique(icogssamples_749_torm)
length(icogssamples_749_torm) #16741
pheno_icogs_749=pheno_icogs_749[!pheno_icogs_749$SG_ID %in% icogssamples_749_torm,]
dim(pheno_icogs_749)
#72653    51
table(pheno_icogs_749$SG_ID %in% pheno_icogs$SG_ID[pheno_icogs$EthnicityGeno=="European"])
# TRUE 
# 72653
table(pheno_icogs_749$SG_ID %in% icogs_gensamples) #all pheno have genotype data
# TRUE 
# 72653
table(pheno_onco_749$Onc_ID %in% pheno_onco$Onc_ID[pheno_onco$EthnicityGeno=="European"])
# TRUE 
# 126306
table(pheno_onco_749$Onc_ID %in% onco_gensamples) #118419 pheno have genotype
# FALSE   TRUE 
# 7887 118419
#keep those with genotype
pheno_onco_749=pheno_onco_749[pheno_onco_749$Onc_ID %in% onco_gensamples,]
dim(pheno_onco_749)
#118419     50

idx=which(olapsamples$ICOGS_ID %in% pheno_icogs$SG_ID)
icogssamples_torm=olapsamples$ICOGS_ID[idx]
length(icogssamples_torm) #18376
sum(duplicated(icogssamples_torm)) #253
icogssamples_torm=unique(icogssamples_torm)
length(icogssamples_torm) #18123
pheno_icogs=pheno_icogs[!pheno_icogs$SG_ID %in% icogssamples_torm,]
dim(pheno_icogs)
#85233    48

table(pheno_icogs$SG_ID %in% icogs_gensamples) #all pheno have genotype data
# TRUE 
# 85233
table(pheno_onco$Onc_ID %in% onco_gensamples) #7887 pheno not have genotype data
# FALSE   TRUE 
# 7887 151615 
pheno_onco=pheno_onco[pheno_onco$Onc_ID %in% onco_gensamples,]
dim(pheno_onco)
#151615     47

#write.table(pheno_icogs,file="../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",row.names = F,sep="\t",quote=F)
#write.table(pheno_onco,file="../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",row.names = F,sep="\t",quote=F)

#remove icogs African
table(pheno_icogs$EthnicityGeno)
# African    Asian European 
# 1758    10822    72653 
unique(pheno_icogs$study[pheno_icogs$EthnicityGeno=="African"])
#"NBHS" "SCCS"
pheno_icogs=pheno_icogs[pheno_icogs$EthnicityGeno!="African",]
dim(pheno_icogs)
#83475    48
table(pheno_onco$EthnicityGeno)
# African    Asian European    other 
# 5569    25214   118419     2413
unique(pheno_onco$study[pheno_onco$EthnicityGeno=="other"])
#"CAMA"    "COLBCCC"
pheno_onco=pheno_onco[pheno_onco$EthnicityGeno!="other",]
dim(pheno_onco)
#149202     47

#select case/control samples and create pheno file. Doesn't remove any specific studies
update_samplefile=function(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_euro_icogs_topmed_1_p1.sample",
                           prefix="zhang_750_euro_icogs_topmed_",
                           pheno=pheno_icogs)
{
  outfolder=dirname(samplefile)
  sampledat=read.table(samplefile,header = T,sep=" ")
  #change sampleID xxx_xxx to xxx
  sampledat2=sampledat[2:nrow(sampledat),]
  tmp=unlist(strsplit(sampledat2$ID_1,"_"))
  sampledat2$ID_1=tmp[seq(1,length(tmp),2)]
  sampledat2$ID_2=sampledat2$ID_1
  sampledat=rbind(sampledat[1,],sampledat2)
  newsamplefile=paste0(outfolder,"/",prefix,".sample")
  write.table(sampledat,file=newsamplefile,sep="\t",row.names = F,quote=F)
  #sampledat$SEX="D" #Sex code ('1' = male, '2' = female, '0' = unknown)
  #sampledat$pheno="D" #	Binary ('0' = control, '1' = case) discrete (categorical, positive integers), or continuous phenotype; missing values represented by 'NA'
  sampledat=sampledat[2:nrow(sampledat),]
  # sampleid=unlist(strsplit(sampledat$ID_1,"_"))
  # sampleid=sampleid[seq(1,length(sampleid),2)]
  sampleid=sampledat$ID_1
  if (sum(colnames(pheno)=="SG_ID")>0)
  {
    pheno$ID=pheno$SG_ID
  }else
  {
    pheno$ID=pheno$Onc_ID
  }
  tmp=sampleid[!sampleid %in% pheno$ID]
  if (length(tmp)>0) warning(paste0(length(tmp)," samples not in phenotype file"))
  sampleid=sampleid[sampleid %in% pheno$ID]
  sampledat=sampledat[match(sampleid,sampledat$ID_1),]
  idx=match(sampleid,pheno$ID)
  sampledat$SEX=2 
  sampledat$pheno=pheno$Behaviour1[idx]
  
  sampledat$pheno[which(is.na(sampledat$pheno))]=0
  sampledat$pheno[which(sampledat$pheno %in% c(2,888))]=NA
  idx=which(!is.na(sampledat$pheno))
  sampledat=sampledat[idx,]
  print(table(sampledat$pheno))
  plinksamplefile=paste0(outfolder,"/",prefix,"plinksamples.txt")
  plinksample=data.frame(FID=sampledat$ID_1,IID=sampledat$ID_2)
  write.table(plinksample,file=plinksamplefile,sep="\t",row.names = F,quote=F)
  phenofile=paste0(outfolder,"/",prefix,"pheno.txt")
  phenodat=sampledat[,-c(3,4)]
  phenodat$pheno=phenodat$pheno+1 #change to 1 and 2
  colnames(phenodat)[1:2]=c("FID","IID")
  write.table(phenodat,file=phenofile,sep=" ",row.names = F,quote=F)
}
# update_samplefile() #34059 (0) 38594 (1)
# update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_asian_icogs_topmed_1.sample",
#                   prefix="zhang_750_asian_icogs_topmed_",
#                   pheno=pheno_icogs) #5940 (0) 4882 (1)
# update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_african_icogs_topmed_1.sample",
#                   prefix="zhang_750_african_icogs_topmed_",
#                   pheno=pheno_icogs) #817 (0) 941 (1)
# 
# update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_topmed_1_p1.sample",
#                   prefix="zhang_750_topmed_",
#                   pheno=pheno_onco)
# #52500 (0) 65919 (1)
# update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_asian_topmed_1.sample",
#                   prefix="zhang_750_asian_topmed_",
#                   pheno=pheno_onco)
# #12344 (0) 12870 (1)
# update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_african_topmed_1.sample",
#                   prefix="zhang_750_african_topmed_",
#                   pheno=pheno_onco)
# #2088 (0) 3481 (1)
# update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_hispanic_topmed_1.sample",
#                   prefix="zhang_750_hispanic_topmed_",
#                   pheno=pheno_onco)
# #1218 (0) 1195 (1)

#compare 749 and 750:
idx=match(pheno_icogs_749$SG_ID,pheno_icogs$SG_ID)
table(pheno_icogs$EthnicityGeno[idx])
sum(pheno_icogs$EthnicityGeno=="European")==nrow(pheno_icogs_749) # 749 are EUR in 750
idx=which(!pheno_icogs$SG_ID %in% pheno_icogs_749$SG_ID)
table(pheno_icogs$EthnicityGeno[idx])
# African   Asian 
# 1758   10822
idx=match(pheno_onco_749$Onc_ID,pheno_onco$Onc_ID)
table(pheno_onco$EthnicityGeno[idx])
sum(pheno_onco$EthnicityGeno=="European")==nrow(pheno_onco_749) # 749 are EUR in 750
idx=which(!pheno_onco$Onc_ID %in% pheno_onco_749$Onc_ID)
table(pheno_onco$EthnicityGeno[idx])
# African   Asian   other 
# 5569   25214    2413 
table(pheno_onco$Onc_ID %in% pheno_icogs$SG_ID)
length(intersect(pheno_icogs$SG_ID,pheno_onco$Onc_ID)) #32679
tmp=intersect(pheno_icogs$SG_ID,pheno_onco$Onc_ID)
idx1=match(tmp,pheno_icogs$SG_ID)
idx2=match(tmp,pheno_onco$Onc_ID)
table(pheno_icogs$BCAC_ID[idx1]==pheno_onco$BCAC_ID[idx2]) #all False
sum(pheno_icogs$BCAC_ID %in% pheno_onco$BCAC_ID) #0, after filtering

oldsup1=as.data.frame(read_excel("/data/BB_Bioinformatics/ProjectData/BCAC/doc/41588_2020_609_MOESM3_ESM.xlsx",sheet = 1,skip=2))
table(unique(pheno_icogs$study) %in% unique(oldsup1$Acronym))
table(unique(pheno_icogs$study) %in% unique(pheno_onco$study))
table(unique(c(pheno_icogs$study,pheno_onco$study)) %in% unique(oldsup1$Acronym))
# FALSE  TRUE 
# 18    79 
unique(c(pheno_icogs$study,pheno_onco$study))[!unique(c(pheno_icogs$study,pheno_onco$study))  %in% unique(oldsup1$Acronym)]
# [1] "ACP"     "GLACIER" "HERPACC" "ICICLE"  "LAABC"   "MYBRCA"  "SBCGS"   "SCCS"   
# [9] "SEBCS"   "SGBCC"   "TWBCS"   "CAMA"    "COLBCCC" "HKBCS"   "KOHBRA"  "NGOBCS" 
# [17] "UBCS"    "WAABCS"
table(unique(oldsup1$Acronym)  %in% unique(c(pheno_icogs$study,pheno_onco$study)))
# FALSE  TRUE 
# 9   79 
unique(oldsup1$Acronym)[!unique(oldsup1$Acronym) %in% unique(c(pheno_icogs$study,pheno_onco$study))]
# [1] "BCFR"                                                                                                                                 
# [2] "BPC3"                                                                                                                                 
# [3] "GC-HBOC"                                                                                                                              
# [4] "HEBCS"                                                                                                                                
# [5] "UK2"                                                                                                                                  
# [6] "VUMC"                                                                                                                                 
# [7] "WHI"   

#749
#Sample size table 1
# European people
# By iCOGs and OncoArray (remove overlapped sample)
# By study
# Only invasive#

#EthnicityClass 1=European (for EthnicitySubClass 1, 2, 3, 4, 5), 
#Status 0=control, 1=invasive case, 2=in-situ case, 3=case unknown invasiveness, 9=excluded sample

# icogs_gensamples0=icogs_gensamples
# icogs_gensamples_749=icogs_gensamples0[!icogs_gensamples0 %in% icogssamples_749_torm] #113547
# table(icogs_gensamples_749 %in% pheno_icogs_749$SG_ID)
# icogs_gensamples=icogs_gensamples0[!icogs_gensamples0 %in% icogssamples_torm] #94638
# tmp=paste0(icogs_gensamples,"_",icogs_gensamples)
# tmp=data.frame(FID=tmp,IID=tmp)
# write.table(tmp,file="../result/icogs_keepsamples.txt",row.names = F,sep=" ",quote=F)

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

sup1$Acronym[!sup1$Acronym %in% oldsup1$Acronym]
# "GLACIER" "ICICLE"  "UBCS" not in the old table
idx=which(!sup1$Acronym %in% oldsup1$Acronym)
sum(as.numeric(sup1$icogs_ctrl[idx]),na.rm=T)+sum(as.numeric(sup1$onco_ctrl[idx]),na.rm=T)#0 added controls
#GLACIER UK (1888 icogs cases)
#ICICLE UK (1 icogs controls, 147 icogs cases)
#UBCS USA (606 onco cases)
sum(as.numeric(sup1$icogs_invasive[idx]),na.rm=T)+sum(as.numeric(sup1$onco_invasive[idx]),na.rm=T) #606 added cases
oldsup1$Acronym[!oldsup1$Acronym %in% sup1$Acronym]
# [1] "BCFR"                                                                                                                                 
# [2] "BPC3"                                                                                                                                 
# [3] "GC-HBOC"                                                                                                                              
# [4] "HEBCS"                                                                                                                                
# [5] "UK2"                                                                                                                                  
# [6] "VUMC"                                                                                                                                 
# [7] "WHI"
idx=which(!oldsup1$Acronym %in% sup1$Acronym)[1:7]
#HEBCS Finland (1059 icogs controls, 177 icogs controls,1515 icogs cases, 281 onco cases), 
# WHI USA (4617 onco controls, onco 4930 cases)
#GC-HBOC Germany (139 icogs controls, 1593 onco controls, 3416 onco cases)

sum(as.numeric(oldsup1$Controls...5[idx]),na.rm=T)+sum(as.numeric(oldsup1$Controls...9[idx]),na.rm=T) #7585
sum(as.numeric(oldsup1$Invasive...6[idx]),na.rm=T)+sum(as.numeric(oldsup1$Invasive...10[idx]),na.rm=T) #10142
sum(as.numeric(oldsup1$Controls...5[1:86]),na.rm=T)+sum(as.numeric(oldsup1$Controls...9[1:86]),na.rm=T) #96201
sup1=sup1[order(sup1$Acronym),]
sum(as.numeric(sup1[,3]),na.rm=T) #33505
sum(as.numeric(sup1[,4]),na.rm=T) #33961
sum(as.numeric(sup1[,5]),na.rm=T) #52500
sum(as.numeric(sup1[,6]),na.rm=T) #64768
write.csv(sup1,file="../result/supplement749_1.csv",row.names = F,quote=F)
sum(as.numeric(unlist(sup1[,3:6])),na.rm=T)
#[1] 184734 total number of samples
sum(as.numeric(unlist(sup1[,3:4])),na.rm=T)
#67466 icogs
sum(as.numeric(unlist(sup1[,5:6])),na.rm=T)
#117268 onco
sum(as.numeric(unlist(sup1[,c(3,5)])),na.rm=T)
#86005 controls
sum(as.numeric(unlist(sup1[,c(4,6)])),na.rm=T)
#98729 cases

#1 We excluded the iCOGS  and OncoArray data from Belarus (the Hannover-Minsk Breast Cancer Study) because it lacked PR and HER2 data; the iCOGs data from Greece (the Triple-Negative Breast Cancer Consortium) because only TN cases were avaialable; and the OncoArray data from Norway (the Norwegian Breast Cancer Study) because there were no controls available from Norway
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
#sup2=sup2[sup2$Acronym %in% oldsup1$Acronym,]
sup2=sup2[order(sup2$Acronym),]
#We excluded the iCOGS  and OncoArray data from Belarus (the Hannover-Minsk Breast Cancer Study) because it lacked PR and HER2 data
# idx=which(sup2$Country=="Belarus")
# sup2=sup2[-idx,]
write.csv(sup2,file="../result/supplement749_2.csv",row.names = F,quote=F)

# source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
# library(ggplot2)
# theme_set(theme_bw(base_size=12))
# p1=ggplot(pheno_icogs, aes(x=pc1, y=pc2, color=EthnicityGeno)) + geom_point() + 
#   labs(title="iCOGS") #+ theme_Publication()
# p2=ggplot(pheno_icogs[idx,], aes(x=pc1, y=pc3, color=EthnicityGeno)) + geom_point(alpha=0.1) + 
#   labs(title="iCOGS") #+ theme_Publication()
# library("cowplot")
# p3=plot_grid(p1, p2, labels=c("A", "B"), ncol = 2, nrow = 1)
# save_plot("../result/icogs_PCA.pdf", p3,ncol=2,
#           base_asp = 1.1 # make room for figure legend
# )
# 
# p1=ggplot(pheno_onco, aes(x=PC_1, y=PC_2, color=EthnicityGeno)) + geom_point(alpha=0.1,size=2) + 
#   labs(title="") +
#   theme(
#     panel.background = element_rect(fill='transparent'),
#     plot.background = element_rect(fill='transparent', color=NA),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill='transparent'),
#     legend.box.background = element_rect(fill='transparent')
#   )
# p2=ggplot(pheno_onco, aes(x=PC_1, y=PC_3, color=EthnicityGeno)) + geom_point(alpha=0.5,size=3) + 
#   labs(title="") +
#   theme(
#     panel.background = element_rect(fill='transparent'),
#     plot.background = element_rect(fill='transparent', color=NA),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill='transparent'),
#     legend.box.background = element_rect(fill='transparent')
#   )
# library("cowplot")
# p3=plot_grid(p1, p2, labels=c("A", "B"), ncol = 2, nrow = 1)
# save_plot("../result/onco_PCA.pdf", p3,ncol=2,
#           base_asp = 1.1 # make room for figure legend
# )


#work on all ancestry, sample size tables 750
#remove samples without age, include control or invasive 
# pheno_icogs1=pheno_icogs
# pheno_icogs1=pheno_icogs1[pheno_icogs1$SG_ID %in% icogs_gensamples,] #94638
# pheno_icogs1=pheno_icogs1[which(pheno_icogs1$Behaviour1==1 |is.na(pheno_icogs1$Behaviour1)),] #89168 control and invasive
# table(is.na(pheno_icogs1$ageInt),is.na(pheno_icogs1$AgeDiagIndex),useNA="ifany")
# idx1=which(is.na(pheno_icogs1$Behaviour1))
# idx2=which(pheno_icogs1$Behaviour1==1)
# table(is.na(pheno_icogs1$age[idx1]))
# # FALSE  TRUE 
# # 41362  3072
# idx3=which(is.na(pheno_icogs1$age[idx1]))
# tmp=unique(pheno_icogs1$study[idx1[idx3]]) #14
# sum(tmp %in% oldsup1$Acronym) #12
# tmp[!tmp %in% oldsup1$Acronym] #ACP, MYBRCA
# table(pheno_icogs1$study[idx1[idx3]])
# table(pheno_icogs1$EthnicityGeno[idx1[idx3]])
# # Asian European 
# # 4     3068 
# unique(pheno_icogs1$StudyCountry[idx1[idx3]])
# table(pheno_icogs1$StudyCountry[idx1[idx3]])
# table(is.na(pheno_icogs1$age[idx2]))
# # FALSE  TRUE 
# # 44512  222
# idx3=which(is.na(pheno_icogs1$age[idx2]))
# tmp=unique(pheno_icogs1$study[idx2[idx3]]) #9
# sum(tmp %in% oldsup1$Acronym) #7
# tmp[!tmp %in% oldsup1$Acronym] #ACP,GLACIER
# table(pheno_icogs1$study[idx2[idx3]])
# table(pheno_icogs1$EthnicityGeno[idx2[idx3]])
# # Asian European 
# # 1     221
# # table(is.na(pheno_icogs1$AgeDiagIndex[idx1]))
# # # TRUE 
# # # 43857 
# # table(is.na(pheno_icogs1$AgeDiagIndex[idx2]))
# # # FALSE  TRUE 
# # # 44512   224 
# idx3=which(is.na(pheno_icogs1$age))
# pheno_icogs1=pheno_icogs1[-idx3,]
# 
# pheno_onco1=pheno_onco
# pheno_onco1=pheno_onco1[pheno_onco1$Onc_ID %in% onco_gensamples,] #169115
# pheno_onco1=pheno_onco1[which(is.na(pheno_onco1$Behaviour1) | pheno_onco1$Behaviour1==1),] #160729
# table(is.na(pheno_onco1$ageInt),is.na(pheno_onco1$AgeDiagIndex),useNA="ifany")
# idx1=which(is.na(pheno_onco1$Behaviour1))
# idx2=which(pheno_onco1$Behaviour1==1)
# table(is.na(pheno_onco1$age[idx1]))
# # FALSE  TRUE 
# # 71639  2094
# table(is.na(pheno_onco1$age[idx2]))
# # FALSE  TRUE 
# # 86865  131
# idx3=which(is.na(pheno_onco1$age[idx1]))
# tmp=unique(pheno_onco1$study[idx1[idx3]]) #21
# tmp[!tmp %in% oldsup1$Acronym] #""ACP"    "HKBCS"  "MYBRCA"
# table(pheno_onco1$EthnicityGeno[idx1[idx3]])
# # Asian European 
# # 405     1689 
# idx3=which(is.na(pheno_onco1$age[idx2]))
# tmp=unique(pheno_onco1$study[idx2[idx3]]) #12
# tmp[!tmp %in% oldsup1$Acronym] #"ACP"    "CAMA"   "KOHBRA" "WAABCS"
# table(pheno_onco1$EthnicityGeno[idx2[idx3]])
# # African    Asian European    other 
# # 2       35       93        1
# 
# # table(is.na(pheno_onco1$AgeDiagIndex[idx1]))
# # # TRUE 
# # # 79202
# # table(is.na(pheno_onco1$AgeDiagIndex[idx2]))
# # # FALSE  TRUE 
# # # 86889   131
# idx3=which(is.na(pheno_onco1$age))
# pheno_onco1=pheno_onco1[-idx3,]
# 
# icogs_europeansamples1=intersect(icogs_europeansamples,pheno_icogs1$SG_ID) #73212
# icogs_asiansamples1=intersect(icogs_asiansamples,pheno_icogs1$SG_ID) #10903
# icogs_africansamples1=intersect(icogs_africansamples,pheno_icogs1$SG_ID) #1759
# onco_europeansamples1=intersect(onco_europeansamples,pheno_onco1$Onc_ID) #124769
# onco_asiansamples1=intersect(onco_asiansamples,pheno_onco1$Onc_ID) #25669
# onco_africansamples1=intersect(onco_africansamples,pheno_onco1$Onc_ID) #5611
# onco_hispanicsamples1=intersect(onco_hispanicsamples,pheno_onco1$Onc_ID) #2455

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
# sup11all=data.frame(Acronym=allstudy,Country=allcountry,icogs_ctrl_eur=0,icogs_ctrl_asn=0,icogs_ctrl_afr=0,
#                     icogs_invasive_eur=0,icogs_invasive_asn=0,icogs_invasive_afr=0,
#                     onco_ctrl_eur=0,onco_ctrl_asn=0,onco_ctrl_afr=0,onco_ctrl_hsp=0,
#                     onco_invasive_eur=0,onco_invasive_asn=0,onco_invasive_afr=0,onco_invasive_hsp=0)
# for (i in 1:nrow(sup11all))
# {
#   idx1=which(pheno_icogs$study==sup11all$Acronym[i] & is.na(pheno_icogs$Behaviour1) & pheno_icogs$SG_ID %in% icogs_europeansamples)
#   sup11all$icogs_ctrl_eur[i]=length(idx1)
#   idx1=which(pheno_icogs$study==sup11all$Acronym[i] & is.na(pheno_icogs$Behaviour1) & pheno_icogs$SG_ID %in% icogs_asiansamples)
#   sup11all$icogs_ctrl_asn[i]=length(idx1)
#   idx1=which(pheno_icogs$study==sup11all$Acronym[i] & is.na(pheno_icogs$Behaviour1) & pheno_icogs$SG_ID %in% icogs_africansamples)
#   sup11all$icogs_ctrl_afr[i]=length(idx1)
#   idx2=which(pheno_icogs$study==sup11all$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$SG_ID %in% icogs_europeansamples)
#   sup11all$icogs_invasive_eur[i]=length(idx2)
#   idx2=which(pheno_icogs$study==sup11all$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$SG_ID %in% icogs_asiansamples)
#   sup11all$icogs_invasive_asn[i]=length(idx2)
#   idx2=which(pheno_icogs$study==sup11all$Acronym[i] & pheno_icogs$Behaviour1==1 & pheno_icogs$SG_ID %in% icogs_africansamples)
#   sup11all$icogs_invasive_afr[i]=length(idx2)
#   
#   if (!sup11all$Acronym[i] %in% "Norway")
#   {
#     idx1=which(pheno_onco$study==sup11all$Acronym[i] & is.na(pheno_onco$Behaviour1) & pheno_onco$Onc_ID %in% onco_europeansamples)
#     sup11all$onco_ctrl_eur[i]=length(idx1)
#     idx1=which(pheno_onco$study==sup11all$Acronym[i] & is.na(pheno_onco$Behaviour1) & pheno_onco$Onc_ID %in% onco_asiansamples)
#     sup11all$onco_ctrl_asn[i]=length(idx1)
#     idx1=which(pheno_onco$study==sup11all$Acronym[i] & is.na(pheno_onco$Behaviour1) & pheno_onco$Onc_ID %in% onco_africansamples)
#     sup11all$onco_ctrl_afr[i]=length(idx1)
#     idx1=which(pheno_onco$study==sup11all$Acronym[i] & is.na(pheno_onco$Behaviour1) & pheno_onco$Onc_ID %in% onco_hispanicsamples)
#     sup11all$onco_ctrl_hsp[i]=length(idx1)
#     idx2=which(pheno_onco$study==sup11all$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Onc_ID %in% onco_europeansamples)
#     sup11all$onco_invasive_eur[i]=length(idx2)
#     idx2=which(pheno_onco$study==sup11all$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Onc_ID %in% onco_asiansamples)
#     sup11all$onco_invasive_asn[i]=length(idx2)
#     idx2=which(pheno_onco$study==sup11all$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Onc_ID %in% onco_africansamples)
#     sup11all$onco_invasive_afr[i]=length(idx2)
#     idx2=which(pheno_onco$study==sup11all$Acronym[i] & pheno_onco$Behaviour1==1 & pheno_onco$Onc_ID %in% onco_hispanicsamples)
#     sup11all$onco_invasive_hsp[i]=length(idx2)
#   }
# }
# sup11all=sup11all[order(sup11all$Acronym),]
# tmp=sup11all
# idx=which(tmp==0,arr.ind = T)
# tmp[idx]=""
# write.csv(tmp,file="../result/supplement750_ancestry.csv",row.names = F,quote=F)


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
# LMBC excluded because inflation was not reduced by the PCs. 
# ICICLE mainly Ductal in situ
# GLACIER mainly lobular invasive

sup11$Acronym[!sup11$Acronym %in% oldsup1$Acronym]
# [1] "ACP"     "GLACIER" "HERPACC" "ICICLE"  "LAABC"   "MYBRCA"  "SBCGS"   "SCCS"    "SEBCS"   "SGBCC"  
# [11] "TWBCS"   "CAMA"    "COLBCCC" "HKBCS"   "KOHBRA"  "NGOBCS"  "UBCS"    "WAABCS" 

idx=which(!sup11$Acronym %in% oldsup1$Acronym)
sum(as.numeric(sup11$icogs_ctrl[idx]),na.rm=T)+sum(as.numeric(sup11$onco_ctrl[idx]),na.rm=T)#20253
sum(as.numeric(sup11$icogs_invasive[idx]),na.rm=T)+sum(as.numeric(sup11$onco_invasive[idx]),na.rm=T) #19843
oldsup1$Acronym[!oldsup1$Acronym %in% sup11$Acronym]
# [1] "BCFR"                                                                                                                                 
# [2] "BPC3"                                                                                                                                 
# [3] "GC-HBOC"                                                                                                                              
# [4] "HEBCS"                                                                                                                                
# [5] "UK2"                                                                                                                                  
# [6] "VUMC"                                                                                                                                 
# [7] "WHI"  
idx=which(!oldsup1$Acronym %in% sup11$Acronym)[1:7]
#HEBCS Finland (1236xx controls, 1796xx cases), WHI USA (4617xx controls, 4930xx cases)
sum(as.numeric(oldsup1$Controls...5[idx]),na.rm=T)+sum(as.numeric(oldsup1$Controls...9[idx]),na.rm=T) #7585
sum(as.numeric(oldsup1$Invasive...6[idx]),na.rm=T)+sum(as.numeric(oldsup1$Invasive...10[idx]),na.rm=T) #10142
sum(as.numeric(oldsup1$Controls...5[1:86]),na.rm=T)+sum(as.numeric(oldsup1$Controls...9[1:86]),na.rm=T) #96201
#sup11=sup11[sup11$Acronym %in% oldsup1$Acronym,]
sup11=sup11[order(sup11$Acronym),]
idx=which(sup11$icogs_ctrl!=""|sup11$icogs_invasive!=""|sup11$onco_ctrl!=""|sup11$onco_invasive!="")
sup11=sup11[idx,]
write.csv(sup11,file="../result/supplement750_1.csv",row.names = F,quote=F)
sum(as.numeric(sup11$icogs_ctrl),na.rm=T)+sum(as.numeric(sup11$icogs_invasive),na.rm=T)
#[1] 78288
sum(as.numeric(sup11$icogs_ctrl),na.rm=T) #39445 control
length(pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym &!pheno_icogs$study %in% study2rm &is.na(pheno_icogs$Behaviour1)]) #39445
sum(as.numeric(sup11$icogs_invasive),na.rm=T) #38843 case
length(pheno_icogs$SG_ID[which(pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm &pheno_icogs$Behaviour1==1)]) #38843
length(pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm]) #78288
#sample used in analysis
icogs_samples=pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm]
write.table(icogs_samples,file="../result/icogs_samples_750.txt",row.names = F,col.names = F,quote=F)
pheno_icogs_final=pheno_icogs[!pheno_icogs$study %in% study2rm,]
dim(pheno_icogs_final) #78288
idx=match(icogs_samples,pheno_icogs$SG_ID)
tmp=pheno_icogs$EthnicityGeno[idx]
icogs_ana_eurosamples=icogs_samples[tmp=="European"]
icogs_ana_asianamples=icogs_samples[tmp=="Asian"]
sum(as.numeric(sup11$onco_ctrl),na.rm=T)+sum(as.numeric(sup11$onco_invasive),na.rm=T) #148051
length(pheno_onco$Onc_ID[pheno_onco$study %in% sup11$Acronym & pheno_onco$StudyCountry!="Norway"])#148051
sum(as.numeric(sup11$onco_ctrl),na.rm=T) #66932 control
sum(as.numeric(sup11$onco_invasive),na.rm=T) #81119
onco_samples=pheno_onco$Onc_ID[pheno_onco$study %in% sup11$Acronym & pheno_onco$StudyCountry!="Norway"]
write.table(onco_samples,file="../result/onco_samples_750.txt",row.names = F,col.names = F,quote=F)
pheno_onco_final=pheno_onco[pheno_onco$StudyCountry!="Norway",]
dim(pheno_onco_final) #148051

#write.table(pheno_icogs_final,file="../data/concept_750_zhang_icogs_pheno_v15_02_age_used.txt",row.names = F,sep="\t",quote=F)
#write.table(pheno_onco_final,file="../data/concept_750_zhang_onco_pheno_v15_02_corrected_age_used.txt",row.names = F,sep="\t",quote=F)
#check correlation among ancestries

allpheno=rbind(pheno_icogs_final[,c("ER_status1","PR_status1","HER2_status1","Grade1","EthnicityGeno")],
               pheno_onco_final[,c("ER_status1","PR_status1","HER2_status1","Grade1","EthnicityGeno")])
idx=which(allpheno==888,arr.ind = T)
allpheno[idx]=NA
idx=which(allpheno$EthnicityGeno=="European")
tmp=allpheno[idx,1:4]
round(cor(tmp,use = "complete.obs"),3)
#              ER_status1 PR_status1 HER2_status1 Grade1
# ER_status1        1.000      0.611       -0.159 -0.390
# PR_status1        0.611      1.000       -0.167 -0.320
# HER2_status1     -0.159     -0.167        1.000  0.202
# Grade1           -0.390     -0.320        0.202  1.000
idx=which(allpheno$EthnicityGeno=="Asian")
tmp=allpheno[idx,1:4]
round(cor(tmp,use = "complete.obs"),3)
#              ER_status1 PR_status1 HER2_status1 Grade1
# ER_status1        1.000      0.634       -0.186 -0.372
# PR_status1        0.634      1.000       -0.169 -0.305
# HER2_status1     -0.186     -0.169        1.000  0.250
# Grade1           -0.372     -0.305        0.250  1.000
idx=which(allpheno$EthnicityGeno=="African")
tmp=allpheno[idx,1:4]
round(cor(tmp,use = "complete.obs"),3)
#              ER_status1 PR_status1 HER2_status1 Grade1
# ER_status1        1.000      0.765       -0.028 -0.526
# PR_status1        0.765      1.000       -0.086 -0.482
# HER2_status1     -0.028     -0.086        1.000  0.126
# Grade1           -0.526     -0.482        0.126  1.000
# idx=which(allpheno$EthnicityGeno=="other")
# tmp=allpheno[idx,1:4]
# round(cor(tmp,use = "complete.obs"),3)
# #              ER_status1 PR_status1 HER2_status1 Grade1
# # ER_status1        1.000      0.795       -0.102 -0.244
# # PR_status1        0.795      1.000       -0.183 -0.274
# # HER2_status1     -0.102     -0.183        1.000  0.052
# # Grade1           -0.244     -0.274        0.052  1.000

#bar plot Within each ancestry, the proportion of ER+ among all cases
pop=c("European","Asian","African")
ancestry=c(rep("EUR",5),rep("ASN",5),rep("AFR",5))
characteristics=c(rep(c("ER+","PR+","HER2+","Grade2","Grade3"),3))
proportion=rep(NA,length(ancestry))
for (i in 1:length(pop))
{
  n=(i-1)*5+1
  idx=which(allpheno$EthnicityGeno==pop[i])
  proportion[n]=sum(allpheno$ER_status1[idx]==1,na.rm = T)/sum(!is.na(allpheno$ER_status1[idx]))
  proportion[n+1]=sum(allpheno$PR_status1[idx]==1,na.rm = T)/sum(!is.na(allpheno$PR_status1[idx]))
  proportion[n+2]=sum(allpheno$HER2_status1[idx]==1,na.rm = T)/sum(!is.na(allpheno$HER2_status1[idx]))
  proportion[n+3]=sum(allpheno$Grade1[idx]==2,na.rm = T)/sum(!is.na(allpheno$Grade1[idx]))
  proportion[n+4]=sum(allpheno$Grade1[idx]==3,na.rm = T)/sum(!is.na(allpheno$Grade1[idx]))
}
data <- data.frame(ancestry,characteristics,proportion)
# Grouped
source("../../PRS_EASLC/code/theme_publication.R")
pdf("../result/characteristics_ancestry_barplot.pdf",width = 12)
ggplot(data, aes(fill=ancestry, y=proportion, x=characteristics)) + 
  geom_bar(position="dodge", stat="identity")+theme_Publication()
dev.off()
pdf("../result/characteristics_ancestry_barplot1.pdf",width = 12)
ggplot(data, aes(fill=characteristics, y=proportion, x=ancestry)) + 
  geom_bar(position="dodge", stat="identity")+theme_Publication()
dev.off()

idx=match(onco_samples,pheno_onco$Onc_ID)
tmp=pheno_onco$EthnicityGeno[idx]
onco_ana_eurosamples=onco_samples[tmp=="European"]
onco_ana_asiansamples=onco_samples[tmp=="Asian"]
# icogs_samples=read.table("../result/icogs_samples_750.txt")$V1
# onco_samples=read.table("../result/onco_samples_750.txt")$V1
# 
# p1=ggplot(pheno_onco[pheno_onco$Onc_ID %in% onco_samples,], aes(x=PC_1, y=PC_2, color=EthnicityGeno)) + geom_point(alpha=0.1,size=2) +
#   labs(title="") +
#   theme(
#     panel.background = element_rect(fill='transparent'),
#     plot.background = element_rect(fill='transparent', color=NA),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill='transparent'),
#     legend.box.background = element_rect(fill='transparent')
#   )

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

sup22=sup22[sup22$Acronym %in% sup11$Acronym,]
sup22=sup22[order(sup22$Acronym),]
write.csv(sup22,file="../result/supplement750_2.csv",row.names = F,quote=F)
sup22=read.csv("../result/supplement750_2.csv")
sup222=data.frame(marker=c(rep(c("ER","PR","HER2","Grade"),each=3),"Grade"),status=c(rep(c("Negative","Positive","Unkown"),3),c("1","2","3","Unknown")),N=NA,prop=NA)
sup222$N[1]=sum(as.numeric(sup22$icogs_ERN),na.rm=T)+sum(as.numeric(sup22$onco_ERN),na.rm=T)
sup222$N[2]=sum(as.numeric(sup22$icogs_ERP),na.rm=T)+sum(as.numeric(sup22$onco_ERP),na.rm=T)
sup222$prop[1]=round(100*(sup222$N[1]/(sup222$N[1]+sup222$N[2])))
sup222$prop[2]=100-sup222$prop[1]
sup222$N[3]=sum(as.numeric(sup22$icogs_ERU),na.rm=T)+sum(as.numeric(sup22$onco_ERU),na.rm=T)
sup222$N[4]=sum(as.numeric(sup22$icogs_PRN),na.rm=T)+sum(as.numeric(sup22$onco_PRN),na.rm=T)
sup222$N[5]=sum(as.numeric(sup22$icogs_PRP),na.rm=T)+sum(as.numeric(sup22$onco_PRP),na.rm=T)
sup222$prop[4]=round(100*(sup222$N[4]/(sup222$N[4]+sup222$N[5])))
sup222$prop[5]=100-sup222$prop[4]
sup222$N[6]=sum(as.numeric(sup22$icogs_PRU),na.rm=T)+sum(as.numeric(sup22$onco_PRU),na.rm=T)
sup222$N[7]=sum(as.numeric(sup22$icogs_HER2N),na.rm=T)+sum(as.numeric(sup22$onco_HER2N),na.rm=T)
sup222$N[8]=sum(as.numeric(sup22$icogs_HER2P),na.rm=T)+sum(as.numeric(sup22$onco_HER2P),na.rm=T)
sup222$prop[7]=round(100*(sup222$N[7]/(sup222$N[7]+sup222$N[8])))
sup222$prop[8]=100-sup222$prop[7]
sup222$N[9]=sum(as.numeric(sup22$icogs_HER2U),na.rm=T)+sum(as.numeric(sup22$onco_HER2U),na.rm=T)
sup222$N[10]=sum(as.numeric(sup22$icogs_grade1),na.rm=T)+sum(as.numeric(sup22$onco_grade1),na.rm=T)
sup222$N[11]=sum(as.numeric(sup22$icogs_grade2),na.rm=T)+sum(as.numeric(sup22$onco_grade2),na.rm=T)
sup222$N[12]=sum(as.numeric(sup22$icogs_grade3),na.rm=T)+sum(as.numeric(sup22$onco_grade3),na.rm=T)
sup222$prop[10]=round(100*(sup222$N[10]/(sup222$N[10]+sup222$N[11]+sup222$N[12])))
sup222$prop[11]=round(100*(sup222$N[11]/(sup222$N[10]+sup222$N[11]+sup222$N[12])))
sup222$prop[12]=100-sup222$prop[10]-sup222$prop[11]
sup222$N[13]=sum(as.numeric(sup22$icogs_gradeU),na.rm=T)+sum(as.numeric(sup22$onco_gradeU),na.rm=T)
write.csv(sup222,file="../result/supplement750_22.csv",row.names = F,quote=F)

#other is Hispanic
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
write.csv(tmp,file="../result/supplement750_4.csv",row.names = F,quote=F)

sup44_euro=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="European",],
                 myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="European",])
write.csv(sup44_euro,file="../result/supplement750_4_euro.csv",row.names = F,quote=F)
sup44_asian=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="Asian",],
                 myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="Asian",])
write.csv(sup44_asian,file="../result/supplement750_4_asian.csv",row.names = F,quote=F)
sup44_african=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="African",],
                 myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="African",])
write.csv(sup44_african,file="../result/supplement750_4_african.csv",row.names = F,quote=F)
# sup44_other=sup44(myphenoicogs = pheno_icogs_final[pheno_icogs_final$EthnicityGeno=="other",],
#                  myphenoonco = pheno_onco_final[pheno_onco_final$EthnicityGeno=="other",])
# write.csv(sup44_other,file="../result/supplement750_4_other.csv",row.names = F,quote=F)
#sup44_all=rbind(sup44_euro,sup44_asian,sup44_african,sup44_other)
sup44_all=rbind(sup44_euro,sup44_asian,sup44_african)
colnames(sup44_all)[which(colnames(sup44_all)=="HRP_HER2N_lowgrade")]="Luminal-A-like"
colnames(sup44_all)[which(colnames(sup44_all)=="HRP_HER2N_highgrade")]="HER2-negative-like"
colnames(sup44_all)[which(colnames(sup44_all)=="HRP_HER2P")]="Luminal-B-like"
colnames(sup44_all)[which(colnames(sup44_all)=="HRN_HER2N")]="Triple-negative"
colnames(sup44_all)[which(colnames(sup44_all)=="HRN_HER2P")]="HER2-enriched-like"
tmp=data.frame(matrix(NA,nrow=4,ncol=ncol(sup44_all)))
colnames(tmp)=colnames(sup44_all)
for (i in 1:nrow(tmp))
{
  tmp[i,]=sup44_all[((i-1)*2+1),]+sup44_all[((i-1)*2+2),]
}
sup44_all=tmp
#bar plot Within each ancestry, the proportion of subtypes among all cases
pop=c("European","Asian","African")
ancestry=c(rep("EUR",5),rep("ASN",5),rep("AFR",5))
subtype=c(rep(colnames(sup44_all),length(pop)))
proportion=rep(NA,length(ancestry))
for (i in 1:length(pop))
{
  n=(i-1)*5+1
  proportion[n]=sup44_all$`Luminal-A-like`[i]/sum(sup44_all[i,])
  proportion[n+1]=sup44_all$`Luminal-B-like`[i]/sum(sup44_all[i,])
  proportion[n+2]=sup44_all$`HER2-negative-like`[i]/sum(sup44_all[i,])
  proportion[n+3]=sup44_all$`HER2-enriched-like`[i]/sum(sup44_all[i,])
  proportion[n+4]=sup44_all$`Triple-negative`[i]/sum(sup44_all[i,])
}
data <- data.frame(ancestry,subtype,proportion)
# Grouped
pdf("../result/subtype_ancestry_barplot.pdf",width=12)
ggplot(data, aes(fill=subtype, y=proportion, x=ancestry)) + 
  geom_bar(position="dodge", stat="identity")+theme_Publication()
dev.off()

pdf("../result/subtype_ancestry_barplot1.pdf",width=12)
ggplot(data, aes(fill=ancestry, y=proportion, x=subtype)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_Publication()+theme(axis.text.x = element_text(size=15,angle = 90, vjust = 0.5, hjust=1))
dev.off()

sup44=read.csv("../result/supplement750_4.csv")
colnames(sup44)[which(colnames(sup44)=="HRP_HER2N_lowgrade")]="Luminal-A-like"
colnames(sup44)[which(colnames(sup44)=="HRP_HER2N_highgrade")]="HER2-negative-like"
colnames(sup44)[which(colnames(sup44)=="HRP_HER2P")]="Luminal-B-like"
colnames(sup44)[which(colnames(sup44)=="HRN_HER2N")]="Triple-negative"
colnames(sup44)[which(colnames(sup44)=="HRN_HER2P")]="HER2-enriched-like"
sup444=data.frame(subtype=c("Luminal-A-like","HER2-negative-like","Luminal-B-like","HER2-enriched-like","Triple-negative"),
                  N=NA,prop=0)
idx=match(sup444$subtype,colnames(sup44))
sup444$N=colSums(sup44[,idx])
for (i in 1:4) sup444$prop[i]=round(sup444$N[i]/sum(sup444$N)*100)
sup444$prop[5]=100-sum(sup444$prop[1:3])
write.csv(sup444,file="../result/supplement750_44.csv",row.names = F,quote=F)

tmp=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_euro_icogs_topmed_plinksamples.txt"))
tmp=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_asian_icogs_topmed_plinksamples.txt"))
tmp=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_topmed_plinksamples.txt"))
tmp=as.data.frame(fread("/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_asian_750_topmed_plinksamples.txt"))

tmp1=unlist(strsplit(tmp$IID,"_"))
tmp1=tmp1[seq(1,length(tmp1),2)]
table(icogs_ana_eurosamples %in% tmp1)
table(icogs_ana_asianamples %in% tmp1)
table(onco_ana_eurosamples %in% tmp1)
table(onco_ana_asiansamples %in% tmp1)

tmp=read.table("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro_mergelist.txt")
allpvar=NULL
tmp1=data.frame(nsample=rep(0,nrow(tmp)),nvar=0,ndupsnp=0)
for(i in 1:nrow(tmp1))
{
  if (i%%10==0) cat(i,'..')
  tmp2=as.data.frame(fread(paste0(tmp$V1[i],".pvar")))
  allpvar=rbind(allpvar,tmp2)
  tmp1$nvar[i]=nrow(tmp2)
  tmp1$ndupsnp[i]=sum(duplicated(tmp2$ID))
  tmp2=read.table(paste0(tmp$V1[i],".psam"))
  tmp1$nsample[i]=nrow(tmp2)
}

tmp=as.data.frame(fread("test1.vcf.gz",nrows=20))
tmp1=unname(unlist(tmp[1,10:ncol(tmp)]))
idx=which(grepl(":$",tmp1))
tmp1[idx]=paste0(tmp1[idx],"0.999")
tmp[1,10:ncol(tmp)]=tmp1
tmp2=unlist(strsplit(tmp1,":"))
tmp3=tmp2[seq(1,length(tmp2),2)]
tmp4=tmp2[seq(2,length(tmp2),2)]
quantile(as.numeric(tmp4))
tmp[1,10]="1/1:."
write.table(tmp[1:1,],file="test.vcf",row.names = F,sep="\t",quote=F)

tmp=as.data.frame(fread("../result/imp_onco/euro/euro_chr.1.1.pvar"))
library(GenomicRanges)
check_ld_block=function(infolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/geno/",
                        pvarfile="../result/imp_icogs/euro/euro.pvar")
{
  ldblock=read.table("/data/BB_Bioinformatics/ProjectData/ld_block/ld_block_EUR",header=T)
  pvar=as.data.frame(fread(pvarfile))
  colnames(pvar)[1]="CHR"
  
  gr_pvar=GRanges(seqnames = pvar$CHR,ranges=IRanges(start=pvar$POS,width = 1))
  gr_ldblock=GRanges(seqnames = gsub("chr","",ldblock$chr),ranges=IRanges(ldblock$start,ldblock$stop-1))
  ldblock$numvar=NA
  for(i in 1:nrow(ldblock))
  {
    ldblock$numvar[i]=length(subsetByOverlaps(gr_pvar,gr_ldblock[i]))
  }
  quantile(ldblock$numvar)
  # 0%   25%   50%   75%  100% 
  # 0  4409  5871  7563 16965
  print(sum(ldblock$numvar==0)) #8
  which.max(ldblock$numvar) #966
  which.max(ldblock$len) #965
  allfiles=list.files(infolder,"*.traw.gz")
  allids=gsub("geno_","",allfiles)
  allids=as.numeric(gsub(".traw.gz","",allids,fixed=T))
  tmp=0:(nrow(ldblock)-1)
  missingids=tmp[!tmp %in% allids]
  print(length(missingids))
  print(ldblock$numvar[missingids+1])
  return(ldblock)
}
check_ld_block_onco=check_ld_block(infolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/geno/",
                        pvarfile="../result/imp_onco/euro/euro.pvar")

check_geno_block=function(prefix="../result/imp_onco/euro/euro")
{
  nvar=2000
  pvar=as.data.frame(fread(paste0(prefix,".pvar")))
  infolder=paste0(dirname(prefix),"/geno/")
  allfiles=list.files(infolder,"\\w+.traw.gz")
  allposfiles=list.files(infolder,"pos_\\d+.txt")
  allidxs=gsub("pos_","",allposfiles)
  allidxs=as.numeric(gsub(".txt","",allidxs,fixed = T))
  allidxs=allidxs[order(allidxs)]
  n=as.integer(nrow(pvar)/nvar)+1
  tmp=1:n
  missingblocks=tmp[!tmp %in% allidxs]
  if (length(missingblocks)>0 ) warning(paste0(paste0(missingblocks,collapse = ",")," are missing"))
  allsnps=NULL
  for (i in 1:length(allidxs))
  {
    if (i %% 500==0) cat(i,'..')
    tmp=as.data.frame(fread(paste0(infolder,"pos_",allidxs[i],".txt"),header=F))
    allsnps=rbind(allsnps,tmp)
  }
  if (sum(duplicated(allsnps[,1]))>0) warning("duplicateds variants found")
  if (nrow(allsnps)!=nrow(pvar)) warning("some variants are missing")
}
check_geno_block(prefix="../result/imp_onco/euro/euro")
check_geno_block(prefix="../result/imp_icogs/euro/euro")

#PLOT PCA on iCOGS
pdf("../result/iCOGS_PCA.pdf")
ggplot(pheno_icogs, aes(x=pc1, y=pc2, color=EthnicityGeno)) + geom_point() + 
  labs(title="iCOGS")+ theme_Publication()
dev.off()
tmp=read.table("../result/imp_icogs/merged1.eigenvec")
colnames(tmp)=c("ID",paste0("plinkPC",1:20))
tmp1=unlist(strsplit(tmp$ID,"_"))
tmp$ID=tmp1[seq(1,length(tmp1),2)]
idx=match(pheno_icogs$SG_ID,tmp$ID)
dat=cbind(pheno_icogs,tmp[idx,2:21])
pdf("../result/iCOGS_plinkPCA.pdf")
ggplot(dat, aes(x=plinkPC1, y=plinkPC2, color=EthnicityGeno)) + geom_point() + 
  labs(title="iCOGS")+ theme_Publication()
dev.off()
#flashpca2
load("../result/icogs_flashpca.RData")
tmp=icogspca$projection
colnames(tmp)=paste0("flashPC",1:10)
tmp1=unlist(strsplit(rownames(tmp),":"))
tmp1=unlist(strsplit(tmp1,"_"))
rownames(tmp)=tmp1[seq(2,length(tmp1),3)]
idx=match(pheno_icogs_final$SG_ID,rownames(tmp))
dat=cbind(pheno_icogs_final,tmp[idx,])
pdf("../result/iCOGS_flashPCA.pdf")
ggplot(dat, aes(x=flashPC1, y=flashPC2, color=EthnicityGeno)) + geom_point() + 
  labs(title="iCOGS")
dev.off()

pdf("../result/Onco_PCA.pdf")
ggplot(pheno_onco_final, aes(x=PC_1, y=PC_2, color=EthnicityGeno)) + geom_point() + 
  labs(title="Onco")
dev.off()
