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
onco_asiansamples=get_sampelnames(paste0(basedir,"onco/",onco_asian[1])) #27515
onco_africansamples=get_sampelnames(paste0(basedir,"onco/",onco_african[1])) #5804
onco_hispanicsamples=get_sampelnames(paste0(basedir,"onco/",onco_hispanic[1])) #2560
onco_europeansamples=get_sampelnames(paste0(basedir,"onco/",onco_european[1])) #133236
onco_gensamples=c(onco_asiansamples,onco_africansamples,onco_hispanicsamples,onco_europeansamples)
onco_gensamples=onco_gensamples[onco_gensamples!=0]
length(onco_gensamples) #169115

#icogs
tmp=list.files(paste0(basedir,"icogs/"),"*.sample")
icogs_asian=tmp[grepl("asian",tmp)]
icogs_african=tmp[grepl("african",tmp)]
icogs_european=tmp[grepl("euro",tmp)]
length(icogs_asian)+length(icogs_african)+length(icogs_european)

icogs_asiansamples=get_sampelnames(paste0(basedir,"icogs/",icogs_asian[1])) #12500, first row is 0 0
icogs_africansamples=get_sampelnames(paste0(basedir,"icogs/",icogs_african[1])) #2047
icogs_europeansamples=get_sampelnames(paste0(basedir,"icogs/",icogs_european[1])) #99000

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
#select case/control samples and create pheno file
update_samplefile=function(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_euro_icogs_topmed_1_p1.sample",
                           prefix="zhang_750_euro_icogs_topmed_",
                           pheno=pheno_icogs)
{
  outfolder=dirname(samplefile)
  sampledat=read.table(samplefile,header = T,sep=" ")
  #sampledat$SEX="D" #Sex code ('1' = male, '2' = female, '0' = unknown)
  #sampledat$pheno="D" #	Binary ('0' = control, '1' = case) discrete (categorical, positive integers), or continuous phenotype; missing values represented by 'NA'
  sampledat=sampledat[2:nrow(sampledat),]
  sampleid=unlist(strsplit(sampledat$ID_1,"_"))
  sampleid=sampleid[seq(1,length(sampleid),2)]
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
  sampledat=sampledat[match(sampleid,sampleid),]
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
update_samplefile()
update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_asian_icogs_topmed_1.sample",
                           prefix="zhang_750_asian_icogs_topmed_",
                           pheno=pheno_icogs)
update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/zhang_750_african_icogs_topmed_1.sample",
                  prefix="zhang_750_african_icogs_topmed_",
                  pheno=pheno_icogs)

update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_topmed_1_p1.sample",
                  prefix="zhang_750_topmed_",
                  pheno=pheno_onco)
#53862 (0) 66012 (1), 5227 extra genotype
update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_asian_topmed_1.sample",
                  prefix="zhang_asian_750_topmed_",
                  pheno=pheno_onco)
#12578 (0) 12905 (1)
update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_african_topmed_1.sample",
                  prefix="zhang_african_750_topmed_",
                  pheno=pheno_onco)
#2088 (0) 3483 (1)
update_samplefile(samplefile="/data/BB_Bioinformatics/ProjectData/BCAC/onco/zhang_750_hispanic_topmed_1.sample",
                  prefix="zhang_hispanic_750_topmed_",
                  pheno=pheno_onco)
#1218 (0) 1196 (1)
#only include controls and invasives Behaviour1(1,888,NA):NAs are controls
pheno_icogs_749=pheno_icogs_749[which(is.na(pheno_icogs_749$Behaviour1)| pheno_icogs_749$Behaviour1==1),]
pheno_onco_749=pheno_onco_749[which(is.na(pheno_onco_749$Behaviour1)| pheno_onco_749$Behaviour1==1),]
pheno_icogs=pheno_icogs[which(is.na(pheno_icogs$Behaviour1)| pheno_icogs$Behaviour1==1),]
pheno_onco=pheno_onco[which(is.na(pheno_onco$Behaviour1)| pheno_onco$Behaviour1==1),]


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
table(pheno_onco$Onc_ID %in% onco_gensamples) #151615 pheno not have genotype data
# FALSE   TRUE 
# 7887 151615 
pheno_onco=pheno_onco[pheno_onco$Onc_ID %in% onco_gensamples,]
dim(pheno_onco)
#151615     47


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
    if(length(idx)>0) sup11$icogs_ancestry[i]=pheno_icogs$EthnicityGeno[idx[1]]
    idx=which(pheno_icogs$study==sup11$Acronym[i] & is.na(pheno_icogs$Behaviour1))
    if (length(idx)>0) sup11$icogs_ctrl[i]=length(idx)
    n1=length(idx)
    idx=which(pheno_icogs$study==sup11$Acronym[i] & pheno_icogs$Behaviour1==1)
    if (length(idx)>0) sup11$icogs_invasive[i]=length(idx)
    n2=length(idx)
    if(n!=n1+n2) stop(i)
    idx=which(pheno_onco$study==sup11$Acronym[i])
    if(length(idx)>0) sup11$onco_ancestry[i]=pheno_onco$EthnicityGeno[idx[1]]
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
#[1] 80046
sum(as.numeric(sup11$icogs_ctrl),na.rm=T) #40262 control
length(pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym &!pheno_icogs$study %in% study2rm &is.na(pheno_icogs$Behaviour1)]) #40262
sum(as.numeric(sup11$icogs_invasive),na.rm=T) #39784 case
length(pheno_icogs$SG_ID[which(pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm &pheno_icogs$Behaviour1==1)]) #39784
length(pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm]) #80046
#sample used in analysis
icogs_samples=pheno_icogs$SG_ID[pheno_icogs$study %in% sup11$Acronym & !pheno_icogs$study %in% study2rm]
write.table(icogs_samples,file="../result/icogs_samples_750.txt",row.names = F,col.names = F,quote=F)
sum(as.numeric(sup11$onco_ctrl),na.rm=T)+sum(as.numeric(sup11$onco_invasive),na.rm=T) #150464
length(pheno_onco$Onc_ID[pheno_onco$study %in% sup11$Acronym & pheno_onco$StudyCountry!="Norway"])#150464
sum(as.numeric(sup11$onco_ctrl),na.rm=T) #68150 control
sum(as.numeric(sup11$onco_invasive),na.rm=T) #82314
onco_samples=pheno_onco$Onc_ID[pheno_onco$study %in% sup11$Acronym & pheno_onco$StudyCountry!="Norway"]
write.table(onco_samples,file="../result/onco_samples_750.txt",row.names = F,col.names = F,quote=F)
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
sup22=sup22[!sup22$Acronym %in% study2rm,]
sup22=sup22[sup22$Acronym %in% sup11$Acronym,]
sup22=sup22[order(sup22$Acronym),]
write.csv(sup22,file="../result/supplement750_2.csv",row.names = F,quote=F)

#other is Hispanic
sup33=data.frame(ancestry=c("European","Asian","African","other"),icogs_ctrl="",icogs_invasive="",onco_ctrl="",onco_invasive="")
for (i in 1:4)
{
  idx=which(sup11$icogs_ancestry==sup33$ancestry[i])
  sup33$icogs_ctrl[i]=sum(as.numeric(sup11$icogs_ctrl[idx]),na.rm=T)
  sup33$icogs_invasive[i]=sum(as.numeric(sup11$icogs_invasive[idx]),na.rm=T)
  idx=which(sup11$onco_ancestry==sup33$ancestry[i])
  sup33$onco_ctrl[i]=sum(as.numeric(sup11$onco_ctrl[idx]),na.rm=T)
  sup33$onco_invasive[i]=sum(as.numeric(sup11$onco_invasive[idx]),na.rm=T)
}
write.csv(sup33,file="../result/supplement750_3.csv",row.names = F,quote=F)

#subtype sample size
# we could define five intrinsic breast cancer subtypes based on the four tumor characteristics: ER, PR, HER2,
# grade. The five intrinsic subtypes are: 1. (ER or PR)+, HER2-, grade 1 or 2; 
# 2. (ER or PR)+, HER2+; # 3. (ER or PR)+, HER2-, grade 3; 
# 4. (ER & PR)-, HER2-; 5. ER-PR-HER2-.


sup44=data.frame(matrix(NA,nrow=2,ncol=5))
rownames(sup44)=c("icogs","onco")
colnames(sup44)=c("HRP_HER2N_lowgrade","HRP_HER2P","HRP_HER2N_highgrade","HRN_HER2P","HRN_HER2N")
idx=which((pheno_icogs1$ER_status1==1|pheno_icogs1$PR_status1==1) & pheno_icogs1$HER2_status1==0 & pheno_icogs1$Grade1 %in% c(1,2))
sup44[1,1]=length(idx)
idx=which((pheno_icogs1$ER_status1==1|pheno_icogs1$PR_status1==1) & pheno_icogs1$HER2_status1==1)
sup44[1,2]=length(idx)
idx=which((pheno_icogs1$ER_status1==1|pheno_icogs1$PR_status1==1) & pheno_icogs1$HER2_status1==0 & pheno_icogs1$Grade1 %in% c(3))
sup44[1,3]=length(idx)
idx=which((pheno_icogs1$ER_status1==0 & pheno_icogs1$PR_status1==0) & pheno_icogs1$HER2_status1==1)
sup44[1,4]=length(idx)
idx=which((pheno_icogs1$ER_status1==0 & pheno_icogs1$PR_status1==0) & pheno_icogs1$HER2_status1==0)
sup44[1,5]=length(idx)

idx=which((pheno_onco1$ER_status1==1|pheno_onco1$PR_status1==1) & pheno_onco1$HER2_status1==0 & pheno_onco1$Grade1 %in% c(1,2))
sup44[2,1]=length(idx)
idx=which((pheno_onco1$ER_status1==1|pheno_onco1$PR_status1==1) & pheno_onco1$HER2_status1==1)
sup44[2,2]=length(idx)
idx=which((pheno_onco1$ER_status1==1|pheno_onco1$PR_status1==1) & pheno_onco1$HER2_status1==0 & pheno_onco1$Grade1 %in% c(3))
sup44[2,3]=length(idx)
idx=which((pheno_onco1$ER_status1==0 & pheno_onco1$PR_status1==0) & pheno_onco1$HER2_status1==1)
sup44[2,4]=length(idx)
idx=which((pheno_onco1$ER_status1==0 & pheno_onco1$PR_status1==0) & pheno_onco1$HER2_status1==0)
sup44[2,5]=length(idx)
write.csv(sup44,file="../result/supplement750_4.csv",row.names = F,quote=F)

#check genotype data
tmp=as.data.frame(fread('test.traw'))
tmp=tmp[,1:1007]
tmp1=unname(unlist(tmp[,7:1007]))
quantile(tmp1,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
checkgeno=function(hardcall=0.1)
{
  idx=which(tmp1>hardcall & tmp1<1-hardcall | tmp1>1+hardcall & tmp1< 2-hardcall)
  res=length(idx)/length(tmp1)
  print(res)
}
checkgeno(hardcall = 0.1) #17.89%
checkgeno(hardcall = 0.2) #10.78%
checkgeno(hardcall = 0.3) #6.43%
checkgeno(hardcall = 0.4) #3.03%
tmp2=as.data.frame(fread('test1.traw',nrows=100))
