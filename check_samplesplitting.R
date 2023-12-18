#!/usr/bin/env Rscript

# Split the EUR and EAS data into training, tuning, validation. iCOGs will be kept as training, OncoArray will be used as training, tuning, validation. The proportion of training, tuning, validation is roughly 70%, 15%, 15% within EUR, EAS. All of AFR and Latino will be used as validation.
#first, in EUR/EAS onco pick tuning and validation
#second, put African into validation
#third, the rest go to training

library(data.table)
#samples with age/casecontrol/genotype
pheno_icogs0=as.data.frame(fread("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T))
pheno_onco0=as.data.frame(fread("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt"))

for (v in c("ER_status1","PR_status1","HER2_status1","Grade1"))
{
  idx=which(pheno_icogs0[,v]==888)
  pheno_icogs0[idx,v]=NA
  idx=which(pheno_onco0[,v]==888)
  pheno_onco0[idx,v]=NA
}
#keep AFR and Hispanic
#only remove Norway data in onco
study2rm=c("ICICLE","GLACIER","LMBC") #remove in icogs
# LMBC excluded because inflation was not reduced by the PCs. 
# ICICLE mainly Ductal in situ
# GLACIER mainly lobular invasive

table(pheno_onco0$StudyCountry)
pheno_onco0=pheno_onco0[pheno_onco0$StudyCountry!="Norway",]
table(pheno_onco0$Behaviour1,pheno_onco0$EthnicityGeno,useNA = "ifany")
#       African Asian European other
# 1       3481 12870    64768  1195
# <NA>    2088 12344    52500  1218
idx=complete.cases(pheno_onco0[,c("ER_status1","PR_status1","HER2_status1","Grade1")])
table(idx)
# FALSE  TRUE 
# 110191  40273
pheno_onco0$complete=idx
table(idx,pheno_onco0$EthnicityGeno)
#        African Asian European other
# FALSE    4192 16870    86982  2147
# TRUE     1377  8344    30286   266

pheno_icogs0=pheno_icogs0[!pheno_icogs0$study %in% study2rm,]
table(pheno_icogs0$Behaviour1,pheno_icogs0$EthnicityGeno,useNA = "ifany")
#      African Asian European
# 1        941  4882    33961
# <NA>     817  5940    33505
idx=complete.cases(pheno_icogs0[,c("ER_status1","PR_status1","HER2_status1","Grade1")])
table(idx)
# FALSE  TRUE 
# 67039 13007
table(idx,pheno_icogs0$EthnicityGeno)
#      African Asian European
# FALSE    1639  9442    55958
# TRUE      119  1380    11508

extract_tuningvalidation=function(pop="European")
{
  set.seed(123)
  ntotalcase=nrow(pheno_icogs0[which(pheno_icogs0$EthnicityGeno==pop & pheno_icogs0$Behaviour1==1),])+nrow(pheno_onco0[which(pheno_onco0$EthnicityGeno==pop & pheno_onco0$Behaviour1==1),])
  ntotalcontr=nrow(pheno_icogs0[which(pheno_icogs0$EthnicityGeno==pop & is.na(pheno_icogs0$Behaviour1)),])+nrow(pheno_onco0[which(pheno_onco0$EthnicityGeno==pop & is.na(pheno_onco0$Behaviour1)),])
  
  #number of cases/controls drawn for tuning/validation 
  ncase=as.integer(ntotalcase*0.3)
  ncontr=as.integer(ntotalcontr*0.3)
  
  pheno_onco1=pheno_onco0[which(pheno_onco0$EthnicityGeno==pop),]
  #complete cases
  idx=complete.cases(pheno_onco1[,c("Behaviour1","ER_status1","PR_status1","HER2_status1","Grade1")])
  pheno_onco2=pheno_onco1[idx,]
  #controls
  pheno_onco3=pheno_onco1[is.na(pheno_onco1$Behaviour1),]
  pheno_onco_tvcases=sample(1:nrow(pheno_onco2),ncase)
  pheno_onco_tvcases=pheno_onco2[pheno_onco_tvcases,]
  pheno_onco_tvcontrs=sample(1:nrow(pheno_onco3),ncontr)
  pheno_onco_tvcontrs=pheno_onco3[pheno_onco_tvcontrs,]
  n=as.integer(nrow(pheno_onco_tvcases)/2)
  pheno_onco_tuningcases=pheno_onco_tvcases[sample(1:nrow(pheno_onco_tvcases),n),]
  pheno_onco_validationcases=pheno_onco_tvcases[!pheno_onco_tvcases$Onc_ID %in% pheno_onco_tuningcases$Onc_ID,]
  n=as.integer(nrow(pheno_onco_tvcontrs)/2)
  pheno_onco_tuningcontrs=pheno_onco_tvcontrs[sample(1:nrow(pheno_onco_tvcontrs),n),]
  pheno_onco_validationcontrs=pheno_onco_tvcontrs[!pheno_onco_tvcontrs$Onc_ID %in% pheno_onco_tuningcontrs$Onc_ID,]
  pheno_onco_tuning=rbind(pheno_onco_tuningcases,pheno_onco_tuningcontrs)
  pheno_onco_validation=rbind(pheno_onco_validationcases,pheno_onco_validationcontrs)
  return(list(pheno_onco_tuning=pheno_onco_tuning,pheno_onco_validation=pheno_onco_validation))
}

EUR_tuningvalidation=extract_tuningvalidation()
EAS_tuningvalidation=extract_tuningvalidation(pop="Asian")

tuning=rbind(EUR_tuningvalidation$pheno_onco_tuning,EAS_tuningvalidation$pheno_onco_tuning)
dim(tuning)#33113
validation=rbind(EUR_tuningvalidation$pheno_onco_validation,EAS_tuningvalidation$pheno_onco_validation)
dim(validation) #33116
#others go to validation
pheno_icogs_other=pheno_icogs0[pheno_icogs0$EthnicityGeno=="African",]
pheno_onco_other=pheno_onco0[pheno_onco0$EthnicityGeno %in% c("African","other"),]
validation_onco=rbind(validation,pheno_onco_other)
dim(validation_onco) #41098
validation_icogs=pheno_icogs_other
dim(validation_icogs) #1758
#the rest go to training
training_onco=pheno_onco0[!pheno_onco0$Onc_ID %in% c(tuning$Onc_ID,validation_onco$Onc_ID),]
dim(training_onco) #76253
table(training_onco$EthnicityGeno)
# Asian European 
# 14404    61849
training_icogs=pheno_icogs0[!pheno_icogs0$SG_ID %in% validation_icogs$SG_ID,]
dim(training_icogs) #78288
table(training_icogs$EthnicityGeno)
# Asian European 
# 10822    67466 
write.table(training_icogs,file="../result/PRS_icogs_trainingpheno.txt",row.names = F,sep="\t",quote = F)
write.table(training_onco,file="../result/PRS_onco_trainingpheno.txt",row.names = F,sep="\t",quote = F)
nrow(training_icogs)+nrow(training_onco)+nrow(tuning)+nrow(validation_icogs)+nrow(validation_onco)==nrow(pheno_icogs0)+nrow(pheno_onco0)

nrow(training_icogs[which(training_icogs$Behaviour1==1),])+nrow(training_onco[which(training_onco$Behaviour1==1),]) #81538
nrow(training_icogs[which(is.na(training_icogs$Behaviour1)),])+nrow(training_onco[which(is.na(training_onco$Behaviour1)),]) #73003

nrow(training_icogs[which(training_icogs$EthnicityGeno=="European" & training_icogs$Behaviour1==1),])+nrow(training_onco[which(training_onco$EthnicityGeno=="European" &training_onco$Behaviour1==1),]) #69111
nrow(training_icogs[which(training_icogs$EthnicityGeno=="European" & is.na(training_icogs$Behaviour1)),])+nrow(training_onco[which(training_onco$EthnicityGeno=="European" &is.na(training_onco$Behaviour1)),]) #60204

nrow(training_icogs[which(training_icogs$EthnicityGeno=="Asian" & training_icogs$Behaviour1==1),])+nrow(training_onco[which(training_onco$EthnicityGeno=="Asian" &training_onco$Behaviour1==1),]) #12427
nrow(training_icogs[which(training_icogs$EthnicityGeno=="Asian" & is.na(training_icogs$Behaviour1)),])+nrow(training_onco[which(training_onco$EthnicityGeno=="Asian" &is.na(training_onco$Behaviour1)),]) #12799

nrow(tuning[which(tuning$Behaviour1==1),]) #17471
nrow(tuning[which(is.na(tuning$Behaviour1)),]) #15642

nrow(tuning[which(tuning$Behaviour1==1 & tuning$EthnicityGeno=="European"),]) #14809
nrow(tuning[which(is.na(tuning$Behaviour1 & tuning$EthnicityGeno=="European")),]) #12900

nrow(tuning[which(tuning$Behaviour1==1 & tuning$EthnicityGeno=="Asian"),]) #2662
nrow(tuning[which(is.na(tuning$Behaviour1 & tuning$EthnicityGeno=="Asian")),]) #2742

table(validation_icogs$EthnicityGeno)
# African 
# 1758
table(validation_onco$EthnicityGeno)
# African    Asian European    other 
# 5569     5406    27710     2413 
nrow(validation_icogs[which(validation_icogs$Behaviour1==1),])+nrow(validation_onco[which(validation_onco$Behaviour1==1),]) #23089
nrow(validation_icogs[which(is.na(validation_icogs$Behaviour1)),])+nrow(validation_onco[which(is.na(validation_onco$Behaviour1)),]) #19767

nrow(validation_icogs[which(validation_icogs$Behaviour1==1),])+nrow(validation_onco[which(validation_onco$EthnicityGeno=="African" & validation_onco$Behaviour1==1),]) #4422
nrow(validation_icogs[which(is.na(validation_icogs$Behaviour1)),])+nrow(validation_onco[which(validation_onco$EthnicityGeno=="African" & is.na(validation_onco$Behaviour1)),]) #2905

nrow(validation_onco[which(validation_onco$EthnicityGeno=="Asian" & validation_onco$Behaviour1==1),]) #2663
nrow(validation_onco[which(validation_onco$EthnicityGeno=="Asian" & is.na(validation_onco$Behaviour1)),]) #2743

nrow(validation_onco[which(validation_onco$EthnicityGeno=="European" & validation_onco$Behaviour1==1),]) #14809
nrow(validation_onco[which(validation_onco$EthnicityGeno=="European" & is.na(validation_onco$Behaviour1)),]) #12901

nrow(validation_onco[which(validation_onco$EthnicityGeno=="other" & validation_onco$Behaviour1==1),]) #1195
nrow(validation_onco[which(validation_onco$EthnicityGeno=="other" & is.na(validation_onco$Behaviour1)),]) #1218

get_table2=function(dat1=training_icogs[training_icogs$EthnicityGeno=="European",],
                    dat2=training_onco[training_onco$EthnicityGeno=="European",])
{
  res=data.frame(nctr=NA,ncase=NA,lumA=NA,lumB=NA,LumB_HN=NA,HER2E=NA,TripN=NA,Unknown=NA)
  res$ncase=sum(dat1$Behaviour1==1,na.rm=T)+sum(dat2$Behaviour1==1,na.rm = T)
  res$nctr=sum(is.na(dat1$Behaviour1))+sum(is.na(dat2$Behaviour1))
  
  y.pheno.mis1 <- dat1[,c("Behaviour1","ER_status1","PR_status1","HER2_status1","Grade1")]
  y.pheno.mis2 <- dat2[,c("Behaviour1","ER_status1","PR_status1","HER2_status1","Grade1")]
  y.pheno.mis1=rbind(y.pheno.mis1,y.pheno.mis2)
  y.pheno.mis1=y.pheno.mis1[!is.na(y.pheno.mis1$Behaviour1),]
  
  
  y.pheno.mis1=y.pheno.mis1[,-which(colnames(y.pheno.mis1)=="Behaviour1")]
  #intrinsic subtypes are defined as follows
  #Luminal-A like: ER or PR +, HER2-, grade is 1 or 2
  #Luminal-B like: ER or PR +, HER2+
  #Luminal B HER2 negative-like: ER or PR+, HER2-, grade 3
  #HER2 enriched-like: both ER and PR-, HER2+
  #Triple negative: ER, PR, HER2-
  #define Luminal-A like
  idx.1 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &(y.pheno.mis1[,4]==1|y.pheno.mis1[,4]==2))
  res$lumA=length(idx.1)
  #define Luminal-B like
  idx.2 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==1)
  res$lumB=length(idx.2)
  #for Luminal B HER2 negative-like
  idx.3 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &y.pheno.mis1[,4]==3)
  res$LumB_HN=length(idx.3)
  #for HER2 enriched-like
  idx.4 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==1)
  res$HER2E=length(idx.4)
  #for Triple negative
  idx.5 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==0)
  res$TripN=length(idx.5)
  res$Unknown=res$ncase-sum(res[1,3:7])
  return(res)
}

get_table1=function(dat1=tuning[tuning$EthnicityGeno=="European",])
{
  res=data.frame(nctr=NA,ncase=NA,lumA=NA,lumB=NA,LumB_HN=NA,HER2E=NA,TripN=NA,Unknown=NA)
  res$ncase=sum(dat1$Behaviour1==1,na.rm=T)
  res$nctr=sum(is.na(dat1$Behaviour1))
  
  y.pheno.mis1 <- dat1[,c("Behaviour1","ER_status1","PR_status1","HER2_status1","Grade1")]
  y.pheno.mis1=y.pheno.mis1[!is.na(y.pheno.mis1$Behaviour1),]
  y.pheno.mis1=y.pheno.mis1[,-which(colnames(y.pheno.mis1)=="Behaviour1")]
  #intrinsic subtypes are defined as follows
  #Luminal-A like: ER or PR +, HER2-, grade is 1 or 2
  #Luminal-B like: ER or PR +, HER2+
  #Luminal B HER2 negative-like: ER or PR+, HER2-, grade 3
  #HER2 enriched-like: both ER and PR-, HER2+
  #Triple negative: ER, PR, HER2-
  #define Luminal-A like
  idx.1 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &(y.pheno.mis1[,4]==1|y.pheno.mis1[,4]==2))
  res$lumA=length(idx.1)
  #define Luminal-B like
  idx.2 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==1)
  res$lumB=length(idx.2)
  #for Luminal B HER2 negative-like
  idx.3 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &y.pheno.mis1[,4]==3)
  res$LumB_HN=length(idx.3)
  #for HER2 enriched-like
  idx.4 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==1)
  res$HER2E=length(idx.4)
  #for Triple negative
  idx.5 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==0)
  res$TripN=length(idx.5)
  res$Unknown=res$ncase-sum(res[1,3:7])
  return(res)
}

table_EURtraining=get_table2()
#  nctr ncase lumA lumB LumB_HN HER2E TripN Unknown
# 60204 69111 8783 2452    2208  1069  3520   51079
table_EURtraining_icogs=get_table1(dat1=training_icogs[training_icogs$EthnicityGeno=="European",])
#  nctr ncase lumA lumB LumB_HN HER2E TripN Unknown
# 33505 33961 7031 1417    1720   722  2393   20678
table_EURtraining_onco=get_table1(dat1=training_onco[training_onco$EthnicityGeno=="European",])
#  nctr ncase lumA lumB LumB_HN HER2E TripN Unknown
# 26699 35150 1752 1035     488   347  1127   30401
table_EURtuning=get_table1()
table_EURvalidation=get_table1(dat1=validation_onco[validation_onco$EthnicityGeno=="European",])

table_EAStraining=get_table2(dat1=training_icogs[training_icogs$EthnicityGeno=="Asian",],
                             dat2=training_onco[training_onco$EthnicityGeno=="Asian",])
table_EAStraining_icogs=get_table1(dat1=training_icogs[training_icogs$EthnicityGeno=="Asian",])
table_EAStraining_onco=get_table1(dat1=training_onco[training_onco$EthnicityGeno=="Asian",])
table_EAStuning=get_table1(dat1=tuning[tuning$EthnicityGeno=="Asian",])
table_EASvalidation=get_table1(dat1=validation_onco[validation_onco$EthnicityGeno=="Asian",])
table_AFRvalidation=get_table2(dat1=validation_icogs[validation_icogs$EthnicityGeno=="African",],
                               dat2=validation_onco[validation_onco$EthnicityGeno=="African",])
#The numbers in S4 (AFRonco)
table_AFRvalidation1=get_table1(dat1=validation_onco[validation_onco$EthnicityGeno=="African",])

table_HISvalidation=get_table1(dat1=validation_onco[validation_onco$EthnicityGeno=="other",])
alltable=rbind(table_EURtraining,table_EURtuning,table_EURvalidation,
               table_EAStraining,table_EAStuning,table_EASvalidation,
               table_AFRvalidation,table_HISvalidation)
rownames(alltable)=c("EUR_training","EUR_tuning","EUR_validation",
                     "EAS_training","EAS_tuning","EAS_validation",
                     "AFR_validation","HIS_validation")
trainingtable=rbind(table_EURtraining_icogs,table_EURtraining_onco,
                    table_EAStraining_icogs,table_EAStraining_onco)
rownames(trainingtable)=c("EUR_training_icogs","EUR_training_onco",
                          "EAS_training_icogs","EAS_training_onco")
write.csv(alltable,file="../result/check_samplespliting.csv",row.names = T,quote=F)
write.csv(trainingtable,file="../result/trainingsamples.csv",row.names = T,quote=F)

