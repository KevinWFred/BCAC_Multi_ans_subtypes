#!/usr/bin/env Rscript
#run on genome-wide for 750 icogs/onco data (without meta)
.libPaths(c("/data/wangx53",.libPaths()))

library(readr)
library(devtools)
library(data.table)
library(bc2)
library(bcutility)

#library(devtools)
#install_github("andrewhaoyu/bc2")
#install_github("andrewhaoyu/bcutility")
#install_github("andrewhaoyu/bc2",branch="version 0.0.2")
#library(bc2,lib.loc ='/Users/zhangh24/Library/R/3.4/library')

###1 represent Icog
###2 represent Onco
###load_all("/Users/zhangh24/GoogleDrive/bc2")

#rm(list=ls())
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
set.seed(123)
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
dataopt=args[1] #icogs or onco
pop=args[2] #euro
i1 = as.numeric(args[[3]]) #ld_block ID
#output
all.sigma.log.odds=list()
all.log.odds=NULL
#genotype file
#genofile=args[2]# "../result/icogs_test.traw" or "../result/onco_test.traw"
if (dataopt=="icogs")
{
  genofile=paste0("../result/imp_icogs/",pop,"/geno/geno_",i1,".traw.gz")
  if (!file.exists(genofile)) stop("no variants in this block!")
  outfolder=paste0("../result/imp_icogs/",pop,"/res")
  if (!dir.exists(outfolder)) dir.create(outfolder)
  #all the samples, some of them will not be used.
  pheno=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
  colnames(pheno)[which(colnames(pheno)=="SG_ID")]="ID"
  samplefile="../result/icogs_samples_750.txt"
}else
{
  genofile=paste0("../result/imp_onco/",pop,"/geno/geno_",i1,".traw.gz")
  if (!file.exists(genofile)) stop("no variants in this block!")
  outfolder=paste0("../result/imp_onco/",pop,"/res")
  if (!dir.exists(outfolder)) dir.create(outfolder)
  pheno=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_onco_pheno_v15_02_corrected.txt",header=T,sep="\t")
  colnames(pheno)[which(colnames(pheno)=="Onc_ID")]="ID"
  samplefile="../result/onco_samples_750.txt"
}

#samples included in the analysis
allsamples=read.table(samplefile)$V1
#read genotype data
genodat=as.data.frame(fread(genofile))
rownames(genodat)=genodat$SNP
genodat=genodat[,7:ncol(genodat)]
#process sample names
tmp=unlist(strsplit(colnames(genodat),"_"))
tmp1=unlist(strsplit(colnames(genodat)[1],"_"))
tmp2=tmp[seq(1,length(tmp),length(tmp1))]
if (length(tmp2) !=ncol(genodat)) stop("Sampel names not correct")
colnames(genodat)=tmp2
##for now, only use EUR samples
#allsamples=intersect(allsamples,pheno$ID[which(pheno$EthnicityGeno=="European")])
#table(colnames(genodat) %in% allsamples)
allsamples=intersect(allsamples,colnames(genodat)) #pick samples from an ancestry
idx=match(allsamples,colnames(genodat))
#if (sum(is.na(idx))>0) stop("Some samples not in genotypedata")
genodat=genodat[,idx]
pheno=pheno[match(allsamples,pheno$ID),]
#number of snps in the current block
nsnps=nrow(genodat) 
data1=pheno
#behaviour1 to (0,1)
data1$Behaviour1[which(is.na(data1$Behaviour1))]=0
#prepare the phenotypes data for iCOGs 
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
#ER, PR, HER2 is binary with negative as 0 and positive as 1
#Grade is ordinal with 1, 2, 3
#controls don't have tumor characteristics data (all NA)
#cases with missing tumor characteristics marked as 888
colnames(y.pheno.mis1) <- c("Behaviour","ER","PR","HER2","Grade")
#idx=which(is.na(y.pheno.mis1[,1]))
#y.pheno.mis1[idx,1]
#generate the z standard matrix
#z standard matrix is link to link subtypes with tumor characteristics
z.standard <- GenerateZstandard(y.pheno.mis1)
#each row of z.standard represent a subtype
#each column represent a tumor marker
#e.g. first row is ER-PR-HER2-Grade 1
#total number of subtypes
#subtypes with less than 10 cases are automatically removed
M <- nrow(z.standard)
#construct the z design matrix for intrinsic subtypes
#intrinsic subtypes are defined as follows
#Luminal-A like: ER or PR +, HER2-, grade is 1 or 2
#Luminal-B like: ER or PR +, HER2+
#Luminal B HER2 negative-like: ER or PR+, HER2-, grade 3
#HER2 enriched-like: both ER and PR-, HER2+
#Triple negative: ER, PR, HER2-
#Prepare z design matrix
z.design <- matrix(0,M,5)
#define Luminal-A like
idx.1 <- which((z.standard[,1]==1|z.standard[,2]==1)
               &z.standard[,3]==0
               &(z.standard[,4]==1|z.standard[,4]==2))
z.design[idx.1,1] <- 1
#define Luminal-B like
idx.2 <- which((z.standard[,1]==1|z.standard[,2]==1)
               &z.standard[,3]==1)
z.design[idx.2,2] <- 1
#for Luminal B HER2 negative-like
idx.3 <- which((z.standard[,1]==1|z.standard[,2]==1)
               &z.standard[,3]==0
               &z.standard[,4]==3)
z.design[idx.3,3] <- 1
#for HER2 enriched-like
idx.4 <- which(z.standard[,1]==0&z.standard[,2]==0
               &z.standard[,3]==1)
z.design[idx.4,4] <- 1
#for Triple negative
idx.5 <- which(z.standard[,1]==0&z.standard[,2]==0
               &z.standard[,3]==0)
z.design[idx.5,5] <- 1

#genotype data for 178 known SNPs (178th SNP only exit on Oncoarray)
x.test.all.mis1 <- t(genodat) #data1[,c(27:203)]
#prepare covariates table: PC1-10
#we only adjusted PC1-10 in known SNPs analyses
#in genome-wide analyses, we will need to adjust age and PC1-10
#need to add age
if ("pc1" %in% colnames(pheno))
{
  allcov=c(paste0("pc",1:10))
}else
{
  allcov=c(paste0("PC_",1:10))
}

x.covar.mis1 =pheno[,allcov]
#x.covar.mis1 <- data1[,5:14]

#################################
#this section is not important
#it's mainly used to convert the effect-size to minor allele
#no need for this in analyses
# idx.control <- which(y.pheno.mis1[,1]==0)
# maf <- sum(x.test.all.mis1[idx.control,i1])/(2*length(idx.control))
# if(maf>=0.5){
#   x.test.all.mis1[,i1] < 2 - x.test.all.mis1[,i1]
# }
#################################


# x.all.mis1 <- as.matrix(cbind(x.test.all.mis1[,i1],x.covar.mis1))
# colnames(x.all.mis1)[1] <- "gene"

gene_value = x.test.all.mis1[,i1,drop = F]
#Fit the two-stage model
#Two stage have several z design matrix structure
#baseline only assume all subtypes have the same effect
#additive assumes all higher order interactions across to be 0
#pair-wise interaction allows main effect and second order interactions across tumor markers
#saturated model allows all interactions
#we can also self design the z matrix
#in this setting, we self define the z matrix to make the definition align with intrinsic subtypes
#we only do this self define z matrix for genetic variant
#meanwhile, we keep the additive model for all other covariates
Sys.time()
Heter.result.Icog = EMmvpolySelfDesign(y.pheno.mis1,
                                       x.self.design = gene_value,
                                       z.design=z.design,
                                       baselineonly = NULL,
                                       additive = x.covar.mis1,
                                       pairwise.interaction = NULL,
                                       saturated = NULL,missingTumorIndicator = 888)
Sys.time()
n.param <- ncol(z.design)
#the log-odds ratio for second-stage parameters are saved in the first elelment
#first M parameters are for intercept
#We don't make any assumptions regarding the intercept
#Therefore, the intercept has the same degree of freedom panelty as the M subtypes
#The next five parameters are for the genetic variants
#They represent the log-odds ratio for the five intrinsic subtypes
log.odds <- Heter.result.Icog[[1]][(M+1):(M+n.param)]
log.odds=as.data.frame(t(log.odds))
rownames(log.odds)=rownames(genodat)[i1]
all.log.odds=rbind(all.log.odds,log.odds)
#nparm <- length(Heter.result.Icog[[1]])  
#variance matrix for the log-odds-ratio are saved in the second component
sigma.log.odds <- Heter.result.Icog[[2]][(M+1):(M+n.param),(M+1):(M+n.param)]
all.sigma.log.odds[[rownames(genodat)[i1]]]=sigma.log.odds
save(all.log.odds,all.sigma.log.odds,file=paste0(outfolder,"/res_",i1,".RData"))

print("Done")
print(Sys.time())




