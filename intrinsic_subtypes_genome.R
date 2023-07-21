#!/usr/bin/env Rscript
#run on genome-wide for 750 icogs/onco data (without meta)

.libPaths(c("/data/wangx53",.libPaths()))
library(readr)
library(devtools)
library(data.table)
library(bc2)
library(bcutility)

#from bc2. set debug=0, no intermiate output to screen
EMStep=function (delta0, y, x.all, z.standard, z.all, missingTumorIndicator) 
{
  tol <- as.numeric(1e-04)
  tolMaxstep <- as.numeric(0.001)
  delta_old <- delta0
  #print(paste0("Begin EM algorithm"))
  #print(paste0("EM round: 1"))
  prob.fit.result <- ProbFitting(delta_old, y, x.all, z.standard, 
                                 z.all, missingTumorIndicator)
  y_em <- prob.fit.result[[1]]
  missing.vec <- as.numeric(as.vector(prob.fit.result[[2]]))
  missing.mat <- prob.fit.result[[3]]
  missing.mat.vec <- as.numeric(as.vector(missing.mat))
  missing.number <- as.integer(length(missing.vec))
  idx.drop = prob.fit.result[[4]]
  if (length(idx.drop) != 0) {
    x.all <- x.all[-idx.drop, , drop = F]
    y_em <- y_em[-idx.drop, , drop = F]
  }
  for (k in 1:length(missing.vec)) {
    missing.vec[k] <- missing.vec[k] - sum(missing.vec[k] >= 
                                             idx.drop)
  }
  N <- as.integer(nrow(x.all))
  M <- as.integer(nrow(z.standard))
  NCOV <- as.integer(ncol(x.all))
  NM <- N * M
  nparm <- as.integer(length(delta0))
  deltai <- as.numeric(delta0)
  NITER <- as.integer(500)
  Y <- as.numeric(as.vector(y_em))
  X <- as.numeric(as.vector(x.all))
  ZallVec = as.numeric(as.vector(z.all))
  Znr = as.integer(nrow(z.all))
  Znc = as.integer(ncol(z.all))
  debug <- as.integer(0) #change here
  ret_rc <- as.integer(1)
  ret_delta <- as.numeric(rep(-9999, nparm))
  ret_info <- as.numeric(rep(-9999, nparm^2))
  ret_p <- as.numeric(rep(0, NM))
  ret_lxx <- as.numeric(rep(0, NM))
  loglikelihood <- as.numeric(-1)
  temp <- .C("EMStep", deltai, nparm, Y = Y, X, ZallVec, Znr, 
             Znc, N, M, NCOV, NITER, tol, tolMaxstep, debug, ret_rc = ret_rc, 
             ret_delta = ret_delta, ret_info = ret_info, ret_p = ret_p, 
             missing.vec, missing.mat.vec, missing.number, loglikelihood = loglikelihood)
  #print(paste0("EM Algorithm Converged"))
  info <- matrix(unlist(temp$ret_info), nparm, nparm)
  result <- list(temp$ret_delta, info, temp$ret_p)
  y_em <- matrix(unlist(temp$Y), N, M)
  delta = result[[1]]
  infor_obs = result[[2]]
  p = result[[3]]
  loglikelihood = temp$loglikelihood
  AIC = 2 * nparm - 2 * loglikelihood
  return(list(delta = delta, infor_obs = infor_obs, p = p, 
              y_em = y_em, M = M, NumberofTumor = ncol(z.standard), 
              loglikelihood = loglikelihood, AIC = AIC))
}

#no change
EMmvpolySelfDesign=function (y, x.self.design, z.design, baselineonly = NULL, additive = NULL, 
                             pairwise.interaction = NULL, saturated = NULL, missingTumorIndicator = 888, 
                             z.all = NULL, delta0 = NULL, cutoff = 10) 
{
  if (is.null(z.all)) {
    missing.data.vec <- GenerateMissingPosition(y, missingTumorIndicator)
    y.pheno.complete <- y[-missing.data.vec, , drop = F]
    initial.set <- InitialSetup(y.pheno.complete, baselineonly, 
                                additive, pairwise.interaction, saturated, x.self.design, 
                                z.design, cutoff = cutoff)
    if (is.null(delta0)) {
      delta0 = initial.set$delta0
    }
    z.all = initial.set$z.all
    z.standard = initial.set$z.standard
    z.deisign.baselineonly = initial.set$z.design.baseline.only
    z.design.additive = initial.set$z.design.additive
    z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
    z.design.saturated = initial.set$z.design.saturated
    x.all <- GenerateSelfXAll(y, x.self.design, baselineonly, 
                              additive, pairwise.interaction, saturated)
    covar.names <- initial.set$covar.names
    tumor.names <- initial.set$tumor.names
    y <- as.matrix(y)
    x.all <- as.matrix(x.all)
    M <- as.integer(nrow(z.standard))
    p.main <- ncol(z.standard) + 1
    model.result = EMStep(delta0, as.matrix(y), x.all, z.standard, 
                          z.all, missingTumorIndicator)
    summary.result <- SummaryResult(model.result, baselineonly, 
                                    additive, pairwise.interaction, saturated, z.standard, 
                                    covar.names, delta, z.design.additive, z.design.pairwise.interaction, 
                                    z.design.saturated, tumor.names, z.all, x.self.design, 
                                    z.design)
    return(summary.result)
  }
  else {
  }
}



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
print(paste0("host:",as.character(Sys.info()["nodename"])))
set.seed(123)
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
print(args)
dataopt=args[1] #icogs or onco
pop=args[2] #euro
i1 = as.numeric(args[[3]]) #block ID, starts with 1
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
  pheno=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
  colnames(pheno)[which(colnames(pheno)=="SG_ID")]="ID"
  samplefile="../result/icogs_samples_750.txt"
}else
{
  genofile=paste0("../result/imp_onco/",pop,"/geno/geno_",i1,".traw.gz")
  if (!file.exists(genofile)) stop("no variants in this block!")
  outfolder=paste0("../result/imp_onco/",pop,"/res")
  if (!dir.exists(outfolder)) dir.create(outfolder)
  pheno=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
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
if (length(tmp2) !=ncol(genodat)) stop("Sampel names are not correct")
colnames(genodat)=tmp2
##for now, only use EUR samples
#allsamples=intersect(allsamples,pheno$ID[which(pheno$EthnicityGeno=="European")])
#table(colnames(genodat) %in% allsamples)
allsamples=intersect(allsamples,colnames(genodat)) #pick samples from an ancestry
#make pheno and genodat consistent in sample orders
idx=match(allsamples,colnames(genodat))
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

#genotype data
x.test.all.mis1 <- t(genodat)
#prepare covariates table: PC1-10
#we only adjusted PC1-10 in known SNPs analyses
#in genome-wide analyses, we will need to adjust age and PC1-10
#need to add age
if ("pc1" %in% colnames(pheno))
{
  allcov=c(paste0("pc",1:10),"age")
}else
{
  allcov=c(paste0("PC_",1:10),"age")
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

n.param <- ncol(z.design)
#output file
outfile=paste0(outfolder,"/res_",i1,".RData")
ilast=1
#if part of the job is done
if (file.exists(outfile))
{
  #load(outfile)
  loaddata = tryCatch(
    expr = {
      load(outfile)
    },
    error = function(e){ 
      return(NULL)
    })
  if (!is.null(loaddata))
  {
    print("Continue the job....")
    if (nrow(all.log.odds)==ncol(x.test.all.mis1))
    {
      stop("All the job was done!")
    }
    ilast=ilast+1 #the next one
  }
}
#for each variant
for (i in ilast: ncol(x.test.all.mis1))
{
  
  gene_value = x.test.all.mis1[,i,drop = F]
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
  #Sys.time()
  Heter.result=tryCatch(
    expr = {
      EMmvpolySelfDesign(y.pheno.mis1,
                         x.self.design = gene_value,
                         z.design=z.design,
                         baselineonly = NULL,
                         additive = x.covar.mis1,
                         pairwise.interaction = NULL,
                         saturated = NULL,missingTumorIndicator = 888)
    },
    error = function(e){ 
      return(NULL)
    })
  # Heter.result = EMmvpolySelfDesign(y.pheno.mis1,
  #                                        x.self.design = gene_value,
  #                                        z.design=z.design,
  #                                        baselineonly = NULL,
  #                                        additive = x.covar.mis1,
  #                                        pairwise.interaction = NULL,
  #                                        saturated = NULL,missingTumorIndicator = 888)
  #Sys.time()
  
  #the log-odds ratio for second-stage parameters are saved in the first elelment
  #first M parameters are for intercept
  #We don't make any assumptions regarding the intercept
  #Therefore, the intercept has the same degree of freedom panelty as the M subtypes
  #The next five parameters are for the genetic variants
  #They represent the log-odds ratio for the five intrinsic subtypes
  if (!is.null(Heter.result))
  {
    log.odds <- Heter.result[[1]][(M+1):(M+n.param)]
    log.odds=as.data.frame(t(log.odds))
    rownames(log.odds)=rownames(genodat)[i]
    all.log.odds=rbind(all.log.odds,log.odds)
    #nparm <- length(Heter.result[[1]])  
    #variance matrix for the log-odds-ratio are saved in the second component
    sigma.log.odds <- Heter.result[[2]][(M+1):(M+n.param),(M+1):(M+n.param)]
    all.sigma.log.odds[[rownames(genodat)[i]]]=sigma.log.odds
  }else
  {
    warning(paste0(i," ",rownames(genodat)[i]," notconverge"))
  }
  
  if (i %% 20==0)
  {
    ilast=i
    cat(i,'..')
    print(Sys.time())
    save(ilast,all.log.odds,all.sigma.log.odds,file=outfile)
  }
}

save(all.log.odds,all.sigma.log.odds,file=outfile)

print("Done")
print(Sys.time())

#create swarm files
create_swarm=function(prefix="../result/imp_icogs/euro/euro",outprefix="euro_icogs",nvar=2000,pop="euro",dataopt="icogs")
{
  pvar=as.data.frame(fread(paste0(prefix,".pvar")))
  infolder=paste0(dirname(prefix),"/geno/")
  n=as.integer(nrow(pvar)/nvar)+1 #total number of blocks
  #each job run two cases 2000*2 (2000 is 1000*2, 2 is p)
  if (n %%4000 ==0 )
  {
    m=n/4000
  }else
  {
    m=as.integer(n/4000)+1
  }
  
  for (i in 1:m)
  {
    if (i<m)
    {
      tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/intrinsic_subtypes_genome.R",2000),
                     dataopt=dataopt,pop=pop,i1=((i-1)*4000+1):(i*4000))
    }else
    {
      tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/intrinsic_subtypes_genome.R"),
                     dataopt=dataopt,pop=pop,i1=((i-1)*4000+1):n)
    }
    
    write.table(tmp,file=paste0(outprefix,i,".swarm"),sep="\t",quote=F,row.names = F,col.names = F)
  }
}
#create_swarm()
#create_swarm(prefix="../result/imp_onco/euro/euro",outprefix="euro_onco",nvar=2000,dataopt="onco")

#cd swarm
#each split into two files
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_onco1.swarm -g 9 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_onco2.swarm -g 9 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2

# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_icogs1.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_icogs2.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2

