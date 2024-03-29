#!/usr/bin/env Rscript
#run conditional analysis on 11 novel snps 

.libPaths(c("/data/wangx53",.libPaths()))
library(readr)
library(devtools)
library(data.table)
library(bc2)
library(bcutility)
library(ACAT)

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

#from bc2. set debug=0, no intermediate output to screen
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

twotests_fun=function(dataopt="icogs",pop="euro",i1=1)
{
  log.odds=sigma.log.odds=score_result=infor_result=NULL
  #intrinsic_subtypes 
  snpid=allnovelsnps1$ID[i1]
  condsnpid=allnovelsnps1$knownvar[i1]
  condsnpid=unlist(strsplit(condsnpid,","))
  condsnpid=c(snpid,condsnpid)
  condsnpfile=paste0("../result/tmp/",snpid,".txt")
  write.table(condsnpid,file=condsnpfile,row.names = F,col.names = F,quote=F)
  if (file.exists(paste0("../result/tmp/",snpid,".traw"))) file.remove(paste0("../result/tmp/",snpid,".traw"))
  cmd=paste0(plink2," --pfile ../result/imp_",dataopt,"/",pop,"/",pop," --extract ",condsnpfile," --memory 16000 --threads 1 --recode A-transpose --out ../result/tmp/",snpid)
  system(cmd)
  if (file.exists(paste0("../result/tmp/",snpid,".traw")))
  {
    rawdat=as.data.frame(fread(paste0("../result/tmp/",snpid,".traw")))
    
    if (dataopt=="icogs")
    {
      # genofile=paste0("../result/imp_icogs/",pop,"/geno/geno_",i1,".traw.gz")
      # if (!file.exists(genofile)) stop("no variants in this block!")
      # outfolder=paste0("../result/imp_icogs/",pop,"/res")
      # if (!dir.exists(outfolder)) dir.create(outfolder)
      #all the samples, some of them will not be used.
      pheno=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
      colnames(pheno)[which(colnames(pheno)=="SG_ID")]="ID"
      samplefile="../result/icogs_samples_750.txt"
    }else
    {
      # genofile=paste0("../result/imp_onco/",pop,"/geno/geno_",i1,".traw.gz")
      # if (!file.exists(genofile)) stop("no variants in this block!")
      # outfolder=paste0("../result/imp_onco/",pop,"/res")
      # if (!dir.exists(outfolder)) dir.create(outfolder)
      pheno=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
      colnames(pheno)[which(colnames(pheno)=="Onc_ID")]="ID"
      samplefile="../result/onco_samples_750.txt"
    }
    
    #samples included in the analysis
    allsamples=read.table(samplefile)$V1
    #read genotype data
    genodat=rawdat
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
    idx=rownames(genodat)==snpid
    condsnpdat=genodat[!idx,]
    genodat=genodat[idx,]
    if (sum(idx)>0 & sum(!idx)>0)
    {
      pheno=pheno[match(allsamples,pheno$ID),]
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
      if ("pc1" %in% colnames(pheno)) #iCOGS
      {
        allcov=c(paste0("plinkPC",1:10),"age")
        tmp=read.table("../result/imp_icogs/merged1.eigenvec")
        colnames(tmp)=c("ID",paste0("plinkPC",1:20))
        tmp1=unlist(strsplit(tmp$ID,"_"))
        tmp$ID=tmp1[seq(1,length(tmp1),2)]
        idx=match(pheno$ID,tmp$ID)
        pheno=cbind(pheno,tmp[idx,2:11])
        
      }else
      {
        allcov=c(paste0("PC_",1:10),"age")
      }
      #add conditional snp
      all(pheno$ID==colnames(condsnpdat))
      ncondsnp=nrow(condsnpdat)
      pheno=cbind(pheno,t(condsnpdat))
      condsnpid1=rownames(condsnpdat)
      allcov=c(allcov,condsnpid1)
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
      # outfile=paste0(outfolder,"/res_",i1,".RData")
      # ilast=1
      
      gene_value = x.test.all.mis1[,1,drop = F]
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
        rownames(log.odds)=rownames(genodat)[1]
        #nparm <- length(Heter.result[[1]])  
        #variance matrix for the log-odds-ratio are saved in the second component
        sigma.log.odds <- Heter.result[[2]][(M+1):(M+n.param),(M+1):(M+n.param)]
      }else
      {
        warning(paste0(rownames(genodat)[1]," notconverge"))
      }
      
      #score test
      score.test.support.ERPRHER2Grade <- ScoreTestSupport(
        y.pheno.mis1,
        baselineonly = NULL,
        additive = x.covar.mis1,
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
      
      ##get the three different z design matrix
      z.design.list = GenerateZDesignCombination(y.pheno.mis1)
      z.additive = z.design.list[[1]]
      z.interaction = z.design.list[[2]]
      z.saturated = z.design.list[[3]]
      #number of second stage parameters
      #if use additive model
      n.second = ncol(z.additive)
      #if use pair-wise interaction model
      #n.second = ncol(z.interaction)
      #if use saturated model
      #n.second = ncol(z.saturated)
      
      #number of subject in the genotype file is n
      n=ncol(genodat)
      #count the number of variants in the file
      num=nrow(genodat)
      all(pheno$ID==colnames(genodat))
      idx.control <- which(y.pheno.mis1[,1]==0)
      #count the number of control in the data
      n.control <- length(idx.control)
      
      #output for each variant
      freq.all=rep(NA,num)
      names(freq.all)=rownames(genodat)
      score_result=data.frame(matrix(NA,nrow=num,ncol=n.second))
      infor_result=data.frame(matrix(NA,nrow=num,ncol=n.second^2))
      rownames(score_result)=rownames(infor_result)=rownames(genodat)
      
      snpvalue = genodat[1,,drop = F]
      snpvalue.control <- snpvalue[idx.control]
      freq <- sum(snpvalue.control,na.rm=T)/(2*n.control)
      freq.all[1] <- freq
      if(freq<0.006|freq>0.994){
        #if the SNP is too rare, just keep as score 0.
        score_result[1,] <- 0
        infor_result[1,] <- 0
      }else
      {
        result=tryCatch(
          expr = {
            score.test<- ScoreTest(y=y.pheno.mis1,
                                   x=unlist(snpvalue),
                                   second.stage.structure="additive",
                                   score.test.support=score.test.support.ERPRHER2Grade,
                                   missingTumorIndicator=888)
          },
          error = function(e){ 
            return(NULL)
          })
        
        if (!is.null(result))
        {
          #the first element is score
          score_result[1,]  <- result[[1]]
          #the second element is the efficient information matrix
          infor_result[1,] <- as.vector(result[[2]])
        }else
        {
          warning(paste0(rownames(genodat)[1]," don't converge"))
        }
      }
      return(list(log.odds=log.odds,sigma.log.odds=sigma.log.odds,score_result=score_result,infor_result=infor_result))
      
    }else
    {
      return(NULL)
    }
 
  }else
  {
    return(NULL)
  }
  
}
LogoddsMetaAnalysis <- function(logoddslist,sigmalist){
  if (length(logoddslist)>1)
  {
    sigma.inv.sum=0
    sigma.inv.logodds.sum=0
    for(i in 1:length(sigmalist))
    {
      if (!is.null(sigmalist[[i]]))
      {
        sigma.inv.sum=sigma.inv.sum+solve(sigmalist[[i]])
        sigma.inv.logodds.sum=sigma.inv.logodds.sum+solve(sigmalist[[i]])%*%unlist(logoddslist[[i]])
        
      }
    }
    sigma.meta=solve(sigma.inv.sum)
    
    logodds.meta <- sigma.meta%*% sigma.inv.logodds.sum
  }else #==1
  {
    logodds.meta=as.numeric(logoddslist[[1]])
    sigma.meta=sigmalist[[1]]
  }
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
  
}
ScoreMetaAnalysis <- function(score=score.meta,infor=infor.meta,second.num=5){
  #get the score
  score=as.matrix(score)
  infor=matrix(unlist(infor),nrow=second.num)
  p.value=tryCatch(
    expr = {
      ScoreGlobalTestForAssoc(score,infor)
    },
    error = function(e){ 
      return(NA)
    })
  
  return(p.value)
}


args <- commandArgs(trailingOnly=T)
i1=as.numeric(args[1])
snpfile=args[2] #"../result/QCallnovelsnps1.txt"
outprefix=args[3] #"QC"

allnovelsnps1=read.table(snpfile,header=T)

#QCed SNP lists:
qc_african_onco=as.data.frame(fread("../result/imp_QC/onco/african_info.snplist",header=F))
qc_asian_onco=as.data.frame(fread("../result/imp_QC/onco/asian_info.snplist",header=F))
qc_euro_onco=as.data.frame(fread("../result/imp_QC/onco/euro_info.snplist",header=F))
qc_asian_icogs=as.data.frame(fread("../result/imp_QC/icogs/asian_info.snplist",header=F))
qc_euro_icogs=as.data.frame(fread("../result/imp_QC/icogs/euro_info.snplist",header=F))
table(allnovelsnps1$ID %in% qc_african_onco$V1) #F1
table(allnovelsnps1$ID %in% qc_asian_onco$V1) #F4
table(allnovelsnps1$ID %in% qc_euro_onco$V1) #F0
table(allnovelsnps1$ID %in% qc_asian_icogs$V1) #F5
table(allnovelsnps1$ID %in% qc_euro_icogs$V1) #F1 

res=data.frame(snp=allnovelsnps1$ID[i1],scoreP=NA,intrinsicP=NA,acatP=NA,euro_icogs_intr=NA,euro_onco_intr=NA,asian_icogs_intr=NA,asian_onco_intr=NA,african_onco_intr=NA,
               euro_icogs_score=NA,euro_onco_score=NA,asian_icogs_score=NA,asian_onco_score=NA,african_onco_score=NA)
outfile=paste0("../result/",outprefix,"conditional_result",i1,".txt")
#this is the result to remove the conditional snp: #allcov=c(allcov,"condsnp")
#outfile=paste0("../result/noconditional_result",i1,".txt")
#for (i1 in 1:nrow(allnovelsnps1))
#{
  print(i1)
  twotests_euro_icogs=NULL
  if (allnovelsnps1$ID[i1] %in% qc_euro_icogs[,1])
    twotests_euro_icogs=twotests_fun(dataopt="icogs",pop="euro",i1=i1)
  twotests_euro_onco=NULL
  if (allnovelsnps1$ID[i1] %in% qc_euro_onco[,1])
    twotests_euro_onco=twotests_fun(dataopt="onco",pop="euro",i1=i1)
  twotests_asian_icogs=NULL
  if (allnovelsnps1$ID[i1] %in% qc_asian_icogs[,1])
    twotests_asian_icogs=twotests_fun(dataopt="icogs",pop="asian",i1=i1)
  twotests_asian_onco=NULL
  if (allnovelsnps1$ID[i1] %in% qc_asian_onco[,1])
    twotests_asian_onco=twotests_fun(dataopt="onco",pop="asian",i1=i1)
  twotests_african_onco=NULL
  if (allnovelsnps1$ID[i1] %in% qc_african_onco[,1])
    twotests_african_onco=twotests_fun(dataopt="onco",pop="african",i1=i1)
  
  logoddslist=sigmalist=list()
  score.meta=data.frame(matrix(0,nrow=1,ncol=5))
  infor.meta=data.frame(matrix(0,nrow=1,ncol=25))
  if (!is.null(twotests_euro_icogs))
  {
    if (!is.null(twotests_euro_icogs$log.odds))
    {
      logoddslist=c(logoddslist,list(twotests_euro_icogs$log.odds))
      sigmalist=c(sigmalist,list(twotests_euro_icogs$sigma.log.odds))
      tmp=LogoddsMetaAnalysis(list(twotests_euro_icogs$log.odds),list(twotests_euro_icogs$sigma.log.odds))
      res$euro_icogs_intr=GlobalTestForAssoc(tmp$logodds.meta,tmp$sigma.meta)
    }
    
    score.meta=score.meta+twotests_euro_icogs$score_result
    infor.meta=infor.meta+twotests_euro_icogs$infor_result
    res$euro_icogs_score=ScoreMetaAnalysis(score=twotests_euro_icogs$score_result,infor=twotests_euro_icogs$infor_result,second.num=5)
  }
  if (!is.null(twotests_euro_onco))
  {
    if (!is.null(twotests_euro_onco$log.odds))
    {
      logoddslist=c(logoddslist,list(twotests_euro_onco$log.odds))
      sigmalist=c(sigmalist,list(twotests_euro_onco$sigma.log.odds))
      tmp=LogoddsMetaAnalysis(list(twotests_euro_onco$log.odds),list(twotests_euro_onco$sigma.log.odds))
      res$euro_onco_intr=GlobalTestForAssoc(tmp$logodds.meta,tmp$sigma.meta)
    }
    
    score.meta=score.meta+twotests_euro_onco$score_result
    infor.meta=infor.meta+twotests_euro_onco$infor_result
    res$euro_onco_score=ScoreMetaAnalysis(score=twotests_euro_onco$score_result,infor=twotests_euro_onco$infor_result,second.num=5)
  }
  if (!is.null(twotests_asian_icogs))
  {
    if (!is.null(twotests_asian_icogs$log.odds))
    {
      logoddslist=c(logoddslist,list(twotests_asian_icogs$log.odds))
      sigmalist=c(sigmalist,list(twotests_asian_icogs$sigma.log.odds))
      tmp=LogoddsMetaAnalysis(list(twotests_asian_icogs$log.odds),list(twotests_asian_icogs$sigma.log.odds))
      res$asian_icogs_intr=GlobalTestForAssoc(tmp$logodds.meta,tmp$sigma.meta)
    }
    
    score.meta=score.meta+twotests_asian_icogs$score_result
    infor.meta=infor.meta+twotests_asian_icogs$infor_result
    res$asian_icogs_score=ScoreMetaAnalysis(score=twotests_asian_icogs$score_result,infor=twotests_asian_icogs$infor_result,second.num=5)
  }
  if (!is.null(twotests_asian_onco))
  {
    if (!is.null(twotests_asian_onco$log.odds))
    {
      logoddslist=c(logoddslist,list(twotests_asian_onco$log.odds))
      sigmalist=c(sigmalist,list(twotests_asian_onco$sigma.log.odds))
      tmp=LogoddsMetaAnalysis(list(twotests_asian_onco$log.odds),list(twotests_asian_onco$sigma.log.odds))
      res$asian_onco_intr=GlobalTestForAssoc(tmp$logodds.meta,tmp$sigma.meta)
    }
    
    score.meta=score.meta+twotests_asian_onco$score_result
    infor.meta=infor.meta+twotests_asian_onco$infor_result
    res$asian_onco_score=ScoreMetaAnalysis(score=twotests_asian_onco$score_result,infor=twotests_asian_onco$infor_result,second.num=5)
  }
  if (!is.null(twotests_african_onco))
  {
    if (!is.null(twotests_african_onco$log.odds))
    {
      logoddslist=c(logoddslist,list(twotests_african_onco$log.odds))
      sigmalist=c(sigmalist,list(twotests_african_onco$sigma.log.odds))
      tmp=LogoddsMetaAnalysis(list(twotests_african_onco$log.odds),list(twotests_african_onco$sigma.log.odds))
      res$african_onco_intr=GlobalTestForAssoc(tmp$logodds.meta,tmp$sigma.meta)
    }
    
    score.meta=score.meta+twotests_african_onco$score_result
    infor.meta=infor.meta+twotests_african_onco$infor_result
    res$african_onco_score=ScoreMetaAnalysis(score=twotests_african_onco$score_result,infor=twotests_african_onco$infor_result,second.num=5)
  }
  metares=LogoddsMetaAnalysis(logoddslist,sigmalist)
  intrinsicP=GlobalTestForAssoc(metares$logodds.meta,metares$sigma.meta)
  scoreP=ScoreMetaAnalysis(score=score.meta,infor=infor.meta,second.num=5)
  res$intrinsicP=intrinsicP
  res$scoreP=scoreP
  if (!is.na(intrinsicP) & !is.na(scoreP))
   res$acatP=ACAT(c(intrinsicP,scoreP))
#}

write.table(res,file=outfile,row.names=F,sep="\t",quote=F)
print(Sys.time())
tmp1=read.table("../result/QCallnovelsnps1.txt",header=T)
tmp=data.frame(code="/data/BB_Bioinformatics/Kevin/BCAC/code/conditional_analysis_QC.R",i1=1:nrow(tmp1),snpfile="../result/QCallnovelsnps1.txt",outprefix="QC")
# write.table(tmp,file="conditional_analysis_QC.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#17689872
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/conditional_analysis_QC.swarm -g 32 --module R/4.3 --time=08:00:00 --gres=lscratch:32

tmp1=read.table("../result/QCfreq01allnovelsnps1.txt",header=T)
tmp=data.frame(code="/data/BB_Bioinformatics/Kevin/BCAC/code/conditional_analysis_QC.R",i1=1:nrow(tmp1),snpfile="../result/QCfreq01allnovelsnps1.txt",outprefix="QCfreq01")
#write.table(tmp,file="conditional_analysis_QCfreq01.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#17689880
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/conditional_analysis_QCfreq01.swarm -g 32 --module R/4.3 --time=08:00:00 --gres=lscratch:32


#read results
read_results=function(snpfile="../result/QCallnovelsnps1.txt",outprefix="QC")
{
  allnovelsnps1=read.table(snpfile,header=T)
  allres=NULL
  for(i in 1:nrow(allnovelsnps1))
  {
    outfile=paste0("../result/",outprefix,"conditional_result",i,".txt")
    tmp=read.table(outfile,header=T)
    allres=rbind(allres,tmp)
  }
  idx=match(allres$snp,allnovelsnps1$ID)
  allres$freq=allnovelsnps1$freq[idx]
  print(table(allres$acatP<1e-6))
  # FALSE  TRUE
  # 2     7
  allres1=allnovelsnps1[,c("ID","rsid","dist2nearestknown1","knownvar_rsid")]
  idx=match(allres1$ID,allres$snp)
  allres1=cbind(allres1,allres[idx,2:4])
  write.csv(allres1,file=paste0("../result/",outprefix,"conditional_analysis_res.csv"),row.names = F,quote=T)
  
}
read_results(snpfile="../result/QCallnovelsnps1.txt",outprefix="QC")
read_results(snpfile="../result/QCfreq01allnovelsnps1.txt",outprefix="QCfreq01")
# FALSE  TRUE 
# 2     4