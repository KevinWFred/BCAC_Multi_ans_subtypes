#!/usr/bin/env Rscript
#run on genome-wide for 750 icogs/onco data (without meta)

.libPaths(c("/data/wangx53",.libPaths()))
library(readr)
library(devtools)
library(data.table)
library(bc2)
library(bcutility)

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
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
print(args)
dataopt=args[1] #icogs or onco
pop=args[2] #euro
i1 = as.numeric(args[3]) #block ID, starts with 1

#genotype file
#genofile=args[2]# "../result/icogs_test.traw" or "../result/onco_test.traw"
if (dataopt=="icogs")
{
  genofile=paste0("../result/imp_icogs/",pop,"/geno/geno_",i1,".traw.gz")
  if (!file.exists(genofile)) stop("no variants in this block!")
  outfolder=paste0("../result/imp_icogs/",pop,"/scoretest")
  if (!dir.exists(outfolder)) dir.create(outfolder)
  #all the samples, some of them will not be used.
  pheno=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
  colnames(pheno)[which(colnames(pheno)=="SG_ID")]="ID"
  samplefile="../result/icogs_samples_750.txt"
}else
{
  genofile=paste0("../result/imp_onco/",pop,"/geno/geno_",i1,".traw.gz")
  if (!file.exists(genofile)) stop("no variants in this block!")
  outfolder=paste0("../result/imp_onco/",pop,"/scoretest")
  if (!dir.exists(outfolder)) dir.create(outfolder)
  pheno=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
  colnames(pheno)[which(colnames(pheno)=="Onc_ID")]="ID"
  samplefile="../result/onco_samples_750.txt"
}

#load the null hypothesis results for other covariates
#this component will be needed in later ScoreTest
load(paste0(outfolder,"/score.test.support.ERPRHER2Grade.Rdata"))

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

#output file
outfile=paste0(outfolder,"/res_",i1,".RData")

#output for each variant
freq.all=rep(NA,num)
names(freq.all)=rownames(genodat)
score_result=data.frame(matrix(NA,nrow=num,ncol=n.second))
infor_result=data.frame(matrix(NA,nrow=num,ncol=n.second^2))
rownames(score_result)=rownames(infor_result)=rownames(genodat)

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
    print(Sys.time())
    if (ilast>=nrow(genodat))
    {
      stop("All the job was done!")
    }
    ilast=ilast+1 #the next one
  }
}

if (ilast<=nrow(genodat))
{
  for (i in ilast: nrow(genodat))
  {
    
    snpvalue = genodat[i,,drop = F]
    snpvalue.control <- snpvalue[idx.control]
    freq <- sum(snpvalue.control,na.rm=T)/(2*n.control)
    freq.all[i] <- freq
    if(freq<0.006|freq>0.994){
      #if the SNP is too rare, just keep as score 0.
      score_result[i,] <- 0
      infor_result[i,] <- 0.1
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
        score_result[i,]  <- result[[1]]
        #the second element is the efficient information matrix
        infor_result[i,] <- as.vector(result[[2]])
      }else
      {
        warning(paste0(i," ",rownames(genodat)[i]," don't converge"))
      }
    }
    
    if (i %% 50==0)
    {
      ilast=i
      cat(i,'..')
      print(Sys.time())
      save(ilast,freq.all,score_result,infor_result,file=outfile)
    }
    gc()
  }
  ilast=i
}  

save(score_result,infor_result,freq.all,ilast,file=outfile)

print("Done")
print(Sys.time())

#create swarm files
create_swarm=function(prefix="../result/imp_icogs/euro/euro",outprefix="scoretest_euro_icogs",nvar=2000,pop="euro",dataopt="icogs")
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
      tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/whole_genome_ERPRHER2Grade_fixed_baseline.R",2000),
                     dataopt=dataopt,pop=pop,i1=((i-1)*4000+1):(i*4000))
    }else
    {
      tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/whole_genome_ERPRHER2Grade_fixed_baseline.R"),
                     dataopt=dataopt,pop=pop,i1=((i-1)*4000+1):n)
    }
    
    write.table(tmp,file=paste0(outprefix,i,".swarm"),sep="\t",quote=F,row.names = F,col.names = F)
  }
}
sum(sapply(ls(),function(x){object.size(get(x))}))/1E6
#create_swarm()
#create_swarm(prefix="../result/imp_onco/euro/euro",outprefix="scoretest_euro_onco",nvar=2000,dataopt="onco")

#cd swarm
#each split into two files
#8383395
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_euro_icogs1.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
#8383396
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_euro_icogs2.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
#8453986,8466937
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_euro_onco1.swarm -g 15 --module R/4.3 --time=5-00:00:00 --gres=lscratch:15 -p 2
#8453993,8466939
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_euro_onco2.swarm -g 15 --module R/4.3 --time=5-00:00:00 --gres=lscratch:15 -p 2


#asian
#create_swarm(prefix="../result/imp_icogs/asian/asian",outprefix="scoretest_asian_icogs",nvar=2000,dataopt="icogs",pop="asian")
#create_swarm(prefix="../result/imp_onco/asian/asian",outprefix="scoretest_asian_onco",nvar=2000,dataopt="onco",pop="asian")
#cd swarm
#each split into two files
#8902629
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_asian_onco1.swarm -g 10 --module R/4.3 --time=5-00:00:00 --gres=lscratch:10 -p 2
#8902631
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_asian_onco2.swarm -g 10 --module R/4.3 --time=5-00:00:00 --gres=lscratch:10 -p 2
#8902832
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_asian_icogs1.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
#8902898
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_asian_icogs2.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2

#african
#create_swarm(prefix="../result/imp_onco/african/african",outprefix="scoretest_african_onco",nvar=2000,dataopt="onco",pop="african")
#8909070
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_african_onco1.swarm -g 6 --module R/4.3 --time=5-00:00:00 --gres=lscratch:6 -p 2
#8909396
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_african_onco2.swarm -g 6 --module R/4.3 --time=5-00:00:00 --gres=lscratch:6 -p 2
#8909407
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/scoretest_african_onco3.swarm -g 6 --module R/4.3 --time=5-00:00:00 --gres=lscratch:6 -p 2
#create_swarm(prefix="../result/imp_icogs/african/african",outprefix="african_icogs",nvar=2000,dataopt="icogs",pop="african")
#5994333
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/african_icogs1.swarm -g 15 --module R/4.3 --time=5-00:00:00 --gres=lscratch:15 -p 2
#5979669
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/african_icogs2.swarm -g 6 --module R/4.3 --time=5-00:00:00 --gres=lscratch:6 -p 2
#not run swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/african_icogs3.swarm -g 6 --module R/4.3 --time=5-00:00:00 --gres=lscratch:6 -p 2
