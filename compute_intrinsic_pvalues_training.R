#!/usr/bin/env Rscript
#intrinsic model p-values
#meta analysis within a population (for intrinsic_subtypes_genome.R) for training samples.
.libPaths(c("/data/wangx53",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(data.table)

#For euro and asian
readres=function(pop="euro",nvar=2000,optreadicogs=T)
{
  print(Sys.time())
  #for pvar
  prefixonco=paste0("../result/imp_onco/",pop,"/",pop)
  prefixicogs=paste0("../result/imp_icogs/",pop,"/",pop)
  
  outprefix=paste0("../result/training/",pop)
  pvar=as.data.frame(fread(paste0(prefixicogs,".pvar")))
  nicogs=as.integer(nrow(pvar)/nvar)+1
  pvar=as.data.frame(fread(paste0(prefixonco,".pvar")))
  nonco=as.integer(nrow(pvar)/nvar)+1
  #beta_icogs=sigma_icogs=data.frame(matrix(NA,nrow=nicogs*nvar,ncol=5))
  #beta_onco=sigma_onco=data.frame(matrix(NA,nrow=nonco*nvar,ncol=5))
  allsigma_icogs=list()
  beta_icogs=sigma_icogs=NULL
  allsigma_onco=list()
  beta_onco=sigma_onco=NULL
  #for training res
  prefixonco1=paste0("../result/training/imp_onco/",pop,"/",pop)
  prefixicogs1=paste0("../result/training/imp_icogs/",pop,"/",pop)
  icogsresfolder=paste0(dirname(prefixicogs1),"/res/")
  oncoresfolder=paste0(dirname(prefixonco1),"/res/")
  for (i in 1:nonco)
  {
    if (i %% 500 ==0) cat(i,'..')
    resfile=paste0(oncoresfolder,"res_",i,".RData")
    if (file.exists(resfile))
    {
      loaddata = tryCatch(
        expr = {
          load(resfile)
        },
        error = function(e){ 
          return(NULL)
        })
      if (!is.null(loaddata))
      {
        idx=match(names(all.sigma.log.odds),rownames(all.log.odds))
        all.log.odds=all.log.odds[idx,]
        allsigma_onco=c(allsigma_onco,all.sigma.log.odds)
        beta_onco=rbind(beta_onco,all.log.odds)
        tmp=do.call("rbind",all.sigma.log.odds)
        colid=rep(1:5,times=nrow(all.log.odds))
        rowid=1:nrow(tmp)
        idx=as.matrix(data.frame(row=rowid,col=colid))
        tmp1=tmp[idx]
        tmp2=matrix(tmp1,ncol=5,byrow = T)
        rownames(tmp2)=rownames(all.log.odds)
        sigma_onco=rbind(sigma_onco,tmp2)
      }else
      {
        warning(paste0(resfile," is missing"))
      }
    }
  }
  print(dim(beta_onco))
  print(Sys.time())
  if (optreadicogs)
  {
    for (i in 1:nicogs)
    {
      if (i %% 500 ==0) cat(i,'..')
      resfile=paste0(icogsresfolder,"res_",i,".RData")
      if (file.exists(resfile))
      {
        loaddata = tryCatch(
          expr = {
            load(resfile)
          },
          error = function(e){ 
            return(NULL)
          })
        if (!is.null(loaddata))
        {
          idx=match(names(all.sigma.log.odds),rownames(all.log.odds))
          all.log.odds=all.log.odds[idx,]
          allsigma_icogs=c(allsigma_icogs,all.sigma.log.odds)
          beta_icogs=rbind(beta_icogs,all.log.odds)
          tmp=do.call("rbind",all.sigma.log.odds)
          colid=rep(1:5,times=nrow(all.log.odds))
          rowid=1:nrow(tmp)
          idx=as.matrix(data.frame(row=rowid,col=colid))
          tmp1=tmp[idx]
          tmp2=matrix(tmp1,ncol=5,byrow = T)
          rownames(tmp2)=rownames(all.log.odds)
          sigma_icogs=rbind(sigma_icogs,tmp2)
        }else
        {
          warning(paste0(resfile," is missing"))
        }
      }
    }
    print(dim(beta_icogs))
    save(beta_icogs,beta_onco,sigma_icogs,sigma_onco,allsigma_icogs,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
  }else #not to save icogs results
  {
    save(beta_onco,sigma_onco,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
  }
  print(Sys.time())
}

#adapted from comput_metapvalues.R (based on QCed data),compute beta/p within training data on EUR and EAS separately.
#the results will be used for subtype PRS
.libPaths(c("/data/wangx53",.libPaths()))
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(data.table)
library(bc2)
startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}

#QCed SNP lists:
qc_asian_onco=as.data.frame(fread("../result/imp_QC/onco/asian_info.snplist",header=F))
qc_euro_onco=as.data.frame(fread("../result/imp_QC/onco/euro_info.snplist",header=F))
qc_asian_icogs=as.data.frame(fread("../result/imp_QC/icogs/asian_info.snplist",header=F))
qc_euro_icogs=as.data.frame(fread("../result/imp_QC/icogs/euro_info.snplist",header=F))

args = commandArgs(trailingOnly = T)
print(Sys.time())
print(args)
i1 = as.numeric(args[1]) #block ID, starts with 1
pop=as.character(args[2]) #euro or asian

#for euro
if(pop=="euro")
{
  load(paste0("../result/training/",pop,"_beta_sigma.RData"))
  euro_allsigma_icogs=allsigma_icogs
  euro_allsigma_icogs=euro_allsigma_icogs[names(euro_allsigma_icogs) %in% qc_euro_icogs[,1]]
  euro_allsigma_onco=allsigma_onco
  euro_allsigma_onco=euro_allsigma_onco[names(euro_allsigma_onco) %in% qc_euro_onco[,1]]
  euro_allbeta_icogs=beta_icogs
  euro_allbeta_icogs=euro_allbeta_icogs[rownames(euro_allbeta_icogs) %in% qc_euro_icogs[,1],]
  all(rownames(euro_allbeta_icogs)==names(euro_allsigma_icogs))
  euro_allbeta_onco=beta_onco
  euro_allbeta_onco=euro_allbeta_onco[rownames(euro_allbeta_onco) %in% qc_euro_onco[,1],]
  all(rownames(euro_allbeta_onco)==names(euro_allsigma_onco))
}else
{
  
  load(paste0("../result/training/",pop,"_beta_sigma.RData"))
  asian_allsigma_icogs=allsigma_icogs
  asian_allsigma_icogs=asian_allsigma_icogs[names(asian_allsigma_icogs) %in% qc_asian_icogs[,1]]
  asian_allsigma_onco=allsigma_onco
  asian_allsigma_onco=asian_allsigma_onco[names(asian_allsigma_onco) %in% qc_asian_onco[,1]]
  asian_allbeta_icogs=beta_icogs
  asian_allbeta_icogs=asian_allbeta_icogs[rownames(asian_allbeta_icogs) %in% qc_asian_icogs[,1],]
  all(rownames(asian_allbeta_icogs)==names(asian_allsigma_icogs))
  asian_allbeta_onco=beta_onco
  asian_allbeta_onco=asian_allbeta_onco[rownames(asian_allbeta_onco) %in% qc_asian_onco[,1],]
  all(rownames(asian_allbeta_onco)==names(asian_allsigma_onco)) 
}

rm(allsigma_icogs,allsigma_onco,beta_icogs,beta_onco,sigma_icogs,sigma_onco)
gc()
#check how much size used
#sum( sapply(ls(),function(x){object.size(get(x))}))/1E6

LogoddsMetaAnalysis <- function(logoddslist,sigmalist){
  if (length(logoddslist)>1)
  {
    sigma.inv.sum=0
    sigma.inv.logodds.sum=0
    for(i in 1:length(sigmalist))
    {
      sigma.inv.sum=sigma.inv.sum+solve(sigmalist[[i]])
      sigma.inv.logodds.sum=sigma.inv.logodds.sum+solve(sigmalist[[i]])%*%unlist(logoddslist[[i]])
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


#Global test for association using logodds and sigma
GlobalTestForAssoc <- function(logodds,sigma){
  sigma <- as.matrix(sigma)
  df <- length(logodds)
  GTA.stat <- t(logodds)%*%solve(sigma)%*%logodds
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  if(p.value.GTA!=0)
  {
    places = 3
    power.number <- floor(-log10(p.value.GTA))+places
    ###format the output with three digits in total
    p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)
  }
  return(p.value.GTA)
}

if (pop=="euro")
{
  allvar=unique(c(rownames(euro_allbeta_icogs),rownames(euro_allbeta_onco)))
}else
{
  allvar=unique(c(rownames(asian_allbeta_icogs),rownames(asian_allbeta_onco)))
}

ntask=1000 #total number of jobs
print(i1)
startendidx=startend(length(allvar),ntask,i1)
print(paste0("startendidx: ",startendidx))
idxstart=startendidx[1]
idxend=startendidx[2]
infor=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=5))
beta=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=5))
beta_se=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=5))
p=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=5))
rownames(infor)=rownames(beta)=rownames(beta_se)=rownames(p)=allvar[idxstart:idxend]
colnames(infor)=c("chr","pos","SNP","a1","a0") #combined P-value
tmp=allvar[idxstart:idxend] #names in the format of chr:pos:alt:ref (pvar)
tmp1=unlist(strsplit(tmp,":"))
infor$chr=tmp1[seq(1,length(tmp1),4)]
infor$chr=gsub("chr","",infor$chr)
infor$pos=as.integer(tmp1[seq(2,length(tmp1),4)])
infor$a1=tmp1[seq(3,length(tmp1),4)]
infor$a0=tmp1[seq(4,length(tmp1),4)]
infor$SNP=tmp
if (pop=="euro")
{
  euro_allbeta_icogs=euro_allbeta_icogs[rownames(euro_allbeta_icogs) %in% rownames(infor),]
  euro_allbeta_onco=euro_allbeta_onco[rownames(euro_allbeta_onco) %in% rownames(infor),]
  euro_allsigma_icogs=euro_allsigma_icogs[names(euro_allsigma_icogs) %in% rownames(infor)]
  euro_allsigma_onco=euro_allsigma_onco[names(euro_allsigma_onco) %in% rownames(infor)]
}else
{
  asian_allbeta_icogs=asian_allbeta_icogs[rownames(asian_allbeta_icogs) %in% rownames(infor),]
  asian_allbeta_onco=asian_allbeta_onco[rownames(asian_allbeta_onco) %in% rownames(infor),]
  asian_allsigma_icogs=asian_allsigma_icogs[names(asian_allsigma_icogs) %in% rownames(infor)]
  asian_allsigma_onco=asian_allsigma_onco[names(asian_allsigma_onco) %in% rownames(infor)]
  
}

gc()

for (i in 1:nrow(infor))
{
  if (i %% 1000 ==0) cat(i,'..')
  varname=rownames(infor)[i]
  logoddslist=sigmalist=list()
  if (pop=="euro")
  {
    if (nrow(euro_allbeta_icogs)>0)
    {
      if (varname %in% rownames(euro_allbeta_icogs))
      {
        logoddslist=c(logoddslist,list(euro_allbeta_icogs[which(rownames(euro_allbeta_icogs)==varname),]))
        sigmalist=c(sigmalist,euro_allsigma_icogs[which(names(euro_allsigma_icogs)==varname)])
      }
    }
    if (nrow(euro_allbeta_onco)>0)
    {
      if (varname %in% rownames(euro_allbeta_onco))
      {
        logoddslist=c(logoddslist,list(euro_allbeta_onco[which(rownames(euro_allbeta_onco)==varname),]))
        sigmalist=c(sigmalist,euro_allsigma_onco[which(names(euro_allsigma_onco)==varname)])
      }
    }
    
  }else
  {
    if (nrow(asian_allbeta_icogs)>0)
    {
      if (varname %in% rownames(asian_allbeta_icogs))
      {
        logoddslist=c(logoddslist,list(asian_allbeta_icogs[which(rownames(asian_allbeta_icogs)==varname),]))
        sigmalist=c(sigmalist,asian_allsigma_icogs[which(names(asian_allsigma_icogs)==varname)])
      }
    }
    if (nrow(asian_allbeta_onco)>0)
    {
      if (varname %in% rownames(asian_allbeta_onco))
      {
        logoddslist=c(logoddslist,list(asian_allbeta_onco[which(rownames(asian_allbeta_onco)==varname),]))
        sigmalist=c(sigmalist,asian_allsigma_onco[which(names(asian_allsigma_onco)==varname)])
      }
    }
    
  }
  
  metares=LogoddsMetaAnalysis(logoddslist,sigmalist)
  beta[i,]=metares$logodds.meta
  beta_se[i,]=sqrt(diag(metares$sigma.meta)) #check this
  z = metares$logodds.meta/sqrt(diag(metares$sigma.meta))
  #function to calculate P-value
  p[i,]=2*pnorm(-abs(z), lower.tail = T)
}
save(infor,beta,beta_se,p,file=paste0("../result/training/metaresult/",pop,"/","QCres_",i1,".RData"))
print(Sys.time())
print("done")

swarmjobs1=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_intrinsic_pvalues_training.R",1000),
                     i1=1:1000,pop="euro")
swarmjobs2=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_intrinsic_pvalues_training.R",1000),
                     i1=1:1000,pop="asian")
swarmjobs=rbind(swarmjobs1,swarmjobs2)
write.table(swarmjobs,file=paste0("compute_intrinsic_pvalues_training.swarm"),row.names=F,col.names=F,sep="\t",quote=F)
#
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_intrinsic_pvalues_training.swarm -g 42 --module R/4.3 --partition=quick --time=4:00:00 --gres=lscratch:40 -p 2
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_intrinsic_pvalues_training.swarm -g 42 --module R/4.3 --time=4:00:00 --gres=lscratch:40 -p 2

#collect the results
subtypes=c("LumA","LumB","LumB_HN","Her2E","TripN")
allinfor=allp=allbeta=allbeta_se=NULL
pop="euro"
for (i in 1:1000)
{
  if (i %% 100 ==0) cat(i,'..')
  load(paste0("../result/training/metaresult/",pop,"/QCres_",i,".RData"))
  allinfor=rbind(allinfor,infor)
  allp=rbind(allp,p)
  allbeta=rbind(allbeta,beta)
  allbeta_se=rbind(allbeta_se,beta_se)
  rm(infor,beta,beta_se,p)
}
colnames(allp)=paste0("p_",subtypes)
colnames(allbeta)=paste0("beta_",subtypes)
colnames(allbeta_se)=paste0("beta_se_",subtypes)
allres=cbind.data.frame(allinfor,allbeta,allbeta_se,allp)
save(allres,file=paste0("../result/PRS_subtype/",pop,"/subtype_sumstats.RData"))

allinfor=allp=allbeta=allbeta_se=NULL
pop="asian"
for (i in 1:1000)
{
  if (i %% 100 ==0) cat(i,'..')
  load(paste0("../result/training/metaresult/",pop,"/QCres_",i,".RData"))
  allinfor=rbind(allinfor,infor)
  allp=rbind(allp,p)
  allbeta=rbind(allbeta,beta)
  allbeta_se=rbind(allbeta_se,beta_se)
  rm(infor,beta,beta_se,p)
}
colnames(allp)=paste0("p_",subtypes)
colnames(allbeta)=paste0("beta_",subtypes)
colnames(allbeta_se)=paste0("beta_se_",subtypes)
allres=cbind.data.frame(allinfor,allbeta,allbeta_se,allp)
save(allres,file=paste0("../result/PRS_subtype/",pop,"/subtype_sumstats.RData"))

load("../result/metascoreinfo4_newQC.RData")

library(data.table)
library("plotrix") #axis.break
library("RColorBrewer")
library("optparse")
library(readr)
library(dplyr)
library(ggplot2)
library(qqman)

source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
#draw QQ plot
qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}
#optbreak=1, break top p-values, used for very low p-value cases
plotqq=function(data,optbreak=1,title="")
{
  dat = data %>%
    mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>%
    select(rsid,CHR,BP,P,MAF) %>%
    rename(SNP = rsid)
  dat=dat[!is.na(dat$P),]
  x = dat$P
  z = qnorm(x / 2)
  lambda = round(median(z^2) / qchisq(0.5,1), 3)
  N.effect = median(data$N)
  lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
  
  yLine <- c(-log10(5E-8))
  colLine <- c("red")
  dat$log10P = -log10(dat$P)
  gwas = as.data.frame(dat)
  # Determine frequency bins and create variable for binned QQ plot
  
  minMAF <- min(gwas$MAF)
  
  freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
  gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
  freqtable <- table(gwas$freqbin)
  freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
  freqtable <- freqtable[freqtable > 0]
  
  ## Generate QQ plot data by frequency bin
  fbin <- character(0)
  fN <- integer(0)
  fx <- numeric(0)
  fy <- numeric(0)
  fcol <- character(0)
  legendcol <- character(0)
  conf <- list()
  allcols <- brewer.pal(4,"Set1")
  ycol <- "log10P"
  for(f in 1:length(freqtable)){
    fbin <- c(fbin,names(freqtable)[f])
    fsnps <- which(gwas$freqbin ==names(freqtable)[f])
    plotdata <- qqplotdata(gwas[[ycol]][fsnps])
    fN <- c(fN,freqtable[f])
    fx <- c(fx,plotdata$e)
    fy <- c(fy,plotdata$o)
    fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
    conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                            'y'=c(plotdata$c975,rev(plotdata$c025)))
    legendcol <- c(legendcol,allcols[f])
  }
  legendtext <- paste0("MAF=",fbin,"; #SNPs=",format(fN,big.mark=",",scientific=FALSE))
  
  opt =  list(break.top = ifelse(optbreak==1,15,ceiling(max(fy))+1),
              top.size = 0.125)
  
  
  #png(filename = paste0(outpath,"/QQ_",eth[i1],"_",trait[i2],".png"), width = 8, height = 8, units = "in",res=300)
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)
  
  par(mar=c(5.1,5.1,4.1,1.1))
  
  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  #top
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  if (any(lab2>max(lab1))) lab2 <- lab2[lab2 > max(lab1)]
  
  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top
  
  # function to rescale the top part
  rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
  
  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
  # Plot confidence intervals
  for(p in 1:length(conf)){
    polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
            col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
            border = NA)
  }
  
  # add points below top
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
  
  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  box()
  if (optbreak==1 & max(fy,na.rm=T)>opt$break.top)
  {
    rescaled.y = rescale(fy[fy>opt$break.top])
    par(las=0)
    points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
    par(las=1)
    axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
    axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
    axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
    lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  }
  
  abline(h=ifelse(yLine<opt$break.top,
                  yLine,
                  rescale(yLine)),
         col=colLine,lwd=1.5,lty=2)
  legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
  text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
  text(5.9,1,paste(lambda_1000),cex = 1.5)
  
  title(title)
  #dev.off()
  
}

plotmanhattan=function(data,title="",filename="../result/man_ILCCOnever.png")
{
  dat = data %>%
    mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>%
    select(rsid,CHR,BP,P,MAF) %>%
    rename(SNP = rsid)
  dat$CHR1=factor(dat$CHR,levels = c(1:22,"X","Y"))
  idx=order(dat$CHR1,dat$BP)
  dat=dat[idx,]
  p.pwas <- 5E-08
  
  nCHR <- length(unique(dat$CHR))
  dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(dat$CHR)){
    nbp[i] <- max(dat[dat$CHR == i,]$BP)
    dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
    s <- s + nbp[i]
  }
  
  axis.set <- dat %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(dat$P)))) + 2
  sig1 <- p.pwas
  
  
  sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
  library(ggplot2)
  manhplot <- ggplot(dat, aes(x = BPcum, y = -log10(P),
                              color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.8, size=0.8) +
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
    guides(color = FALSE) +
    labs(x = NULL,
         y = "-log10(p)",
         linetype = "",
         title = title)+
    #subtitle = "A2: Critically ill COVID19+ vs. population controls;\nB1: Hospitalized COVID19+ vs non-hospitalized COVID19+;\nB2: Hospitalized COVID19+ vs. population controls;\nC2: Reported SARS-CoV-2 infection vs. population controls") +
    theme_Publication()+
    theme(
      legend.position = "top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  ggsave(filename=filename,
         plot=manhplot, device="png",
         width=12, height=4, units="in", dpi=300)
  #outpath <-"../result"
  # myknownvar=knownvar
  # myallnovelsnps=allnovelsnps
  # if (!is.null(myknownvar))
  # {
  #   gr_myknownvar=GRanges(seqnames=myknownvar$CHR,ranges = IRanges(myknownvar$hg38position,width=1))
  #   gr_dat=GRanges(seqnames = dat$CHR,ranges = IRanges(start=dat$BP,width = 1))
  #   tmp=distanceToNearest(gr_myknownvar,gr_dat)
  #   manhplot1=manhplot+
  #     geom_vline(data = dat[tmp@to,], aes(xintercept = BPcum), color= "green", linetype="solid",size=0.3,alpha=0.3)
  #   filename1=gsub(".png","",filename)
  #   filename1=paste0(filename1,"_knownvar.png")
  #   ggsave(filename=filename1,
  #          plot=manhplot1, device="png",
  #          width=12, height=4, units="in", dpi=300)
  #   gr_myallnovelsnps=GRanges(seqnames = myallnovelsnps$chr,ranges = IRanges(start=myallnovelsnps$pos,width=1))
  #   tmp=distanceToNearest(gr_myallnovelsnps,gr_dat)
  #   manhplot2=manhplot+
  #     geom_vline(data = dat[tmp@to,], aes(xintercept = BPcum), color= "red", linetype="solid",size=0.3,alpha=0.3)
  #   filename1=gsub(".png","",filename)
  #   filename2=paste0(filename1,"_novelvar.png")
  #   ggsave(filename=filename2,
  #          plot=manhplot2, device="png",
  #          width=12, height=4, units="in", dpi=300)
  # }
  # 
  
}
#table(rownames(allpvalues0) %in% names(metascoreinfo4$freq))

myplot=function(myallpvalues=allpvalues,outprefix="../result/PRS_subtype/euro_LumA",N=1000)
{
  idx0=match(names(myallpvalues),names(metascoreinfo4$freq))
  all(names(metascoreinfo4$freq[idx0])==names(myallpvalues)) #T
  
  tmp=unlist(strsplit(names(myallpvalues),":"))
  data=data.frame(SNP=names(myallpvalues),CHR=tmp[seq(1,length(tmp),4)],BP=tmp[seq(2,length(tmp),4)],EAF=metascoreinfo4$freq[idx0],P=myallpvalues,N=0)
  data$CHR=gsub("chr","",data$CHR)
  data$BP=as.integer(data$BP)
  data$CHR=as.integer(data$CHR)
  data$N=N
  data$EAF[is.na(data$EAF)]=0.01
  
  data1 = data %>%
    mutate(CHR = as.integer(CHR),
           BP = as.integer(BP),
           FREQ_A1 = as.numeric(EAF),
           P = as.numeric(P),
           rsid = SNP,
           N = as.integer(N)) %>%
    mutate(P = ifelse(P==0,1E-300,P)) %>%
    select(rsid,CHR,BP,FREQ_A1,P,N)
  data1=data1[data1$FREQ_A1>0.01 & data1$FREQ_A1<0.99,]
  png(filename = paste0(outprefix,"_qqplot.png"), width = 8, height = 8, units = "in",res=300)
  plotqq(data=data1,optbreak=1,title="")
  dev.off()
  plotmanhattan(data=data1,title="",filename=paste0(outprefix,"_man.png"))
}

#compute effective sample sizes for each pop/subtype
#training samples
trainingsamples_icogs=read.table("../result/PRS_icogs_trainingpheno.txt",sep="\t",header=T)
trainingsamples_onco=read.table("../result/PRS_onco_trainingpheno.txt",sep="\t",header=T)
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
compute_train_effsize_subtype=function(pop="euro")
{
  phenoicogs1=phenoicogs[phenoicogs$ID %in% trainingsamples_icogs$SG_ID,]
  phenoonco1=phenoonco[phenoonco$ID %in% trainingsamples_onco$Onc_ID,]
  if (pop=="euro")
  {
    phenoicogs1=phenoicogs1[phenoicogs1$EthnicityGeno %in% "European",]
    phenoonco1=phenoonco1[phenoonco1$EthnicityGeno %in% "European",]
  }
  if (pop=="asian")
  {
    phenoicogs1=phenoicogs1[phenoicogs1$EthnicityGeno %in% "Asian",]
    phenoonco1=phenoonco1[phenoonco1$EthnicityGeno %in% "Asian",]
  }
  effsize=rep(NA,length(subtypes))
  names(effsize)=subtypes
  for (i in 1:length(subtypes))
  {
    subtype=subtypes[i]
    effsize[i]=1/(1/sum(phenoicogs1[,subtype]==0,na.rm = T)+1/(sum(phenoicogs1[,subtype]==1,na.rm = T)))+
      1/(1/sum(phenoonco1[,subtype]==0,na.rm = T)+1/(sum(phenoonco1[,subtype]==1,na.rm = T)))
  }
  return(effsize)
}
euro_effsize=compute_train_effsize_subtype(pop="euro")
asian_effsize=compute_train_effsize_subtype(pop="asian")
for (pop in c("euro","asian"))
{
  load(paste0("../result/PRS_subtype/",pop,"/subtype_sumstats.RData"))
  for (i in 1:length(subtypes))
  {
    subtype=subtypes[i]
    allpvalues=allres[,paste0("p_",subtype)]
    names(allpvalues)=allres$SNP
    if (pop=="euro")
    {
      myplot(myallpvalues=allpvalues,outprefix=paste0("../result/PRS_subtype/",pop,"_",subtype),N=euro_effsize[i])
    }else
    {
      myplot(myallpvalues=allpvalues,outprefix=paste0("../result/PRS_subtype/",pop,"_",subtype),N=asian_effsize[i])
    }
    
  }
}
