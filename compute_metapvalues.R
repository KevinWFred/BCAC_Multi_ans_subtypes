#!/usr/bin/env Rscript
#aggregate logodds and sigma across populations
.libPaths(c("/data/wangx53",.libPaths()))
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
# read_metap=function(pop="euro")
# {
#   allbeta=NULL
#   allsigma=list()
#   for (i in 1:1000)
#   {
#     if (i %% 100 ==0) cat(i,'..')
#     load(paste0("../result/metaresult/",pop,"/res_",i,".RData"))
#     allbeta=rbind(allbeta,allmetares$metares)
#     allsigma=c(allsigma,allmetares$metaallsigma)
#     
#   }
#   res=list(allbeta=allbeta,allsigma=allsigma)
#   save(res,file=paste0("../result/",pop,"_meta_beta_sigma.RData"))
# }
#for euro
pop="euro"
load(paste0("../result/",pop,"_beta_sigma.RData"))
euro_allsigma_icogs=allsigma_icogs
euro_allsigma_onco=allsigma_onco
euro_allbeta_icogs=beta_icogs
euro_allbeta_onco=beta_onco

pop="asian"
load(paste0("../result/",pop,"_beta_sigma.RData"))
asian_allsigma_icogs=allsigma_icogs
asian_allsigma_onco=allsigma_onco
asian_allbeta_icogs=beta_icogs
asian_allbeta_onco=beta_onco
#for african
pop="african"
#save(beta_onco,sigma_onco,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
load(paste0("../result/",pop,"_beta_sigma.RData"))
african_allsigma_onco=allsigma_onco
african_allbeta_onco=beta_onco

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



# LogoddsMetaAnalysis2 <- function(logodds1,sigma1,logodds2,sigma2){
#   sigma1.inv <- solve(sigma1)
#   sigma2.inv <- solve(sigma2)
#   sigma.meta <- solve(sigma1.inv+sigma2.inv)
#   logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2 )
#   
#   return(list(logodds.meta = logodds.meta,
#               sigma.meta = sigma.meta))
# }

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


allvar=unique(c(rownames(euro_allbeta_icogs),rownames(euro_allbeta_onco),
                rownames(asian_allbeta_icogs),rownames(asian_allbeta_onco),
                rownames(african_allbeta_onco)))

args = commandArgs(trailingOnly = T)
print(Sys.time())
print(args)

i1 = as.numeric(args[1]) #block ID, starts with 1
ntask=1000 #total number of jobs
print(i1)
startendidx=startend(length(allvar),ntask,i1)
print(paste0("startendidx: ",startendidx))
idxstart=startendidx[1]
idxend=startendidx[2]
allres=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=1))
rownames(allres)=allvar[idxstart:idxend]
colnames(allres)="P" #combined P-value
euro_allbeta_icogs=euro_allbeta_icogs[rownames(euro_allbeta_icogs) %in% rownames(allres),]
euro_allbeta_onco=euro_allbeta_onco[rownames(euro_allbeta_onco) %in% rownames(allres),]
asian_allbeta_icogs=asian_allbeta_icogs[rownames(asian_allbeta_icogs) %in% rownames(allres),]
asian_allbeta_onco=asian_allbeta_onco[rownames(asian_allbeta_onco) %in% rownames(allres),]
african_allbeta_onco=african_allbeta_onco[rownames(african_allbeta_onco) %in% rownames(allres),]
euro_allsigma_icogs=euro_allsigma_icogs[names(euro_allsigma_icogs) %in% rownames(allres)]
euro_allsigma_onco=euro_allsigma_onco[names(euro_allsigma_onco) %in% rownames(allres)]
asian_allsigma_icogs=asian_allsigma_icogs[names(asian_allsigma_icogs) %in% rownames(allres)]
asian_allsigma_onco=asian_allsigma_onco[names(asian_allsigma_onco) %in% rownames(allres)]
african_allsigma_onco=african_allsigma_onco[names(african_allsigma_onco) %in% rownames(allres)]
gc()

for (i in 1:nrow(allres))
{
  if (i %% 1000 ==0) cat(i,'..')
  varname=rownames(allres)[i]
  logoddslist=sigmalist=list()
  if (varname %in% rownames(euro_allbeta_icogs))
  {
    logoddslist=c(logoddslist,list(euro_allbeta_icogs[which(rownames(euro_allbeta_icogs)==varname),]))
    sigmalist=c(sigmalist,euro_allsigma_icogs[which(names(euro_allsigma_icogs)==varname)])
  }
  if (varname %in% rownames(euro_allbeta_onco))
  {
    logoddslist=c(logoddslist,list(euro_allbeta_onco[which(rownames(euro_allbeta_onco)==varname),]))
    sigmalist=c(sigmalist,euro_allsigma_onco[which(names(euro_allsigma_onco)==varname)])
  }

  if (varname %in% rownames(asian_allbeta_icogs))
  {
    logoddslist=c(logoddslist,list(asian_allbeta_icogs[which(rownames(asian_allbeta_icogs)==varname),]))
    sigmalist=c(sigmalist,asian_allsigma_icogs[which(names(asian_allsigma_icogs)==varname)])
  }
  if (varname %in% rownames(asian_allbeta_onco))
  {
    logoddslist=c(logoddslist,list(asian_allbeta_onco[which(rownames(asian_allbeta_onco)==varname),]))
    sigmalist=c(sigmalist,asian_allsigma_onco[which(names(asian_allsigma_onco)==varname)])
  }
  if (varname %in% rownames(african_allbeta_onco))
  {
    logoddslist=c(logoddslist,list(african_allbeta_onco[which(rownames(african_allbeta_onco)==varname),]))
    sigmalist=c(sigmalist,african_allsigma_onco[which(names(african_allsigma_onco)==varname)])
  }
  
  metares=LogoddsMetaAnalysis(logoddslist,sigmalist)
  allres$P[i]=GlobalTestForAssoc(metares$logodds.meta,metares$sigma.meta)
}
save(allres,file=paste0("../result/metaresult/all/","res_",i1,".RData"))
print(Sys.time())
print("done")

# swarmjobs=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_metapvalues.R",1000),
#                      i1=1:1000)
# write.table(swarmjobs,file=paste0("compute_metapvalues.swarm"),row.names=F,col.names=F,sep="\t",quote=F)
#9757248
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_metapvalues.swarm -g 56 --module R/4.3 --partition=quick --time=4:00:00 --gres=lscratch:40 -p 2

#collect the results
allpvalues=NULL
for (i in 1:1000)
{
  if (i %% 100 ==0) cat(i,'..')
  load(paste0("../result/metaresult/all","/res_",i,".RData"))
  allpvalues=rbind(allpvalues,allres)
}
intrinsicpvalues=allpvalues

load("../result/scoretestresult/allpvalues.Rdata")
scoretestallpvaues=data.frame(P=allpvalues)
rownames(scoretestallpvaues)=names(allpvalues)
tmp=intersect(rownames(scoretestallpvaues),rownames(intrinsicpvalues))
library(ACAT)
allpvalues=cbind(scoretestallpvaues[match(tmp,rownames(scoretestallpvaues)),],intrinsicpvalues[match(tmp,rownames(intrinsicpvalues)),])
idx=complete.cases(allpvalues)
allpvalues=allpvalues[idx,]
allpvalues=ACAT(t(allpvalues))
names(allpvalues)=tmp[idx]

#qqplot,manhattan plot
load("../result/metascoreinfo4.RData")

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
  if (optbreak==1)
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

  #outpath <-"../result"


  ggsave(filename=filename,
         plot=manhplot, device="png",
         width=9, height=4, units="in", dpi=300)
}
table(rownames(allpvalues0) %in% names(metascoreinfo4$freq))

myplot=function(myallpvalues=allpvalues,outprefix="waldtest_meta")
{
  idx0=match(names(myallpvalues),names(metascoreinfo4$freq))
  all(names(metascoreinfo4$freq[idx0])==names(myallpvalues)) #T

  tmp=unlist(strsplit(names(myallpvalues),":"))
  data=data.frame(SNP=names(myallpvalues),CHR=tmp[seq(1,length(tmp),4)],BP=tmp[seq(2,length(tmp),4)],EAF=metascoreinfo4$freq[idx0],P=myallpvalues,N=0)
  data$CHR=gsub("chr","",data$CHR)
  data$BP=as.integer(data$BP)
  data$CHR=as.integer(data$CHR)
  #euronco
  idx=which(metascoreinfo4$avail[idx0,1]==1)
  data$N[idx]=data$N[idx]+1/(1/34367+1/52500)
  #euricogs
  idx=which(metascoreinfo4$avail[idx0,2]==1)
  data$N[idx]=data$N[idx]+1/(1/13283+1/33505)
  #asianicogs
  idx=which(metascoreinfo4$avail[idx0,3]==1)
  data$N[idx]=data$N[idx]+1/(1/1682+1/5940)
  idx=which(metascoreinfo4$avail[idx0,4]==1)
  data$N[idx]=data$N[idx]+1/(1/9064+1/12344)
  idx=which(metascoreinfo4$avail[idx0,5]==1)
  data$N[idx]=data$N[idx]+1/(1/1815+1/2088)
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
  png(filename = paste0("../result/QQplot_",outprefix,".png"), width = 8, height = 8, units = "in",res=300)
  plotqq(data=data1,optbreak=1,title="")
  dev.off()
  plotmanhattan(data=data1,title="",filename=paste0("../result/man_",outprefix,".png"))
}

myplot(myallpvalues=allpvalues,outprefix="ACAT_combined")
tmp=intrinsicpvalues$P
names(tmp)=rownames(intrinsicpvalues)
myplot(myallpvalues=tmp,outprefix="Intrinsic_meta")
