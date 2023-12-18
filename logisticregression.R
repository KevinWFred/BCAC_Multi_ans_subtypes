#!/usr/bin/env Rscript

.libPaths(c("/data/wangx53",.libPaths()))

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

training_icogs=read.table("../result/PRS_icogs_trainingpheno.txt",sep="\t",header = T)
training_icogs$ID=training_icogs$SG_ID
icogs_samplefile="../result/icogs_samples_750.txt"
icogs_samples=read.table(icogs_samplefile)
all(training_icogs$SG_ID %in% icogs_samples$V1)
tmp=read.table("../result/imp_icogs/merged1.eigenvec")
colnames(tmp)=c("ID",paste0("PC_",1:20))
tmp1=unlist(strsplit(tmp$ID,"_"))
tmp$ID=tmp1[seq(1,length(tmp1),2)]
idx=match(training_icogs$ID,tmp$ID)
training_icogs=cbind(training_icogs,tmp[idx,2:11])

training_onco=read.table("../result/PRS_onco_trainingpheno.txt",header = T,sep="\t")
training_onco$ID=training_onco$Onc_ID
onco_samplefile="../result/onco_samples_750.txt"
onco_samples=read.table(onco_samplefile)
all(training_onco$Onc_ID %in% onco_samples$V1)

rungwas=function(phenodat=training_icogs[training_icogs$EthnicityGeno=="European",],
                 genoprefix="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro",
                 outprefix="../result/GWAS_icogs_euro")
{
  
  #generate sample files
  tmp=data.frame(FID=phenodat$ID,IID=phenodat$ID)
  samplefile=paste0(outprefix,".samples")
  write.table(tmp,file=samplefile,row.names = F,quote=F,sep=" ")
  
  #generate covariate files----------
  tmp=data.frame(FID=phenodat$ID,IID=phenodat$ID)
  tmp=cbind(tmp,phenodat[,c("age",paste0("PC_",1:10))])
  covfile=paste0(outprefix,".covar")
  write.table(tmp,file=covfile,row.names = F,col.names=T,sep=" ",quote=F)
  
  #generate phenotype file
  tmp=data.frame(FID=phenodat$ID,IID=phenodat$ID,case=as.integer(!is.na(phenodat$Behaviour1))+1)
  phenofile=paste0(outprefix,".pheno")
  write.table(tmp,file=phenofile,row.names = F,quote=F,sep=" ")
  
  cmd=paste0(plink2," --pfile ",genoprefix," --keep ",samplefile," --pheno ",phenofile," --covar ",covfile," --logistic hide-covar --ci 0.95 --out ",outprefix," --threads 8 --memory 64000 --covar-variance-standardize")
  system(cmd)
}
rungwas()
rungwas(phenodat=training_icogs[training_icogs$EthnicityGeno=="Asian",],
                 genoprefix="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian",
                 outprefix="../result/GWAS_icogs_asian")
rungwas(phenodat=training_onco[training_onco$EthnicityGeno=="European",],
                 genoprefix="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro",
                 outprefix="../result/GWAS_onco_euro")
rungwas(phenodat=training_onco[training_onco$EthnicityGeno=="Asian",],
        genoprefix="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian",
        outprefix="../result/GWAS_onco_asian")

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
addrsid=function(sumfile="../result/GWAS_icogs_euro.pheno.glm.logistic.hybrid")
{
  print(sumfile)
  dat=as.data.frame(fread(sumfile))
  colnames(dat)[which(colnames(dat)=="#CHROM")]="CHR"
  dat=dat[!is.na(dat$P),]
  positions <- GPos(seqnames = dat$CHR, pos = dat$POS)
  all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
  ## query the genome with positions
  my_snps <- snpsByOverlaps(all_snps, positions, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  ## this gives us a GPos object
  my_snps = as.data.frame(my_snps)
  
  idx = match(paste0(dat$CHR,":",dat$POS),paste0(my_snps$seqnames,":",my_snps$pos))
  dat$rsid=my_snps$RefSNP_id[idx]
  print(paste0("there are ",sum(is.na(idx))/nrow(dat)," variants without rsid"))
  tmp=dat[order(dat$P),]
  tmp1=sum(is.na(tmp$rsid[1:10]))
  if (tmp1>0) print(paste0("there are ",tmp1," variant without rsid in top 10"))
  #dat=as.data.frame(fread(paste0(sumfile,".rsid")))
  idx=which(colnames(dat)=="POS")
  if (length(idx)==0) stop("no pos")
  colnames(dat)[idx]="BP"
  dat$A2=dat$REF
  idx=which(dat$A1==dat$REF)
  dat$A2[idx]=dat$ALT[idx]
  idx=which(colnames(dat)=="LOG(OR)_SE")
  if (length(idx)==0) stop("no LOG(OR)_SE")
  colnames(dat)[idx]="SE"
  #dat$SNP=dat$rsid
  dat$SNP=dat$ID
  write.table(dat,file=paste0(sumfile,".rsid"),row.names=F,sep="\t",quote=F)
  #return(dat)
}

addrsid(sumfile="../result/GWAS_icogs_euro.pheno.glm.logistic.hybrid")
addrsid(sumfile="../result/GWAS_onco_euro.pheno.glm.logistic.hybrid")
addrsid(sumfile="../result/GWAS_icogs_asian.pheno.glm.logistic.hybrid")
addrsid(sumfile="../result/GWAS_onco_asian.pheno.glm.logistic.hybrid")

runmeta=function(sumfiles=c("../result/GWAS_icogs_euro.pheno.glm.logistic.hybrid.rsid",
                            "../result/GWAS_icogs_asian.pheno.glm.logistic.hybrid.rsid",
                            "../result/GWAS_onco_euro.pheno.glm.logistic.hybrid.rsid",
                            "../result/GWAS_onco_asian.pheno.glm.logistic.hybrid.rsid"),
                 outprefix="../result/meta_lr")
{
  cmd=paste0(plink," --meta-analysis ",paste0(sumfiles,collapse = " ")," + report-all --threads 8 --memory 100000 --out ",outprefix)
  system(cmd)
  metares=as.data.frame(fread(paste0(outprefix,".meta"))) 
  metares=metares[!is.na(metares$SNP),]
  idx=which(is.na(metares$P))
  if (length(idx)>0)
  {
    onestudysnps=metares$SNP[idx]
    nonestudysnps=rep(0,length(idx))
    for (i in 1:length(sumfiles))
    {
      if (length(onestudysnps)>0)
      {
        gwas=as.data.frame(fread(sumfiles[i]))
        tmp=intersect(gwas$SNP,onestudysnps)
        if (length(tmp)>0)
        {
          idx1=match(tmp,gwas$SNP)
          idx2=match(tmp,metares$SNP)
          metares$P[idx2]=gwas$P[idx1]
          metares$OR[idx2]=gwas$OR[idx1]
          onestudysnps=onestudysnps[!onestudysnps %in% tmp]
        }
      }
    }
  }
  idx=which(is.na(metares$P))
  print(length(idx))
  metares=metares[!is.na(metares$P),]
  write.table(metares,file=paste0(outprefix,".meta"),sep="\t",quote=F,row.names = F)
}
runmeta()
metares=as.data.frame(fread("../result/meta_lr.meta"))
#previous results
load("../result/compute_metapvalues_new.RData")
tmp=intersect(names(allpvalues),metares$SNP)
idx1=match(tmp,names(allpvalues))
idx2=match(tmp,metares$SNP)
plot(-log10(allpvalues[idx1]),-log10(metares$P[idx2]),xlab="-log10(ACAT p-value)",ylab="-log10(Logistic reg p-value")
abline(0,1,col="red")


metapvalues=metares$P
names(metapvalues)=metares$SNP

load("../result/metascoreinfo4_new.RData")

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
#table(rownames(allpvalues0) %in% names(metascoreinfo4$freq))

myplot=function(myallpvalues=allpvalues,outprefix="waldtest_meta")
{
  myallpvalues=myallpvalues[names(myallpvalues) %in% names(metascoreinfo4$freq)]
  idx0=match(names(myallpvalues),names(metascoreinfo4$freq))
  all(names(metascoreinfo4$freq[idx0])==names(myallpvalues)) #T
  
  tmp=unlist(strsplit(names(myallpvalues),":"))
  data=data.frame(SNP=names(myallpvalues),CHR=tmp[seq(1,length(tmp),4)],BP=tmp[seq(2,length(tmp),4)],EAF=metascoreinfo4$freq[idx0],P=myallpvalues,N=0)
  data$CHR=gsub("chr","",data$CHR)
  data$BP=as.integer(data$BP)
  data$CHR=as.integer(data$CHR)
  tmp1=1/(1/12427+1/12799)
  tmp2=1/(1/12799+1/60204)
  data$N=1/(1/tmp1+1/tmp2)
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

myplot(myallpvalues=metapvalues,outprefix="logisticregression")
