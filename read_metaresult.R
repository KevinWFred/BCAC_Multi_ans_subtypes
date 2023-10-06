#!/usr/bin/env Rscript
#after running compute_pvalues, collect the results
library(data.table)
read_metap=function(pop="euro")
{
  load(paste0("../result/",pop,"_beta_sigma_pvalue.RData")) #allres
  for (i in 1:1000)
  {
    if (i %% 100 ==0) cat(i,'..')
    load(paste0("../result/metaresult/",pop,"/res_",i,".RData"))
    res=allmetares$metares
    idx=match(rownames(res),rownames(allres))
    if (sum(is.na(idx))>0) warning(paste0(i," has NA SNPs"))
    allres[idx,31:45]=res
  }
  save(allres,file=paste0("../result/",pop,"_beta_sigma_pvalue.RData"))
#EURO
  # sum(!is.na(allres$icogs_beta2))
  # [1] 10338723
  # sum(!is.na(allres$onco_beta1))
  # [1] 10347049
  # sum(!is.na(allres$meta_beta1))
  # [1] 10149705
}


#draw plots
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

#For European
tmp=unlist(strsplit(rownames(allres),":"))
allres$CHR=tmp[seq(1,length(tmp),4)]
allres$CHR=gsub("chr","",allres$CHR)
allres$BP=tmp[seq(2,length(tmp),4)]
allres$SNP=rownames(allres)
N=1/(1/33505+1/33961)+1/(1/52500+1/64768)
tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro.afreq"))
allres$EAF=0.05
tmp1=intersect(allres$SNP,tmp$ID)
idx1=match(tmp1,tmp$ID)
idx2=match(tmp1,allres$SNP)
allres$EAF[idx2]=tmp$ALT_FREQS[idx1]
data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p1),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_EURO_meta1.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_EURO_meta1.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p2),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_EURO_meta2.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_EURO_meta2.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p3),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_EURO_meta3.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_EURO_meta3.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p4),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_EURO_meta4.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_EURO_meta4.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p5),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_EURO_meta5.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_EURO_meta5.png")

#For Asian
pop="asian"
load(paste0("../result/",pop,"_beta_sigma_pvalue.RData"))
tmp=unlist(strsplit(rownames(allres),":"))
allres$CHR=tmp[seq(1,length(tmp),4)]
allres$CHR=gsub("chr","",allres$CHR)
allres$BP=tmp[seq(2,length(tmp),4)]
allres$SNP=rownames(allres)
N=1/(1/5940+1/4882)+1/(1/12344+1/12870)
tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian.afreq"))
allres$EAF=0.05
tmp1=intersect(allres$SNP,tmp$ID)
idx1=match(tmp1,tmp$ID)
idx2=match(tmp1,allres$SNP)
allres$EAF[idx2]=tmp$ALT_FREQS[idx1]
data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p1),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_ASN_meta1.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_ASN_meta1.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p2),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_ASN_meta2.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_ASN_meta2.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p3),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_ASN_meta3.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_ASN_meta3.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p4),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_ASN_meta4.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_ASN_meta4.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(meta_p5),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.01 & data$FREQ_A1<0.99,]
png(filename = "../result/QQplot_ASN_meta5.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_ASN_meta5.png")

#African, only using onco data
pop="african"
load(paste0("../result/",pop,"_beta_sigma_pvalue.RData"))
tmp=unlist(strsplit(rownames(allres),":"))
allres$CHR=tmp[seq(1,length(tmp),4)]
allres$CHR=gsub("chr","",allres$CHR)
allres$BP=tmp[seq(2,length(tmp),4)]
allres$SNP=rownames(allres)
N=1/(1/2088+1/3481)
tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african.afreq"))
allres$EAF=0.05
tmp1=intersect(allres$SNP,tmp$ID)
idx1=match(tmp1,tmp$ID)
idx2=match(tmp1,allres$SNP)
allres$EAF[idx2]=tmp$ALT_FREQS[idx1]
data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(onco_p1),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.05 & data$FREQ_A1<0.95,]
png(filename = "../result/QQplot_AFR_meta1.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_AFR_meta1.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(onco_p2),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.05 & data$FREQ_A1<0.95,]
png(filename = "../result/QQplot_AFR_meta2.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_AFR_meta2.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(onco_p3),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.05 & data$FREQ_A1<0.95,]
png(filename = "../result/QQplot_AFR_meta3.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_AFR_meta3.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(onco_p4),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.05 & data$FREQ_A1<0.95,]
png(filename = "../result/QQplot_AFR_meta4.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_AFR_meta4.png")

data = allres %>% 
  mutate(CHR = as.integer(CHR),
         BP = as.integer(BP),
         FREQ_A1 = as.numeric(EAF),
         P = as.numeric(onco_p5),
         rsid = SNP,
         N = as.integer(N)) %>% 
  mutate(P = ifelse(P==0,1E-300,P)) %>% 
  select(rsid,CHR,BP,FREQ_A1,P,N)
data=data[data$FREQ_A1>0.05 & data$FREQ_A1<0.95,]
png(filename = "../result/QQplot_AFR_meta5.png", width = 8, height = 8, units = "in",res=300)
plotqq(data=data,optbreak=0,title="")
dev.off()
plotmanhattan(data=data,title="",filename="../result/man_AFR_meta5.png")



euro=as.data.frame(fread("../result/imp_icogs/euro/euro.pvar"))
asian=as.data.frame(fread("../result/imp_icogs/asian/asian.pvar"))
african=as.data.frame(fread("../result/imp_icogs/african/african.pvar"))

idx=which(!euro$ID %in% african$ID)
tmp=paste0(euro$`#CHROM`[idx],":",euro$POS[idx]) 
tmp1=paste0(african[,1],":",african$POS)
table(tmp %in% tmp1)
# FALSE    TRUE 
# 1971001   10949

table(euro$ID %in% asian$ID)
# FALSE    TRUE 
# 2028961 8310957
idx=which(!euro$ID %in% asian$ID)
tmp=paste0(euro$`#CHROM`[idx],":",euro$POS[idx]) 
tmp1=paste0(asian[,1],":",asian$POS)
table(tmp %in% tmp1)
# FALSE    TRUE 
# 2024497    4464

#generate files for meta-analysis
for (pop in c("euro","asian","african"))
{
  print(pop)
  load(paste0("../result/",pop,"_beta_sigma_pvalue.RData"))
  tmp=unlist(strsplit(rownames(allres),":"))
  #allres$CHR=tmp[seq(1,length(tmp),4)]
  #allres$CHR=gsub("chr","",allres$CHR)
  #allres$BP=tmp[seq(2,length(tmp),4)]
  allres$SNP=rownames(allres)
  allres$A1=tmp[seq(3,length(tmp),4)]
  allres$A2=tmp[seq(4,length(tmp),4)]
  if (pop=="euro")
  {
    N=1/(1/33505+1/33961)+1/(1/52500+1/64768)
    tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro.afreq"))
    allres$EAF=0.05
    tmp1=intersect(allres$SNP,tmp$ID)
    idx1=match(tmp1,tmp$ID)
    idx2=match(tmp1,allres$SNP)
    allres$EAF[idx2]=tmp$ALT_FREQS[idx1]
  }
  if (pop=="asian")
  {
    N=1/(1/5940+1/4882)+1/(1/12344+1/12870)
    tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian.afreq"))
    allres$EAF=0.05
    tmp1=intersect(allres$SNP,tmp$ID)
    idx1=match(tmp1,tmp$ID)
    idx2=match(tmp1,allres$SNP)
    allres$EAF[idx2]=tmp$ALT_FREQS[idx1]
  }
  if (pop=="african")
  {
    N=1/(1/2088+1/3481)
    tmp=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african.afreq"))
    allres$EAF=0.05
    tmp1=intersect(allres$SNP,tmp$ID)
    idx1=match(tmp1,tmp$ID)
    idx2=match(tmp1,allres$SNP)
    allres$EAF[idx2]=tmp$ALT_FREQS[idx1]
  }
  allres$N=N
  for (i in 1:5)
  {
    if (pop!="african")
    {
      res=data.frame(SNP=allres$SNP,A1=allres$A1,A2=allres$A2,EAF=allres$EAF,N=allres$N,SE=allres[,paste0("meta_se",i)],Beta=allres[,paste0("meta_beta",i)],P=allres[,paste0("meta_p",i)])
    }else
    {
      res=data.frame(SNP=allres$SNP,A1=allres$A1,A2=allres$A2,EAF=allres$EAF,N=allres$N,SE=allres[,paste0("onco_se",i)],Beta=allres[,paste0("onco_beta",i)],P=allres[,paste0("onco_p",i)])
    }
    
    fwrite(res,file=paste0("../result/",pop,"_type",i,"_metalinput.txt"),row.names = F,sep="\t",quote=F)
  }
}

metainput_euro=as.data.frame(fread(paste0("../result/","euro","_type",1,"_metalinput.txt")))
metainput_asian=as.data.frame(fread(paste0("../result/","asian","_type",1,"_metalinput.txt")))
metainput_african=as.data.frame(fread(paste0("../result/","african","_type",1,"_metalinput.txt")))

plot_meta=function(type=1)
{
  data=as.data.frame(fread(paste0("/data/BB_Bioinformatics/Kevin/BCAC/result/meta_type",type,"1.tbl")))
  tmp=unlist(strsplit(data$MarkerName,":"))
  #HetDf = the number of degrees of freedom (i.e. non-missing studies minus one) in the meta-analysis for each variant.
  data$CHR=gsub("chr","",tmp[seq(1,length(tmp),4)])
  data$BP=tmp[seq(2,length(tmp),4)]
  data$N=0
  idx=which(data$MarkerName %in% metainput_euro$SNP)
  data$N[idx]=data$N[idx]+metainput_euro$N[1]
  idx=which(data$MarkerName %in% metainput_asian$SNP)
  data$N[idx]=data$N[idx]+metainput_asian$N[1]
  idx=which(data$MarkerName %in% metainput_african$SNP)
  data$N[idx]=data$N[idx]+metainput_african$N[1]
  data = data %>% 
    mutate(CHR = as.integer(CHR),
           BP = as.integer(BP),
           FREQ_A1 = as.numeric(Freq1),
           P = as.numeric(`P-value`),
           rsid = MarkerName,
           N = as.integer(N)) %>% 
    mutate(P = ifelse(P==0,1E-300,P)) %>% 
    select(rsid,CHR,BP,FREQ_A1,P,N)
  png(filename = paste0("../result/QQplot_meta_type",type,".png"), width = 8, height = 8, units = "in",res=300)
  plotqq(data=data,optbreak=0,title="")
  dev.off()
  plotmanhattan(data=data,title="",filename=paste0("../result/man_meta_type",type,".png"))
}
plot_meta()
plot_meta(type=2)
plot_meta(type=3)
plot_meta(type=4)
plot_meta(type=5)
