#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(data.table)
library(bc2)
library(bcutility)

startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}

readres=function(pop="euro",dataopt="onco",nvar=2000)
{
  scorefolder=paste0("../result/imp_",dataopt,"/",pop,"/scoretest")
  outfile=paste0("../result/imp_",dataopt,"/",pop,"/scoretest/scoreinfor.RData")
  pvar=as.data.frame(fread(paste0("../result/imp_",dataopt,"/",pop,"/",pop,".pvar")))
  nblock=as.integer(nrow(pvar)/nvar)+1
  score=infor=freq=NULL
  for (i in 1:nblock)
  {
    if (i %% 500 ==0) cat(i,'..')
    resfile=paste0(scorefolder,"/res_",i,".RData")
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
        score=rbind(score,score_result)
        infor=rbind(infor,infor_result)
        freq=c(freq,freq.all)
        rm(freq.all,score_result,infor_result)
      }else
      {
        warnings(paste0(i," block can't load!"))
      }
    }else
    {
      warnings(paste0(i," block is missing!"))
    }
  }
  res=list(freq=freq,score=score,infor=infor)
  save(res,file=outfile)
  return(res)
}
# euro_onco=readres()
# euro_icogs=readres(pop="euro",dataopt="icogs",nvar=2000)
# asian_onco=readres(pop="asian",dataopt="onco",nvar=2000)
# asian_icogs=readres(pop="asian",dataopt="icogs",nvar=2000)
# african_onco=readres(pop="african",dataopt="onco",nvar=2000)

#read and combine the score/infor
load(paste0("../result/imp_","onco","/","euro","/scoretest/scoreinfor.RData"))
euro_onco=res
rm(res)
load(paste0("../result/imp_","icogs","/","euro","/scoretest/scoreinfor.RData"))
euro_icogs=res
rm(res)
load(paste0("../result/imp_","onco","/","asian","/scoretest/scoreinfor.RData"))
asian_onco=res
rm(res)
load(paste0("../result/imp_","icogs","/","asian","/scoretest/scoreinfor.RData"))
asian_icogs=res
rm(res)
load(paste0("../result/imp_","onco","/","african","/scoretest/scoreinfor.RData"))
african_onco=res
rm(res)
#first time opt="new", then opt="add". keep availability of variants in each data
aggregate_scoreinfor=function(res1=euro_onco,
                              res2=euro_icogs,n.second=5,opt="new",oldavail=NULL)
{
  #restrict to SNPs with allele frequency more than 0.006
  idx=which(res1$freq>0.006 & res1$freq < 0.994)
  res1$score=res1$score[idx,]
  res1$infor=res1$infor[idx,]
  res1$freq=res1$freq[idx]
  idx=which(res2$freq>0.006 & res2$freq < 0.994)
  res2$score=res2$score[idx,]
  res2$infor=res2$infor[idx,]
  res2$freq=res2$freq[idx]

  allvar=unique(c(rownames(res1$score),rownames(res2$score)))
  print(length(allvar))
  var1=rownames(res1$score)
  var2=rownames(res2$score)
  score=data.frame(matrix(0,nrow=length(allvar),ncol=n.second))
  infor=data.frame(matrix(0,nrow=length(allvar),ncol=n.second^2))
  rownames(score)=rownames(infor)=allvar
  idx1=match(var1,allvar)
  score[idx1,]=score[idx1,]+res1$score
  infor[idx1,]=infor[idx1,]+res1$infor
  idx2=match(var2,allvar)
  score[idx2,]=score[idx2,]+res2$score
  infor[idx2,]=infor[idx2,]+res2$infor
  #freq is minimum freq
  freq=maxfreq=rep(0,length(allvar))
  names(freq)=names(maxfreq)=allvar
  var1only=var1[!var1 %in% var2]
  var2only=var2[!var2 %in% var1]
  idx3=match(var1only,allvar)
  idx4=match(var1only,names(res1$freq))
  freq[idx3]=maxfreq[idx3]=res1$freq[idx4]
  idx3=match(var2only,allvar)
  idx4=match(var2only,names(res2$freq))
  freq[idx3]=maxfreq[idx3]=res2$freq[idx4]
  var12=intersect(var1,var2)
  idx3=match(var12,allvar)
  idx4=match(var12,names(res1$freq))
  idx5=match(var12,names(res2$freq))
  freq[idx3]=pmin(res1$freq[idx4],res2$freq[idx5])
  maxfreq[idx3]=pmax(res1$freq[idx4],res2$freq[idx5])
  if (opt=="new")
  {
    avail=data.frame(matrix(0,nrow=length(allvar),ncol=2))
    rownames(avail)=allvar
    avail[idx1,1]=1
    avail[idx2,2]=1
  }else
  {
    avail=data.frame(matrix(0,nrow=length(allvar),ncol=ncol(oldavail)+1))
    rownames(avail)=allvar
    avail[idx1,1:ncol(oldavail)]=oldavail
    avail[idx2,ncol(avail)]=1
  }
  return(list(score=score,infor=infor,freq=freq,maxfreq=maxfreq,avail=avail))
}
#
metascoreinfo1=aggregate_scoreinfor()
metascoreinfo2=aggregate_scoreinfor(res1=metascoreinfo1,res2=asian_icogs,
                                    opt="add",oldavail=metascoreinfo1$avail)
metascoreinfo3=aggregate_scoreinfor(res1=metascoreinfo2,res2=asian_onco,
                                    opt="add",oldavail=metascoreinfo2$avail)
metascoreinfo4=aggregate_scoreinfor(res1=metascoreinfo3,res2=african_onco,
                                    opt="add",oldavail=metascoreinfo3$avail)
# calculate the p-value for based on score statistics and information matrix
# to perform fixed-effect test using chi-square test
# we can use function ScoreGlobalTestForAssoc(score, infor)
# to perform random-effect test using mixture chi-square test
# we can use function ScoreMixedGlobalTestForHeter(score, infor)
# table(metascoreinfo4$freq>0.01 & metascoreinfo4$freq<0.99)
# FALSE     TRUE
# 3870641 14994842
save(metascoreinfo4,file="../result/metascoreinfo4_new.RData") #new is the one remove MAF<0.006
load("../result/metascoreinfo4_new.RData")

#compute p-values based on score and infor
Pfunction <- function(score=metascoreinfo4$score[1,],infor=metascoreinfo4$infor[1,],second.num=5){
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


print(paste0("host:",as.character(Sys.info()["nodename"])))
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
print(Sys.time())
print(args)
i1 = as.numeric(args[1]) #block ID, starts with 1
ntask=1000 #total number of jobs
#inputfile="../result/euro_beta_sigma.RData"

allvar=rownames(metascoreinfo4$score)
print(length(allvar))
nvar=length(allvar)
startendidx=startend(nvar,ntask,i1)
print(paste0("startendidx: ",startendidx))
tmp=sum( sapply(ls(),function(x){object.size(get(x))}))/1E6
print(paste0("memory size: ",tmp))

start <- startendidx[1]
end <- startendidx[2]
#count the p_value
pvalue_sub <- rep(0,end-start+1)
names(pvalue_sub)=allvar[start:end]
temp = 1

for(j in start:end){
  if (temp %% 1000 ==0) print(j)

  pvalue_sub[temp] <- Pfunction(score=metascoreinfo4$score[j,],infor=metascoreinfo4$infor[j,])
  temp = temp+1
}


#save the pvalue_sub file to a folder
save(pvalue_sub,file=paste0("../result/scoretestresult/p_value_sub",i1,".RData"))

print(Sys.time())
print("done")

createswarm=function()
{
  swarmjobs=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_pvalues_scoretest.R",1000),
                       i1=1:1000)
  write.table(swarmjobs,file=paste0("compute_scoretest_pvalues.swam"),row.names=F,col.names=F,sep="\t",quote=F)
}
#createswarm()
#9226087
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_scoretest_pvalues.swam -g 30 --module R/4.3 --partition=quick --time=4:00:00 --gres=lscratch:30 -p 2

#collect all the results
allpvalues=NULL
for (i1 in 1:1000)
{
  if (i1 %% 100==0) cat(i1,'..')
#   load(paste0("../result/scoretestresult/p_value_sub",i1,".Rdata"))
  load(paste0("../result/scoretestresult/p_value_sub",i1,".RData"))
  allpvalues=c(allpvalues,pvalue_sub)
}
save(allpvalues,file="../result/scoretestresult/allpvalues.Rdata")
# 
# #qqplot and manhattan plot
# #draw plots
# library(data.table)
# library("plotrix") #axis.break
# library("RColorBrewer")
# library("optparse")
# library(readr)
# library(dplyr)
# library(ggplot2)
# library(qqman)
# 
# source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
# #draw QQ plot
# qqplotdata <- function(logpvector){
#   o = sort(logpvector,decreasing=T)
#   e = -log10(ppoints(length(o)))       
#   qqdata <- data.frame(o,e)
#   qqdata$o <- round(qqdata$o,3)
#   qqdata$e <- round(qqdata$e,3)
#   keepU <- which(!duplicated(qqdata))
#   qqdata <- qqdata[keepU,]
#   
#   N <- length(logpvector) ## number of p-values
#   ## create the confidence intervals
#   qqdata$c975 <- NA
#   qqdata$c025 <- NA
#   
#   ## the jth order statistic from a
#   ## uniform(0,1) sample
#   ## has a beta(j,n-j+1) distribution
#   ## (Casella & Berger, 2002,
#   ## 2nd edition, pg 230, Duxbury)
#   
#   for(i in 1:length(keepU)){
#     j <- keepU[i]
#     qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
#     qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
#   }
#   return(qqdata)
# }
# #optbreak=1, break top p-values, used for very low p-value cases
# plotqq=function(data,optbreak=1,title="")
# {
#   dat = data %>% 
#     mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>% 
#     select(rsid,CHR,BP,P,MAF) %>% 
#     rename(SNP = rsid)
#   dat=dat[!is.na(dat$P),]
#   x = dat$P
#   z = qnorm(x / 2)
#   lambda = round(median(z^2) / qchisq(0.5,1), 3)
#   N.effect = median(data$N)
#   lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)
#   
#   yLine <- c(-log10(5E-8))
#   colLine <- c("red")
#   dat$log10P = -log10(dat$P)
#   gwas = as.data.frame(dat)
#   # Determine frequency bins and create variable for binned QQ plot
#   
#   minMAF <- min(gwas$MAF)
#   
#   freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
#   gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
#   freqtable <- table(gwas$freqbin)
#   freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
#   freqtable <- freqtable[freqtable > 0]
#   
#   ## Generate QQ plot data by frequency bin
#   fbin <- character(0)
#   fN <- integer(0)
#   fx <- numeric(0)
#   fy <- numeric(0)
#   fcol <- character(0)
#   legendcol <- character(0)
#   conf <- list()
#   allcols <- brewer.pal(4,"Set1")
#   ycol <- "log10P"
#   for(f in 1:length(freqtable)){
#     fbin <- c(fbin,names(freqtable)[f])
#     fsnps <- which(gwas$freqbin ==names(freqtable)[f])
#     plotdata <- qqplotdata(gwas[[ycol]][fsnps])
#     fN <- c(fN,freqtable[f])
#     fx <- c(fx,plotdata$e)
#     fy <- c(fy,plotdata$o)
#     fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
#     conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
#                             'y'=c(plotdata$c975,rev(plotdata$c025)))
#     legendcol <- c(legendcol,allcols[f])
#   }
#   legendtext <- paste0("MAF=",fbin,"; #SNPs=",format(fN,big.mark=",",scientific=FALSE))
#   
#   opt =  list(break.top = ifelse(optbreak==1,15,ceiling(max(fy))+1),
#               top.size = 0.125)
#   
#   
#   #png(filename = paste0(outpath,"/QQ_",eth[i1],"_",trait[i2],".png"), width = 8, height = 8, units = "in",res=300)
#   xlim <- c(0,max(fx,na.rm=T))
#   ylim <- c(0,max(fy,na.rm=T))
#   maxY <- max(fy,na.rm=T)
#   
#   par(mar=c(5.1,5.1,4.1,1.1))
#   
#   lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
#   lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
#   #top
#   lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
#   if (any(lab2>max(lab1))) lab2 <- lab2[lab2 > max(lab1)]
#   
#   # resulting range of top scale in bottom scale units
#   top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
#   top.data = max(lab2)-opt$break.top
#   
#   # function to rescale the top part
#   rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
#   
#   plot(0,0,
#        ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
#        xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
#        ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
#        cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
#        main=opt$maintitle,pch=19)
#   
#   # Plot confidence intervals	
#   for(p in 1:length(conf)){
#     polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
#             col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
#             border = NA)
#   }
#   
#   # add points below top
#   points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
#   
#   # identify line & add axis break
#   lines(xlim,xlim,col="black",lty = 2)
#   axis(1,cex.axis=1.5,cex.lab=1.5)
#   par(las=1)
#   axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
#   box()
#   if (optbreak==1)
#   {
#     rescaled.y = rescale(fy[fy>opt$break.top])
#     par(las=0)
#     points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
#     par(las=1)
#     axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
#     axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
#     axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
#     lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
#   }
#   
#   abline(h=ifelse(yLine<opt$break.top,
#                   yLine,
#                   rescale(yLine)),
#          col=colLine,lwd=1.5,lty=2)
#   legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
#   text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
#   text(5.9,1,paste(lambda_1000),cex = 1.5)
#   
#   title(title)
#   #dev.off()
#   
# }
# 
# 
# plotmanhattan=function(data,title="",filename="../result/man_ILCCOnever.png")
# {
#   dat = data %>% 
#     mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>% 
#     select(rsid,CHR,BP,P,MAF) %>% 
#     rename(SNP = rsid)
#   dat$CHR1=factor(dat$CHR,levels = c(1:22,"X","Y"))
#   idx=order(dat$CHR1,dat$BP)
#   dat=dat[idx,]
#   p.pwas <- 5E-08
#   
#   nCHR <- length(unique(dat$CHR))
#   dat$BPcum <- NA
#   s <- 0
#   nbp <- c()
#   for (i in unique(dat$CHR)){
#     nbp[i] <- max(dat[dat$CHR == i,]$BP)
#     dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
#     s <- s + nbp[i]
#   }
#   
#   axis.set <- dat %>% 
#     group_by(CHR) %>% 
#     summarize(center = (max(BPcum) + min(BPcum)) / 2)
#   ylim <- abs(floor(log10(min(dat$P)))) + 2 
#   sig1 <- p.pwas
#   
#   
#   sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
#   library(ggplot2)
#   manhplot <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
#                               color = as.factor(CHR), size = -log10(P))) +
#     geom_point(alpha = 0.8, size=0.8) + 
#     scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
#     scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
#     scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
#     scale_size_continuous(range = c(0.5,3)) +
#     geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
#     guides(color = FALSE) + 
#     labs(x = NULL, 
#          y = "-log10(p)", 
#          linetype = "",
#          title = title)+
#     #subtitle = "A2: Critically ill COVID19+ vs. population controls;\nB1: Hospitalized COVID19+ vs non-hospitalized COVID19+;\nB2: Hospitalized COVID19+ vs. population controls;\nC2: Reported SARS-CoV-2 infection vs. population controls") + 
#     theme_Publication()+
#     theme(
#       legend.position = "top",
#       panel.border = element_blank(),
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
#       plot.title = element_text(size = 12, face = "bold"),
#       plot.subtitle = element_text(size = 8)
#     )
#   
#   #outpath <-"../result"  
#   
#   
#   ggsave(filename=filename,
#          plot=manhplot, device="png",
#          width=9, height=4, units="in", dpi=300)
# }
# load("../result/scoretestresult/allpvalues.Rdata")
idx=match(names(allpvalues),rownames(metascoreinfo4$score))
metascoreinfo4$score=metascoreinfo4$score[idx,]
metascoreinfo4$infor=metascoreinfo4$infor[idx,]
metascoreinfo4$freq=metascoreinfo4$freq[idx]
metascoreinfo4$avail=metascoreinfo4$avail[idx,]
# all(names(metascoreinfo4$freq)==names(allpvalues)) #T
# idx=which(is.na(metascoreinfo4$freq))
# quantile(allpvalues[idx])
# tmp=unlist(strsplit(names(allpvalues),":"))
# data=data.frame(SNP=names(allpvalues),CHR=tmp[seq(1,length(tmp),4)],BP=tmp[seq(2,length(tmp),4)],EAF=metascoreinfo4$freq,P=allpvalues,N=0)
# data$CHR=gsub("chr","",data$CHR)
# data$BP=as.integer(data$BP)
# data$CHR=as.integer(data$CHR)
# #euronco
# idx=which(metascoreinfo4$avail[,1]==1)
# data$N[idx]=data$N[idx]+1/(1/34367+1/52500)
# #euricogs
# idx=which(metascoreinfo4$avail[,2]==1)
# data$N[idx]=data$N[idx]+1/(1/13283+1/33505)
# #asianicogs
# idx=which(metascoreinfo4$avail[,3]==1)
# data$N[idx]=data$N[idx]+1/(1/1682+1/5940)
# idx=which(metascoreinfo4$avail[,4]==1)
# data$N[idx]=data$N[idx]+1/(1/9064+1/12344)
# idx=which(metascoreinfo4$avail[,5]==1)
# data$N[idx]=data$N[idx]+1/(1/1815+1/2088)
# data$EAF[is.na(data$EAF)]=0.01
# 
# data1 = data %>% 
#   mutate(CHR = as.integer(CHR),
#          BP = as.integer(BP),
#          FREQ_A1 = as.numeric(EAF),
#          P = as.numeric(P),
#          rsid = SNP,
#          N = as.integer(N)) %>% 
#   mutate(P = ifelse(P==0,1E-300,P)) %>% 
#   select(rsid,CHR,BP,FREQ_A1,P,N)
# data1=data1[data1$FREQ_A1>0.01 & data1$FREQ_A1<0.99,]
# png(filename = "../result/QQplot_scoretest_meta.png", width = 8, height = 8, units = "in",res=300)
# plotqq(data=data1,optbreak=1,title="")
# dev.off()
# plotmanhattan(data=data1,title="",filename="../result/man_scoretest_meta.png")

